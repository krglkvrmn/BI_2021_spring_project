import argparse
import re
import os
from collections import defaultdict
from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_files", nargs="+")
parser.add_argument("-o", "--output_dir", default="HOGs")
parser.add_argument("-g", "--group", default="HOG", choices=["HOG", "OMA"])
parser.add_argument("--oma-groups")
parser.add_argument("--full_HOG", default=False, action="store_true")
parser.add_argument("--incl-nonPP", default=False, action="store_true", help="Include groups withoud polyprolines in output")
args = parser.parse_args()


X1_strength = {"A": "S", "D": "S", "E": "S", "G": "S", "H": "S", "I": "S",
               "K": "S", "P": "S", "Q": "S", "R": "S", "V": "S", "W": "S",
               "Y": "S", "F": "M", "M": "M", "N": "M", "C": "W", "L": "W",
               "S": "W", "T": "W", "X": "W"}

X2_strength = {"A": "S", "D": "S", "G": "S", "P": "S", "S": "M", "E": "M",
               "C": "W", "F": "W", "H": "W", "I": "W", "K": "W", "L": "W",
               "M": "W", "N": "W", "Q": "W", "R": "W", "T": "W", "V": "W",
               "W": "W", "Y": "W", "X": "W"}

X3_strength = {"D": "S", "E": "S", "G": "S", "N": "S", "P": "S", "S": "S",
               "W": "S", "A": "M", "K": "M", "Q": "M", "T": "M", "C": "W",
               "F": "W", "H": "W", "I": "W", "L": "W", "M": "W", "R": "W",
               "V": "W", "Y": "W", "X": "W"}

# N3_RULE[x1][x3]
N3_RULE = {"S": {"S": "S", "M": "S", "W": "S"},
           "M": {"S": "S", "M": "M", "W": "M"},
           "W": {"S": "S", "M": "M", "W": "W"}}

# N2_RULE[x3][x2][x1]
N2_RULE = {"W": {"S": {"S": "S", "M": "M", "W": "W"},
                 "M": {"S": "M", "M": "M", "W": "W"},
                 "W": {"S": "W", "M": "W", "W": "W"}},
           "M": {"S": {"S": "S", "M": "M", "W": "W"},
                 "M": {"S": "M", "M": "M", "W": "M"},
                 "W": {"S": "M", "M": "M", "W": "W"}},
           "S": {"S": {"S": "S", "M": "S", "W": "S"},
                 "M": {"S": "M", "M": "M", "W": "M"},
                 "W": {"S": "W", "M": "W", "W": "W"}}}

OUTPUT_DIR = args.output_dir


def setup():
    if OUTPUT_DIR in os.listdir():
        for filename in os.listdir(OUTPUT_DIR):
            os.remove(os.path.join(OUTPUT_DIR, filename))
    else:
        os.mkdir(OUTPUT_DIR)


def find_polyproline_motifs(sequence):
    matches = []
    for match in re.finditer(r"(?=(\w\wPP[^P]))", sequence):
        matches.append((match.group(1), match.start(1), match.end(1)))
    return matches


def classify_motif(motif):
    x1, x2, x3 = motif[0], motif[1], motif[-1]
    x1s = X1_strength[x1]
    x2s = X2_strength[x2]
    x3s = X3_strength[x3]
    if x1 == x2 == "P":
        return "S"
    elif x1 != "P" and x2 == "P":
        return N3_RULE[x1s][x3s]
    else:
        return N2_RULE[x3s][x2s][x1s]


def write_into_hog(sequence, hog, taxid, motifs, prot_id):
    hog_filename = hog.replace(":", "")
    writable_motifs = map(lambda x: f"{x[0]}:{x[1]}:{x[2]}:{x[3]}", motifs)
    with open(f"{OUTPUT_DIR}/{hog_filename}.fasta", "a") as out_file:
        description = f">{prot_id} | {taxid} | {hog} | {' '.join(writable_motifs)}\n"
        out_file.write(description)
        out_file.write(f"{sequence}\n")


def make_oma_index(path):
    # 640 MB size
    oma_index = {}
    with open(path) as file_in:
        for line in file_in:
            if not line.startswith("#"):
                oma_group, *oma_entries = line.strip().split()
                for entry in oma_entries:
                    oma_index[entry] = oma_group
    return oma_index


def extract_hogs(path, args):
    if args.group == "HOG" and args.full_HOG:
        group_pattern = re.compile(r"(HOG:.+?)\s")
    elif args.group == "HOG" and not args.full_HOG:
        group_pattern = re.compile(r"(HOG:\d+)")
    elif args.group == "OMA":
        oma_index = make_oma_index(args.oma_groups)
        group_pattern = re.compile(r"(HOG:\d+)")
    hog_motif_status = defaultdict(bool)
    for idx, record in enumerate(SeqIO.parse(path, "fasta")):
        # r"(HOG:.+?)\s" for same IDs as in initial files
        # r"(HOG:\d+)" for short IDs 
        if hog_search := re.search(group_pattern, record.description):
            sequence = str(record.seq)
            strain = path.split("/")[-1].split(".")[0]
            prot_id = record.description.split(" | ")[0]
            if args.group == "HOG":
                hog = hog_search.group(1)
            elif args.group == "OMA":
                try:
                    hog = oma_index[prot_id]
                except KeyError:
                    print(f"Protein {prot_id} has no OMA group")
            matches = find_polyproline_motifs(sequence)
            motifs = []
            for motif, start, end in matches:
                motifs.append((motif, start, end, classify_motif(motif)))
            hog_motif_status[hog] |= bool(matches)
            write_into_hog(sequence=sequence, hog=hog,
                           taxid=strain, motifs=motifs,
                           prot_id=prot_id)
    if not args.incl_nonPP:
        for hog, status in hog_motif_status.items():
            if not status:
                hog_filename = hog.replace(":", "") + ".fasta"
                os.remove(os.path.join(OUTPUT_DIR, hog_filename))


if __name__ == "__main__":
    setup()
    for filename in args.input_files:
        extract_hogs(filename, args)
