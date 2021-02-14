import argparse
import re
import os
from collections import defaultdict
from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_files", nargs="+")
parser.add_argument("-o", "--output_dir", default="HOGs")
args = parser.parse_args()


X1_strength = {"A": "S", "D": "S", "E": "S", "G": "S", "H": "S", "I": "S",
               "K": "S", "P": "S", "Q": "S", "R": "S", "V": "S", "W": "S",
               "Y": "S", "F": "M", "M": "M", "N": "M", "C": "W", "L": "W",
               "S": "W", "T": "W"}

X2_strength = {"A": "S", "D": "S", "G": "S", "P": "S", "S": "M", "E": "M",
               "C": "W", "F": "W", "H": "W", "I": "W", "K": "W", "L": "W",
               "M": "W", "N": "W", "Q": "W", "R": "W", "T": "W", "V": "W",
               "W": "W", "Y": "W"}

X3_strength = {"D": "S", "E": "S", "G": "S", "N": "S", "P": "S", "S": "S",
               "W": "S", "A": "M", "K": "M", "Q": "M", "T": "M", "C": "W",
               "F": "W", "H": "W", "I": "W", "L": "W", "M": "W", "R": "W",
               "V": "W", "Y": "W"}

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


def write_into_hog(sequence, hog, taxid, motifs):
    hog_filename = hog.replace(":", "")
    writable_motifs = map(lambda x: f"{x[0]}:{x[1]}:{x[2]}:{x[3]}", motifs)
    with open(f"{OUTPUT_DIR}/{hog_filename}.fasta", "a") as out_file:
        description = f">{taxid} | {hog} | {' '.join(writable_motifs)}\n"
        out_file.write(description)
        out_file.write(f"{sequence}\n")


def extract_hogs(path):
    hog_motif_status = defaultdict(bool)
    for idx, record in enumerate(SeqIO.parse(path, "fasta")):
        if hog_search := re.search(r"HOG:\d+", record.description):
            sequence = str(record.seq)
            strain = record.description.split(" | ")[0]
            hog = hog_search.group()
            matches = find_polyproline_motifs(sequence)
            motifs = []
            for motif, start, end in matches:
                motifs.append((motif, start, end, classify_motif(motif)))
            hog_motif_status[hog] |= bool(matches)
            write_into_hog(sequence, hog, strain, motifs)
    for hog, status in hog_motif_status.items():
        if not status:
            hog_filename = hog.replace(":", "") + ".fasta"
            os.remove(os.path.join(OUTPUT_DIR, hog_filename))

if __name__ == "__main__":
    setup()
    for filename in args.input_files:
        extract_hogs(filename)
