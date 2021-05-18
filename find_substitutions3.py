import argparse
from multiprocessing import Process
import math
import os
import re
import shutil
from itertools import combinations
from functools import reduce
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import numpy as np
import pandas as pd


def all_best_paiwise_alignments(sequences):
    alignments = [max(pairwise2.align.globalds(*pair, matlist.blosum62, -10, -0.5), key=lambda x: x[2]) for pair in combinations(sequences, 2)]
    return alignments, list(combinations(range(len(sequences)), 2))


def find_aligned_pp(aligned_seq):
    motif_candidates = []
    for i, _ in enumerate(aligned_seq):
        match = re.match(r"(P-*P)[^P-]?", aligned_seq[i:])
        if match:
            motif_candidates.append((match.group(1), i, match.start(), match.end() - 1))
    return motif_candidates

def complex_slice(string, coords, exclude=False):
    fragment = string
    if exclude:
        for coord in coords:
            if exclude:
                fragment = fragment[:coord[0]] + "@" * (coord[1] - coord[0]) + fragment[coord[1]:]
    else:
        fragment = "".join(map(lambda x: fragment[slice(*x)], coords))
    return fragment.replace("@", "")

def substitutions(seq1, seq2):
    subst = pd.Series(index=list("ACDEFGHIKLMNPQRSTVWYX-"), dtype=int)
    for s1, s2 in zip(seq1, seq2):
        if s1 == "P":
            subst[s2] += 1
        elif s2 == "P":
            subst[s1] += 1
    return subst

def substitutions2(seq1, seq2):
    motif_subst = pd.Series(index=list("ACDEFGHIKLMNPQRSTVWYX-"), dtype=int)
    single_subst = pd.Series(index=list("ACDEFGHIKLMNPQRSTVWYX-"), dtype=int)
    
    try:
        seq1_motifs = reduce(lambda x, y: x | y, [set(range(start, start+end+1)) for _, start, _, end in find_aligned_pp(seq1)])
    except TypeError:
        seq1_motifs = set()
    try:
        seq2_motifs = reduce(lambda x, y: x | y, [set(range(start, start+end+1)) for _, start, _, end in find_aligned_pp(seq2)])
    except TypeError:
        seq2_motifs = set()
    for idx, aa in enumerate(seq1):
        if aa == "P" and idx in seq1_motifs:
            motif_subst[seq2[idx]] += 1
        elif aa == "P" and idx not in seq1_motifs:
            single_subst[seq2[idx]] += 1
    
    for idx, aa in enumerate(seq2):
        if aa == "P" and idx in seq2_motifs:
            motif_subst[seq1[idx]] += 1
        elif aa == "P" and idx not in seq2_motifs:
            single_subst[seq1[idx]] += 1
            
    return motif_subst, single_subst

def parse_description_and_motifs(description):
    desc_splitted = description.split(" | ")
    if len(desc_splitted) > 3:
        all_motifs_joined = desc_splitted[3]
        all_motifs_splitted = all_motifs_joined.split(" ")
        motifs = list(map(lambda k: (k[0], int(k[1]), int(k[2]), k[3]), map(lambda x: x.split(":"), all_motifs_splitted)))
        return desc_splitted[0], motifs
    else:
        return desc_splitted[0], None

def main(base_dir, file_set, worker_id, file_out):
    for file_num, oma_group_file in enumerate(file_set):
        print(f"Worker {worker_id}: {file_num+1}/{len(file_set)}")
        oma_group_id = oma_group_file.split(".")[0]
        oma_group_records = list(SeqIO.parse(os.path.join(base_dir, oma_group_file), "fasta"))
        descriptions, sequences = list(map(lambda x: x.description, oma_group_records)), list(map(lambda x: x.seq, oma_group_records))
        alignments, idc = all_best_paiwise_alignments(sequences)

        aligned_motifs = []
        aligned_no_motifs = []

        for alignment, seqs in zip(alignments, idc):
            if seqs[0] == seqs[1]:
                continue
            seqA_id, motifs1 = parse_description_and_motifs(descriptions[seqs[0]])
            seqB_id, motifs2 = parse_description_and_motifs(descriptions[seqs[1]])
            motif_pairs = []
            motif_idxs = []

            seqA_pp = find_aligned_pp(alignment[0])
            for ppa_match in seqA_pp:
                motif_seq, pos, _, end = ppa_match
                start, end = pos, pos + end
                coords = (start, end)
                if coords not in motif_idxs:
                    motif_pairs.append({seqA_id: alignment[0][start:end], seqB_id: alignment[1][start:end]})
                    motif_idxs.append(coords)

            seqB_pp = find_aligned_pp(alignment[1])
            for ppb_match in seqB_pp:
                motif_seq, pos, _, end = ppb_match
                start, end = pos, pos + end
                coords = (start, end)
                if coords not in motif_idxs:
                    motif_pairs.append({seqA_id: alignment[0][start:end], seqB_id: alignment[1][start:end]})
                    motif_idxs.append(coords)

            if motif_idxs:
                seqA_no_motifs = complex_slice(alignment[0], motif_idxs, exclude=True)
                seqB_no_motifs = complex_slice(alignment[1], motif_idxs, exclude=True)
            else:
                seqA_no_motifs = alignment[0]
                seqB_no_motifs = alignment[1]
            aligned_motifs.append(motif_pairs)
            aligned_no_motifs.append({seqA_id: seqA_no_motifs, seqB_id: seqB_no_motifs})
        for seq_related_motifs in aligned_motifs:
            for motif_id, pair in enumerate(seq_related_motifs):
                motif_seqs_no_gaps = list(map(lambda x: x.replace("-", ""), pair.values()))
                motif_seqs = list(pair.values())
                if "PP" not in motif_seqs_no_gaps[0] or "PP" not in motif_seqs_no_gaps[1]:
                    aligned_no_motifs.append(pair)
                if len(motif_seqs[0]) == 2 and len(motif_seqs[1]) == 2:   # No indels only subst
                    a, _ = substitutions2(*pair.values())
                    subst = ",".join((map(lambda x: f"{x[0]}|{x[1]}", filter(lambda y: y[1] > 0, zip(a.index, a.values)))))
                    print(oma_group_id, *pair.keys(), 2, motif_id + 1, "m", "s", subst, sep="\t", file=file_out)
                    
                elif len(motif_seqs[0]) > 2 and len(motif_seqs[1]) > 2 and motif_seqs[0][0] == motif_seqs[1][0] and motif_seqs[0][-1] == motif_seqs[1][-1]: 
                    if "PP" in motif_seqs_no_gaps[0]:
                        print(oma_group_id, *pair.keys(), 3, motif_id + 1, "m", "i", motif_seqs[1][1:-1], sep="\t", file=file_out)
                    elif "PP" in motif_seqs_no_gaps[1]:
                        print(oma_group_id, *pair.keys(), 3, motif_id + 1, "m", "i", motif_seqs[0][1:-1], sep="\t", file=file_out)
                        
                elif (len(motif_seqs[0]) > 2 and len(motif_seqs[1]) > 2) and (motif_seqs[0][0] != motif_seqs[1][0] or motif_seqs[0][-1] != motif_seqs[1][-1]): 
                    a, _ = substitutions2(*pair.values())
                    subst = ",".join((map(lambda x: f"{x[0]}|{x[1]}", filter(lambda y: y[1] > 0, zip(a.index, a.values)))))
                    print(oma_group_id, *pair.keys(), 4, motif_id + 1, "m", "s", subst, sep="\t", file=file_out)
                    if "PP" in motif_seqs_no_gaps[0]:
                        print(oma_group_id, *pair.keys(), 4, motif_id + 1, "m", "i", motif_seqs[1][1:-1], sep="\t", file=file_out)
                    elif "PP" in motif_seqs_no_gaps[1]:
                        print(oma_group_id, *pair.keys(), 4, motif_id + 1, "m", "i", motif_seqs[0][1:-1], sep="\t", file=file_out)
        for pair in aligned_no_motifs:
            _, a = substitutions2(*pair.values())
            subst = ",".join((map(lambda x: f"{x[0]}|{x[1]}", filter(lambda y: y[1] > 0, zip(a.index, a.values)))))
            print(oma_group_id, *pair.keys(), 1, 0, "s", "s", subst, sep="\t", file=file_out)
        
        
def worker(worker_id, base_dir, file_set, out_dir):
    with open(os.path.join(out_dir, f"output_{worker_id}.tsv"), "a") as file_out:
        main(base_dir, file_set, worker_id, file_out)
    
def set_up(args):
    if args.output not in os.listdir():
        os.mkdir(args.output)
        all_files = os.listdir(args.oma_groups_dir)
        batch_size = math.ceil(len(all_files) / args.threads)
        file_sets = [all_files[batch_size*i:batch_size*(i+1)] for i in range(args.threads)]
        return file_sets
    else:
        shutil.rmtree(args.output)
        return set_up(args)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("oma_groups_dir", help="Path to directory containing fasta files which represent OMA groups")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=1, help="# of threads to use")
    args = parser.parse_args()
    
    file_sets = set_up(args)
    processes = [Process(target=worker, args=(worker_id, args.oma_groups_dir, file_sets[worker_id], args.output)) for worker_id in range(args.threads)]
    for pid, _ in enumerate(processes):
        processes[pid].start()
    for pid, _ in enumerate(processes):
        processes[pid].join()

