import argparse
from multiprocessing import Process
import math
import os
import re
import shutil
from itertools import combinations
from typing import List, Union, Tuple, Match, TextIO, Generator

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
from pandas import Series


class ProlineSubstitutionRecord:
    """Instances of this class represent groups of proline substitution events and store
    following information:
        1. Number of proline substitutions for each amino-acid + gap (subst_counts)
        2. Number of motif in the sequence (motif_num)
        3. Sequence info that shows location of reference prolines (seq_num) 0 for 1st seq, 1 for 2nd seq
        4. Type of event. One of the following: single proline substitution, substitution of proline in
        motif without insertions and substitution of proline in motif with insertion (type)"""
    def __init__(self, substitution_counts: Series, motif_num: int, seq_num: int, type_: str):
        self.motif_num = motif_num
        self.subst_counts = substitution_counts
        self.seq_num = seq_num  # 0 - seqA, 1 - seqB
        if type_ == "single":
            self.type_id = 1
        elif type_ == "motif_subst_only":
            self.type_id = 2
        elif type_ == "motif_subst_mixed":
            self.type_id = 4
        else:
            raise RuntimeError("{type_} is invalid type for {type(self)} instance")
        self.type = type_
        self.type_short = type_[0]

    def create_substitution_encoding_string(self) -> str:
        """Format substitution data for convenient storage.
        Output example: A|1,P|4,G|2. This means proline was substituted by A and G 1 time
        and saved 4 times. Substitutions that took place 0 times are not saved."""
        real_substitutions = filter(lambda y: y[1] > 0, zip(self.subst_counts.index, self.subst_counts.values))
        formatted_substitutions = map(lambda x: f"{x[0]}|{x[1]}", real_substitutions)
        subst_string = ",".join(formatted_substitutions)
        return subst_string

    def print_formatted(self, oma_group_id: str, seqA_id: str = "seqA", seqB_id: str = "seqB", **kwargs):
        """Print all information about mutation event in the following format:
        GROUP_ID SEQ_A_ID SEQ_B_ID TYPE_ID MOTIF_NUMBER LOCATION_TYPE EVENT_TYPE ENCODED_SUBSTITUTION_STRING"""
        if sum(self.subst_counts):
            subst_string = self.create_substitution_encoding_string()
            seq_id_pair = (seqA_id, seqB_id) if not self.seq_num else (seqB_id, seqA_id)
            print(oma_group_id, *seq_id_pair, self.type_id, self.motif_num,
                  self.type_short, "s", subst_string, **kwargs)


class ProlineMotifInsertionRecord:
    """Instances of this class represent  insertion events and store
        following information:
            1. Insertion sequence (insertion_sequence)
            2. Number of motif in the sequence (motif_num)
            3. Sequence info that shows location of reference prolines (seq_num) 0 for 1st seq, 1 for 2nd seq
            4. Type of event. One of the following: insertion in motif without substitutions and
            insertion in motif with substitutions (type)"""
    def __init__(self, insertion_sequence: Union[str, None], motif_num: int, seq_num: int, type_: str):
        self.insertion_sequence = insertion_sequence
        self.motif_num = motif_num
        self.seq_num = seq_num  # 0 - seqA, 1 - seqB
        if type_ == "motif_gaps_only":
            self.type_id = 3
        elif type_ == "motif_gaps_mixed":
            self.type_id = 4
        else:
            raise RuntimeError("{type_} is invalid type for {type(self)} instance")
        self.type = type_
        self.type_short = type_[0]

    def print_formatted(self, oma_group_id: str, seqA_id: str = "seqA", seqB_id: str = "seqB", **kwargs):
        """Print all information about mutation event in the following format:
        GROUP_ID SEQ_A_ID SEQ_B_ID TYPE_ID MOTIF_NUMBER LOCATION_TYPE EVENT_TYPE ENCODED_SUBSTITUTION_STRING"""
        if self.insertion_sequence:
            seq_id_pair = (seqA_id, seqB_id) if not self.seq_num else (seqB_id, seqA_id)
            print(oma_group_id, *seq_id_pair, self.type_id, self.motif_num,
                  self.type_short, "i", self.insertion_sequence, **kwargs)


def parse_description_and_motifs(description: str) -> Tuple[str, Union[List[Tuple[str, int, int, str]], None]]:
    """Extracts PROTEIN_ID and MOTIFS from a fasta description having following format:
    >PROTEIN_ID | TAXONOMIC_ID | ORTHOLOGY_GROUP | MOTIF1_SEQ:MOTIF1_START:MOTIF1_END:MOTIF1_STRENGTH MOTIF2..."""
    desc_splitted = description.split(" | ")
    if len(desc_splitted) > 3:    # Polyproline motifs are present
        all_motifs_joined = desc_splitted[3]
        all_motifs_splitted = all_motifs_joined.split(" ")   # Motifs are space-separated
        individual_motifs_joined = map(lambda x: x.split(":"), all_motifs_splitted)
        # List of tuples each representing single motif (sequence, start, end, strength)
        individual_motifs_splitted = list(map(lambda k: (k[0], int(k[1]), int(k[2]), k[3]), individual_motifs_joined))
        return desc_splitted[0], individual_motifs_splitted
    else:
        return desc_splitted[0], None


def generate_all_best_pairwise_alignments(sequences: List[Seq]) -> Generator[Tuple[Tuple[str, str, float, int, int],
                                                                                   Tuple[int, int]], None, None]:
    """Returns best pairwise alignments for given set of sequences and
    pairs of indices of sequences in these alignments. For each pair of sequence
    an alignment with the highest score are returned"""
    for pair, idx_pair in zip(combinations(sequences, 2), combinations(range(len(sequences)), 2)):
        best_alignment = max(pairwise2.align.globalds(*pair, matlist.blosum62, -10, -0.5), key=lambda x: x[2])
        yield best_alignment, idx_pair


def find_aligned_polyproline_motifs(aligned_seq: str) -> List[Tuple[int, Match]]:
    """Find polyproline motifs in aligned sequence (sequence with gaps taken from alignment).
    It is assumed that polyproline motif is two adjacent prolines followed by non-proline.
    The function returns list of tuples containing position of polyproline motif and its Match object.
    Motif sequence is available in group 1, gap sequence (only '-') is available in group 2"""
    motif_candidates = []
    for i, _ in enumerate(aligned_seq):
        match = re.match(r"(P(-*)P)[^P-]?", aligned_seq[i:])
        if match:
            motif_candidates.append((i, match))
    return motif_candidates


def find_single_prolines(aligned_seq: str) -> List[int]:
    """Find single prolines in aligned sequence (sequence with gaps taken from alignment).
        It is assumed that single proline has no adjacent prolines even over gaps.
        E.g. AYP-----PKYT. Proline at index 8 is not single because it is adjacent to proline at index 2.
        The function returns list of integers which are indices of single prolines in sequence"""
    motif_candidates = []
    for idx, aa in enumerate(aligned_seq):
        if aa == 'P':
            if is_single(aligned_seq, idx):
                motif_candidates.append(idx)
    return motif_candidates


def is_single(seq: str, position: int) -> bool:
    """Check if proline at given position is single by seeking adjacent prolines on the left and on the right"""
    if position == 0:
        directions = (1,)
    elif position == len(seq) - 1:
        directions = (-1,)
    else:
        directions = (-1, 1)
    for direction in directions:
        pos = position + direction
        while pos != -1 and pos != len(seq):
            if seq[pos] == "-":
                pos += direction
            elif seq[pos] == "P":
                return False
            else:
                break
    return True


def describe_prolines_in_alignment(seq1: str, seq2: str) -> Tuple[List[ProlineSubstitutionRecord],
                                                                  List[ProlineSubstitutionRecord],
                                                                  List[ProlineMotifInsertionRecord],
                                                                  List[ProlineSubstitutionRecord],
                                                                  List[ProlineMotifInsertionRecord]]:
    """Describes proline-related mutations in a pairwise alignment. This includes single prolines substitutions,
    substitutions in polyproline motifs and insertions in polyproline motifs.
    Output consists of 5 lists. Each containing ProlineSubstitutionRecord instance or ProlineMotifInsertionRecord
    depending on event (substitution or insertion):
        1. Single proline substitutions
        2. Substitutions of prolines in polyproline motifs without insertions
        3. Insertions in polyproline motifs without substitutions
        4. Substitutions of prolines in polyproline motifs with insertions (rare case)
        5. Insertions in polyproline motifs with substitutions (rare case)"""
    # Final lists (plural)
    single_substs = []
    motif_substs_only = []
    motif_substs_mixed = []
    motif_gaps_only = []
    motif_gaps_mixed = []

    # Describe one seq than another (swap them)
    for seq_num, (seq, inv_seq) in enumerate(((seq1, seq2), (seq2, seq1))):
        for motif_num, (pos, match) in enumerate(find_aligned_polyproline_motifs(seq)):
            # Temporary containers for each motif (singular)
            motif_subst = pd.Series(index=list("ACDEFGHIKLMNPQRSTVWYX-"), dtype=int)
            motif_subst_mixed = pd.Series(index=list("ACDEFGHIKLMNPQRSTVWYX-"), dtype=int)
            gaps_only = [None]
            gaps_mixed = [None]

            start, end = pos, pos + match.end(1) - 1
            gap_seq = match.group(2)
            if gap_seq:
                if inv_seq[start] != seq[start] or inv_seq[end] != seq[end]:
                    gaps_mixed.append(inv_seq[start + 1:end])
                    motif_subst_mixed[inv_seq[start]] += 1
                    motif_subst_mixed[inv_seq[end]] += 1
                else:
                    gaps_only.append(inv_seq[start + 1:end])
                    motif_subst[inv_seq[start]] += 1
                    motif_subst[inv_seq[end]] += 1
            else:
                motif_subst[inv_seq[start]] += 1
                motif_subst[inv_seq[end]] += 1
            motif_substs_only.append(ProlineSubstitutionRecord(motif_subst, motif_num + 1,
                                                               seq_num, "motif_subst_only"))
            motif_substs_mixed.append(ProlineSubstitutionRecord(motif_subst_mixed, motif_num + 1,
                                                                seq_num, "motif_subst_mixed"))
            motif_gaps_only.append(ProlineMotifInsertionRecord(gaps_only[-1], motif_num + 1,
                                                               seq_num, "motif_gaps_only"))
            motif_gaps_mixed.append(ProlineMotifInsertionRecord(gaps_mixed[-1], motif_num + 1,
                                                                seq_num, "motif_gaps_mixed"))
    for seq_num, (seq, inv_seq) in enumerate(((seq1, seq2), (seq2, seq1))):
        single_subst = pd.Series(index=list("ACDEFGHIKLMNPQRSTVWYX-"), dtype=int)
        for pos in find_single_prolines(seq):
            single_subst[inv_seq[pos]] += 1
        single_substs.append(ProlineSubstitutionRecord(single_subst, 0, seq_num, "single"))
    return single_substs, motif_substs_only, motif_gaps_only, motif_substs_mixed, motif_gaps_mixed


def get_substitution_data(base_dir: str, file_set: List[str], worker_id: int, file_out: TextIO) -> None:
    """Describes proline substitutions in a set of files representing orthology groups. Output has
    following format:
        GROUP_ID SEQ_A_ID SEQ_B_ID TYPE_ID MOTIF_NUMBER LOCATION_TYPE EVENT_TYPE ENCODED_SUBSTITUTION_STRING
    Function is designed for working in multiprocessing system, thus file_set for each worker and worker_id are
    required"""
    for file_num, oma_group_file in enumerate(file_set):
        # Progress visualization
        print(f"Worker {worker_id}: {file_num + 1}/{len(file_set)}", oma_group_file)
        oma_group_id = oma_group_file.split(".")[0]
        oma_group_records = list(SeqIO.parse(os.path.join(base_dir, oma_group_file), "fasta"))
        descriptions, sequences = \
            list(map(lambda x: x.description, oma_group_records)), list(map(lambda x: x.seq, oma_group_records))

        for algn_idx, (alignment, seqs) in enumerate(generate_all_best_pairwise_alignments(sequences)):
            if seqs[0] == seqs[1]:
                continue
            seqA_id, _ = parse_description_and_motifs(descriptions[seqs[0]])
            seqB_id, _ = parse_description_and_motifs(descriptions[seqs[1]])
            single_substs, motif_substs_only, motif_insertions_only, motif_substs_mixed, motif_insertions_mixed = \
                describe_prolines_in_alignment(alignment[0], alignment[1])
            for single_subst in single_substs:
                single_subst.print_formatted(oma_group_id, seqA_id, seqB_id, sep="\t", file=file_out)
            for motif_subst_only in motif_substs_only:
                motif_subst_only.print_formatted(oma_group_id, seqA_id, seqB_id, sep="\t", file=file_out)
            for motif_subst_mixed in motif_substs_mixed:
                motif_subst_mixed.print_formatted(oma_group_id, seqA_id, seqB_id, sep="\t", file=file_out)
            for motif_insertion_only in motif_insertions_only:
                motif_insertion_only.print_formatted(oma_group_id, seqA_id, seqB_id, sep="\t", file=file_out)
            for motif_insertion_mixed in motif_insertions_mixed:
                motif_insertion_mixed.print_formatted(oma_group_id, seqA_id, seqB_id, sep="\t", file=file_out)


def worker(worker_id: int, base_dir: str, file_set: List[str], out_dir: str):
    """Wrapper function that sets file handler for get_substitution_data.
        base_dir is folder with orthology groups
        file_set is set of files in this directory
        out_dir directory to write output into
    Files are named 'output_X.tsv' where X is worker_id """
    with open(os.path.join(out_dir, f"output_{worker_id}.tsv"), "a") as file_out:
        get_substitution_data(base_dir, file_set, worker_id, file_out)


def set_up(arguments: argparse.Namespace):
    """Preparation for program run. Creates output directory if it does not exist otherwise remove its contents.
    Creates file sets of nearly equal size for each worker depending on processes number and returns it.
        'data_subset' specifies number of files to use (for testing)"""
    if arguments.output not in os.listdir():
        os.mkdir(arguments.output)
        all_files = os.listdir(arguments.oma_groups_dir)[slice(0, arguments.subset, 1)]
        batch_size = math.ceil(len(all_files) / arguments.threads)
        file_sets_for_workers = [all_files[batch_size * i:batch_size * (i + 1)] for i in range(arguments.threads)]
        return file_sets_for_workers
    else:
        shutil.rmtree(arguments.output)
        return set_up(arguments)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("oma_groups_dir", help="Path to directory containing fasta files which represent OMA groups")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=1, help="# of threads to use")
    parser.add_argument("-s", "--subset", type=int, default=None, help="Number of files to process (for testing)")
    args = parser.parse_args()

    file_sets = set_up(args)
    processes = [Process(target=worker, args=(worker_id, args.oma_groups_dir, file_sets[worker_id], args.output)) for
                 worker_id in range(args.threads)]
    for pid, _ in enumerate(processes):
        processes[pid].start()
    for pid, _ in enumerate(processes):
        processes[pid].join()
