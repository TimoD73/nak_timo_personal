#!/usr/bin/env python

"""
Perform multiple sequence alignments for the produced amplicons of the ingroup as well as the outgroup if applicable.
"""

__author__ = "Timo Dijkstra BSc."
__version__ = "1.0"
__title__ = "MSA of ingroup and outgroup amplicons"
__author_email__ = "T.Dijkstra@naktuinbouw.nl"
__date__ = "03-12-2021"

import ast
import argparse
from typing import List, DefaultDict, Dict, Tuple, Any
import subprocess
from pathlib import Path
from collections import defaultdict
from colorama import Fore

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from utils.utils import is_valid_path, is_valid_directory

FASTA_FILE = "Amplicon_sequences"
SEQUENCE_EXTENSION = ".fasta"
DELAY_IDENTITY = 0


def parse_args() -> argparse:
    """ Parse command line arguments """
    parser = argparse.ArgumentParser(description="Perform multiple sequence alignments for the produced amplicons of "
                                                 "the ingroup as well as the outgroup if applicable.")
    parser.add_argument("-i", "--input", required=True, type=lambda x: is_valid_path(parser, x),
                        help="Input .tsv file of the amplicons produced by Combine_pipeline_results.py")
    parser.add_argument("-id", required=True, type=lambda x: is_valid_path(parser, x),
                        help="Path to the .tsv file containing the genome identifiers.")
    parser.add_argument("-cw", "--clustalw", required=True, type=lambda x: is_valid_path(parser, x),
                        help="Path to the ClustalW executable: clustalw2")
    parser.add_argument("-o", "--output", required=False, default="./",
                        type=lambda x: is_valid_directory(parser, x), help="Output directory.")
    args = parser.parse_args()

    return args


def read_target_amplicons(target_column: pd.Series, id_df: pd.DataFrame) -> \
        Tuple[DefaultDict[int, List[SeqRecord]], DefaultDict[str, List[Any]]]:
    """
    The function will convert the input column with target amplicons together with the IDs to a list of seqrecords.
    We will assume that due to the nature of the pipeline, every assay has a valid target amplicon for every genome.
    :param target_column: The column of the dataframe that contains all the target amplicons.
    :param id_df: Dataframe containing the genome ids.
    :return: Default dictionary with as key the assay and as value a list of target amplicon seqrecords.
    id_dictionary to match the IDs with the file names.
    """
    dict_amplicon_seq_records = defaultdict(list)
    id_dict = defaultdict(list)
    for assay_index, assay in enumerate(target_column):
        assay = assay.split("\n")
        for genome_index, amplicon in enumerate(assay):
            if amplicon != "NA":
                genome = id_df["Ingroup genomes raw"][genome_index]
                id = ["Target", assay_index + 1, genome_index + 1]
                dict_amplicon_seq_records[assay_index].append(SeqRecord(
                    Seq(amplicon),
                    id=".".join([str(string_element) for string_element in id]),
                    description=genome
                ))
                id_dict["File name"].append(genome)
                id_dict["Group"].append("Ingroup")
                id_dict["Target"].append(True)
                id_dict["Assay index"].append(assay_index + 1)
                id_dict["Genome index"].append(genome_index + 1)
                id_dict["Amplicon index"].append(1)
    return dict_amplicon_seq_records, id_dict


def read_off_target_amplicons(amplicon_list_as_string: pd.Series, id_df: pd.DataFrame, group: str) -> \
        Tuple[DefaultDict[int, List[SeqRecord]], DefaultDict[str, List[Any]]]:
    """
    The function will convert the input column with off-target amplicons together with the IDs to a list of seqrecords.
    The list is stored as string and therefore it needs to be converted to a list with the literal_eval function.
    Amplicons that are not available are skipped.
    :param amplicon_list_as_string: The column of the dataframe that contains the off-target amplicons.
    :param id_df: Dataframe containing the genome ids.
    :param group: Either <Ingroup> or <Outgroup>.
    :return: Default dictionary with as key the assay and as value a list of off-target amplicon seqrecords.
    """
    dict_amplicon_seq_records = defaultdict(list)
    id_dict = defaultdict(list)
    for assay_index, assay in enumerate(amplicon_list_as_string):
        assay = np.asarray(ast.literal_eval(assay), dtype=object)
        for genome_index, genome_amplicons in enumerate(assay):
            genome = id_df[group + " genomes raw"][genome_index]
            for amplicon_index, amplicon in enumerate(genome_amplicons):
                if amplicon != "NA":
                    id = ["Off-target", group, assay_index + 1, genome_index + 1, amplicon_index + 1]
                    dict_amplicon_seq_records[assay_index].append(SeqRecord(
                        Seq(amplicon),
                        id=".".join([str(string_element) for string_element in id]),
                        description=genome
                    ))
                    id_dict["File name"].append(genome)
                    id_dict["Group"].append(group)
                    id_dict["Target"].append(False)
                    id_dict["Assay index"].append(assay_index + 1)
                    id_dict["Genome index"].append(genome_index + 1)
                    id_dict["Amplicon index"].append(amplicon_index + 1)
    return dict_amplicon_seq_records, id_dict


def get_amplicon_sequences(mfe_tsv_file: Path, id_tsv_file: Path) -> \
        Tuple[Dict[str, Dict[int, List[SeqRecord]]], Dict[str, Dict[int, List[SeqRecord]]], Dict[str, List[Any]]]:
    """
    Generate a list off seqrecords for all amplicons predicted by MFEprimer.
    :param mfe_tsv_file: Path to the .tsv file of the amplicons produced by Combine_pipeline_results.py
    :param id_tsv_file: Path to the .tsv file containing the genome identifiers.
    :return: Nested dictionary of seqrecords for all amplicons predicted by MFEprimer for every assay both for regular
    amplicons and 50bp extended amplicons on both sides. In addition, return a dictionary to match the fasta headers
    to genomes and other information.
    """
    amplicons_regular = {}
    amplicons_extended = {}

    amplicon_df = pd.read_csv(mfe_tsv_file, sep="\t")
    id_df = pd.read_csv(id_tsv_file, sep="\t")

    targets = amplicon_df["amplicon_target_ingroup"]
    off_targets_ingroup = amplicon_df["qpcr_detectable_off_target_amplicon_sequences_ingroup"]
    amplicons_outgroup = amplicon_df["qpcr_detectable_amplicon_sequences_outgroup"]

    amplicons_regular["targets"], id_dict_target = read_target_amplicons(targets, id_df)
    amplicons_regular["off_targets_ingroup"], id_dict_off_target_ingroup = read_off_target_amplicons(off_targets_ingroup, id_df, "Ingroup")
    amplicons_regular["amplicons_outgroup"], id_dict_off_target_outgroup = read_off_target_amplicons(amplicons_outgroup, id_df, "Outgroup")

    targets_extended = amplicon_df["amplicon_target_ingroup_extended"]
    off_targets_ingroup_extended = amplicon_df["qpcr_detectable_off_target_amplicon_sequences_ingroup_extended"]
    amplicons_outgroup_extended = amplicon_df["qpcr_detectable_amplicon_sequences_outgroup_extended"]

    amplicons_extended["targets_extended"], id_dict_target = read_target_amplicons(targets_extended, id_df)
    amplicons_extended["off_targets_ingroup_extended"], id_dict_off_target_ingroup = read_off_target_amplicons(off_targets_ingroup_extended, id_df, "Ingroup")
    amplicons_extended["amplicons_outgroup_extended"], id_dict_off_target_outgroup = read_off_target_amplicons(amplicons_outgroup_extended, id_df, "Outgroup")

    # Transpose dictionary
    amplicons_dict_regular = pd.DataFrame(amplicons_regular).transpose().to_dict()
    amplicons_dict_extended = pd.DataFrame(amplicons_extended).transpose().to_dict()

    # Merge defaultdicts with ID information
    id_dict = {key: value + id_dict_off_target_ingroup[key] + id_dict_off_target_outgroup[key] for
               key, value in id_dict_target.items()}

    return amplicons_dict_regular, amplicons_dict_extended, id_dict


def generate_fasta_files(amplicons_dict: Dict[str, Dict[int, List[SeqRecord]]], output_dir: Path):
    """
    Write the seqrecords to a separate output fasta file for every assay.
    :param amplicons_dict: Nested dictionary of seqrecords for all amplicons predicted by MFEprimer for every assay.
    :param output_dir: Path to the output directory.
    """
    for assay in list(amplicons_dict.keys()):
        amplicons = []
        file_name = FASTA_FILE + "_assay" + str(assay + 1)
        output_path = output_dir / file_name
        output_path = output_path.with_suffix(SEQUENCE_EXTENSION)
        for amplicon_type in list(amplicons_dict[assay].keys()):
            if type(amplicons_dict[assay][amplicon_type]) == list:
                amplicons += amplicons_dict[assay][amplicon_type]
        with open(output_path, "w") as output_handle:
            SeqIO.write(amplicons, output_handle, "fasta")


def run_clustalw(clustal_path: str, in_file: Path):
    """
    Command to run ClustalW for all input fasta files.
    :param clustal_path: Path to the ClustalW executable: clustalw2
    :param in_file: Path to the fasta file of N-terminus effector sequences for every Signal family.
    """
    amplicon_fasta_files = sorted(list(in_file.glob("*" + SEQUENCE_EXTENSION)))
    print(amplicon_fasta_files)
    for assay in amplicon_fasta_files:
        clustalw_cline = ClustalwCommandline(clustal_path,
                                             infile=assay,
                                             type="DNA",
                                             maxdiv=DELAY_IDENTITY,
                                             OUTORDER="INPUT")
        subprocess.run(F"{str(clustalw_cline)}", shell=True)
        print(F"{Fore.BLUE}DONE: multiple sequence alignment for {assay.stem}")


def main():
    args = parse_args()

    output = Path(args.output)
    output_regular_amplicons = output / "regular"
    output_extended_amplicons = output / "extended"
    Path.mkdir(output_regular_amplicons, exist_ok=True)
    Path.mkdir(output_extended_amplicons, exist_ok=True)

    amplicons_dict_regular, amplicons_dict_extended, id_dict = get_amplicon_sequences(Path(args.input), args.id)

    generate_fasta_files(amplicons_dict_regular, output_regular_amplicons)
    generate_fasta_files(amplicons_dict_extended, output_extended_amplicons)

    run_clustalw(args.clustalw, output_regular_amplicons)
    run_clustalw(args.clustalw, output_extended_amplicons)

    id_df = pd.DataFrame.from_dict(id_dict)
    # noinspection PyTypeChecker
    id_df.to_csv(output / "MSA.tsv", sep="\t")

    print(F"{Fore.GREEN}FINISHED")


if __name__ == '__main__':
    main()

