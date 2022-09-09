#!/usr/bin/env python

"""
Scan annotated proteins by Prokka for T3S effectors. In addition, extract the gene from the corresponding genome.
"""

__author__ = "Timo Dijkstra BSc."
__version__ = "1.0"
__title__ = "Scan for effectors"
__author_email__ = "T.Dijkstra@naktuinbouw.nl"
__date__ = "07-10-2021"

import argparse
import pandas as pd
from utils.utils import directory_exists, make_dataframe_effectors
from typing import TextIO, Tuple, List
from Bio import SeqIO, SeqRecord
from pathlib import Path
from colorama import Fore


def parse_args() -> argparse:
    """ Parse command line arguments """
    parser = argparse.ArgumentParser(description='Scan annotated proteins for T3S effectors and extract the genes.')
    parser.add_argument('-ffn', '--transcripts', required=True, type=argparse.FileType("r"),
                        help='Path to the .ffn output file from Prokka; the nucleotide FASTA file of all the prediction'
                             ' transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA).')
    parser.add_argument('-faa', '--proteins', help='Prokka .faa output file with protein sequences.', required=True,
                        type=argparse.FileType("r"))
    parser.add_argument('-t', '--t3sepp', help='The T3SEpp.out.txt file.', required=True,
                        type=argparse.FileType("r"))
    parser.add_argument('-og', '--output_genes', help='The output path of a fasta file containing effector genes.',
                        required=True, type=Path)
    parser.add_argument('-op', '--output_proteins', help='The output path of a fasta file containing effector proteins.',
                        required=True, type=Path)

    return parser.parse_args()


def parse_ffn(ffn: TextIO) -> List[SeqRecord.SeqRecord]:
    """
    Compute a list of seqrecords.
    :param ffn: Opened .ffn output file from Prokka.
    :retur: list of seqrecords.
    """
    transcript_records = []
    for record in SeqIO.parse(ffn, "fasta"):
        transcript_records.append(record)
    return transcript_records


def extract_effector_genes(df_effectors: pd.DataFrame, transcripts: List[SeqRecord.SeqRecord]) -> \
        Tuple[List[SeqRecord.SeqRecord], List[str]]:
    """
    Retrieve the effector genes using the IDs from the T3SEpp output.
    :param df_effectors: Dataframe containing only proteins predicted as effector.
    :param transcripts: A list of seqrecords of all genes predicted by Prokka.
    :return: List of effector genes.
    """
    list_effector_ids = df_effectors["prot"].tolist()
    effector_genes = [protein for protein in transcripts if protein.id in list_effector_ids]
    return effector_genes, list_effector_ids


def process_extract_effector_input(ffn: TextIO, t3sepp_out: TextIO) -> Tuple[List[SeqRecord.SeqRecord], List[str]]:
    """
    Master function that processes the two input files and returns a dataframe with solely effectors, the corresponding
    genome, and a list of all annotated proteins as class SeqFeature. The function also checks if T3SEpp produced
    valid output.
    :param ffn: Opened .ffn output file from Prokka.
    :param t3sepp_out: Direct output from T3SEpp (T3SEpp.out.txt file).
    :return effector_genes: List of effector genes.
    :return list_effector_ids: List with all the IDs corresponding to a putative effector.
    """
    try:
        df_effectors = make_dataframe_effectors(t3sepp_out)
        transcripts = parse_ffn(ffn)
        effector_genes, list_effector_ids = extract_effector_genes(df_effectors, transcripts)
        return effector_genes, list_effector_ids
    except pd.errors.ParserError:
        print(F"{Fore.LIGHTRED_EX}{'#'*59}\n{Fore.RED}T3SEpp did not produce valid output. Check the row lengths.\n"
              F"{Fore.LIGHTRED_EX}{'#'*59}")
        quit()
    except pd.errors.EmptyDataError:
        print(F"{Fore.LIGHTRED_EX}{'#'*59}\n{Fore.RED}T3SEpp did not produce output.\n"
              F"{Fore.LIGHTRED_EX}{'#'*59}")
        quit()


def generate_fasta_genes(effector_genes: List[SeqRecord.SeqRecord], output: Path):
    """
    Function that makes an output fasta file of the effector gene sequences.
    Function that makes an output fasta file of the effector gene sequences.
    :param effector_genes: List of effector genes.
    :param output: Output path of the fasta file.
    """
    with open(output, "w") as output_handle:
        fasta_out = SeqIO.FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_out.write_file(effector_genes)


def generate_fasta_proteins(list_effector_ids: List[str], faa: TextIO, output: Path):
    """
    Function that makes an output fasta file of the effector gene sequences.
    :param list_effector_ids: List with all the IDs corresponding to a putative effector.
    :param faa: Opened Prokka .faa output file with protein sequences.
    :param output: Output path of the fasta file.
    """
    protein_seq_records = [record for record in SeqIO.parse(faa, "fasta") if
                           record.id in list_effector_ids]
    with open(output, "w") as output_handle:
        fasta_out = SeqIO.FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_out.write_file(protein_seq_records)


def main():
    args = parse_args()
    directory_exists(args.output_genes.parent)
    directory_exists(args.output_proteins.parent)
    effector_genes, list_effector_ids = process_extract_effector_input(args.transcripts, args.t3sepp)
    generate_fasta_genes(effector_genes, args.output_genes)
    generate_fasta_proteins(list_effector_ids, args.proteins, args.output_proteins)
    print(F"{Fore.GREEN}FINISHED")


if __name__ == '__main__':
    main()
