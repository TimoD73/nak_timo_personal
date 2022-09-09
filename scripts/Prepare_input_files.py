#!/usr/bin/env python

'''
Make a ready to run input primer file for MFEprimer.
'''

__author__ = "Timo Dijkstra BSc."
__version__ = "1.0"
__title__ = "MFEprimer preparation"
__author_email__ = "T.Dijkstra@naktuinbouw.nl"
__date__ = "25-10-2021"


import argparse
import pandas as pd
from typing import TextIO
from utils.utils import is_valid_directory


def parse_args() -> argparse:
    """ Parse command line arguments """
    parser = argparse.ArgumentParser(description='Compute a ready to run .fa file for MFEprimer.')
    parser.add_argument('-i', '--input', type=argparse.FileType("r"), required=True,
                        help='Table  generated by Primer_Design_script.py.'
                             'with all the primers needed for quality control.')
    parser.add_argument('-n', '--row', type=int, required=True,
                        help='Row number to extract.')
    parser.add_argument('-o', '--output', type=lambda x: is_valid_directory(parser, x), required=False, default="./",
                        help='Output directory')

    args = parser.parse_args()
    return args


def process_input_df(table_path: TextIO, row_number: int) -> pd.DataFrame:
    """
    Take the output dataframe from the Primer_Design_script.py script and
    extract the columns that identify the primer set and the corresponding sequences.
    :param table_path: A .tsv file computed by Primer_Design_script.py.
    :param row_number: The row number that determines which primerset to extract.
    :return: A dataframe with 1 row and 5 columns specific for one primer set.
    """
    primer_table = pd.read_csv(table_path, sep="\t")
    # Column 0 and 1 identify the primer. Column 14, 15, 16 contain the fwd, rev, and probe sequences respectively.
    primer_table_subset = primer_table.iloc[:, [0, 1, 14, 15, 16]].loc[row_number]

    return primer_table_subset


def check_succesful_primer_design(primer_df: pd.DataFrame):
    """
    If primer design was not successful, terminate the script.
    :param primer_df: A dataframe with 1 row and 5 columns specific for one primer set.
    """
    if primer_df.isnull().values.any():
        quit()


def generate_fasta_primer(primer_df: pd.DataFrame, output_path: str):
    """
    Format a fasta file for the left primer, the right primer and the internal oligo.
    :param primer_df: A dataframe with 1 row and 5 columns specific for one primer set.
    :param output_path: String that defines the output path.
    """
    with open(output_path + "Primer_sequences.fa", "w") as fa:
        fa.write(">" + primer_df.iloc[0].replace("s ", "_") + list(primer_df.keys())[2].replace("PRIMER", "") +
                 "|" + str(primer_df.loc["Fasta header"]) + "\n")
        fa.write(primer_df.loc["PRIMER_LEFT_SEQUENCE"] + "\n")
        fa.write(">" + primer_df.iloc[0].replace("s ", "_") + list(primer_df.keys())[3].replace("PRIMER", "") +
                 "|" + str(primer_df.loc["Fasta header"]) + "\n")
        fa.write(primer_df.loc["PRIMER_RIGHT_SEQUENCE"] + "\n")
        fa.write(">" + primer_df.iloc[0].replace("s ", "_") + list(primer_df.keys())[4].replace("PRIMER", "") +
                 "|" + str(primer_df.loc["Fasta header"]) + "\n")
        fa.write(primer_df.loc["PRIMER_INTERNAL_SEQUENCE"])


def main():
    args = parse_args()
    primer_df = process_input_df(args.input, args.row - 1)
    check_succesful_primer_design(primer_df)
    generate_fasta_primer(primer_df, args.output)

    print(primer_df["Fasta header"])


if __name__ == '__main__':
    main()