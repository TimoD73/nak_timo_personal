#!/usr/bin/env python

'''
Generate primers for effector gene sequences with Primer3
'''

__author__ = "Timo Dijkstra BSc."
__version__ = "1.2"
__title__ = "Generate primers"
__author_email__ = "T.Dijkstra@naktuinbouw.nl"
__date__ = "22-09-2021"

import argparse
import re
import primer3
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from itertools import repeat
from typing import TextIO, List, Dict, Any, Optional
from colorama import Fore

PRIMER_PRODUCT_OPTIMUM_MARGIN = 20
NUMBER_OF_PRIMERS = 5


def parse_args() -> argparse:
    """ Parse command line arguments """
    parser = argparse.ArgumentParser(description='Runs Primer3 on specified nucleotide sequences to generate primer'
                                                 'pairs. Writes output in tabular .tsv format. The fasta header'
                                                 'will be used as identifier.')
    parser.add_argument('-i', '--input', type=argparse.FileType("r"), required=True,
                        help='Sequences for which primers will be generated')
    parser.add_argument('-p', '--parameters', type=str, required=True,
                        help='Config file in .tsv format that stores the parameters for primer design')
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output file (.tsv format)')
    parser.add_argument('-cl', '--cut_off_left', type=int, required=False, default=None,
                        help='Slice the input sequences starting from base number.')
    parser.add_argument('-cr', '--cut_off_right', type=int, required=False, default=None,
                        help='Slice the input sequences until base number.')

    args = parser.parse_args()
    return args


def read_effectors(input_database: TextIO) -> List[SeqIO.SeqRecord]:
    """
    Process the fasta file using the parse function of SeqIO.
    :param input_database: A Fasta file downloaded from T3Enc containing protein sequences.
    The database can be found using the following link: http://www.szu-bioinf.org/T3Enc/index.html.
    :return: list_seq_objects: A list with fasta records.
    """
    list_seq_objects = []
    for record in SeqIO.parse(input_database, "fasta"):
        list_seq_objects.append(record)

    return list_seq_objects


def read_config_file(parameter_tsv_file: str) -> Dict[str, int]:
    """
    Read the input config .tsv file and save the Primer3 parameters in a dictionary.
    :param parameter_tsv_file: .tsv file containing the parameters.
    :return: A dictionary with as key the parameter name as described by Primer3 and as value the parameter.
    """
    parameters_table = pd.read_csv(parameter_tsv_file, sep="\t", names=["Name_parameter", "Parameter"])
    parameters = dict(zip(parameters_table["Name_parameter"].tolist(), parameters_table["Parameter"].tolist()))

    return parameters


def slice_input_sequence(sequence: str, left_boundary: Optional[int] = None, right_boundary: Optional[int] = None) \
        -> str:
    """
    A function that gives the user the option to search for primers in only a region of the sequence.
    The function extracts the specified region.
    :param sequence: The target sequence.
    :param left_boundary: Starting base number.
    :param right_boundary: Ending base number.
    :return: The (extracted region of the) sequence.
    """
    if left_boundary is None:
        return sequence[0:right_boundary]
    elif right_boundary is None:
        return sequence[left_boundary - 1:len(sequence)]
    elif left_boundary is not None and right_boundary is not None:
        return sequence[left_boundary - 1:right_boundary]
    else:
        return sequence


def generate_primers(sequence_id: str, sequence: str, primer_design_parameters: Dict[str, int]) -> Dict[str, Any]:
    """
    Compute primers for the given sequences.
    :param sequence_id: String that identifies the target sequence.
    :param sequence: The target sequence.
    :param primer_design_parameters: Dictionary that contains the Primer3 parameters.
    :return: Dictionary containing the primer design results for five possible primers.
    """
    primer_dict = primer3.designPrimers(
        {
            'SEQUENCE_ID': sequence_id,
            'SEQUENCE_TEMPLATE': sequence
        },
        {
            'PRIMER_OPT_SIZE': primer_design_parameters['PRIMER_OPT_SIZE'],
            'PRIMER_PICK_INTERNAL_OLIGO': primer_design_parameters['PRIMER_PICK_INTERNAL_OLIGO'],
            'PRIMER_INTERNAL_MAX_SELF_END': primer_design_parameters['PRIMER_INTERNAL_MAX_SELF_END'],
            'PRIMER_MIN_SIZE': primer_design_parameters['PRIMER_MIN_SIZE'],
            'PRIMER_MAX_SIZE': primer_design_parameters['PRIMER_MAX_SIZE'],
            'PRIMER_MIN_TM': primer_design_parameters['PRIMER_MIN_TM'],
            'PRIMER_MAX_TM': primer_design_parameters['PRIMER_MAX_TM'],
            'PRIMER_OPT_TM': primer_design_parameters['PRIMER_OPT_TM'],
            'PRIMER_MIN_GC': primer_design_parameters['PRIMER_MIN_GC'],
            'PRIMER_MAX_GC': primer_design_parameters['PRIMER_MAX_GC'],
            'PRIMER_INTERNAL_MIN_SIZE': primer_design_parameters['PRIMER_INTERNAL_MIN_SIZE'],
            'PRIMER_INTERNAL_MAX_SIZE': primer_design_parameters['PRIMER_INTERNAL_MAX_SIZE'],
            'PRIMER_INTERNAL_MIN_TM': primer_design_parameters['PRIMER_INTERNAL_MIN_TM'],
            'PRIMER_INTERNAL_OPT_TM': primer_design_parameters['PRIMER_INTERNAL_OPT_TM'],
            'PRIMER_INTERNAL_MAX_TM': primer_design_parameters['PRIMER_INTERNAL_MAX_TM'],
            'PRIMER_INTERNAL_MIN_GC': primer_design_parameters['PRIMER_INTERNAL_MIN_GC'],
            'PRIMER_INTERNAL_MAX_GC': primer_design_parameters['PRIMER_INTERNAL_MAX_GC'],
            'PRIMER_MAX_POLY_X': primer_design_parameters['PRIMER_MAX_POLY_X'],
            'PRIMER_INTERNAL_MAX_POLY_X': primer_design_parameters['PRIMER_INTERNAL_MAX_POLY_X'],
            'PRIMER_SALT_MONOVALENT': primer_design_parameters['PRIMER_SALT_MONOVALENT'],
            'PRIMER_DNA_CONC': primer_design_parameters['PRIMER_DNA_CONC'],
            'PRIMER_MAX_NS_ACCEPTED': primer_design_parameters['PRIMER_MAX_NS_ACCEPTED'],
            'PRIMER_MAX_SELF_ANY': primer_design_parameters['PRIMER_MAX_SELF_ANY'],
            'PRIMER_MAX_SELF_END': primer_design_parameters['PRIMER_MAX_SELF_END'],
            'PRIMER_PAIR_MAX_COMPL_ANY': primer_design_parameters['PRIMER_PAIR_MAX_COMPL_ANY'],
            'PRIMER_PAIR_MAX_COMPL_END': primer_design_parameters['PRIMER_PAIR_MAX_COMPL_END'],
            'PRIMER_PRODUCT_SIZE_RANGE': [[primer_design_parameters['PRIMER_PRODUCT_SIZE_MINIMUM'],
                                           primer_design_parameters['PRIMER_PRODUCT_SIZE_MINIMUM'] +
                                           PRIMER_PRODUCT_OPTIMUM_MARGIN],
                                          [primer_design_parameters['PRIMER_PRODUCT_SIZE_MINIMUM'],
                                           primer_design_parameters['PRIMER_PRODUCT_SIZE_MAXIMUM']]]
        }
    )
    print(F"{Fore.BLUE}Primer design completed for {sequence_id}.")
    return primer_dict


def generate_nested_dictionary(primers: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Make a nested dictionary for all the five primers. First, define the keys for the nested dictionary.
    To do this, replace the _0, _1, _2, _3, or _4 with a general whitespace, resulting in one underscore.
    Subsequently, fill in the nested dictionary.
    :param primers: Dictionary containing the primer design results for five possible primers.
    :return: Nested dictionary storing the primer design results.
    Access the filled in nested dictionary using: primer_dictionary[<Primers> + <number>][<Primer3 output information>]
    """
    primer_dictionary = {}

    replacements = {f"_{i}": "" for i in list(range(0, NUMBER_OF_PRIMERS))}
    replacements = dict((re.escape(k), v) for k, v in replacements.items())
    pattern = re.compile("|".join(replacements.keys()))

    for i in range(1, NUMBER_OF_PRIMERS + 1):
        primer_dictionary["Primers " + str(i)] = {}
        for k, v in primers.items():
            k = pattern.sub(lambda m: replacements[re.escape(m.group(0))], k)
            primer_dictionary["Primers " + str(i)][k] = v

    return primer_dictionary


def generate_dataframe_primer_set(primer_dictionary: Dict[str, Dict[str, Any]], sequence_header: str) -> pd.DataFrame:
    """
    Compute a pandas dataframe as output for the primers.
    Add an extra column with the protein accession to the dataframe.
    :param primer_dictionary: Nested dictionary storing the primer design results.
    :param sequence_header: String that identifies the target sequence.
    :return: A dataframe in which each row is a primer set and each column a feature of the primer set.
    """
    sequence_header_column = []
    sequence_header_column.extend(repeat(sequence_header, NUMBER_OF_PRIMERS))

    primer_df = pd.DataFrame(primer_dictionary).transpose()
    primer_df.insert(0, "Fasta header", sequence_header_column)

    return primer_df


def output_concatenated_dataframe(list_dataframes: List[pd.DataFrame], output: str):
    """
    Concatenate a list of dataframes for each target sequence to one dataframe.
    :param list_dataframes: List of dataframes for each target sequence.
    :param output: String that defines that output path.
    :return: A concatenated dataframe. Each row is an instance of a sequence - primer set combination.
    Each column describes a feature of the primer set.
    """
    concat_dataframe = pd.concat(list_dataframes).drop_duplicates().dropna()
    # noinspection PyTypeChecker
    concat_dataframe.to_csv(output, sep="\t")


def main():
    args = parse_args()
    list_primer_dataframes = []

    list_seq_objects = read_effectors(args.input)
    parameters = read_config_file(args.parameters)

    for i in range(len(list_seq_objects)):
        sequence = slice_input_sequence(str(list_seq_objects[i].seq), args.cut_off_left, args.cut_off_right)
        primers = generate_primers(list_seq_objects[i].description, sequence, parameters)
        primer_dictionary = generate_nested_dictionary(primers)
        list_primer_dataframes.append(generate_dataframe_primer_set(primer_dictionary, list_seq_objects[i].description))
    output_concatenated_dataframe(list_primer_dataframes, args.output)

    print(F"{Fore.GREEN}FINISHED")


if __name__ == '__main__':
    main()
