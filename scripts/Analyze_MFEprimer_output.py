#!/usr/bin/env python

"""
Analyze MFEprimer-3.2.3 .json output for multiple files simultaneously.
"""

__author__ = "Timo Dijkstra BSc."
__version__ = "1.0"
__title__ = "Analyze MFEprimer-3.2.3 output"
__author_email__ = "T.Dijkstra@naktuinbouw.nl"
__date__ = "30-09-2021"

import argparse
import json

import pandas as pd
from typing import List, Dict, Optional, Tuple, Union
from pathlib import Path
from colorama import Fore
from utils.utils import is_valid_directory, is_valid_path
from Bio.Seq import Seq
from Bio import SeqIO

MFE_PRIMER_EXTENSION = ".json"
HAIRPIN = None
DIMER = None
QPCR_MIN_AMPLICON_LENGTH = 50
QPCR_MAX_AMPLICON_LENGTH = 300
FORWARD_PRIMER_ASSAY = "LEFT"
REVERSE_PRIMER_ASSAY = "RIGHT"
PROBE_ASSAY = "INTERNAL"
BP_AMPLICON_EXTENSION = 50


def parse_args() -> argparse:
    """ Parse command line arguments """
    parser = argparse.ArgumentParser(description='Analyze the results generated from MFEprimer-3.2.3.')
    parser.add_argument('-i', '--input_dir', type=lambda x: is_valid_directory(parser, x), metavar="DIR", required=True,
                        help='Path of the directory with the MFEprimer-3.2.3 .json files.')
    parser.add_argument('-g', '--genome', type=Path, metavar="FILE", required=True,
                        help="Path to the genome analyzed by MFEprimer corresponding to the input directory "
                             "with json files.")
    parser.add_argument("-c", "--control", required=False, default="", type=lambda x: is_valid_path(parser, x),
                        help="Path to a fasta file containing a control assay.")
    parser.add_argument('-cg', '--control_genomes', required=False, nargs="+", default="",
                        help="List of genome paths that corresponds to the control assay. This parameter must be"
                             "given in combination with the control parameter.")
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Output file (.tsv format)')

    args = parser.parse_args()
    return args


def get_files_from_directory(input_path: str) -> Optional[List[Path]]:
    """
    List the .mfe.json files in the given directory. Check if there is at least one file in the given directory.
    :param input_path: String that defines the path of the directory.
    :return: List of file paths.
    """
    mfe_primer_output = list(Path(input_path).glob("*" + MFE_PRIMER_EXTENSION))

    if len(mfe_primer_output) >= 1:
        return mfe_primer_output
    else:
        print(F"{Fore.RED}No .json files in the given directory.")
        quit()


def get_control_primer_names(control_primer_file: str) -> List[str]:
    """
    Read the fasta file of the control primers and save the primer names.
    :param control_primer_file: Path to the fasta file of the control primers.
    :return: List of control primer names.
    """
    control_primers = open(control_primer_file)

    lines_control_primers = control_primers.read().split("\n")
    control_primer_names = [control_primer_name.strip(">") for control_primer_name in lines_control_primers if
                            control_primer_name.startswith(">")]
    control_primers.close()
    return control_primer_names


def make_assay_control_combinations(control_primer_names: List[str]) -> Dict[str, Dict[str, str]]:
    """
    Formulate the possible primer-probe combinations to generate qPCR amplicons.
    :param control_primer_names: List of names for the control primers.
    :return: Nested Dictionary with primer-probe combinations.
    """
    primer_combinations_dict = {"combination_1": {"forward": FORWARD_PRIMER_ASSAY,
                                                  "reverse": control_primer_names[1],
                                                  "probe": PROBE_ASSAY},
                                "combination_2": {"forward": FORWARD_PRIMER_ASSAY,
                                                  "reverse": control_primer_names[1],
                                                  "probe": control_primer_names[2]},
                                "combination_3": {"forward": control_primer_names[0],
                                                  "reverse": REVERSE_PRIMER_ASSAY,
                                                  "probe": PROBE_ASSAY},
                                "combination_4": {"forward": control_primer_names[0],
                                                  "reverse": REVERSE_PRIMER_ASSAY,
                                                  "probe": control_primer_names[2]}
                                }
    return primer_combinations_dict


def find_forward_reverse_primer_amplicons(primer_quality_control_dict: dict, forward_name=FORWARD_PRIMER_ASSAY,
                                          reverse_name=REVERSE_PRIMER_ASSAY) -> Dict[int, dict]:
    """
    Search for amplicons that originate from a forward and reverse primer.
    :param primer_quality_control_dict: .json file converted to a dictionary.
    :param forward_name: Part of name of the forward primer.
    :param reverse_name: Part of name of the reverse primer.
    :return: Nested dictionary with amplicons that are build from a forward and reverse primer.
    The dictionary contains the sequence and genome coordinates.
    """
    fwd_rev_amplicons = {}

    for count, amplicon in enumerate(primer_quality_control_dict["AmpList"]):
        if forward_name in amplicon['F']['Seq']['ID'] and reverse_name in amplicon['R']['Seq']['ID']:
            fwd_rev_amplicons[count] = {
                "seq": amplicon["P"]["Seq"]["Seq"].upper(),
                "end_fwd": amplicon['F']['End'],
                "end_rev": amplicon['R']['End'],
                "dg_fwd": amplicon['F']['Dg'],
                "dg_rev": amplicon['R']['Dg'],
                "seq_true": amplicon["P"]["Seq"]["Seq"].upper(),
                "amplicon_start": amplicon['F']['Start'],
                "amplicon_end": amplicon['R']['Start'],
                "contig": amplicon['Hid'],
                "strand": "+",
                "database": Path(amplicon["DB"])
            }
        elif forward_name in amplicon['R']['Seq']['ID'] and reverse_name in amplicon['F']['Seq']['ID']:
            fwd_rev_amplicons[count] = {
                "seq": amplicon["P"]["Seq"]["Seq"].upper(),
                "end_fwd": amplicon['F']['End'],
                "end_rev": amplicon['R']['End'],
                "dg_fwd": amplicon['F']['Dg'],
                "dg_rev": amplicon['R']['Dg'],
                "seq_true": str(Seq(amplicon["P"]["Seq"]["Seq"].upper()).reverse_complement()),
                "amplicon_start": amplicon['F']['Start'],
                "amplicon_end": amplicon['R']['Start'],
                "contig": amplicon['Hid'],
                "strand": "-",
                "database": Path(amplicon["DB"])
            }
    return fwd_rev_amplicons


def find_primer_probe_amplicons(primer_quality_control_dict: dict, probe_name=PROBE_ASSAY) -> Dict[int, dict]:
    """
    Search for amplicons that originate from a probe and a forward or reverse primer.
    :param primer_quality_control_dict: .json file converted to a dictionary.
    :return: Nested dictionary with amplicons that are build from a forward primer and probe
    or a reverse primer and probe. The dictionary contains the sequence and genome coordinates.
    """
    primer_probe_amplicons = {}

    for count, amplicon in enumerate(primer_quality_control_dict["AmpList"]):
        if probe_name in amplicon['F']['Seq']['ID'] and probe_name not in amplicon['R']['Seq']['ID']:
            primer_probe_amplicons[count] = {
                "seq": amplicon["P"]["Seq"]["Seq"].upper(),
                "probe_start": amplicon['F']['Start'],
                "probe_end": amplicon['F']['End'],
                "dg_probe": amplicon['F']['Dg']
            }
        elif probe_name in amplicon['R']['Seq']['ID'] and probe_name not in amplicon['F']['Seq']['ID']:
            primer_probe_amplicons[count] = {
                "seq": amplicon["R"]["Seq"]["Seq"].upper(),
                "probe_start": amplicon['R']['End'],
                "probe_end": amplicon['R']['Start'],
                "dg_probe": amplicon['R']['Dg']
            }
    return primer_probe_amplicons


def find_qpcr_amplicons(primer_quality_control_dict: dict, analyze_control: dict=None) -> \
        Union[Tuple[int, str, str, dict, Dict[int, dict], List[Path]], None]:
    """
    Determine whether an amplicon is build from a forward and reverse primer and check if the amplicon
    contains a probe. Meaning, is the amplicon detectable with qPCR. Make sure that the probe is within the region
    between the forward and reverse primers. In addition, filter amplicons on qPCR length criteria.
    :param primer_quality_control_dict: .json file converted to a dictionary.
    :param analyze_control: A dictionary with control primer names. If given. Run function for control instead of
    new assay.
    :return: Statistics of amplicons that could be detected by qPCR and the amplicons constructing by solely a
    forward primer and reverse primer.
    Return None if there are no amplicons or no amplicons with a probe.
    """
    amplicons_with_probe = []
    dg_amplicons_fwd, dg_amplicons_rev, dg_amplicons_probe, database = [], [], [], []

    if primer_quality_control_dict["AmpList"] is None:
        return None

    if analyze_control is not None:
        fwd_rev_amplicons = find_forward_reverse_primer_amplicons(primer_quality_control_dict,
                                                                  forward_name=analyze_control["forward"],
                                                                  reverse_name=analyze_control["reverse"])
        primer_probe_amplicons = find_primer_probe_amplicons(primer_quality_control_dict,
                                                             probe_name=analyze_control["probe"])
    else:
        fwd_rev_amplicons = find_forward_reverse_primer_amplicons(primer_quality_control_dict)
        primer_probe_amplicons = find_primer_probe_amplicons(primer_quality_control_dict)

    if fwd_rev_amplicons.__len__() == 0 and primer_probe_amplicons.__len__() == 0:
        return None
    for key1, probe in primer_probe_amplicons.items():
        for key2, amplicon in fwd_rev_amplicons.items():
            if probe["seq"] in amplicon["seq"] and amplicon["seq"] not in amplicons_with_probe and \
                    QPCR_MIN_AMPLICON_LENGTH <= len(amplicon["seq"]) <= QPCR_MAX_AMPLICON_LENGTH and \
                    amplicon["end_fwd"] < probe["probe_start"] < probe["probe_end"] < amplicon["end_rev"]:
                amplicons_with_probe.append(amplicon["seq_true"])
                dg_amplicons_fwd.append(amplicon["dg_fwd"])
                dg_amplicons_rev.append(amplicon["dg_rev"])
                dg_amplicons_probe.append(probe["dg_probe"])
                database.append(amplicon["database"])

    if len(amplicons_with_probe) != 0:
        length_amplicons_qpcr_amplicons = ",".join([str(len(amplicon)) for amplicon in amplicons_with_probe])
        qpcr_amplicon_sequences = ",".join([str(amplicon) for amplicon in amplicons_with_probe])
        number_of_qpcr_amplicons = len(length_amplicons_qpcr_amplicons.split(","))
        dg_amplicons_fwd = ",".join([str(dg) for dg in dg_amplicons_fwd])
        dg_amplicons_rev = ",".join([str(dg) for dg in dg_amplicons_rev])
        dg_amplicons_probe = ",".join([str(dg) for dg in dg_amplicons_probe])
        dg_amplicons = {
            "forward": dg_amplicons_fwd,
            "reverse": dg_amplicons_rev,
            "probe": dg_amplicons_probe
        }

        return number_of_qpcr_amplicons, length_amplicons_qpcr_amplicons, qpcr_amplicon_sequences, dg_amplicons, \
               fwd_rev_amplicons, database
    else:
        return None


def extend_amplicon_upstream_downstream(fwd_rev_amplicon_dict, genome_path):
    """
    Retrieve the amplicon with an additional 50bp upstream and downstream.
    :param fwd_rev_amplicon_dict: Nested dictionary with amplicons that are build from a forward and reverse primer.
    :param genome_path: Path to the genome analyzed by MFEprimer corresponding to the input directory with json files.
    :return: String of extended amplicons separated by a comma.
    """
    extended_amplicons = []
    with open(genome_path) as handle:
        genome = list(SeqIO.parse(handle, "fasta"))
        contig_dictionary = {record.id: record for record in genome}
        for i, amplicon in fwd_rev_amplicon_dict.items():
            contig = contig_dictionary[amplicon["contig"]]
            extended_amplicon = contig.seq[amplicon['amplicon_start'] - BP_AMPLICON_EXTENSION:
                                           amplicon['amplicon_end'] + BP_AMPLICON_EXTENSION]
            if amplicon["strand"] == "-":
                extended_amplicon = extended_amplicon.reverse_complement()
            extended_amplicons.append(extended_amplicon)
    return ",".join(map(str, extended_amplicons))


def process_json(list_files: List[Path], genome_path, file_name: str, control: str = "",
                 control_genomes: List[str] = "") -> Dict[str, dict]:
    """
    Extract valuable information from the MFEprimer-3.2.3 .json output for the phylogenetic group type.
    :param list_files: List of file paths.
    :param genome_path: Path to the genome analyzed by MFEprimer corresponding to the input directory with json files.
    :param file_name: Name of the file to determine to which file the results belong to.
    :param control: Path to a fasta file containing a control assay.
    :param control_genomes: List of genome paths for the control assay.
    :return: Nested dictionary in the format primer_quality_control_dict[<primer_id>][<MFEprimer_output_information>]
    """
    nested_primer_quality_control_dict = {}
    control_is_specific = True

    for file in list_files:
        primer_quality_control_dict = json.load(file.open())

        key = file.name.replace(MFE_PRIMER_EXTENSION, "")
        nested_primer_quality_control_dict[key] = {
            "Left primer sequence": primer_quality_control_dict["PrimerList"][0]["Seq"]["Seq"],
            "Right primer sequence": primer_quality_control_dict["PrimerList"][1]["Seq"]["Seq"],
            "Probe sequence": primer_quality_control_dict["PrimerList"][2]["Seq"]["Seq"],
            "Left primer GC": primer_quality_control_dict["PrimerList"][0]["GC"],
            "Right primer GC": primer_quality_control_dict["PrimerList"][1]["GC"],
            "Probe primer GC": primer_quality_control_dict["PrimerList"][2]["GC"],
            "Left primer Tm": primer_quality_control_dict["PrimerList"][0]["Tm"],
            "Right primer Tm": primer_quality_control_dict["PrimerList"][1]["Tm"],
            "Probe primer Tm": primer_quality_control_dict["PrimerList"][2]["Tm"],
            "No hairpin": primer_quality_control_dict["HairpinList"] == HAIRPIN,
            "No dimer": primer_quality_control_dict["DimerList"] == DIMER,
            "Amplicons": len(find_forward_reverse_primer_amplicons(primer_quality_control_dict))
            if primer_quality_control_dict['AmpList'] is not None else 0,
            "Genome name": file_name
        }
        amplicon_data = find_qpcr_amplicons(primer_quality_control_dict)
        if control and isinstance(control_genomes, list):
            specific = True
            genome_is_control = True

            control_genomes = [Path(name).stem for name in control_genomes]
            control_primer_names = get_control_primer_names(control)
            control_primer_names_dict = {
                "forward": control_primer_names[0],
                "reverse": control_primer_names[1],
                "probe": control_primer_names[2]
            }

            # Determine presence of amplicon for control assay.
            control_amplicon_data = find_qpcr_amplicons(primer_quality_control_dict, analyze_control=control_primer_names_dict)
            if control_amplicon_data:
                number_of_control_amplicons_qpcr = control_amplicon_data[0]
                database_of_control_amplicons_qpcr = control_amplicon_data[5]
                if database_of_control_amplicons_qpcr[0].stem in control_genomes:
                    genome_is_control = number_of_control_amplicons_qpcr == 1
            # Determine presence of amplicons mixed between assay and control.
            assay_control_primer_combinations = make_assay_control_combinations(control_primer_names)
            for primer_combination in assay_control_primer_combinations.values():
                qpcr_amplicon_results = find_qpcr_amplicons(primer_quality_control_dict, analyze_control=primer_combination)
                if qpcr_amplicon_results:
                    specific = False
            control_is_specific = specific and genome_is_control

        if amplicon_data is not None:
            number_of_amplicons_qpcr = amplicon_data[0]
            length_amplicons_qpcr = amplicon_data[1]
            amplicon_sequences_qpcr = amplicon_data[2]
            dg_amplicons = amplicon_data[3]
            fwd_rev_amplicon_dict = amplicon_data[4]
            extended_amplicons = extend_amplicon_upstream_downstream(fwd_rev_amplicon_dict, genome_path)
        else:
            number_of_amplicons_qpcr = None
            length_amplicons_qpcr = None
            amplicon_sequences_qpcr = None
            dg_amplicons = None
            extended_amplicons = None

        nested_primer_quality_control_dict[key]["Number of amplicons qpcr"] = number_of_amplicons_qpcr
        nested_primer_quality_control_dict[key]["Length amplicons qpcr"] = length_amplicons_qpcr
        nested_primer_quality_control_dict[key]["Amplicon sequences qpcr"] = amplicon_sequences_qpcr
        nested_primer_quality_control_dict[key]["Extended amplicon sequences qpcr"] = extended_amplicons
        nested_primer_quality_control_dict[key]["Delta G forward primer"] = None if dg_amplicons is None \
            else dg_amplicons["forward"]
        nested_primer_quality_control_dict[key]["Delta G reverse primer"] = None if dg_amplicons is None \
            else dg_amplicons["reverse"]
        nested_primer_quality_control_dict[key]["Delta G probe primer"] = None if dg_amplicons is None \
            else dg_amplicons["probe"]
        nested_primer_quality_control_dict[key]["Control is specific"] = control_is_specific

        print(F"{Fore.BLUE}DONE: {file.stem}")

    return nested_primer_quality_control_dict


def compute_output_dataframe(nested_dictionary: Dict[str, dict], output: Path):
    """
    Compute a dataframe and write it to an output .tsv file.
    :param nested_dictionary: Nested dictionary in the format
    primer_quality_control_dict[<primer_id>][<MFEprimer_output_information>]
    :param output: String that defines that output path.
    """
    primer_df = pd.DataFrame.from_dict(nested_dictionary).transpose()
    # noinspection PyTypeChecker
    primer_df.to_csv(output, sep="\t")


def main():
    args = parse_args()
    args.output = Path(args.output)
    mfe_primer_output = get_files_from_directory(args.input_dir)
    nested_primer_quality_control_dict = process_json(mfe_primer_output, args.genome, args.output.stem,
                                                      control=args.control, control_genomes=args.control_genomes)
    compute_output_dataframe(nested_primer_quality_control_dict, args.output)
    print(F"{Fore.GREEN}FINISHED")


if __name__ == '__main__':
    main()
