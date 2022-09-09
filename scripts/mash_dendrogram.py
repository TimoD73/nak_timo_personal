#!/usr/bin/env python

"""
compute a dendrogram.
"""

__author__ = "Ronald de Jongh MSc, Timo Dijkstra BSc."
__version__ = "1.0"
__title__ = "Compute a dendrogram for the ingroup and outgroup."
__author_email__ = "r.d.jongh@naktuinbouw.nl, T.Dijkstra@naktuinbouw.nl"
__date__ = "23-11-2021"


import json
import argparse
import subprocess
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
from typing import Tuple, List
from pathlib import Path
from utils.utils import is_valid_directory
from scipy.cluster import hierarchy

GENOME_EXTENSIONS = ["*.fna", "*.fasta", "*.fa"]
THREADS = 16
MASH_SIZE = 10000
INPUT_FILENAME = "Input_genome_list"
INPUT_EXTENSION = ".txt"
MASH_SKETCH_FILENAME = "Mash_sketch"
MASH_SKETCHFILE_EXTENSION = ".msh"
MASH_DISTANCE_FILENAME = "Mash_distance"
MASH_DIST_EXTENSION = ".mashdist"
REFSEQ_ASSEMBLY_ACCESSION = "GCF"
LABELS_FILE = "Labels.csv"
LABEL_JSON = "Labels.json"
LABELS_INGROUP_FILE = "Labels_ingroup.csv"
DENDROGRAM_FILE_NAME = "Mash_tree"
INCHES_PER_LABEL = 0.2


def parse_args():
    """ Parse command line arguments """

    parser = argparse.ArgumentParser(description='Create a hierarchical tree plot based on output of mash triangle.')
    parser.add_argument('-ig', '--ingroup', type=lambda x: is_valid_directory(parser, x), required=True,
                        help='Directory with ingroup genomes.')
    parser.add_argument('-og', '--outgroup', type=lambda x: is_valid_directory(parser, x), required=True,
                        help='Directory with outgroup genomes.')
    parser.add_argument('-g', '--genus', type=str, required=False, default="",
                        help="The genus of the ingroup and outgroup, which captures the dendrogram's title.")
    parser.add_argument('-m', '--metadata', required=False, type=str, default="",
                        help="Metadata file for the genomes used as input for the dendrogram. The file is an in-house"
                             "file for Naktuinbouw and can be found at /0_db_hdd/primer_tool/genus/data.tsv.")
    parser.add_argument('-o', '--output', type=lambda x: is_valid_directory(parser, x), required=False, default="./",
                        help='Output directory of the Mash tree.')

    return parser.parse_args()


def list_input_files(ingroup_dir: Path, outgroup_dir: Path, input_genomes_path: Path) -> List[Path]:
    """
    Write all input genome paths to an output text file.
    :param ingroup_dir: Directory containing all ingroup genomes.
    :param outgroup_dir: Directory containing all outgroup genomes.
    :param input_genomes_path: Output path of the text file.
    :return: The genome paths for the ingroup.
    """
    genome_paths_ingroup = [file for ext in GENOME_EXTENSIONS for file in ingroup_dir.rglob(ext)]
    genome_paths_outgroup = [file for ext in GENOME_EXTENSIONS for file in outgroup_dir.rglob(ext)]
    genome_paths = genome_paths_ingroup + genome_paths_outgroup

    with open(input_genomes_path, "w") as txt:
        for genome_path in genome_paths:
            txt.write(str(genome_path) + "\n")

    return genome_paths_ingroup


def make_db_class(ingroup_dir: Path, outgroup_dir: Path, output: Path) -> Tuple[List[Path], Path]:
    """
    Run Mash to get a distance matrix for all ingroup and outgroup genomes.
    :param ingroup_dir: Directory with ingroup genomes.
    :param outgroup_dir: Directory with outgroup genomes.
    :param output: Output directory of the Mash tree.
    :return: The genome paths for the ingroup and the path to the file with distances for the genomes.
    """
    input_genomes_path = (output / INPUT_FILENAME).with_suffix(INPUT_EXTENSION)
    output_mash_path = (output / MASH_SKETCH_FILENAME).with_suffix(MASH_SKETCHFILE_EXTENSION)
    output_dist_path = (output / MASH_DISTANCE_FILENAME).with_suffix(MASH_DIST_EXTENSION)

    ingroup_files = list_input_files(ingroup_dir, outgroup_dir, input_genomes_path)
    subprocess.run(F"mash sketch -p {THREADS} -o {output_mash_path} -l {input_genomes_path}", shell=True)
    subprocess.run(F"mash triangle -p {THREADS} {output_mash_path}  > {output_dist_path}", shell=True)
    return ingroup_files, output_dist_path


def retrieve_organism_name(metadata: dict, genome_names: List[str], outpath: Path = "", save_metadata: bool = False,
                           save_ingroup: bool = False) -> List[str]:
    """
    Function that finds the organism name from the metadata file if available.
    :param metadata: Metadata dictionary for the genomes used as input for the dendrogram. The file is an in-house file
    for Naktuinbouw and can be found at /0_db_hdd/primer_tool/genus/data.tsv.
    :param genome_names: List of labels derived from the filename of the genome.
    :param outpath: Output directory of the files.
    :param save_metadata: Boolean whether the metadata should be saved.
    :param save_ingroup: Boolean whether the ingroup names should be saved.
    :return: Newly formatted labels in which the names are replaced by the organism name if available.
    """
    labels = []
    label_dict = {}
    assembly_accession_dict = metadata.get("assembly_accession")
    # If the column name does not exist, the assembly accession number is the first column.
    assembly_accession_dict = metadata.get("Unnamed: 0") if assembly_accession_dict is None else assembly_accession_dict
    # If the dictionary is still None, return the genome names.
    if assembly_accession_dict is None:
        return genome_names

    organism_name_dict = metadata.get("organism_name")
    for genome_name in genome_names:
        list_genome_name_elements = genome_name.split("_")
        if REFSEQ_ASSEMBLY_ACCESSION in list_genome_name_elements:
            presumed_assembly_accession = list_genome_name_elements.index(REFSEQ_ASSEMBLY_ACCESSION) + 1
            assembly_accession = f"{REFSEQ_ASSEMBLY_ACCESSION}_{list_genome_name_elements[presumed_assembly_accession]}"
            try:
                column_index = list(assembly_accession_dict.keys())[list(assembly_accession_dict.values()).index(assembly_accession)]
                organism_name = organism_name_dict.get(column_index)
            except ValueError:      # Catch error if assembly accession is not in list.
                organism_name = None
            if organism_name is not None:
                labels.append(f"{organism_name} ({assembly_accession})")
                label_dict[assembly_accession] = organism_name
            else:
                labels.append(genome_name)
        else:
            labels.append(genome_name)

    if outpath and save_metadata:
        with open(outpath / LABELS_FILE, "w") as csv:
            csv.write(",".join(labels))
        with open(outpath / LABEL_JSON, "w") as js:
            json.dump(label_dict, js)
    if outpath and save_ingroup:
        with open(outpath / LABELS_INGROUP_FILE, "w") as csv:
            csv.write(",".join(labels))

    return labels


def make_tree_plot(lt_distmatrix_path: Path, ingroup: List[Path], outpath: str, title: str, metadata: pd.DataFrame = "",
                   img_fmt: str = ".png"):
    """
    Function that reads the file produced by `mash triangle` and compiles the results to a dendrogram.
    :param lt_distmatrix_path: The path to the file with distances for the genomes.
    :param ingroup: The genome paths for the ingroup.
    :param outpath: Output directory of the Mash tree.
    :param title: Genus that classifies the ingroup and outgroup.
    :param metadata: Metadata dataframe for the genomes used as input for the dendrogram. The file is an in-house file
    for Naktuinbouw and can be found at /0_db_hdd/primer_tool/genus/data.tsv.
    :param img_fmt: Extension that specifies the figure output.
    """
    ingroup_file_names = [file.stem for file in ingroup]
    output = (Path(outpath) / DENDROGRAM_FILE_NAME).with_suffix(img_fmt)
    distances = []
    seqnames = []

    with open(lt_distmatrix_path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue

            # Manipulate the file names to usable xlabels.
            l = line.strip().split()
            file_name = Path(l[0])
            seqnames.append(file_name.name.strip(file_name.suffix))

            # Get the distances from the input file.
            distances.append([float(_) for _ in l[1:]] + [0.0, ])
            for j in range(i - 1):
                distances[j].append(distances[i - 1][j])

    if isinstance(metadata, pd.DataFrame):
        seqnames = retrieve_organism_name(metadata.to_dict(), seqnames, outpath=Path(outpath), save_metadata=True)
        ingroup_file_names = retrieve_organism_name(metadata.to_dict(), ingroup_file_names, outpath=Path(outpath),
                                                    save_ingroup=True)

    arr = np.array(distances)
    # convert the redundant n*n square matrix form into a condensed nC2 array
    Z = hierarchy.linkage(ssd.squareform(arr))
    # Style the figure.
    hierarchy.set_link_color_palette(['b', 'g', 'y', 'm'])
    fig, ax = plt.subplots(1, 1, figsize=(16 if INCHES_PER_LABEL*len(seqnames) < 16
                                          else INCHES_PER_LABEL*len(seqnames), 8))
    dn1 = hierarchy.dendrogram(Z, ax=ax, above_threshold_color='#bcbddc', orientation='top',
                               labels=np.array(seqnames), color_threshold=1)
    hierarchy.set_link_color_palette(None)  # reset to default after use

    ax.set_ylabel('Distance (mash)', fontsize=20)
    ax.set_xlabel('Genome', fontsize=20)

    # Calculate the distances between labels and find the center coordinate.
    labels = list(ax.get_xticklabels())
    x_ticks_coordinates = ax.get_xticks()
    label_distance = x_ticks_coordinates[1] - x_ticks_coordinates[0]
    label_center = math.floor(label_distance / 2)

    for index, label in enumerate(labels):
        label.set_rotation(45)
        label.set_ha('right')
        if label.get_text() in ingroup_file_names:
            plt.axvspan(x_ticks_coordinates[index] - label_center,
                        x_ticks_coordinates[index] - label_center + label_distance,
                        color="#20d60f", alpha=0.5)    # Add vertical span for the ingroup.
            label.set_color("#084203")                 # Change xlabel color for ingroup.
    if title:
        plt.title(title, fontsize=25)
    fig.tight_layout()
    plt.savefig(output, dpi=200)


def main():
    args = parse_args()
    ingroup_files, output_dist_path = make_db_class(Path(args.ingroup), Path(args.outgroup), Path(args.output))
    if args.metadata:
        metadata = pd.read_csv(args.metadata, sep="\t")
        make_tree_plot(output_dist_path, ingroup_files, args.output, args.genus, metadata=metadata)
    else:
        make_tree_plot(output_dist_path, ingroup_files, args.output, args.genus)


if __name__ == '__main__':
    main()

