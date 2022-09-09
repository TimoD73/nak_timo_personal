#!/usr/bin/env python

"""
Find genes by CheckM and select them for primer targets.
"""

__author__ = "Timo Dijkstra BSc."
__version__ = "1.0"
__title__ = "Filter genes by CheckM"
__author_email__ = "T.Dijkstra@naktuinbouw.nl"
__date__ = "09-02-2021"

import re
import json
import argparse
import subprocess
import shutil
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import Dict, Any, List
from utils.utils import is_valid_path

TAXON_LIST_FILE = "checkm_markers.txt"
TAXON_MARKER_FILE = "taxon_markers.txt"
MARKER_GENES_FILE = "marker_genes.fasta"
FASTAS_PER_FILE = 50


def parse_args():
    """ Parse command line arguments """

    parser = argparse.ArgumentParser(description='Create a hierarchical tree plot based on output of mash triangle.')
    parser.add_argument('-g', '--genome', type=lambda x: is_valid_path(parser, x), required=True,
                        help='Path to the reference genome.')
    parser.add_argument('-gn', '--genus', type=str, required=True, help='Genus of interest')
    parser.add_argument('-c', '--checkm', type=str, required=True, help='Path to the checkm executable.')
    parser.add_argument('-o', '--output_dir', help='Output directory.', required=True, type=Path)
    parser.add_argument('-t', '--threads', required=False, default="12", type=str,
                        help='Number of threads for running CheckM')

    return parser.parse_args()


def is_marker_set_available(marker_file: Path, genus: str):
    """
    Read the table of CheckM with available markers for each taxa. Check if the requested genus is among them.
    :param marker_file: Path to the marker file produced by the checkm taxon_list command.
    :param genus: Genus of interest.
    :return: True if there is one hit.
    """
    with open(marker_file) as txt:
        data = [line.split("\n") for line in txt.readlines()][7:-2]
        data = [re.split(r"\s{2,}", line[0])[1:-1] for line in data]
    marker_df = pd.DataFrame(data, columns=["Rank", "Taxon", "# genomes", "# marker gene", "# marker sets"])
    available_genera = marker_df[marker_df.Rank == "genus"]["Taxon"]
    number_matching_hits = sum(available_genera.isin([genus]))
    return number_matching_hits == 1


def run_checkm(checkm: str, genus: str, genome: Path, outdir: Path, threads: str):
    """
    Series of commands to run CheckM.
    :param checkm: Path to the checkm executable.
    :param genus: Genus of interest.
    :param genome: Path to the reference genome.
    :param outdir: Output directory.
    :param threads: Number of threads for running CheckM
    """
    Path.mkdir(outdir, exist_ok=True)
    Path.mkdir(outdir / "genome", exist_ok=True)
    shutil.copy(genome, outdir / "genome")
    subprocess.run(f"{checkm} taxon_list > {outdir / TAXON_LIST_FILE}", shell=True)
    genus = genus.capitalize()
    taxon_level = f"genus {genus}" if is_marker_set_available(outdir / TAXON_LIST_FILE, genus) else "domain Bacteria"
    subprocess.run(f"{checkm} taxon_set {taxon_level} {outdir / TAXON_MARKER_FILE}", shell=True)
    subprocess.run(f"{checkm} analyze {outdir / TAXON_MARKER_FILE} {outdir / 'genome'} -x {genome.suffix} -t {threads} {outdir}", shell=True)
    subprocess.run(f"{checkm} qa {outdir / TAXON_MARKER_FILE} {outdir}", shell=True)


def parse_gff(gff: Path) -> Dict[str, Any]:
    """
    Read the gff file from CheckM as dataframe. Convert this dataframe to a dictionary with as key the hit ID and the
    gff information as values.
    :param gff: Path to the gff file from CheckM.
    :return: A dictionary with as key the hit ID and the gff information as values.
    """
    gff_df = pd.read_csv(gff, sep="\t", comment="#",
                         names=["Sequence_ID", "Source", "Feature_type", "Feature_start",
                                "Feature_end", "Score", "Strand", "Phase", "Attributes"])
    col_names = [hit_id.split(";")[0].split("=")[1] for hit_id in gff_df["Attributes"]]
    gff_df = gff_df.loc[:, gff_df.columns != "Attributes"].transpose()
    gff_df.columns = col_names
    gff_dict = gff_df.to_dict()
    return gff_dict


def parse_marker_gene_stats(marker_gene_stats: Path) -> List[str]:
    """
    Parse the marker_gene_stats.tsv file produced by CheckM.
    :param marker_gene_stats: Path to the marker_gene_stats.tsv file.
    :return: List with marker IDs present in the genome of interest found by CheckM.
    """
    marker_gene_stats = open(marker_gene_stats)
    literal_json = marker_gene_stats.read().split("\t")[1].replace("'", '"')
    marker_dict = json.loads(literal_json)
    marker_gene_stats.close()
    marker_gene_list = list(marker_dict.keys())
    return marker_gene_list


def extract_marker_genes(genome: Path, outdir: Path, genome_name):
    """
    Create fasta files per 50 marker genes for the genome of interest.
    :param genome: Path to the reference genome.
    :param outdir: Output directory.
    :param genome_name: String of the name of the genome.
    """
    marker_genes = []
    gff = outdir / "bins" / genome_name / "genes.gff"
    marker_gene_stats = outdir / "storage" / "marker_gene_stats.tsv"
    gff_dict = parse_gff(gff)
    marker_hits = parse_marker_gene_stats(marker_gene_stats)
    marker_hits = [hit.split(".")[1] for hit in marker_hits]
    marker_dict = {hit_id: gff_dict.get(hit_id) for hit_id in marker_hits}
    contigs = {}
    with open(genome) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contigs[record.id] = record
    for marker in marker_dict:
        if marker_dict[marker] is not None:
            contig = contigs.get(marker_dict[marker]["Sequence_ID"]).seq
            # GFF is 1-based. Therefore, convert to 0-based.
            marker_gene = contig[marker_dict[marker]["Feature_start"] - 1:marker_dict[marker]["Feature_end"]]
            marker_gene = marker_gene if marker_dict[marker]["Strand"] == "+" else marker_gene.reverse_complement()
            record = SeqRecord(seq=marker_gene, id=marker.split("_")[-1], description="", name="")
            marker_genes.append(record)
    marker_genes.sort(key=lambda x: int(x.id)) # Sort on ID.

    for i in range(1, len(marker_genes), FASTAS_PER_FILE):
        i_end = i + FASTAS_PER_FILE - 1 if i + FASTAS_PER_FILE - 1 < len(marker_genes) else len(marker_genes)
        with open(outdir / f"{str(i)}_{str(i_end)}_{MARKER_GENES_FILE}", "w") as output_handle:
            fasta_out = SeqIO.FastaIO.FastaWriter(output_handle, wrap=None)
            fasta_out.write_file(marker_genes[i - 1:i + FASTAS_PER_FILE - 1])
    return marker_genes



def main():
    args = parse_args()
    genome = Path(args.genome)
    run_checkm(args.checkm, args.genus, genome, args.output_dir, args.threads)
    extract_marker_genes(args.genome, args.output_dir, genome.stem)


if __name__ == '__main__':
    main()
