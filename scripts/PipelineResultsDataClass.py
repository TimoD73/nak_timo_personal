from dataclasses import dataclass
from typing import List, Union
import numpy as np

GUANINE = "G"
QPCR_N_AMPLICONS_INGROUP = 1
QPCR_N_AMPLICONS_OUTGROUP = 0


@dataclass
class GeneralResults:
    prokka_protein: str
    t3sepp_prediction_score: Union[float, str]
    forward_primer_sequence: str
    reverse_primer_sequence: str
    probe_sequence: str
    forward_primer_length: int
    reverse_primer_length: int
    probe_primer_length: int
    gc_forward_primer: float
    gc_reverse_primer: float
    gc_probe_primer: float
    product_size: int
    probe_does_not_start_with_g: bool = False
    no_dimers_and_no_hairpins: Union[bool, str] = "NA"
    specific_with_control: bool = False
    great_primers_ingroup: float = 0
    great_primers_outgroup: float = 0

    def check_five_prime_probe(self):
        """
        A guanine has a mild quenching effect.
        So, check if the 5' site does not contain a guanine.
        """
        self.probe_does_not_start_with_g = self.probe_sequence.upper()[0] != GUANINE


@dataclass
class Primer3:
    forward_tm: float
    reverse_tm: float
    probe_tm: float
    forward_self_complementarity: float
    reverse_self_complementarity: float
    probe_self_complementarity: float
    forward_self_homodimer: float
    reverse_self_homodimer: float
    probe_self_homodimer: float
    forward_hairpin: float
    reverse_hairpin: float
    probe_hairpin: float
    forward_end_stability: float
    reverse_end_stability: float
    pair_complementarity: float
    pair_complementarity_end: float
    good_delta_tm_primers: bool = False
    good_delta_tm_probe_primer: bool = False

    def check_tm_primers(self, delta_tm_primer: float):
        """
        Check whether the primer tm difference between two primers is within the given range.
        :param delta_tm_primer: A maximum temperature range between two primers.
        :return: True if the tm difference is within the maximum temperature range.
        """
        self.good_delta_tm_primers = abs(self.forward_tm - self.reverse_tm) <= delta_tm_primer

    def check_tm_probe(self, delta_tm_probe_primer_min: float, delta_tm_probe_primer_max: float):
        """
        Check whether the probe tm is within a given temperature range calculated from the lowest primer tm.
        :param delta_tm_probe_primer_min: Minimum value of the temperature range.
        :param delta_tm_probe_primer_max: Maximum value of the temperature range.
        :return: True if the probe's tm is not below or over the temperature range
        calculated from the primer with the lowest tm.
        """
        self.good_delta_tm_primers = min(self.forward_tm, self.reverse_tm) + delta_tm_probe_primer_min <= self.probe_tm\
                                     <= min(self.forward_tm, self.reverse_tm) + delta_tm_probe_primer_max


@dataclass
class MfePrimer:
    forward_tm: float
    reverse_tm: float
    probe_tm: float
    no_hairpin: bool
    no_dimer: bool
    control_specific_ingroup: Union[bool, str]
    control_specific_outgroup: Union[bool, str]
    number_of_amplicons_ingroup: List[int]
    number_of_qpcr_detectable_amplicons_ingroup: list
    lengths_of_qpcr_detectable_amplicons_ingroup: List[list]
    number_of_amplicons_outgroup: List[int]
    number_of_qpcr_detectable_amplicons_outgroup: list
    lengths_of_qpcr_detectable_amplicons_outgroup: Union[List[list], None]
    delta_g_forward_ingroup_qpcr_detectable_amplicons: list
    delta_g_reverse_ingroup_qpcr_detectable_amplicons: list
    delta_g_probe_ingroup_qpcr_detectable_amplicons: list
    delta_g_forward_outgroup_qpcr_detectable_amplicons: list
    delta_g_reverse_outgroup_qpcr_detectable_amplicons: list
    delta_g_probe_outgroup_qpcr_detectable_amplicons: list
    good_delta_tm_primers: bool = False
    good_delta_tm_probe_primer: bool = False

    def check_tm_primers(self, delta_tm_primer: float):
        """
        Check whether the tm values are within the given range.
        :param delta_tm_primer: A maximum temperature range between two primers.
        :return: True if the temperature range is within the maximum temperature range.
        """
        self.good_delta_tm_primers = abs(self.forward_tm - self.reverse_tm) <= delta_tm_primer

    def check_tm_probe(self, delta_tm_probe_primer_min: float, delta_tm_probe_primer_max: float):
        """
        Check whether the probe tm is within a given temperature range calculated from the lowest primer tm.
        :param delta_tm_probe_primer_min: Minimum value of the temperature range.
        :param delta_tm_probe_primer_max: Maximum value of the temperature range.
        :return: True if the probe's tm is not below or over the temperature range
        calculated from the primer with the lowest tm.
        """
        self.good_delta_tm_probe_primer = min(self.forward_tm, self.reverse_tm) + delta_tm_probe_primer_min <= \
                                          self.probe_tm <= min(self.forward_tm, self.reverse_tm) + \
                                          delta_tm_probe_primer_max


@dataclass
class GenomeNames:
    Genome_names_ingroup: List[str]
    Genome_names_outgroup: List[str]


@dataclass
class Amplicons:
    amplicon_target_ingroup: str
    amplicon_target_ingroup_extended: str
    qpcr_detectable_off_target_amplicon_sequences_ingroup: Union[list, None]
    qpcr_detectable_off_target_amplicon_sequences_ingroup_extended: Union[list, None]
    qpcr_detectable_amplicon_sequences_outgroup: Union[list, None]
    qpcr_detectable_amplicon_sequences_outgroup_extended: Union[list, None]


@dataclass
class PipelineResults:
    general_results: GeneralResults
    primer3: Primer3
    mfe_primer: MfePrimer
    genome_names: GenomeNames
    amplicons: Amplicons

    def check_great_primers(self):
        """
        Determine whether there are no dimers and hairpins found by MFEprimer. In addition return the percentage of
        genomes in the ingroup that produces only one qPCR detectable amplicon and the percentage of outgroup genomes
        that produces 0 qPCR detectable amplicons.
        """
        self.general_results.no_dimers_and_no_hairpins = self.mfe_primer.no_hairpin and self.mfe_primer.no_dimer
        self.general_results.great_primers_ingroup = sum([length == QPCR_N_AMPLICONS_INGROUP for length in
                                                        self.mfe_primer.number_of_qpcr_detectable_amplicons_ingroup])/ \
                                                     len(self.mfe_primer.number_of_qpcr_detectable_amplicons_ingroup) \
                                                     * 100
        self.general_results.great_primers_outgroup = sum([length == QPCR_N_AMPLICONS_OUTGROUP for length in
                                                        self.mfe_primer.number_of_qpcr_detectable_amplicons_outgroup]) \
                                                   / len(self.mfe_primer.number_of_qpcr_detectable_amplicons_outgroup) \
                                                      * 100

    def check_specific_with_control(self):
        """
        Check if the assay is specific with the control.
        """
        if isinstance(self.mfe_primer.control_specific_ingroup, np.bool_) and \
                isinstance(self.mfe_primer.control_specific_outgroup, np.bool_):
            self.general_results.specific_with_control = self.mfe_primer.control_specific_ingroup and \
                                                         self.mfe_primer.control_specific_outgroup
        else:
            self.general_results.specific_with_control = "NA"
