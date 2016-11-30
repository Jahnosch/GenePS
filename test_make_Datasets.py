#!/usr/bin/env python3
import unittest
import os
from unittest.mock import patch
import make_Datasets
import tempfile as tmp
import compute_msa
from run_command import tempdir
from compute_msa import msa_operations, MsaObject
from unittest.mock import mock_open

script_path = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(script_path, "test_data/compile_script")
file_name = "eef_3.5_no_cea.fa"


class TestHashFasta(unittest.TestCase):

    fa_hash = make_Datasets.hash_fasta(os.path.join(test_data, file_name))

    @patch("make_Datasets.hash_fasta", return_value={">a": 1, ">b": 2})
    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_check_and_hash_fasta(self, logger, fasta_hash):
        fasta_hash = make_Datasets.check_and_hash_fasta(os.path.join(test_data, file_name), "filename")
        self.assertEqual(fasta_hash, None)

    @patch("make_Datasets.hash_fasta", return_value={})
    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_check_and_hash_fasta(self, logger, fasta_hash):
        fasta_hash = make_Datasets.check_and_hash_fasta(os.path.join(test_data, file_name), "filename")
        self.assertEqual(fasta_hash, None)

    @patch("make_Datasets.hash_fasta", return_value={">a": 1, ">b": 2, ">c": 3})
    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_check_and_hash_fasta(self, logger, fasta_hash):
        fasta_hash = make_Datasets.check_and_hash_fasta(os.path.join(test_data, file_name), "filename")
        self.assertEqual(fasta_hash, {">a": 1, ">b": 2, ">c": 3})

    def test_hash_fasta_equal_number_header_seq(self):
        number_seq = len(self.fa_hash.keys())
        number_header = len(self.fa_hash.values())
        self.assertEqual(number_header, number_seq)

    def test_hash_fasta_number_of_entries(self):
        number_entries = len(self.fa_hash.keys())
        self.assertEqual(number_entries, 190)

    def test_hash_fasta_correct_seq_length(self):
        single_seq = len(self.fa_hash[">AMELL.GB42352-PA"][0])
        self.assertEqual(single_seq, 942)


class TestCleanFastaHash(unittest.TestCase):

    def test_remove_two_keys(self):
        clean_hash = make_Datasets.clean_fasta_hash({">a": 1, ">b": 2, ">c": 3}, [">c"])
        self.assertEqual(clean_hash, {">c": 3})

    def test_key_not_in_hash(self):
        clean_hash = make_Datasets.clean_fasta_hash({">a": 1, ">b": 2, ">c": 3}, [">t", ">p"])
        self.assertEqual(clean_hash, {})

    def test_empty_dictionary(self):
        clean_hash = make_Datasets.clean_fasta_hash({}, [">t", ">p"])
        self.assertEqual(clean_hash, {})


class TestWriteHashToFasta(unittest.TestCase):

    def test_fasta_dictionary(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"})
            self.assertSetEqual(set(open(file_path).readlines()), {'>2\n', 'ADFAT\n', '>1\n', 'ATGSAD\n'})

    def test_empty_hash(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {})
            self.assertEqual(file_path, None)

    def test_rigth_optional_formatting(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"}, line_style="{}\t{}\t")
            self.assertSetEqual(set(open(file_path).readlines()[0].split("\t")), {'>2', 'ADFAT', '>1', 'ATGSAD', ''})

    def test_formatting_without_dictionary_value(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"}, line_style="{}")
            self.assertSetEqual(set(open(file_path).readlines()[0].split(">")), {'', '2', '1'})

    def test_formatting_None_if_more_then_two_brackets(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"}, line_style="{}\n{}\n{}\n")
            self.assertEqual(file_path, None)

    def test_value_is_a_list(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1": ["ATGSAD"], ">2": "ADFAT"})
            self.assertSetEqual(set(open(file_path).readlines()), {'>2\n', 'ADFAT\n', '>1\n', 'ATGSAD\n'})


class TestIntersectionPoint(unittest.TestCase):

    def test_no_point_smaller_tp_mean(self):
        interpoint = make_Datasets.find_intersection(0.585692592593, 0.397771428571, 0.324498504734, 0.21228085651)
        self.assertEqual(interpoint, 0.729)

    def test_extremly_high_point_in_list(self):
        interpoint = make_Datasets.find_intersection(1.10492361111, 0.449, 0.232010418823, 0.327981237737)
        self.assertEqual(interpoint, 0.874)

    def test_best_point_pick_between_too_large_and_perfect(self):
        interpoint = make_Datasets.find_intersection(1.31034831461, 0.620403669725, 0.261504825592, 0.268836948615)
        self.assertEqual(interpoint, 1.038)

    def test_identical_distributions(self):
        interpoint = make_Datasets.find_intersection(500, 500, 300, 300)
        self.assertEqual(interpoint, None)


class TestLengthBinnedFasta(unittest.TestCase):

    def test_estimate_bin_number_never_below_amount_data(self):
        self.assertEqual(make_Datasets.estimate_bin_number(4, 3, 2, forced_bins=3), 2)

    def test_estimate_bin_force_over_guess(self):
        self.assertEqual(make_Datasets.estimate_bin_number(5, 3, 10, forced_bins=4), 4)

    def test_estimate_bin_emit_guess_if_possible(self):
        self.assertEqual(make_Datasets.estimate_bin_number(5, 3, 10), 5)

    def test_estimate_bin_minimum_before_guess(self):
        self.assertEqual(make_Datasets.estimate_bin_number(5, 6, 10), 6)

    def test_bin_sequence_length_into_three(self):
        bin_dict = make_Datasets.bin_sequence_lengths({"0": 0, "4": 4, "12": 12, "16a": 16, "16b": 16, "18": 18, "24": 24, "26": 26, "28": 28}, 3)
        self.assertEqual(bin_dict, {1: '4', 2: '18', 3: '28'})

    def test_huge_window_size_reduces_bin_number(self):
        self.assertEqual(make_Datasets.bin_sequence_lengths({"0": 0, "4": 4, "12": 12, "16a": 16, "16b": 16, "18": 18, "24": 24, "26": 26, "100": 100}, 3), {1: '26', 3: '100'})

    def test_write_length_binned_fasta_correct_length_hash(self):
        with tmp.NamedTemporaryFile() as tempf:
            length_dict = make_Datasets.write_length_binned_fasta({">first" : ["ADFADF"], ">second" : ["DFAF"]}, "test", tempf.name)
            self.assertEqual(length_dict, {'>second': 4, '>first': 6})

    def test_write_length_binned_fasta_raise_error_if_seq_not_list(self):
        with tmp.NamedTemporaryFile() as tempf:
            self.assertRaises(AttributeError, make_Datasets.write_length_binned_fasta, {">first" : "ADFADF", ">second" : "DFAF"}, "test", tempf.name)

    def test_write_length_binned_fasta_converts_from_bin_to_header_into_bin_to_sequence_correctly(self):
        with tmp.NamedTemporaryFile() as tempf:
            length_dict = make_Datasets.write_length_binned_fasta({">first" : ["ADFADF"], ">second" : ["DFAF"]}, "test", tempf.name)
            self.assertSetEqual(set(open(tempf.name).readlines()), {'>1_test\n', 'DFAF\n', '>2_test\n', 'ADFADF\n'})


class TestGetpHmmScores(unittest.TestCase):

    hmm_search_score_file = open(os.path.join(test_data, "eef_hmmer_search.txt"))
    hmm_search_score_list = hmm_search_score_file.readlines()
    hmm_search_score_file.close()

    @patch("make_Datasets.run_cmd", return_value=hmm_search_score_list)
    def test_hmmsearch_probing_score_dict_length(self, run_command):
        score_dict = make_Datasets.get_phmm_score("adf", "dafs")
        self.assertEqual(len(score_dict), 190)
        self.assertEqual(1581, score_dict[">OVOLV.OVOC1360"])

    @patch("make_Datasets.run_cmd", return_value=hmm_search_score_list)
    def test_hmmsearch_length_normalization(self, run_command):
        fasta_hash = make_Datasets.hash_fasta(os.path.join(test_data, file_name))
        len_hash = {header: len(seq[0]) for header, seq in fasta_hash.items()}
        len_hash[">OVOLV.OVOC1360"] = 10
        score_dict = make_Datasets.get_phmm_score("adf", "dafs", header_to_length=len_hash)
        self.assertEqual(score_dict[">OVOLV.OVOC1360"], 1581 / 10)


class TestScoreObject(unittest.TestCase):

    def setUp(self, size=5):
        self.fasta_dict = {">ACRAS.cds.Contig10403m.5077": ["QTFRRIVENINVIIATYGDDDGPMGPIMVDPALGNVGFGSGLHGWAFTLKQFAEMYASKFGVQVDKLMKNLWGDRFFNMKTKKRSTSQEDGAVRGFTQFVLDPIFKVF*"],
                  ">ALUMB.ALUE_0000951001-mRNA-1": ["XVCVQTETVLRQAIAERIKPVLFMNKMDRALLELQLGQEELYQTFQRIVENTNVIIATYGDDDGPMGQIMVDPAIGNVGFGSGLHGWAFTLKQFAEMYSEKFGVQ"],
                  ">ACRAS.cds.Contig3658m.1561": ["MEMLYEGPHDDEVAVAIKNCDPNGPLMMYVSKMVPTSDKGRFYAFGRVFSGKVATGMKARIQGPNYVPGKKEDLYEKTIQRTILMMGRYVEPIEDIPSGNIAGLVGVDQYLVKGGTITTFKDAHNLRVMKFSVSPVVRVAVEPKNAGDLPKLVEGLKRLSKSDPMV"]}
        self.score_obj = make_Datasets.ScoreObject(self.fasta_dict, os.path.join(test_data, "eef.hmm"))

    def test_iterative_scoring_hand_checked_score(self):
        score_dict = self.score_obj.iterative_score_computation()
        self.assertEqual(score_dict[">ALUMB.ALUE_0000951001-mRNA-1"], 37)
        self.assertNotIn(">ACRAS.cds.Contig3658m.1561", score_dict)

    def test_iterative_scoring_length_normalized(self):
        length_dict = {">ACRAS.cds.Contig10403m.5077": 10, ">ALUMB.ALUE_0000951001-mRNA-1": 10, ">ACRAS.cds.Contig3658m.1561": 10}
        self.score_obj.length_dict = length_dict
        score_dict = self.score_obj.iterative_score_computation(length_normalized=True)
        self.assertEqual(score_dict[">ALUMB.ALUE_0000951001-mRNA-1"], 37 / 10)

    def test_bulk_scoring_hand_checked(self):
        score_dict = self.score_obj.bulk_score_computation()
        self.assertEqual(score_dict[">ALUMB.ALUE_0000951001-mRNA-1"], 193)

    def test_bulk_scoring_hand_checked_length_normalized(self):
        length_dict = {">ACRAS.cds.Contig10403m.5077": 10, ">ALUMB.ALUE_0000951001-mRNA-1": 10, ">ACRAS.cds.Contig3658m.1561": 10}
        self.score_obj.length_dict = length_dict
        score_dict = self.score_obj.bulk_score_computation(length_normalized=True)
        self.assertEqual(score_dict[">ALUMB.ALUE_0000951001-mRNA-1"], 193 / 10)

    def test_length_distribution_parameters(self):
        self.score_obj.length_dict = {">ACRAS.cds.Contig10403m.5077": 10, ">ALUMB.ALUE_0000951001-mRNA-1": 20, ">ACRAS.cds.Contig3658m.1561": 30}
        self.assertListEqual(list(self.score_obj.calculate_length_distribution_parameters()), [10.0, 30.0])

    def test_score_distribution_return_minimum(self):
        self.score_obj.score_dict = {"4": 4, "12": 12, "16": 16, "18": 18, "24": 24, "26": 26}
        self.assertEqual(4, self.score_obj.calculate_score_distribution_parameters())

    def test_score_distribution_return_intersection(self):
        self.score_obj.score_dict = {"4": 4, "12": 12, "16": 16, "18": 18, "24": 24, "26": 26}
        self.assertEqual(13.35, self.score_obj.calculate_score_distribution_parameters(true_negative_scores=[2, 3, 9, 10, 5, 11]))

    def test_score_distribution_too_less_values(self):
        self.score_obj.score_dict = {"4": 4, "12": 12, "16": 16, "18": 18, "24": 24, "26": 26}
        self.assertEqual(4, self.score_obj.calculate_score_distribution_parameters(true_negative_scores=[2, 3, 9]))


class TestMsaObject(unittest.TestCase):

    msa_path = open(os.path.join(test_data, "eef.aln"))
    msa_list = msa_path.readlines()
    msa_path.close()

    @patch("compute_msa.run_cmd", return_value=msa_list)
    def test_msa_operations(self, run_cmd_mock):
        msa_list = compute_msa.msa_operations("blaaa")
        self.assertEqual(len(msa_list)/2, 190)




if __name__ == '__main__':
    unittest.main()



