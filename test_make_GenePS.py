#!/usr/bin/env python3
import unittest
import os
import re
from make_GenePS import walk_through_input, check_size_len_after_trimal, InputFileError, MsaLengthError, MsaSizeError, \
     hash_fasta, ScoreObject, generate_hmm, get_consensus, get_phmm_score, check_for_sufficient_taxa, get_outdir
from run_command import tempdir
from compute_msa import msa_operations, MsaObject

script_path = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(script_path, "test_data")
single_file = os.path.join(test_data, "group1/eef_test.fa")
file_name = "eef_test"

msa_test_2 = open(os.path.join(test_data, "group1/msa_test_2.fa"))
msa_list = [line.strip() for line in msa_test_2]


class TestGetOutDir(unittest.TestCase):

    def test_get_outdir_if_exits(self):
        self.assertTrue(get_outdir(test_data) == test_data + "/")

    def test_get_outdir_file_error(self):
        self.assertRaises(SystemExit, lambda: get_outdir(single_file))

    def test_get_outdir_added_dir(self):
        folder_str = "test_folder"
        new_outdir = get_outdir(test_data, folder_str)
        self.assertEqual(os.path.join(test_data, folder_str), new_outdir)
        os.rmdir(new_outdir)


class TestWalkThroughInput(unittest.TestCase):

    dir_tree = walk_through_input(test_data)
    file_tree = walk_through_input(single_file)
    single_folder = walk_through_input(os.path.join(test_data, "group1"))

    def test_input_is_file(self):
        for path, file in self.file_tree.items():
            file_path = os.path.join(path, file[0])
            self.assertTrue(file_path == single_file)

    def test_number_of_groups(self):
        self.assertTrue(len(self.dir_tree), 2)

    def test_number_of_file(self):
        num_files = 0
        for folder, file_list in self.dir_tree.items():
            num_files += len(file_list)
        self.assertEqual(num_files, 14)

    def test_folder_name(self):
        for folder in self.dir_tree.keys():
            folder = folder.split("/")[-1]
            self.assertIn(folder, ["group1", "group2", "test_data"])

    def test_file_name_true(self):
        kog_file = self.dir_tree[os.path.join(test_data, "group2")][0]
        self.assertRegex(kog_file, "KOG0018.fas")


class TestHashFasta(unittest.TestCase):

    fa_hash = hash_fasta(single_file)

    def test_hash_fasta_correct_header(self):
        self.assertIn('>ASUUM1.ASU_00947', self.fa_hash.keys())

    def test_check_for_sufficient_taxa_enough_taxa(self):
        self.assertWarns(None, check_for_sufficient_taxa(self.fa_hash))

    def test_check_for_sufficient_taxa_should_raise_error(self):
        false_hash = {"1": "size", "2": "less_than", "3": "four"}
        self.assertRaises(InputFileError, lambda: check_for_sufficient_taxa(false_hash))

    def test_hash_fasta_correct_sequence(self):
        seq_list = [x[0] for x in self.fa_hash.values()]
        for seq in seq_list:
            self.assertNotIn(">", seq)

    def test_hash_fasta_equal_number_header_seq(self):
        number_seq = len(self.fa_hash.keys())
        number_header = len(self.fa_hash.values())
        self.assertEqual(number_header, number_seq)

    def test_hash_fasta_number_of_entries(self):
        number_entries = len(self.fa_hash.keys())
        self.assertEqual(number_entries, 11)

    def test_hash_fasta_correct_seq_length(self):
        single_seq = len(self.fa_hash[">AMELL.GB42352-PA"][0])
        self.assertEqual(single_seq, 942)


class TestMsaObject(unittest.TestCase):

    out_dir = os.path.join(test_data, "group1")

    def test_generate_msa(self): # check own test input
        self.assertEqual((len(msa_list)/2), 11)

    def test_msa_obj_type(self):
        with tempdir() as tmp_dir:
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            self.assertTrue(type(msa_obj) == MsaObject)

    def test_msa_path(self):
        with tempdir() as tmp_dir:
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            self.assertTrue(msa_obj.path == os.path.join(tmp_dir, "eef_test.msaGenePS"))

    def test_msa_operations(self):
        msa_test_file = os.path.join(test_data, "group1/msa_eef.fa")
        should_be_msa_list = msa_operations("less " + msa_test_file)
        self.assertListEqual(should_be_msa_list, msa_list)

    def test_check_size_len_after_trimal_no_error(self):
        with tempdir() as tmp_dir:
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            msa_obj.size.append(6)
            msa_obj.lengths.append(746)
            self.assertRaises(None, check_size_len_after_trimal(msa_obj))

    def test_check_size_len_after_trimal_expect_size_error(self):
        with tempdir() as tmp_dir:
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            msa_obj.size.append(3)
            msa_obj.lengths.append(746)
            self.assertRaises(MsaSizeError, lambda: check_size_len_after_trimal(msa_obj))

    def test_check_size_len_after_trimal_expect_length_error(self):
        with tempdir() as tmp_dir:
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            msa_obj.size.append(6)
            msa_obj.lengths.append(19)
            self.assertRaises(MsaLengthError, lambda: check_size_len_after_trimal(msa_obj))

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_trim_remove(self):
        with tempdir() as tmp_dir:
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            msa_obj.msa_to_fasta()
            msa_obj.trim_remove()
            org_size, trim_size = msa_obj.size[0], msa_obj.size[1]
            self.assertEqual(int((org_size - trim_size)), 5)

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_trim_length(self):
        with tempdir() as tmp_dir:
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            msa_obj.msa_to_fasta()
            msa_obj.trim_length()
            org_len, trim_len = msa_obj.lengths[0], msa_obj.lengths[1]
            self.assertEqual(int((org_len - trim_len)), 372)

    def test_all_header(self):
        with tempdir() as tmp_dir:
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            msa_obj.msa_to_fasta()
            header_list = msa_obj.all_header()
            for header in header_list:
                self.assertIn(">", header)

    def test_all_aln(self):
        with tempdir() as tmp_dir:
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            msa_obj.msa_to_fasta()
            seq_list = msa_obj.all_aln()
            regex = re.compile('[^a-zA-Z]')
            for seq in seq_list:
                seq = seq.replace("-", "")
                self.assertFalse(re.match(regex, seq))


class TestScoreObject(unittest.TestCase):

    fa_hash = hash_fasta(single_file)
    eef_scores = [1981, 1808, 1893, 1981, 1981, 1889]
    left_taxa = ['>MINCO1.Minc01141 850 bp', '>PREDI.g17259.t1 850 bp',
                 '>PJUL7.cds.JU765comp10941_c0_seq1m.15227 850 bp', '>MINCO2.g11856.t1 850 bp',
                 '>MJAVA.g3839.t1 850 bp', '>ASUUM1.ASU_00947 850 bp']
    test_consensus = "MVNFTVDEIRALMDKKKNIRNMSVIAHVDHGKSTLTDSLVSKAGIIAGAKAGETRFTDTRKDEQDRCITIKSTAISLFFELD" \
                     "EKDLDFVKGDNQIDIVDGAKKKYNGFLINLIDSPGHVDFSSEVTAALRVTDGALVVVDCVSGVCVQTETVLRQAIAERIKPV" \
                     "LFMNKMDRALLELQLGQEELYQTFQRIVENINVIIATYGDDDGPMGAIQVDPAIGNVGFGSGLHGWAFTLKQFAEMYADKFG" \
                     "VQVDKLMKNLWGDRFFNLKTKKWTSTQEDDTKRGFVQFVLDPIFKVFDAVMNVKKEETTKLVEKLNVKLAAEEKDLEGKALL" \
                     "KVLMRKWLPAGDTMLQMICIHLPSPVTAQKYRMEMLYEGPHDDEAAIAIKACDPNGPLMMYISKMVPTSDKGRFYAFGRVFS" \
                     "GKVATGMKARIQGPNYVVGKKEDLYEKTIQRTILMMGRYIEPIEDIPAGNIAGLVGVDQYLVKGGTITTFKDAHNLRVMKFS" \
                     "VSPVVRVAVEPKNAGDLPKLVEGLKRLAKSDPMVQCIFEESGEHIIAGAGELHLEICLKDLEEDHACIPIKKSDPVVSYRET" \
                     "VTEESDIMCLSKSPNKHNRLFCKAKPLADGLAEAIEKGEVSARDEAKNRAKILAEKYEFDATDARKIWCFGPDGTGANLLVD" \
                     "VTKGVQYLNEIKDSVVAGFQWATKEGVLCDENLRGVRFDIHDVTLHADAIHRGGGQIIPTARRVLYASVLTAKPRLLEPVYL" \
                     "VEIQCPEAAVGGIYGVLNRRRGVVFEESQIAGTPMFIVKAYLPVNESFGFTADLRSNTGGQAFPQCVFDHWQILPGDPLEST" \
                     "SKPAQVVAETRKRKGLKEGIPALDNFLDKL"

    def test_obj_name(self):
        with tempdir() as tmp_dir:
            scoring_obj = ScoreObject(self.fa_hash, self.left_taxa, file_name, tmp_dir)
            self.assertTrue(scoring_obj.name == file_name + ".hmmGenePS")

    def test_get_phmm_score(self):
        set_score = 1777.0
        test_hmm_search = os.path.join(test_data, "group1/hmmsearch_eef.test")
        command = ["less", test_hmm_search]
        score_to_test = get_phmm_score("no", "no", command)
        self.assertEqual(set_score, score_to_test)

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_phmm_path_not_none(self):
        with tempdir() as tmp_dir:
            scoring_obj = ScoreObject(self.fa_hash, self.left_taxa, file_name, tmp_dir)
            self.assertFalse(scoring_obj.hmm_path)
            scoring_obj.compute_full_phmm()
            self.assertTrue(scoring_obj.hmm_path)

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_correct_score_list(self):
        with tempdir() as tmp_dir:
            scoring_obj = ScoreObject(self.fa_hash, self.left_taxa, file_name, tmp_dir)
            score_list = scoring_obj.compute_scores()
            self.assertListEqual(score_list, self.eef_scores)

    def test_query_to_fasta(self):
        with tempdir() as tmp_dir:
            scoring_obj = ScoreObject(self.fa_hash, self.left_taxa, file_name, tmp_dir)
            fasta_string = scoring_obj.query_for_fasta(">MINCO1.Minc01141 850 bp")
            control = ">MINCO1.Minc01141" + "\n" + "".join(self.fa_hash[">MINCO1.Minc01141"])
            self.assertTrue(fasta_string == control)

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_get_consensus(self):
        with tempdir() as tmp_dir:
            cons_hmm = os.path.join(tmp_dir, file_name + ".chmm")
            msa_obj = MsaObject(msa_list, file_name, tmp_dir)
            msa_obj.msa_to_fasta()
            msa_obj.trim_remove()
            msa_obj.trim_length()
            generate_hmm(cons_hmm, msa_obj.path)
            consensus_seq = get_consensus(cons_hmm)
            self.assertTrue(consensus_seq == self.test_consensus)


if __name__ == '__main__':
    unittest.main()



