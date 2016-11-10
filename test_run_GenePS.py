#!/usr/bin/env python3
import unittest
import os
from unittest.mock import patch
from unittest.mock import MagicMock
from run_command import tempdir
from run_GenePS import Overseer, DataProviderObject, check_arguments, coverage_filter
from find_regions import read_blast_output, parse_blastdb, make_blast_db, BlastObject, HspListObject
from exonerate_parser import remove_non_letter_signs, clear_hashed_bases, grap_values, aacode_3to1, ExonerateObject

script_path = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(script_path, "test_data")


class TestBlastObject(unittest.TestCase):

    blast_file = os.path.join(test_data, "remanei_OG1685.blast")
    blast = read_blast_output(blast_file, os.path.join(test_data,"c_brenneri.PRJNA20035.WS249.genomic.fa"))
    blast.infer_regions()

    def test_inferred_regions_single_hsp_in_hsp_list(self):
        single_hsp_region = self.blast.inferred_regions["Cbre_Contig419"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region.s_start, 60160)
        self.assertEqual(single_hsp_region.s_end, 70339)
        self.assertEqual(single_hsp_region.strand, "+")
        self.assertEqual(single_hsp_region.chunk_cov, 100)
        self.assertEqual(single_hsp_region.query_cov, 11)
        self.assertEqual(single_hsp_region.q_len, 466)

    def test_inferred_regions_merging_reverse_compliment(self):
        single_hsp_region = self.blast.inferred_regions["Cbre_Contig6"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region.s_start, 1479061)
        self.assertEqual(single_hsp_region.s_end, 1492025)
        self.assertEqual(single_hsp_region.strand, "-")
        self.assertEqual(single_hsp_region.chunk_cov, 82)
        self.assertEqual(single_hsp_region.query_cov, 82)
        self.assertEqual(single_hsp_region.q_len, 466)

    def test_inferred_regions_merging_plus_strand(self):
        single_hsp_region = self.blast.inferred_regions["Cbre_Contig231"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region.s_start, 10349)
        self.assertEqual(single_hsp_region.s_end, 24214)
        self.assertEqual(single_hsp_region.strand, "+")
        self.assertEqual(single_hsp_region.chunk_cov, 82)
        self.assertEqual(single_hsp_region.query_cov, 82)
        self.assertEqual(single_hsp_region.q_len, 466)

    def test_unmerged_blast_output(self):
        single_hsp_region = self.blast.blast_out["Cbre_Contig419"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region['s_start'], 65160)
        self.assertEqual(single_hsp_region['s_end'], 65339)

    def test_unmerged_blast_reverse_strand(self):
        blast = read_blast_output(self.blast_file, "fake_db")
        single_hsp_region = blast.blast_out["Cbre_Contig6"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertLess(single_hsp_region['s_end'], single_hsp_region['s_start'])

    def test_compare_merged_unmerged_number_regions(self):
        blast = read_blast_output(self.blast_file, "fake_db")
        number_unmerged_regions = len(blast.blast_out["Cbre_Contig231"]["OrthologousGroups_I3.5.OGoverlapp.txt"])
        self.assertEqual(number_unmerged_regions, 4)
        number_merged_regions = len(self.blast.inferred_regions["Cbre_Contig231"]["OrthologousGroups_I3.5.OGoverlapp.txt"])
        self.assertEqual(number_merged_regions, 1)

    def test_compute_coverage(self):
        test_obj = HspListObject("fake", "fake")
        chunk, query = test_obj.compute_coverage([35, 40, 80], [50, 60, 100], 100)
        self.assertEqual(chunk, 69)
        self.assertEqual(query, 45)
        chunk_single, query_single = test_obj.compute_coverage([35], [50], 100)
        self.assertEqual(chunk_single, 100)
        self.assertEqual(query_single, 15)


class TestExonerateObject(unittest.TestCase):

    # with hashes in sequence
    exonerate_file = os.path.join(test_data, "elegans_eef_true.exonerate")
    exonerate_obj = ExonerateObject(exonerate_file)
    target_hash = "I;9859005-9873867"
    target_range_hash = ('4791', '9857')
    target_prot_hash = exonerate_obj.target_prot[target_hash][target_range_hash]
    target_dna_hash = exonerate_obj.target_dna[target_hash][target_range_hash]

    # no hashes in sequence
    target2 = 'I;9159307-9171689'
    target_range2 = ("819", "10142")
    target_prot2 = exonerate_obj.target_prot[target2][target_range2]
    target_dna2 = exonerate_obj.target_dna[target2][target_range2]

    def test_target_target_range(self):
        self.assertIn(self.target_hash, self.exonerate_obj.query_prot)
        self.assertIn(self.target_range_hash, list(self.exonerate_obj.query_prot[self.target_hash].keys()))

    def test_sequences_dividability(self):
        self.assertEqual(len(self.target_dna_hash) / 3, len(self.target_prot_hash))
        self.assertEqual(len(self.target_dna2) / 3, len(self.target_prot2))

    def test_protein_length(self):
        self.assertEqual(827, len(self.target_prot_hash))
        self.assertEqual(961, len(self.target_prot2))

    def test_dna_length(self):
        self.assertEqual(2481, len(self.target_dna_hash))
        self.assertEqual(2883, len(self.target_dna2))

    def test_grap_values(self):
        self.assertEqual({2481, 2883}, set([len(x) for x in grap_values(self.exonerate_obj.target_dna)]))

    def test_gff_before_rewriting(self):
        first_element = 'I;9859005-9873867\texonerate:protein2genome:bestfit\tgene\t4792\t9857\t485\t+\t.\tgene_id 0 ; sequence eef_3.5 ; gene_orientation +'
        self.assertEqual(self.exonerate_obj.gff[self.target_hash][self.target_range_hash][0], first_element)

    def test_gff_cds_phase(self):
        self.assertEqual(self.exonerate_obj.gff[self.target_hash][self.target_range_hash][11].split("\t")[-2], '1')

    def test_gff_mRNA(self):
        self.assertEqual(self.exonerate_obj.gff[self.target_hash][self.target_range_hash][1].split("\t")[2], 'mRNA')

    def test_gff_no_splice_no_similarity(self):
        self.assertNotIn("similarity", self.exonerate_obj.gff[self.target_hash][self.target_range_hash][-1].split("\t")[2])
        self.assertIn("exon", self.exonerate_obj.gff[self.target_hash][self.target_range_hash][-1].split("\t")[2])
        self.assertNotIn("splice", self.exonerate_obj.gff[self.target_hash][self.target_range_hash][3].split("\t")[2])

    def test_aa3_to_1_coding_dict(self):
        test_string = 'Val***XaaUnkSer'
        out_string = 'VXXXS'
        failed_string = "ValUnk*"
        self.assertEqual(out_string, aacode_3to1(test_string))
        self.assertRaises(Exception, aacode_3to1, failed_string)

    def test_clear_hashed_bases(self):
        testprot = "        sHisValGluLeuLeu##LysSer------------------SerGlyIleSerLeuLeuCy"
        testdna = " 5030 : ACACGTGGAATTACTATGAAATCC------------------AGTGGAATCAGTTTGCTCTG : 5073"
        out_dna = " 5030 : ACACGTGGAATTACTAAAATCC------------------AGTGGAATCAGTTTGCTCTG : 5073"
        self.assertEqual(out_dna, clear_hashed_bases(testprot, testdna))

    def test_remove_non_letter_signs(self):
        test_string = "AdfA%y*8l#+-w"
        out_string = "AdfAy*lw"
        self.assertEqual(remove_non_letter_signs(test_string), out_string)


class TestDataProviderObject(unittest.TestCase):

    data_base = DataProviderObject(test_data)

    def test_all_attributes(self):
        self.assertEqual(self.data_base.cluster_scope, 5)
        self.assertEqual(self.data_base.group_names[0], "next_best_blast_eef")
        self.assertEqual(self.data_base.group_to_single_copy_ortholog["next_best_blast_eef"], "True")
        self.assertEqual(self.data_base.group_to_group_size["next_best_blast_eef"], 5)
        self.assertEqual(self.data_base.group_by_cluster_to_consensus_length["next_best_blast_eef"]["eef_3.5"], 840)
        self.assertTrue(self.data_base.group_by_cluster_to_hmm["next_best_blast_eef"]["eef_3.5"].split(".")[-1] == "hmmGenePS")
        self.assertEqual(self.data_base.group_by_cluster_to_score_cutoff["next_best_blast_eef"]["eef_3.5"], 1.124957886666315)
        self.assertEqual(self.data_base.group_by_cluster_to_len_confidence["next_best_blast_eef"]["eef_3.5"], ['282.6490986677823', '953.9017487898448'])


class TestOverseer(unittest.TestCase):
    elegans_db = os.path.join(test_data, "c_elegans.PRJNA13758.WS254.genomic.fa")
    overseer = Overseer("c_elegans", elegans_db, "prediction_location", DataProviderObject(test_data))
    blast_file = os.path.join(test_data, "elegans_blast_eef_nextBest.txt")
    blast = read_blast_output(blast_file, elegans_db)
    blast.infer_regions()
    print(blast.inferred_regions.values())
    overseer.group_to_blast_obj["next_best_blast_eef"] = blast

    def test_attributes_after_blast(self):
        with tempdir() as tmp:
            self.overseer.predict_on_all_regions(tmp)








if __name__ == '__main__':
    unittest.main()
