#!/usr/bin/env python3

'''
Usage: run_GenePS.py                          -m <DIR> -g <FILE> [-c <INT>] [-f <INT>]

    Options:
        -h, --help                            show this screen.

        General
        -m, --GenePS_result_dir <DIR>         folder with consensus/score results from make_GenePS.py
        -g, --genome <FILE>                   Target genome
        Optional
        -c, --coverage_filer <INT>            Minimal aligned length of a Blast query to the target genome (used to filter Blast hits)
        -f, --HMM_filter <INT>                Factor to multiply standard deviation of the HMM score distribution with (validate predictions)

'''

import os
from docopt import docopt
import tempfile as tmp
from collections import defaultdict
from statistics import mean, stdev

from run_command import tempdir, check_programs
from exonerate_parser import run_exonerate, grap_values
from find_regions import run_tblastn, make_blast_db, parse_blastdb
from check_accuracy import true_coordinates_hash, RegionObj
from make_GenePS import get_phmm_score, write_to_tempfile, get_outdir


class ExonerateError(Exception):
    pass


def judge_score(phmm_score):
    adj_average = round((score_mean + phmm_score) / 2, 2)
    print("\nPrediction Score: {}; Adjusted score {}".format(str(phmm_score), str(adj_average)))
    if confidence_inter[0] < adj_average:
        return True, {"prediction": phmm_score, "confidence": confidence_inter, "adj_mean": adj_average}
    else:
        return False, {"prediction": phmm_score, "confidence": confidence_inter, "adj_mean": adj_average}


def coverage_filter(area):
    if area.chunk_cov == 100 and area.query_cov < 40:
        return False
    elif area.chunk_cov > coverage_min:
        return True


def number_blast_regions(query_blast_out):
    sum_region = 0
    for contig_summary, regions_of_contig in query_blast_out.items():
        sum_region += len(regions_of_contig)
    return sum_region


def score_prediction(ex_obj, hmm):
    ex_seq = grap_values(ex_obj.target_prot)[0]
    q_name = grap_values(ex_obj.header)[0]["query"]
    with tmp.NamedTemporaryFile() as ex_file:
        write_to_tempfile(ex_file.name, ">{}\n{}".format(q_name, ex_seq))
        final_score = get_phmm_score(hmm, ex_file.name)[0]
    return final_score


def make_prediction(q_name, cons_file, dir_path, area, blast_db):
    region_str = parse_blastdb(blast_db, area.contig, area.s_start, area.s_end)
    with tmp.NamedTemporaryFile() as reg_file:
        write_to_tempfile(reg_file.name, region_str)
        ex_name = q_name + ".exon"
        ex_obj = run_exonerate(ex_name, dir_path, reg_file.name, cons_file)
        if bool(ex_obj.header) is False:
            raise ExonerateError
        return ex_obj


class ResultsObject:
    def __init__(self, make_output):
        self.path = make_output
        self.group_name = ""
        self.group_size = 0
        self.seq_length = {}
        self.consensus = {}
        self.phmm = {}
        self.score_list = {}
        self.exonerate_out = defaultdict(list)

    def read_gene_ps_consensus_file(self):
        """reads the results file.makeGenePS form the first "make_GenePS" script, containing: consensus seq and scores
        for each cluster. There is one results file for each Group-folder the user in "make_GenePS specified as input.
        :return: init the attributes (all clusters per group-folder in one instance)
        """
        def mod_next(): return next(mg).strip().split(":")
        with open(self.path) as mg:
            self.group_name = mg.readline().split(":")[1].strip()
            self.group_size = int(mg.readline().split(":")[1].strip())
            for line in mg:
                if line.startswith(">name:"):
                    cluster = line.split(":")[1].strip()
                    self.phmm[cluster] = mod_next()[1].strip()
                    scores = mod_next()[1].strip().strip("[|]")
                    cons_seq = mod_next()[0]
                    self.seq_length[cluster] = (len(cons_seq))
                    self.consensus[cluster] = cons_seq
                    self.score_list[cluster] = [float(x) for x in scores.split(",")]
                else:
                    pass

    def consensus_to_fa(self, cons_header, file_path):
        """ writes a fasta file for consensus sequences
        :param cons_header: all consensus seq header to take the corresponding seq from self.consensus[header]
        :param file_path: where to write the file
        :return: fasta file
        """
        directory = file_path + ".consensus"
        c_file = open(directory, "w")
        if type(cons_header) == str:
            cons_header = [cons_header]
        for header_key in cons_header:
            fasta_str = ">{}\n{}\n".format(header_key, self.consensus[header_key])
            c_file.write(fasta_str)
        c_file.close()
        return directory


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    gene_ps_results = os.path.abspath(args['--GenePS_result_dir'])
    genome = args['--genome']
    coverage_min = args['--coverage_filer']
    std_multip = args['--HMM_filter']
    if coverage_min is None:
        coverage_min = 30
    else:
        coverage_min = int(coverage_min)
    if std_multip is None:
        std_multip = 2
    else:
        std_multip = int(std_multip)

    check_programs("tblastn", "makeblastdb", "exonerate")

    out_dir = get_outdir(gene_ps_results, add_dir="Predictions")
    print("#" * 32)
    print("# writing run_GenePS results to: {}".format(out_dir))
    print("#" * 32 + "\n")

    # make database
    db_path = make_blast_db(genome, out_dir)
    with tempdir() as tmp_dir:
        for subdir, dirs, files in os.walk(gene_ps_results):
            for file_path_str in files:
                # file = group of many cluster
                if file_path_str.split(".")[-1] == "makeGenePS":
                    group_file = os.path.join(subdir, file_path_str)
                    group_result = ResultsObject(group_file)
                    group_result.read_gene_ps_consensus_file()
                    print("\n# Analyzing {} - containing {} files"
                          .format(group_result.group_name, group_result.group_size))
                    # all consensus of a folder
                    header_cons = group_result.consensus.keys()
                    group_cons = os.path.join(tmp_dir, group_result.group_name)
                    group_cons = group_result.consensus_to_fa(header_cons, group_cons)
                    # run t-blastn
                    blast_obj = run_tblastn(db_path, group_cons)
                    blast_obj.infer_regions()
                    for query, contig_regions in blast_obj.inferred_regions.items():
                        print("\n### {} - {} potential regions identified\n".format(query, number_blast_regions(contig_regions)))
                        single_cons = os.path.join(tmp_dir, query)
                        single_cons = group_result.consensus_to_fa(query, single_cons)
                        score_list = group_result.score_list[query]
                        score_mean, score_std = mean(score_list), stdev(score_list)
                        confidence_inter = (score_mean - (std_multip * score_std), score_mean + (std_multip * score_std))
                        print("### Confidence Interval {} - {} \n".format(round(confidence_inter[0]), round(confidence_inter[1])))

                        # writing to cluster specific output
                        cluster_spec_dir = os.path.join(out_dir, query)
                        cluster_passed_fasta = open(cluster_spec_dir + "_PASSED.fa", "w")
                        cluster_passed_fasta.write("# confidence interval {}-{}\n".format(round(confidence_inter[0]), round(confidence_inter[1])))
                        cluster_filtered_fasta = open(cluster_spec_dir + "_FILTERED.fa", "w")
                        cluster_filtered_fasta.write("# confidence interval {}-{}\n".format(round(confidence_inter[0]), round(confidence_inter[1])))
                        cluster_passed_gff = open(cluster_spec_dir + "_PASSED.gff", "w")
                        cluster_filtered_gff = open(cluster_spec_dir + "_FILTERED.gff", "w")

                        for contig_regions in contig_regions.values():
                            for region in contig_regions:
                                if coverage_filter(region) is True:
                                    try:
                                        exo_obj = make_prediction(query, single_cons, tmp_dir, region, db_path)
                                        score = score_prediction(exo_obj, group_result.phmm[query])
                                    except ExonerateError:
                                        print("[!] {}, {}, {}, {}\t\t NO EXONERATE PREDICTION"
                                              .format(query, region.contig, region.s_start, region.s_end))
                                        continue
                                    # HMM filter
                                    valid, adj_score = judge_score(score)
                                    if valid is False:
                                        print("[!] {}, {}, {}, {}\t\t Filtered by HMM-score"
                                              .format(query, region.contig, region.s_start, region.s_end))
                                        cluster_filtered_fasta.write(">{} {};{}-{} HMM_score: {} Adjusted_Score: {}\n"
                                                                     .format(query, region.contig, region.s_start, region.s_end, score, adj_score["adj_mean"]))
                                        cluster_filtered_fasta.write(grap_values(exo_obj.target_prot)[0] + "\n")
                                        cluster_filtered_gff.write("\n".join(grap_values(exo_obj.gff)[0]) + "\n")
                                        continue
                                    # Passed
                                    print("[+] {}, {}, {}, {}\t\t PASSED".format(query, region.contig, region.s_start, region.s_end))
                                    cluster_passed_fasta.write(">{} {};{}-{} HMM_score: {} Adjusted_Score: {}\n"
                                                               .format(query, region.contig, region.s_start, region.s_end, score, adj_score["adj_mean"]))
                                    cluster_passed_fasta.write(grap_values(exo_obj.target_prot)[0] + "\n")
                                    cluster_passed_gff.write("\n".join(grap_values(exo_obj.gff)[0]) + "\n")
                        cluster_passed_fasta.close()
                        cluster_passed_gff.close()
                        cluster_filtered_fasta.close()
                        cluster_filtered_gff.close()

                    '''
                    test_Regobj = RegionObj(region.contig, region.s_start, region.s_end)
                    for reg in true_coordinates_hash[query]:
                        if reg.contains(test_Regobj):
                            print("\tfits existing model: {}, {}, {}".format(reg.chrom, reg.start, reg.end))
                            print("\toverlap: ", reg.get_overlap_length(test_Regobj))
                    group_result.exonerate_out[query].append(exo_obj)
                    #print(query)
                    #print(grap_values(exo_obj.target_prot)[0])
    # print(group_result.group_name)
    # print((len(group_result.exonerate_out) / group_result.group_size) * 100)'''
            print("\n")


# to add later
'''
#q_len = group_result.seq_length[query]
#prediction_len = len(exo_obj.grap_values(exo_obj.target_prot))
#fragmentation = (prediction_len / q_len) * 100
#print(q_len, prediction_len)
#print(fragmentation)

#print to result file
#print missing genes to doc file'''
