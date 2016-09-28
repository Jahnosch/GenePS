#!/usr/bin/env python3

'''
Usage: run_GenePS.py                          -m <DIR> -g <FILE>

    Options:
        -h, --help                            show this screen.

        General
        -m, --GenePS_result_dir <DIR>         folder with consensus/score results from make_GenePS.py
        -g, --genome <FILE>                   Target genome

'''

import os
from docopt import docopt
import tempfile as tmp
from collections import defaultdict
from statistics import mean, stdev

from run_command import tempdir, check_programs
from exonerate_parser import run_exonerate
from find_regions import run_tblastn, make_blast_db, parse_blastdb
from make_GenePS import get_phmm_score, write_to_tempfile, get_outdir


class ExonerateError(Exception):
    pass


class ScoreError(Exception):
    pass


def judge_score(list_scores, phmm_score):
    average = mean(list_scores)
    conf_inter = (average - stdev(list_scores),
                  average + stdev(list_scores))
    adj_average = round((average + phmm_score) / 2, 2)
    if conf_inter[0] < adj_average < conf_inter[1]:
        return {"aver": average, "conf_int": conf_inter, "adj_aver": adj_average}
    else:
        raise ScoreError


def coverage_filter(area):
    if area.chunk_cov == 100 and area.query_cov < 40:
        return False
    elif area.chunk_cov > 41:
        return True


def score_prediction(ex_obj, hmm):
    ex_seq = ex_obj.grap_values(ex_obj.target_prot)
    q_name = ex_obj.grap_values(ex_obj.header)["query"]
    with tmp.NamedTemporaryFile() as ex_file:
        write_to_tempfile(ex_file.name, ">{}\n{}".format(q_name, ex_seq))
        final_score = get_phmm_score(hmm, ex_file.name)
    return final_score


def make_prediction(q_name, cons_file, dir_path, area):
    region_str = parse_blastdb(db_path, area.contig, area.s_start, area.s_end)
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
    gene_ps_results = args['--GenePS_result_dir']
    genome = args['--genome']

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
                    print("\n# Analyzing {} - containing {} files\n"
                          .format(group_result.group_name, group_result.group_size))
                    # all consensus of a folder
                    header_cons = group_result.consensus.keys()
                    group_cons = os.path.join(out_dir, group_result.group_name)
                    group_cons = group_result.consensus_to_fa(header_cons, group_cons)
                    # run t-blastn
                    blast_obj = run_tblastn(db_path, group_cons)
                    blast_obj.infer_regions()
                    for query, contig_regions in blast_obj.inferred_regions.items():
                        single_cons = os.path.join(tmp_dir, query)
                        single_cons = group_result.consensus_to_fa(query, single_cons)
                        hmm_file = group_result.phmm[query]
                        for contig_regions in contig_regions.values():
                            for region in contig_regions:
                                if coverage_filter(region) is True:
                                    try:
                                        exo_obj = make_prediction(query, single_cons, tmp_dir, region)
                                        score = score_prediction(exo_obj, hmm_file)
                                        score_valid = judge_score(group_result.score_list[query], score)
                                    except ExonerateError:
                                        print("[!] NO EXONERATE PREDICTION {}, {}, {}\t"
                                              .format(query, region.contig, region.s_start))
                                        continue
                                    except ScoreError:
                                        print("[!] {}, {}, {}\t got filtered"
                                              .format(query, region.contig, region.s_start))
                                        continue
                                    print("[+] {}, {}, {}\t FOUND".format(query, region.contig, region.s_start))
                                    group_result.exonerate_out[query].append(exo_obj)
                    # print(group_result.group_name)
                    # print((len(group_result.exonerate_out) / group_result.group_size) * 100)
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
