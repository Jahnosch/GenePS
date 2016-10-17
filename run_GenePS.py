#!/usr/bin/env python3

'''
Usage: run_GenePS.py                          -m <DIR> -g <FILE> [-c <INT>] [-s <INT>] [-o <DIR>] [-k <STR>]

    Options:
        -h, --help                            show this screen.

        General
        -m, --GenePS_result_dir <DIR>         folder with consensus/score results from make_GenePS.py
        -g, --genome <FILE|LIST>              Target genome or a ".genomes" file in style of: path/to/genome <TAB> name_prefix
        Optional
        -c, --coverage_filer <INT>            Minimal aligned length of a Blast query to the target genome (used to filter Blast hits)
        -o, --out_dir <DIR>                   Directory for the output files and folders (default: same as input directory)
        -s, --HMM_filter <INT>                Factor to multiply standard deviation of the HMM score distribution with (validate predictions)
        -k, --keep <STR>                      If "Yes" intermediate files will be stored (default is "No")

'''

import os
import sys
from docopt import docopt
import tempfile as tmp
from collections import defaultdict
from statistics import mean, stdev

from run_command import tempdir, check_programs
from exonerate_parser import run_exonerate, grap_values
from find_regions import run_tblastn, make_blast_db, parse_blastdb
#from check_accuracy import true_coordinates_hash, RegionObj
from make_GenePS import get_phmm_score, write_to_tempfile, get_outdir


class ExonerateError(Exception):
    pass


def judge_score(phmm_score, score_mean, confidence_inter):
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


def write_merged_region_to_intermediate(blast_ob):
    results_list = [
        "# Merging Distance: {}, Flanking Distance {}".format(blast_ob.merging_distance, blast_ob.flanking_distance),
        "# Fields: name, contig, subject_start, subject_end, strand, chunk_coverage, total_coverage, query_length"]
    for inferred_region in blast_ob.inferred_regions:
        for contig in blast_ob.inferred_regions[inferred_region]:
            for region in blast_ob.inferred_regions[inferred_region][contig]:
                results_list.append("\t".join([inferred_region, region.contig, str(region.s_start), str(region.s_end), str(region.strand), str(region.chunk_cov), str(region.query_cov), str(region.q_len)]))
    return "\n".join(results_list)


def write_blast_raw_to_intermediate(blast_ob):
    results_list = ["# Fields: query acc., subject acc., evalue, q. start, q. end, s. start, s. end, query length, strand"]
    for query in blast_ob.blast_out:
        for hsp_list in blast_ob.blast_out[query]:
            for hsp in blast_ob.blast_out[query][hsp_list]:
                results_list.append("\t".join([query, hsp["contig"], str(hsp["evalue"]), str(hsp["q_start"]), str(hsp["q_end"]), str(hsp["s_start"]), str(hsp["s_end"]), str(hsp["q_len"]), str(hsp["strand"])]))
    return "\n".join(results_list)


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


def run_geneps_on_genome():
    with tempdir() as tmp_dir:
        for subdir, dirs, files in os.walk(gene_ps_results):
            for file_path_str in files:
                # file = group of many cluster
                if file_path_str.split(".")[-1] == "makeGenePS":
                    group_file = os.path.join(subdir, file_path_str)
                    group_result = ResultsObject(group_file)
                    group_result.read_gene_ps_consensus_file()
                    print("\n# Analyzing {} - containing {} files".format(group_result.group_name, group_result.group_size))

                    # writing to group specific output
                    # if no name because single genome set genome name to ""
                    # add the genome name to 3 things 1) the folder name here, 2) the filename below 3) in the fasta-header
                    if g_prefix is None:
                        genome_group_name = group_result.group_name
                    else:
                        genome_group_name = g_prefix + "_" + group_result.group_name

                    group_specific_dir = os.path.join(prediction_dir, genome_group_name)
                    if not os.path.exists(group_specific_dir):
                        os.mkdir(group_specific_dir)

                    group_passed_fasta_protein = open(os.path.join(group_specific_dir, genome_group_name) + "_PASSED_protein.fa", "w")
                    group_passed_fasta_dna = open(os.path.join(group_specific_dir, genome_group_name) + "_PASSED_dna.fa", "w")
                    group_passed_gff = open(os.path.join(group_specific_dir, genome_group_name) + "_PASSED.gff", "w")

                    group_filtered_fasta_protein = open(os.path.join(group_specific_dir, genome_group_name) + "_FILTERED_protein.fa", "w")
                    group_filtered_fasta_dna = open(os.path.join(group_specific_dir, genome_group_name) + "_FILTERED_dna.fa", "w")
                    group_filtered_gff = open(os.path.join(group_specific_dir, genome_group_name) + "_FILTERED.gff", "w")

                    # all consensus of a folder/group
                    header_cons = group_result.consensus.keys()
                    group_cons = os.path.join(tmp_dir, group_result.group_name)
                    group_cons = group_result.consensus_to_fa(header_cons, group_cons)

                    # run t-blastn
                    blast_obj = run_tblastn(db_path, group_cons)
                    blast_obj.infer_regions()

                    if keep is True:
                        with open(os.path.join(group_specific_dir, genome_group_name) + "_intermediate_merged_blast_regions.txt", "w") as merged_regions:
                            merged_regions.write(write_merged_region_to_intermediate(blast_obj))
                        with open(os.path.join(group_specific_dir, genome_group_name) + "_intermediate_merged_blast_regions.txt", "w") as blast_file:
                            blast_file.write(write_blast_raw_to_intermediate(blast_obj))
                     #   exonerate_file_list = []

                    for query, contig_regions in blast_obj.inferred_regions.items():
                        print("\n### {} - {} potential regions identified\n".format(query, number_blast_regions(contig_regions)))

                        # single consensus per cluster for exonerate
                        cluster_cons = os.path.join(tmp_dir, query)
                        cluster_cons = group_result.consensus_to_fa(query, cluster_cons)

                        # cluster specific score distribution attributes
                        score_list = group_result.score_list[query]
                        score_mean, score_std = mean(score_list), stdev(score_list)
                        confidence_inter = (score_mean - (std_multip * score_std), score_mean + (std_multip * score_std))
                        print("### Confidence Interval {} - {} \n".format(round(confidence_inter[0]), round(confidence_inter[1])))

                        # processing all regions of a cluster
                        for contig_regions in contig_regions.values():
                            for region in contig_regions:
                                if coverage_filter(region) is True:
                                    try:
                                        exo_obj = make_prediction(query, cluster_cons, tmp_dir, region, db_path)
                                        #if keep is True:
                                         #   exonerate_file_list.append(exo_obj.path)
                                        score = score_prediction(exo_obj, group_result.phmm[query])
                                    except ExonerateError:
                                        print("[!] {}, {}, {}, {}\t\t NO EXONERATE PREDICTION".format(query, region.contig, region.s_start, region.s_end))
                                        continue
                                    # HMM filter
                                    valid, adj_score = judge_score(score, score_mean, confidence_inter)
                                    fasta_header = ">{} Cluster:{} Location:{};{}-{} HMM_score:{} Adjusted_Score:{}\n".format(g_prefix, query, region.contig, region.s_start, region.s_end, score, adj_score["adj_mean"])
                                    if valid is False:
                                        print("[!] {}, {}, {}, {}\t\t Filtered by HMM-score".format(query, region.contig, region.s_start, region.s_end))
                                        group_filtered_fasta_protein.write(fasta_header + grap_values(exo_obj.target_prot)[0] + "\n")
                                        group_filtered_fasta_dna.write(fasta_header + grap_values(exo_obj.target_dna)[0] + "\n")
                                        group_filtered_gff.write("\n".join(grap_values(exo_obj.gff)[0]) + "\n")
                                        continue
                                    # Passed
                                    print("[+] {}, {}, {}, {}\t\t PASSED".format(query, region.contig, region.s_start, region.s_end))
                                    group_passed_fasta_protein.write(fasta_header + grap_values(exo_obj.target_prot)[0] + "\n")
                                    group_passed_fasta_dna.write(fasta_header + grap_values(exo_obj.target_dna)[0] + "\n")
                                    group_passed_gff.write("\n".join(grap_values(exo_obj.gff)[0]) + "\n")
               #     with open(os.path.join(group_specific_dir, genome_group_name) + "_intermediate_exonerate.txt", "wb") as exo_file:
               #         for exof in exonerate_file_list:
               #             with open(exof, "rb") as infile:
               #                 exo_file.write(infile.read())
                    group_passed_fasta_protein.close()
                    group_passed_gff.close()
                    group_passed_fasta_dna.close()
                    group_filtered_fasta_protein.close()
                    group_filtered_gff.close()
                    group_filtered_fasta_dna.close()
            print("\n")


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    gene_ps_results = os.path.abspath(args['--GenePS_result_dir'])
    genome = args['--genome']
    coverage_min = args['--coverage_filer']
    std_multip = args['--HMM_filter'] #--HMM_filter
    out_dir = args["--out_dir"]
    keep = args["--keep"]

    if coverage_min is None:
        coverage_min = 30
    else:
        coverage_min = int(coverage_min)
    if std_multip is None:
        std_multip = 2
    else:
        std_multip = int(std_multip)
    if out_dir is None:
        out_dir = "/".join(gene_ps_results.split("/")[:-1])
    if keep is not None:
        if keep == "Yes":
            keep = True
    if genome.split(".")[-1] == "txt":
        genome_dict = {}
        with open(genome) as g_file:
            for genome_line in g_file:
                try:
                    genome_line = genome_line.strip("\n").split("\t")
                    genome_dict[genome_line[0]] = genome_line[1]
                except IndexError:
                    print("ERROR: genome file not in correct format")
                    sys.exit()
    else:
        genome_dict = None

    check_programs("tblastn", "makeblastdb", "exonerate")

    print("#" * 32)
    print("# writing run_GenePS results to: {}".format(out_dir))
    print("#" * 32 + "\n")


    blast_db_dir = get_outdir(out_dir, add_dir="Blast_dbs")
    prediction_dir = get_outdir(out_dir, add_dir="Predictions")
    if genome_dict is None:
        db_path = make_blast_db(genome, blast_db_dir)
        g_prefix = None
        run_geneps_on_genome()
    else:
        for genome_path, g_prefix in genome_dict.items():
            blast_db_dir = get_outdir(blast_db_dir, add_dir=g_prefix)
            prediction_dir = get_outdir(prediction_dir, add_dir=g_prefix)
            db_path = make_blast_db(genome_path, blast_db_dir)
            print("starting with genome: {}".format(g_prefix))
            run_geneps_on_genome()
            # go back in file order
            blast_db_dir = get_outdir(out_dir, add_dir="Blast_dbs")
            prediction_dir = get_outdir(out_dir, add_dir="Predictions")


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

# to add later
'''
#q_len = group_result.seq_length[query]
#prediction_len = len(exo_obj.grap_values(exo_obj.target_prot))
#fragmentation = (prediction_len / q_len) * 100
#print(q_len, prediction_len)
#print(fragmentation)

#print to result file
#print missing genes to doc file'''
