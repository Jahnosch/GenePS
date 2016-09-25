#!/usr/bin/env python3
from statistics import mean, stdev
from docopt import docopt, DocoptExit
from exonerate_parser import * # try not to use * ... import only what you want to import
from find_regions import * # try not to use * ... import only what you want to import
from make_GenePS import get_phmm_score, write_to_tempfile

"""
usage: run_GenePS.py    -g <FILE> -m <FILE> [-o <PREFIX>]
                        [-h|--help]

    Options:
        -h --help                           show this

        General
            -g, --genome <FILE>                 Target genome
            -m, --make_file <FILE>              make file
            -o, --out_prefix <STR>              Prefix
"""


class ExonerateError(Exception):
    pass


class ScoreError(Exception):
    pass


def get_out_folder(makeGenePS_dir):
    if os.path.isfile(makeGenePS_dir):
        print(makeGenePS_dir, " is NOT a directory!")
        print("Please specify an output directory")
        sys.exit()
    else:
        if not os.path.exists(os.path.join(makeGenePS_dir, "Predictions")):
            os.mkdir(os.path.join(makeGenePS_dir, "Predictions"))
        return os.path.join(makeGenePS_dir, "Predictions")


def check_score(list_scores, phmm_score):
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

        # initialize automatically
        self.process_make_output() # do this outside of init

    def process_make_output(self):
        def mod_next(): return next(mg).strip().split(":")
        with open(self.path) as mg: # mg?
            self.group_name = mg.readline().split(":")[1].strip()
            self.group_size = int(mg.readline().split(":")[1].strip())
            for line in mg:
                if line.startswith(">name:"):
                    cluster = line.split(":")[1].strip()
                    self.phmm[cluster] = mod_next()[1].strip()
                    scores = mod_next()[1].strip()
                    cons_seq = mod_next()[0]
                    self.seq_length[cluster] = (len(cons_seq))
                    self.consensus[cluster] = cons_seq
                    self.score_list[cluster] = [float(x) for x in scores.split(",")]
                else:
                    pass

    def consensus_to_fa(self, cons_header, file_path):
        directory = file_path + ".consensus"
        c_file = open(directory, "w")
        if type(cons_header) == str:
            cons_header = [cons_header]
        for key in cons_header:
            fasta_str = ">{}\n{}\n".format(key, self.consensus[key])
            c_file.write(fasta_str)
        c_file.close()
        return directory


if __name__ == "__main__":
        __version__ = 0.1
    args = docopt(__doc__)
    try:
        make_file = args['--make_file'] # better variable name?
        genome = args['--genome'] # target_genome ?
    except docopt.DocoptExit:
        print __doc__.strip()

    get_out_folder(make_file) # return gets not saved ?

    # make database
    db_path = make_blast_db(genome, make_file)
    print("\ngenerating blast db")

    with tempdir() as tmp_dir:
        for subdir, dirs, files in os.walk(make_file):
            for file in files: # you are overwriting "file" built-in, use tmp_file
                # file = group of many cluster
                if file.split(".")[-1] == "makeGenePS":
                    group_file = os.path.join(subdir, file)
                    group_result = ResultsObject(group_file)
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
                                        score_valid = check_score(group_result.score_list[query], score)
                                    except ExonerateError:
                                        print("NO EXONERATE PREDICTION")
                                        continue
                                    except ScoreError:
                                        print("{}, {}, {} was filtered out"
                                              .format(query, region.contig, region.s_start))
                                        continue
                                    print("{}, {}, {} TRUE".format(query, region.contig, region.s_start))
                                    group_result.exonerate_out[query].append(exo_obj)
                    print(group_result.group_name)
                    print((len(group_result.exonerate_out) / group_result.group_size) * 100)

# to add later
'''
#q_len = group_result.seq_length[query]
#prediction_len = len(exo_obj.grap_values(exo_obj.target_prot))
#fragmentation = (prediction_len / q_len) * 100
#print(q_len, prediction_len)
#print(fragmentation)

#print to result file
#print missing genes to doc file'''
