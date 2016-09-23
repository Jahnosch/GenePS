#!/usr/bin/env python3
import subprocess
import os
import sys
import tempfile
import shutil
from statistics import mean, stdev
import exonerate_parser as ex
import find_regions as fr

def run_cmd(**kwargs):
    cmd = kwargs['command']
    if type(cmd) != list:
        cmd = kwargs['command'].split() # sanitation
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True,
                         bufsize=-1) # buffersize of system
    wait = kwargs.get('wait', False)
    if wait:
        p.wait()
        if p.returncode == 0:
            pass
    else:
        return iter(p.stdout.readline, '')

def which(program):
    if shutil.which(program):
        return program
    else:
        print(program + " not Found\n")
        sys.exit(1)

def check_programs():
    which("mafft")
    which("trimal")
    which("exonerate")
    which("makeblastdb")
    which("hmmsearch")
    which("hmmemit")
    which("hmmbuild")
    which("blastdbcmd")

def hash_fasta(fasta_file):
    fasta = {}
    active_sequence_name = ""
    with open(fasta_file) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line.split(" ")[0]
                if active_sequence_name not in fasta:
                    fasta[active_sequence_name] = []
                continue
            sequence = line
            fasta[active_sequence_name].append(sequence)
    file_one.close()
    return fasta


def get_out_dir(outfile):
    #add later more ...
    if os.path.split(outfile)[1] == "":
        outfile += "still_no_name"
    try:
        if not os.path.exists(outfile):
            os.mkdir(outfile)
            os.mkdir(outfile + "/temp")
    except:
        raise OSError

    return outfile, outfile + "/temp"


def msa_operations(infile, program ="mafft"):
    if program == "trimal":
        command = infile
    else:
        command = "mafft " + infile
    read_flag = 0
    seq = []
    msa_list = []
    for line in run_cmd(command=command, wait=False):
        if read_flag == 0:
            if line.startswith(">"):
                read_flag = 1
        if read_flag == 1:
            if line.startswith(">"):
                if len(msa_list) > 0:
                    msa_list.append("".join(seq))
                msa_list.append(line.rstrip("\n"))
                seq = []
            else:
                seq.append(line.rstrip("\n"))
    msa_list.append("".join(seq))
    # try and error suit: last entry seq, first entry >, length/2 == even
    return msa_list

class MsaObject:

    def __init__(self, msa_list, cluster_name, output_dir):
        self.name = cluster_name + ".msa"
        self.path = os.path.join(output_dir, self.name)
        self.size = [len(msa_list)/2]
        self.msa_string = "\n".join(msa_list)
        self.lengths = [len(msa_list[1])]
        self.trim_cmd = ("trimal -in " + self.path + " -resoverlap 0.6 -seqoverlap 80",
                          "trimal -in " + self.path + " -automated1")

    def update_msa(self, msa_list):
        self.size.append(len(msa_list)/2)
        self.lengths.append(len(msa_list[1]))
        self.msa_string = "\n".join(msa_list)

    # flag to ensure that not to many get removed and minimum amount of taxa
    # (omit first trimming/remove cluster)
    def trim_msa (self):
        with open(msa.path, "w") as mf:
            mf.write(msa.msa_string + "\n")
        for trimmer in self.trim_cmd:
            self.update_msa(msa_operations(trimmer, program ="trimal"))
            with open(self.path, "w") as trmf:
                trmf.write(self.msa_string + "\n")


class HmmFileObject:
    # always a file
    def __init__(self,cluster_name, output_dir, msa):
        self.name = cluster_name + ".hmm"
        self.path = os.path.join(output_dir, self.name)
        self.source = msa
        if type(msa) == MsaObject:
            self.hmm_build(msa.path)
        else:
            try:
                self.hmm_build(msa)
            except ValueError:
                print("can not build HMM from" + msa)

    def hmm_build(self, msa_file):
        command = "hmmbuild " + self.path + " " + msa_file
        run_cmd(command=command, wait=True)

    def hmm_search(self, query_file):
        command = "hmmsearch " + "--noali " + self.path + " " + query_file
        read_count = 0
        for line in run_cmd(command=command, wait=False):
            if "E-value" in line or read_count ==1:
                read_count += 1
            elif read_count == 2:
                line = line.strip("\n").split()
                return float(line[1])

    def get_consensus(self):
        command = "hmmemit " + " -c " + self.path
        cons_list = []
        header = ">" + hmm_trimmed_msa.name.split(".hmm")[0]
        for line in run_cmd(command=command, wait=False):
            if not line.startswith(">"):
                cons_list.append(line.strip("\n"))
        return header + "\n" + "".join(cons_list)


def handle_temps(query, remainder, fa_dict, hmm_name, hmm_dir):
    temp_query = tempfile.NamedTemporaryFile(mode="r")
    tq = open(temp_query.name, "w")
    temp_msa = tempfile.NamedTemporaryFile(mode="r")
    tm = open(temp_msa.name, "w")
    temp_remainder = tempfile.NamedTemporaryFile(mode="r")
    tr = open(temp_remainder.name, "w")
    # write to temp files
    try:
        if query == "all":
            tq.write(query + "\n")
        else:
            tq.write(query + "\n" + "".join(fa_dict[query]))
        seq_list =[]
        for key in remainder:
            seq_list.append(key)
            seq_list.append("".join(fa_dict[key]))
        tr.write("\n".join(seq_list))
    finally:
        tq.close()
        tr.close()
    try:
        tm.write("\n".join(msa_operations(temp_remainder.name)))
    finally:
        tm.close()
    profile = HmmFileObject(hmm_name, hmm_dir, temp_msa.name)
    score = profile.hmm_search(temp_query.name)
    temp_query.close()
    temp_remainder.close()
    temp_msa.close()
    return profile, score

def scoring(fasta, msa, hmm_name, hmm_dir):
    score_list = []
    fa_dict = hash_fasta(fasta)
    msa_dict = hash_fasta(msa.path)
    msa_key_list = list(msa_dict.keys())
    for idx in range(0, len(msa_key_list)):
        if idx == 0:
            query = msa_key_list[idx]
            remainder = msa_key_list[idx+1:]
        elif idx == len(msa_key_list)-1:
            query = msa_key_list[idx]
            remainder = msa_key_list[:idx]
        else:
            query = msa_key_list[idx]
            remainder = msa_key_list[:idx]+msa_key_list[idx+1:]
        hmm, score = handle_temps(query, remainder,fa_dict, hmm_name, hmm_dir)
        score_list.append(score)
    # full msa for prediction check, no query
    hmm, score = handle_temps("all", msa_key_list,fa_dict, hmm_name, hmm_dir)
    return score_list, hmm


def benchmark(list_scores, phmm_score):
    average = mean(list_scores)
    conf_inter = (average - stdev(list_scores),
                  average + stdev(list_scores))
    adj_average = round((average + phmm_score) / 2, 2)
    print("exonerate score = {}".format(str(phmm_score)))
    if conf_inter[0] < adj_average < conf_inter[1]:
        print("adjusted average: ", adj_average, " WITHIN confidence interval ", conf_inter)
        return True
    else:
        print("adjusted average: ", adj_average, " OUTSIDE confidence interval ", conf_inter)
        return False


def predict_protein(cluster, blast_ob, temp_dir, hmm, score_list, flank=2000):
    def get_prot_seq():
        for y in ex_ob.target_prot.values():
            for prot in y.values():
                return prot
    consensus = temp_dir+"/consensus_hmm.fasta"
    dir_string = os.path.join(temp_dir,cluster)
    for idx in range(0, len(blast_ob.hsp_idx[cluster])):
        dict = blast_ob.get_column(cluster, idx)
        if dict["chunk_coverage"] > 70:
            region = fr.parse_blastdb(
                blast_ob.db_path, dict["hsp_contig"],
                dict["hsp_region"][0], dict["hsp_region"][1], flank)
            with open(dir_string+".blast", "w") as blast_file:
                blast_file.write("\n".join(region))
            ex_ob = ex.run_exonerate(cluster+".exn", temp_dir, dir_string +".blast", consensus)
            prediction = get_prot_seq()
            with open(dir_string+".exn", "w") as ex_file:
                ex_file.write(">" + cluster + "_" + region[0] + "\n")
                ex_file.write(prediction)
            phmm_score = hmm.hmm_search(dir_string+".exn")
            print(region[0] + ";" + cluster)
            benchmark(score_list, phmm_score)
            # append exon.object to list or dict when passed and return
            # every seq is an exonerate object


def progress(iteration, steps, max_value):
    if int(iteration) == int(max_value):
        sys.stdout.write('\r')
        print ("[PROGRESS]\t: %d%%" % (100))
    elif int(iteration) % int(steps+1) == 0:
        sys.stdout.write('\r')
        print ("[PROGRESS]\t: %d%%" % (float(int(iteration)/int(max_value))*100))
        sys.stdout.flush()
    else:
        pass

if __name__ == "__main__":
    check_programs()
    infile = sys.argv[1]
    cluster_name = os.path.splitext(infile)[0].split("/")[-1]
    folder_name = "candidate_regions.fa"

    genome = sys.argv[2]

    outfile = sys.argv[3] # will create a directory
    out_dir, temp_dir = get_out_dir(outfile)

    consensus = temp_dir+"/consensus_hmm.fasta"
    consensus_file = open(consensus,"w")
    cutoff_dict = {}

    print("\nanalyzing:\t {}\n".format(cluster_name))
    # first MSA for pHMM consensus
    msa = MsaObject(msa_operations(infile), cluster_name, temp_dir)
    msa.trim_msa()
    progress(1,0,4)
    # pHMM consensus
    hmm_trimmed_msa = HmmFileObject(cluster_name, temp_dir, msa)
    consensus_file.write(hmm_trimmed_msa.get_consensus())
    consensus_file.close()
    progress(2,0,4)
    #compute score
    score_list, hmm2evaluate = scoring(infile, msa, "score_" + str(cluster_name), temp_dir)
    cutoff_dict[cluster_name] = score_list
    progress(3,0,4)
    # exonerate_parser
    db_path = fr.make_blast_db(genome, temp_dir)
    blast = fr.run_tblastn(db_path, consensus)
    progress(4,0,4)
    print()
    predict_protein(cluster_name, blast, temp_dir, hmm2evaluate, score_list)
    print()

