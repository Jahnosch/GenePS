#!/usr/bin/env python3
import re
import os
from collections import defaultdict
from run_command import run_cmd

# if is going to be a parser, deal with softmasking

aa3_to_1_coding_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', '***': 'X',
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M','Xaa': 'X', 'Unk': 'X'}


def aacode_3to1(seq):
    """Turn a three letter protein into a one letter protein.
The 3 letter code can be upper, lower, or any mix of cases
The seq input length should be a factor of 3 or else results
in an error
"""
    if len(seq) % 3 == 0:
        single_seq = []
        for i in range(0, len(seq), 3):
            single_seq.append(aa3_to_1_coding_dict.get(seq[i:i+3]))
        return "".join(single_seq)
    else:
        raise False


def kill_char(string, n):
    if n > len(string):
        raise IndexError
    begin = string[:n]
    end = string[n+1:]
    return begin + end


def find_hashes(protein_string):
    for idx, sign in enumerate(protein_string):
        if sign == "#":
            yield idx


def clear_hashed_bases(protein_string, dna_string):
    """removes nucleotides which are opposite to #'s in the protein sequence (non-coding-nucleotides)"""
    if "#" in protein_string:
        for idx in sorted(find_hashes(protein_string), reverse=True):
            dna_string = kill_char(dna_string, idx)
    return dna_string


def remove_non_letter_signs(string):
    """removes all non-alphabetic letters from string"""
    regex = re.compile('[^a-zA-Z*]')
    return regex.sub("", string)


def remove_lower(text_string):
    return re.sub('[a-z]', '', text_string)


def extract_geneps_data(exonerate_obj):
    """returns a fasta string containing all predicted protein sequences, with query name as header. Additionally
    a hash will be returned with query as key and sequence length as value"""
    len_dict = {}
    fasta_list = []
    for query in exonerate_obj.query_prot:
        for trange in exonerate_obj.query_prot[query]:
            fasta_list.append(">" + query)
            fasta_list.append(exonerate_obj.target_prot[query][trange])
            len_dict[">" + query] = len(exonerate_obj.target_prot[query][trange])
    return "\n".join(fasta_list), len_dict


def correct_gene_position(strand, blast_location, exonerate_location):
    if strand == "-":
        gene_start = max(blast_location) - min(exonerate_location)
        gene_end = max(blast_location) - max(exonerate_location)
    elif strand == "+":
        gene_start = min(blast_location) + min(exonerate_location)
        gene_end = min(blast_location) + max(exonerate_location)
    else:
        print("[!]\t ERROR: {} can not be identified as + or - strand".format(strand))
        return None
    return str(gene_start), str(gene_end)


def write_exonerate_gff(gff_list, off_set_tuple, strand):
    last_phase = 0
    new_gff_list = []
    mrna_copy = None
    for line in gff_list:
        line = line.split("\t")
        if strand != line[6]:
            print("blast strand != exonerate strand")
        line[3], line[4] = correct_gene_position(line[6], off_set_tuple, (int(line[3]), int(line[4])))
        if "cds" in line[2]:
            current_phase = ((abs(int(line[4]) - int(line[3])) + 1) + last_phase) % 3
            line[7] = str(last_phase)
            last_phase = current_phase
        elif "gene" in line[2]:
            mrna_copy = line[:]
            mrna_copy[2] = "mRNA"
        elif "similarity" in line[2] or "splice" in line[2]:
            continue
        new_gff_list.append("\t".join(line))
    new_gff_list.insert(1, "\t".join(mrna_copy))
    return new_gff_list


class ExonerateObject:

    def __init__(self, exonerate_file):
        self.path = os.path.realpath(exonerate_file)
        self.header = defaultdict(dict)
        self.query_prot = defaultdict(lambda: defaultdict(list)) #three leafs dict
        self.target_dna = defaultdict(lambda: defaultdict(list))
        self.target_prot = defaultdict(lambda: defaultdict(list))
        self.gff = defaultdict(lambda: defaultdict(list))
        self.exonerate_processor(exonerate_file)

    def exonerate_processor(self, exonerate_file):
        ''' protein against 1 or more genome regions (targets),
        file with aln and/or gff output '''
        del_intron = lambda seq: seq.replace("TargetIntron","")                  # no lambda
        prep_line = lambda line: line.strip("\n").split(": ")[1].strip(" ")
        next_line = lambda: prep_line(next(ex))
        next_block = lambda size: [next_line() for x in range(0, size)]

        read_flag = 0
        with open(exonerate_file, "r") as ex:
            for line in ex:
                if line.strip():
                    if read_flag == 0:
                        if line.startswith("C4 Alignment"):
                            read_flag = 1
                    elif read_flag == 1:
                        query, target, model, score = next_block(4)
                        qrange, trange = [tuple(x.split(" -> ")) for x in next_block(2)]
                        self.header[query][trange] = {"model": model, "score": score,
                                                       "qrange": qrange, "query": query}
                        read_flag = 2
                    elif read_flag == 2:
                        if not line.startswith("#"):
                            self.query_prot[query][trange].append(del_intron(remove_non_letter_signs(line)))
                            next(ex)
                            target_prot = next(ex)
                            target_dna = clear_hashed_bases(target_prot, next(ex))
                            self.target_prot[query][trange].append(remove_non_letter_signs(target_prot))
                            self.target_dna[query][trange].append(remove_lower(remove_non_letter_signs(target_dna)))
                        else:
                            self.query_prot[query][trange] = "".join(self.query_prot[query][trange])
                            self.target_dna[query][trange] = "".join(self.target_dna[query][trange])
                            self.target_prot[query][trange] = aacode_3to1("".join(self.target_prot[query][trange]))
                            if len(self.target_dna[query][trange]) / 3 != len(self.target_prot[query][trange]):
                                print(self.header[query][trange]["query"], "exonerate length error")
                            if "GFF" in line:
                                read_flag = 3
                            else:
                                read_flag = 0
                    elif read_flag == 3:
                        if not line.startswith("#"):
                            line = line.strip("\n")
                            self.gff[query][trange].append(line)
                        elif "END OF GFF DUMP" in line:
                            read_flag = 0
                    else:
                        pass
            ex.seek(0)


def make_exonerate_command(model, query_file, region_file):
    cmd = "exonerate {} --softmaskquery no -Q protein -T dna -n 1 " \
      "--softmasktarget no --showvulgar no --minintron 20 --maxintron 50000 " \
      "--showalignment yes --showtargetgff yes -q {} -t {}".format(model, query_file, region_file)
    return cmd


#### getting query from fasta hash ans writing new output file within the tmp file function
def exonerate_best_raw_score(model, query_file, region_file):
    cmd = make_exonerate_command(model, query_file, region_file)
    read_flag = 0
    hit_dict, query = {}, None
    for line in run_cmd(command=cmd, wait=False):
        line = line.strip("\n")
        if line.startswith("C4 Alignment"):
            read_flag = 1
        elif read_flag == 1 and "Query" in line:
            query = line.split(": ")[1]
        elif read_flag == 1 and "Raw score" in line:
            hit_dict[query] = line.split(": ")[1]
            read_flag = 0
    if hit_dict:
        return max(hit_dict, key=lambda key: hit_dict[key])
    else:
        return None


def get_exonerate_object(output_path, command):
    line_count = 0
    with open(output_path, "w") as ex:
        for line in run_cmd(command=command, wait=False):
            ex.write(line)
            line_count += 1
    if line_count < 10:
        return None
    else:
        return ExonerateObject(output_path)


# sometimes the exhaustive command fails due to an exonerate buck
# if so, try non-exhaustive
def run_exonerate(mode_string, name, directory, region, query):
    cmd = make_exonerate_command(mode_string, query, region)     # "-m p2g:b -E yes" is the best exhaustive command
    out_file = os.path.join(directory, name)
    exonerate_obj = get_exonerate_object(out_file, cmd)
    if get_exonerate_object(out_file, cmd) is None and "-E yes" in mode_string:
        cmd = make_exonerate_command("-m p2g -E no", query, region)
        exonerate_obj = get_exonerate_object(out_file, cmd)
    return exonerate_obj


if __name__ == "__main__":
    test = run_exonerate("-m p2g -E no", "test_exonerate.out", "/home/jgravemeyer/Dropbox/MSc_project/data",
                  "/home/jgravemeyer/Desktop/experiment_data/blast_region.fasta",
                  "/home/jgravemeyer/Desktop/experiment_data/eef_consensus.fa")

    test_raw_score = exonerate_best_raw_score("-m p2g -E no", "/home/jgravemeyer/Desktop/experiment_data/eef_consensus.fa", "/home/jgravemeyer/Desktop/experiment_data/blast_region.fasta")
    print(test_raw_score)
#    print("\n".join(test.gff["eef_3.5"][('5000', '7747')]))
#    print(write_exonerate_gff(test.gff["eef_3.5"][('5000', '7747')], (9158949, 9171695), "+"))
#    print(write_exonerate_gff(test.gff["eef_3.5"][('5000', '7747')], (9158949, 9171695), "+")[0].split("\t")[3:5])

