#!/usr/bin/env python3
import os
from operator import itemgetter
from collections import defaultdict, namedtuple
from run_command import run_cmd


def make_blast_db(genome, temp_dir):
    name = os.path.splitext(genome)[0].split("/")[-1]
    out_dir = os.path.join(temp_dir,name)
    if os.path.exists(out_dir + ".nhr"):
        return out_dir
    command = "makeblastdb -in " + genome + " -dbtype nucl -parse_seqids -out " + out_dir
    run_cmd(command=command, wait=True)
    return out_dir


def parse_blastdb(db_path, contig, start, end):
    region = str(start) + "-" + str(end)
    command = "blastdbcmd -db " + db_path + " -dbtype nucl -entry " + contig + " -range " + region
    results = []
    for line in run_cmd(command=command, wait=False):
        line = line.strip("\n")
        if line.startswith(">"):
            results.append(line.strip() + ";" + region + "\n")
        else:
            results.append(line)
    return "".join(results)


class HspListObject:
    def __init__(self, subject_hsp_list, merging_dist):
        self.merge_dist = merging_dist
        self.hsp_list = subject_hsp_list
        self.hsp_sorted = []
        self.s_start = None
        self.s_end = None
        self.q_len = None
        self.strand = None
        self.q_start = None
        self.q_end = None

    def get_sorted_attributes(self):
        self.s_start = [x["s_start"] for x in self.hsp_sorted]
        self.s_end = [x["s_end"] for x in self.hsp_sorted]
        self.q_len = [x["q_len"] for x in self.hsp_sorted]
        self.strand = [x["strand"] for x in self.hsp_sorted]
        self.q_start = [x["q_start"] for x in self.hsp_sorted]
        self.q_end = [x["q_end"] for x in self.hsp_sorted]

    def sort_hsp_list(self):
        self.hsp_sorted = sorted(self.hsp_list, key=itemgetter("s_start"))
        self.get_sorted_attributes()

    def merge_to_region(self):
        if self.hsp_sorted is []:
            self.sort_hsp_list()
        last_idx = 0
        single_region = None
        all_regions = []
        for idx in range(1, len(self.hsp_sorted)):
            if (self.s_end[idx] - self.s_start[last_idx]) < self.merge_dist \
                    and self.strand[idx] == self.strand[last_idx]:
                if single_region is None:
                    single_region = [last_idx]
                single_region.append(idx)
            else:
                if single_region is not None:
                    all_regions.append(single_region)
                    single_region = None
                else:
                    all_regions.append([idx])
            last_idx = idx
        all_regions.append(single_region)
        return all_regions

    def compute_cov(self, pos_list, denominator):
        """ coverage of aligned sequences """
        aligned_pos = 0
        pos_list = sorted(pos_list)
        for x in range(0, len(pos_list)-1,2):
            aligned_pos += pos_list[x+1] - pos_list[x]
        return round((aligned_pos / denominator) * 100)


class BlastObject:

    def __init__(self, blast_out, db_path, merging_dist=10000, flanking_dist=5000):
        self.blast_out = blast_out
        self.db_path = db_path
        self.merging_distance = merging_dist
        self.flanking_distance = flanking_dist
        self.inferred_regions = None

    def infer_regions(self):
        inferred_regions = {}
        Region = namedtuple('Region', 'contig, s_start, s_end, strand, chunk_cov, query_cov, q_len')
        for query in self.blast_out:
            if query not in inferred_regions:
                inferred_regions[query] = defaultdict(list)
            for subject in self.blast_out[query]:
                hsp_list = self.blast_out[query][subject]
                if len(hsp_list) == 1:
                    chunk_cov = 100  # -> just one start and end pos
                    query_cov = round(((hsp_list[0]["s_end"] - hsp_list[0]["s_start"]) * 100) / hsp_list[0]["q_len"])
                    strand = hsp_list[0]["strand"]
                    s_start = hsp_list[0]["s_start"] - self.flanking_distance
                    s_end = hsp_list[0]["s_end"] + self.flanking_distance
                    region = Region(contig=subject, s_start=s_start, s_end=s_end, strand=strand,
                                    chunk_cov=chunk_cov, query_cov=query_cov, q_len = hsp_list[0]["q_len"])
                    inferred_regions[query][subject].append(region)
                else:
                    # test in test suit: if subject I 2*pop(0) == [[1],[3,4]]
                    hits = HspListObject(hsp_list, self.merging_distance)
                    hits.sort_hsp_list()
                    all_regions = hits.merge_to_region()
                    for region in all_regions:
                        begin, stop = region[0], region[-1]
                        strand = hits.strand[begin]
                        s_start = hits.s_start[begin] - self.flanking_distance
                        s_end = hits.s_end[stop] + self.flanking_distance
                        q_pos = hits.q_start[begin:stop+1] + hits.q_end[begin:stop+1]
                        query_cov = hits.compute_cov(q_pos, (hits.q_len[0]))
                        chunk_cov = hits.compute_cov(q_pos, (max(q_pos) - min(q_pos)))
                        region = Region(contig=subject, s_start=s_start, s_end=s_end, strand=strand,
                                        chunk_cov=chunk_cov, query_cov=query_cov, q_len= hits.q_len[0])
                        inferred_regions[query][subject].append(region)
        self.inferred_regions = inferred_regions


# since every consensus is know in one file, tblastn should always give a hit
# just if no consensus hits --> quit entiry analysis
def run_tblastn(db_path, q_file):
    command =["tblastn", "-query", q_file,"-db", db_path, "-outfmt",
              """7 qacc sacc evalue qstart qend sstart send qlen sframe""", "-evalue", "1e-5"]

    blast_dict = {}
    results_flag = False
    for line in run_cmd(command=command, wait=False):
        if not line.startswith("#"):
            results_flag = True
            line = line.strip("\n").split("\t")
            if "-" in line[8]:
                strand = "-"
            else:
                strand = "+"
            query = line[0]
            subject = line[1]
            row = { "contig": subject,          # can be removed
                    "evalue": float(line[2]),
                    "q_start": int(line[3]),
                    "q_end": int(line[4]),
                    "s_start": int(line[5]),
                    "s_end": int(line[6]),
                    "q_len": int(line[7]),
                    "strand": strand}
            if query not in blast_dict:
                blast_dict[query] = defaultdict(list)
            blast_dict[query][subject].append(row)

    if results_flag:
        blastObject = BlastObject(blast_dict, db_path)
        return blastObject
    return None


if __name__ == "__main__":
    db = "/home/jgravemeyer/Dropbox/MSc_project/data/test_out/Predictions/c_elegans.PRJNA13758.WS254.genomic"
    query = "/home/jgravemeyer/Desktop/kog0018.consensus"
    blast = run_tblastn(db, query)
    blast.infer_regions()
    print(blast.inferred_regions)
    #for query, subject in blast.inferred_regions.items():
    #    print(query)
    #    for subject in subject.values():
    #        print(subject)
    #        for region in subject:
    #            print(region.s_start)
    print(blast.inferred_regions["KOG0018"]["I"][0].contig)
    print(blast.inferred_regions["KOG0018"]["I"][0].q_len)
    print(blast.inferred_regions["KOG0018"]["I"][0].s_start)
    print(blast.inferred_regions["KOG0018"]["I"][0].s_end)


