#!/usr/bin/env python3
import os
import re
from operator import itemgetter
from collections import defaultdict, namedtuple
from run_command import run_cmd, tempdir, check_programs


def run_genblastg(query, genome, outdir, prefix):
    print(os.getcwd())
    cmd = "genblast_v138_linux_x86_64 -p genblastg -q {} -t {}  -gff -pro -cdna -r 1 -e 1e-5 -pid -o {}".format(query, genome, os.path.join(outdir, prefix))
    run_cmd(command=cmd, wait=True)
    output_files = [filename for filename in os.listdir(outdir) if filename.startswith(prefix) or "blast" in filename]
    return output_files


if __name__ == "__main__":
    check_programs("genblast_v138_linux_x86_64")
    query = "/home/jgravemeyer/Dropbox/MSc_project/data/geblast_test/B0205.10_Inf5.0_OG0011445.fasta"
    genome = "/home/jgravemeyer/Dropbox/MSc_project/res/c_elegans.PRJNA13758.WS254.genomic.fa"
    outdir = "/home/jgravemeyer/Dropbox/MSc_project/data/geblast_test"
    prefix = "python_genblast_test"
    with tempdir() as tmp:
        print(run_genblastg(query, genome, outdir, prefix))
    pass
