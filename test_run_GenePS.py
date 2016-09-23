#!/usr/bin/env python3
import unittest
from run_GenePS import *


script_path = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(script_path, "group1.makeGenePS")
test_db = os.path.join(script_path, "c_elegans.PRJNA13758.WS254.genomic")

print(test_db)
print(script_path)
