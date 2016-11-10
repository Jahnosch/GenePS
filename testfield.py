#!/usr/bin/env python3
import unittest
import os
from unittest import mock
import run_GenePS


def simple_urandom(length):
    return 'f' * length


class TestRandom(unittest.TestCase):
    @mock.patch('run_GenePS.get_phmm_score', return_value={"abc" : 0.9})
    def test_summary(self, clusterwise_prediction_scoring_function):
        assert run_GenePS.clusterwise_prediction_scoring("asdf", "dfa", "dafsd") == {"abc" : 0.9}



class TestDataProvider(unittest.TestCase):
    def setUp(self):
        self.mock_data_base = mock.patch(run_GenePS, "data_base", return_value="hallo")
        self.overseer = run_GenePS.Overseer("a", "d", "h")
        self.letter = "dfa"

    def test_1(self):
        with self.mock_data_base:
            print(self.overseer.blast_all_consensus())






if __name__ == '__main__':
    unittest.main()

