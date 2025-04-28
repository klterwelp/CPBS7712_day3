import unittest
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from alignment import build_index
from alignment import scan_read_for_seeds
from alignment import extending_seed
from alignment import extend_seed
from alignment import identify_most_common_diagonal
from alignment import prioritize_seeds_by_diagonal
from alignment import align_seq_to_contig
from alignment import align_read_to_contig
from alignment import align_reads_to_contig
from alignment import align_reads_to_contig_parallel
from alignment import read_fasta_with_sequence_name
from alignment import reverse_complement


class TestBuildIndex(unittest.TestCase):
    def test_build_index_basic(self):
        seq = "TTTTAAAACCCCGGGGTTTT"
        k = 4
        expected_index = {
                "TTTT": [0, 16],
                "TTTA": [1],
                "TTAA": [2],
                "TAAA": [3],
                "AAAA": [4],
                "AAAC": [5],
                "AACC": [6],
                "ACCC": [7],
                "CCCC": [8],
                "CCCG": [9],
                "CCGG": [10],
                "CGGG": [11],
                "GGGG": [12],
                "GGGT": [13],
                "GGTT": [14],
                "GTTT": [15]
            }
        index = build_index(seq, k)
        self.assertEqual(index, expected_index)

    def test_build_index_empty_sequence(self):
        seq = ""
        k = 4
        # expect error
        with self.assertRaises(ValueError):
            build_index(seq, k)

    def test_build_index_k_larger_than_sequence(self):
        seq = "ACGT"
        k = 5
        # expect error
        with self.assertRaises(ValueError):
            build_index(seq, k)

    def test_build_index_repeated_kmers(self):
            seq = "AAAAAA"
            k = 2
            expected_index = {
                "AA": [0, 1, 2, 3, 4]
            }
            index = build_index(seq, k)
            self.assertEqual(index, expected_index)

    def test_build_index_single_kmer(self):
            seq = "ACGT"
            k = 4
            expected_index = {
                "ACGT": [0]
            }
            index = build_index(seq, k)
            self.assertEqual(index, expected_index)

class TestScanReadForSeeds(unittest.TestCase):
    def test_scan_read_for_seeds_basic(self):
        read_seq = "TTAATTTT"
        read_len = len(read_seq)
        # 0 TTAA, 2 AATT 4 TTTT
        contig_index = {
            "CCCC": [0],
            "CCCT": [1],
            "CCTT": [2],
            "CTTA": [3],
            "TTAA": [4],
            "TAAT": [5],
            "AATT": [6],
            "ATTT": [7],
            "TTTT": [8, 9]
        }
        k = 4
        I = 2
        expected_seeds = [(4, 0), (6, 2), (8, 4), (9, 4)]
        seeds = scan_read_for_seeds(read_seq, read_len, contig_index, k, I)
        self.assertEqual(seeds, expected_seeds)
    
    def test_scan_read_for_seeds_no_seeds(self):
        read_seq = "ACGT"
        read_len = len(read_seq)
        contig_index = {
            "TTAA": [0],
            "TAAT": [1],
            "AATT": [2],
            "ATTT": [3],
            "TTTT": [4]
        }
        k = 4
        I = 2
        expected_seeds = []
        seeds = scan_read_for_seeds(read_seq, read_len, contig_index, k, I)
        self.assertEqual(seeds, expected_seeds)

class TestExtendingSeed(unittest.TestCase):
    def test_extending_seed_basic(self):
        query = "TTAAGGGG"
        read = "TTAATTTT"
        q_pos = 1
        r_pos = 1
        q_len = len(query)
        r_len = len(read)
        max_mismatch = 1
        expected_right_index = 2
        expected_left_index = 1
        best_index_right = extending_seed(query, read, q_pos, r_pos, q_len, r_len, 1, max_mismatch)
        best_index_left = extending_seed(query, read, q_pos, r_pos, q_len, r_len, -1, max_mismatch)
        self.assertEqual(best_index_right, expected_right_index)
        self.assertEqual(best_index_left, expected_left_index)

class TestExtendSeed(unittest.TestCase):
    def test_extend_seed_basic(self):
        seed = (0, 0)
        query = "TTAAGGGG"
        q_len = len(query)
        read = "TTAATTTT"
        r_len = len(read)
        max_mismatch = 1
        expected_extended_seed = (0, 3, 0, 3, 4)
        extended_seed = extend_seed(seed, query, q_len, read, r_len, max_mismatch, length_threshold = 4)
        self.assertEqual(extended_seed, expected_extended_seed)

class TestIdentifyMostCommonDiagonal(unittest.TestCase):
    def test_identify_most_common_diagonal(self):
        seeds = [(0, 0), (1, 1), (2, 2)]
        expected_diagonal = 0 
        actual_diagonal = identify_most_common_diagonal(seeds) 
        self.assertEqual(actual_diagonal, expected_diagonal)

class TestPrioritizeSeedsByDiagonal(unittest.TestCase):
    def test_prioritize_seeds_by_diagonal(self):
        seeds = [(0, 0), (1, 2), (2, 3)]
        expected_seed = (1, 2)
        prioritized_seed = prioritize_seeds_by_diagonal(seeds)
        self.assertEqual(prioritized_seed, expected_seed)

class TestAlignSeqToContig(unittest.TestCase):
    def test_align_seq_to_contig(self):
        read = "TTAAGCCCTTTAAAA"
        read_len = len(read)
        contig = "TTAAGCCCCGTAAACCCC"
        contig_len = len(contig)
        contig_index = build_index(contig, k=4)
        expected_alignment = (0, 7, 0, 7, 8)
        alignment = align_seq_to_contig(read, read_len, contig, contig_len, contig_index, k=4, I=1)
        self.assertEqual(alignment, expected_alignment)

class TestAlignReadToContig(unittest.TestCase):
    def test_align_read_to_contig_f(self):
        read = "CGTTAAGG"
        contig = "ATGCAACGTTAAGGTTA"
        contig_len = len(contig)
        contig_index = build_index(contig, k=4)
        alignment = align_read_to_contig(read=read, 
                                         contig=contig, 
                                         contig_len=contig_len, 
                                         contig_index = contig_index, 
                                         k=4, 
                                         I=1)
        expected_alignment = (0,7, 6, 13, 8, 1)
        self.assertEqual(alignment, expected_alignment)
    
    def test_align_read_to_contig_r(self):
        read = "CCTTAACG"
        contig = "ATGCAACGTTAAGGTTA"
        contig_len = len(contig)
        contig_index = build_index(contig, k=4)
        alignment = align_read_to_contig(read=read, 
                                         contig=contig, 
                                         contig_len=contig_len, 
                                         contig_index = contig_index, 
                                         k=4, 
                                         I=1)
        expected_alignment = (7,0, 6, 13, 8, -1)
        self.assertEqual(alignment, expected_alignment)
    
    def test_align_read_to_contig_no_alignment(self):
        read = "TCGGATT"
        contig = "ATGCAACGTTAAGGTTA"
        contig_len = len(contig)
        contig_index = build_index(contig, k=4)
        alignment = align_read_to_contig(read=read, 
                                         contig=contig, 
                                         contig_len=contig_len, 
                                         contig_index = contig_index, 
                                         k=4, 
                                         I=1)
        expected_alignment = None
        self.assertEqual(alignment, expected_alignment)
    
class TestAlignReadsToContig(unittest.TestCase):
    def test_align_reads_to_contig(self):
        reads = {
            "read1": "CGTTAAGG",
            "read2": "ATGCAACG",
            "read3": "CCTTAAGG",
            "read4": "AAAAAAA" 
            }
        
        contig = "GGGGCCTTAACGATGCAACGAAAAAA"
        contig_name = "contig1"
        expected_alignments = [
            ("read1", "contig1", 7, 0, 4, 11, 8, -1),
            ("read2", "contig1", 0, 7, 12, 19, 8, 1), 
            ("read4", "contig1", 0, 5, 20, 25, 6, 1)
        ]
        k = 4
        I = 1
        alignments = align_reads_to_contig(reads, contig, contig_name, k, I)
        print(alignments)
        self.assertEqual(alignments, expected_alignments)

class TestAlignReadsToContigParallel(unittest.TestCase):
    def test_align_reads_to_contig(self):
        reads = {
            "read1": "CGTTAAGG",
            "read2": "ATGCAACG",
            "read3": "CCTTAAGG",
            "read4": "AAAAAAA" 
            }
        
        contig = "GGGGCCTTAACGATGCAACGAAAAAA"
        contig_name = "contig1"
        expected_alignments = [
            ("read1", "contig1", 7, 0, 4, 11, 8, -1),
            ("read2", "contig1", 0, 7, 12, 19, 8, 1), 
            ("read4", "contig1", 0, 5, 20, 25, 6, 1)
        ]
        k = 4
        I = 1
        alignments = align_reads_to_contig_parallel(reads, contig, contig_name, k, I)
        print(alignments)
        self.assertEqual(alignments, expected_alignments)

        
if __name__ == "__main__":
    unittest.main()
