import unittest
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
try:
    from deBruijnGraph import deBruijnGraph
except ModuleNotFoundError:
    raise ModuleNotFoundError("Ensure 'directedGraph.py' exists in the '../src' directory relative to this test file.")

class TestDeBruijnGraph(unittest.TestCase):
    def test_initialization(self):
        k = 5
        dbg = deBruijnGraph(k)
        self.assertEqual(dbg.graph, {'k': 5, 'compacted': False})

    def test_inheritance(self):
        k = 3
        dbg = deBruijnGraph(k)
        self.assertTrue(isinstance(dbg, deBruijnGraph))
    
    def test_reverse_complement(self):
        k = 3
        dbg = deBruijnGraph(k)
        seq = 'TACCG'
        self.assertEqual(dbg.reverse_complement(seq), 'CGGTA')
    
    def test_find_or_add_node_existing(self):
        k = 4
        dbg = deBruijnGraph(k)
        seq = "ATGC"
        dbg.add_node(seq)
        dbg.set_node_attributes(seq, coverage=1, sequence=seq, id=0, visited=False, compacted=False, sum=0, dist=0, tip=False, isolated=False, low_coverage=False)
        node = dbg.find_or_add_node(seq)
        self.assertEqual(node['sequence'], seq)
        self.assertEqual(node['id'], 0)
        self.assertEqual(node['coverage'], 1)

    def test_find_or_add_node_new(self):
        k = 4
        dbg = deBruijnGraph(k)
        seq = "ATGC"
        node = dbg.find_or_add_node(seq)
        self.assertEqual(node['sequence'], seq)
        self.assertEqual(node['id'], 0)
        self.assertEqual(node['coverage'], 1)
        self.assertFalse(node['visited'])
        self.assertFalse(node['compacted'])
        self.assertEqual(node['sum'], 0)
        self.assertEqual(node['dist'], 0)
        self.assertFalse(node['tip'])
        self.assertFalse(node['isolated'])
        self.assertFalse(node['low_coverage'])
    
    def test_add_nodes_with_edge_new_edge(self):
        k = 4
        dbg = deBruijnGraph(k)
        seq1 = "ATGC"
        seq2 = "TGCA"
        node1 = dbg.find_or_add_node(seq1)
        node2 = dbg.find_or_add_node(seq2)
        dbg.add_nodes_with_edge(node1['sequence'], node2['sequence'])
        self.assertTrue(dbg.has_successor(seq1, seq2))
        edge_attrs = dbg.get_edge_attributes(seq1, seq2)
        self.assertEqual(edge_attrs['coverage'], 1)

    def test_add_nodes_with_edge_new_nodes(self):
        k = 4
        dbg = deBruijnGraph(k)
        seq1 = "ATGC"
        seq2 = "TGCA"
        node1 = dbg.find_or_add_node(seq1)
        node2 = dbg.find_or_add_node(seq2)
        dbg.add_nodes_with_edge(node1['sequence'], node2['sequence'])
        node1_attrs = dbg.get_node_attributes(seq1)
        node2_attrs = dbg.get_node_attributes(seq2)
        self.assertEqual(node1_attrs['sequence'], seq1)
        self.assertEqual(node2_attrs['sequence'], seq2)
        self.assertTrue(dbg.has_successor(seq1, seq2))

    def test_add_nodes_with_edge_one_new_one_old_node(self):
        k = 4
        dbg = deBruijnGraph(k)
        seq1 = "ATGC"
        seq2 = "TGCA"
        node1 = dbg.find_or_add_node(seq1)
        dbg.set_node_attributes(seq1, coverage=1, sequence=seq1, id=0, visited=False, compacted=False, sum=0, dist=0, tip=False, isolated=False, low_coverage=False)
        node2 = dbg.find_or_add_node(seq2)
        dbg.add_nodes_with_edge(node1['sequence'], node2['sequence'])
        node1_attrs = dbg.get_node_attributes(seq1)
        node2_attrs = dbg.get_node_attributes(seq2)
        self.assertEqual(node1_attrs['sequence'], seq1)
        self.assertEqual(node2_attrs['sequence'], seq2)
        self.assertTrue(dbg.has_successor(seq1, seq2))

    def test_add_nodes_with_edge_old_edge(self):
        k = 4
        dbg = deBruijnGraph(k)
        seq1 = "ATGC"
        seq2 = "TGCA"
        node1 = dbg.find_or_add_node(seq1)
        node2 = dbg.find_or_add_node(seq2)
        dbg.add_nodes_with_edge(node1['sequence'], node2['sequence'])
        dbg.add_nodes_with_edge(node1['sequence'], node2['sequence'])
        edge_attrs = dbg.get_edge_attributes(seq1, seq2)
        self.assertEqual(edge_attrs['coverage'], 2)

    def test_add_nodes_with_edge_compacted_graph(self):
        k = 4
        dbg = deBruijnGraph(k)
        dbg.graph['compacted'] = True
        seq1 = "ATGC"
        seq2 = "TGCA"
        with self.assertRaises(ValueError) as context:
            dbg.add_nodes_with_edge(seq1, seq2)
        self.assertEqual(str(context.exception), "Cannot add edge to a compacted graph.")

    def test_add_nodes_with_edge_invalid_base(self):
        k = 4
        dbg = deBruijnGraph(k)
        seq1 = "ATGC"
        seq2 = "TGCZ"
        with self.assertRaises(AssertionError) as context:
            dbg.add_nodes_with_edge(seq1, seq2)
        self.assertEqual(str(context.exception), "Invalid base in sequence 2: TGCZ")
    
    def test_add_kmer(self):
        k = 4
        dbg = deBruijnGraph(k)
        kmer = "ATGC"
        dbg.add_kmer(kmer)
        self.assertTrue(dbg.has_successor("ATG", "TGC"))
    
    def test_add_kmer_invalid_length(self):
        k = 4
        dbg = deBruijnGraph(k)
        kmer = "ATGCA"
        with self.assertRaises(AssertionError) as context:
            dbg.add_kmer(kmer)
        self.assertEqual(str(context.exception), "k-mer length does not match graph k value")

    def test_add_kmer_and_reverse_complement(self):
        k = 4
        dbg = deBruijnGraph(k)
        kmer = "ATGC"
        rc_kmer = "GCAT"
        dbg.add_kmer(kmer)
            
        # Check original k-mer nodes and edges
        self.assertTrue(dbg.has_successor("ATG", "TGC"))
        self.assertTrue(dbg.has_predecessor("TGC", "ATG"))
            
        # Check reverse complement k-mer nodes and edges
        self.assertTrue(dbg.has_successor(rc_kmer[:-1], rc_kmer[1:]))
        self.assertTrue(dbg.has_predecessor(rc_kmer[1:], rc_kmer[:-1]))

    def test_add_multiple_kmers_and_reverse_complements(self):
        k = 4
        dbg = deBruijnGraph(k)
        kmers = ["ATGC", "GCTA", "TACG"]
        for kmer in kmers:
            dbg.add_kmer(kmer)
            rc_kmer = dbg.reverse_complement(kmer)
            
            # Check original k-mer nodes and edges
            self.assertTrue(dbg.has_successor(kmer[:-1], kmer[1:]))
            
            # Check reverse complement k-mer nodes and edges
            self.assertTrue(dbg.has_successor(rc_kmer[:-1], rc_kmer[1:]))

    def test_add_kmer_reverse_complement_edge_coverage(self):
        k = 4
        dbg = deBruijnGraph(k)
        kmer = "ATGC"
        dbg.add_kmer(kmer)
        dbg.add_kmer(kmer)  # Add the same k-mer again
        rc_kmer = dbg.reverse_complement(kmer)
            
        # Check edge coverage for original k-mer
        edge_attrs = dbg.get_edge_attributes("ATG", "TGC")
        self.assertEqual(edge_attrs['coverage'], 2)
            
        # Check edge coverage for reverse complement k-mer
        rc_edge_attrs = dbg.get_edge_attributes(rc_kmer[:-1], rc_kmer[1:])
        self.assertEqual(rc_edge_attrs['coverage'], 2)

    def test_add_kmer_reverse_complement_with_overlap(self):
        k = 4
        dbg = deBruijnGraph(k)
        kmer1 = "ATGC"
        kmer2 = "GCAT"  # Overlaps with reverse complement of kmer1
        dbg.add_kmer(kmer1)
        dbg.add_kmer(kmer2)
            
        # Check original k-mer nodes and edges
        self.assertTrue(dbg.has_successor("ATG", "TGC"))
        self.assertTrue(dbg.has_successor("GCA", "CAT"))
            
        # Check reverse complement k-mer nodes and edges
        rc_kmer1 = dbg.reverse_complement(kmer1)
        rc_kmer2 = dbg.reverse_complement(kmer2)
        self.assertTrue(dbg.has_successor(rc_kmer1[:-1], rc_kmer1[1:]))
        self.assertTrue(dbg.has_successor(rc_kmer2[:-1], rc_kmer2[1:]))

    def test_add_kmers_from_fasta(self):

        k = 5
        dbg = deBruijnGraph(k)

        # small fasta file 
        small_fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/test_add_kmers.fa'))

        # Add kmers from small fasta file
        dbg.add_kmers_from_fasta(small_fasta)

        # Check that all k-mers and their reverse complements were added
        kmers = ["ATGCA", "TGCAG", "GCTAT", "CTATT"]
        for kmer in kmers:
            rc_kmer = dbg.reverse_complement(kmer)

            # Check original k-mer nodes and edges
            self.assertTrue(dbg.has_successor(kmer[:-1], kmer[1:]))

            # Check reverse complement k-mer nodes and edges
            self.assertTrue(dbg.has_successor(rc_kmer[:-1], rc_kmer[1:]))

    def test_add_kmers_from_fasta_invalid_file(self):
        k = 5
        dbg = deBruijnGraph(k)
        invalid_fasta = "invalid_file.fa"
        with self.assertRaises(FileNotFoundError) as context:
            dbg.add_kmers_from_fasta(invalid_fasta)
        self.assertEqual(str(context.exception), f"[Errno 2] No such file or directory: '{invalid_fasta}'")
        
    def test_add_degrees_attributes(self):
        k = 4
        dbg = deBruijnGraph(k)
        
        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)
        
        # Add degrees to nodes
        dbg.add_degrees_attributes()
        
        # Hard-code checks for in and out degrees
        node1_attrs = dbg.get_node_attributes("ATG")
        self.assertEqual(node1_attrs['out_degree'], 1)
        self.assertEqual(node1_attrs['in_degree'], 0)
        
        node2_attrs = dbg.get_node_attributes("TGC")
        self.assertEqual(node2_attrs['out_degree'], 1)
        self.assertEqual(node2_attrs['in_degree'], 1)
        
        node3_attrs = dbg.get_node_attributes("GCA")
        self.assertEqual(node3_attrs['out_degree'], 1)
        self.assertEqual(node3_attrs['in_degree'], 1)
        
        node4_attrs = dbg.get_node_attributes("CAT")
        self.assertEqual(node4_attrs['out_degree'], 0)
        self.assertEqual(node4_attrs['in_degree'], 1)
        
    def test_query_kmers_from_fasta(self):
        k = 5
        dbg = deBruijnGraph(k)

        # Load small FASTA file for testing
        test_fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/test_query_kmers.fa'))

        # Query k-mers from the FASTA file
        query_kmers = list(dbg.query_kmers_from_fasta(test_fasta))

        # Check that the correct k-1mers are returned
        expected_kmers = ["ATGC", "TGCA", "GCAT", "CATG"]
        for kmer in expected_kmers:
            self.assertIn(kmer, query_kmers)

        # Check reverse complement k-1mers
        for kmer in expected_kmers:
            rc_kmer = dbg.reverse_complement(kmer)
            self.assertIn(rc_kmer, query_kmers)


    def test_identify_query_connected_nodes(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT", "CATG"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Query connected nodes
        query_kmers = ["ATG"]
        subgraph_dist = dbg.identify_query_connected_nodes(query_kmers)

        # Check distances
        self.assertEqual(subgraph_dist["ATG"], 0)
        self.assertEqual(subgraph_dist["TGC"], 1)
        self.assertEqual(subgraph_dist["GCA"], 2)
        self.assertEqual(subgraph_dist["CAT"], 1)

    def test_filter_graph_to_subgraph(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT", "CATG", "GATT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Identify query-connected nodes
        query_kmers = {"ATG": 0}
        subgraph_dist = dbg.identify_query_connected_nodes(query_kmers)

        # Filter graph to subgraph
        dbg.filter_graph_to_subgraph(subgraph_dist)

        # Check that only query-connected nodes remain
        remaining_nodes = list(dbg.get_nodes())
        self.assertIn("ATG", remaining_nodes)
        self.assertIn("TGC", remaining_nodes)
        self.assertIn("GCA", remaining_nodes)
        self.assertIn("CAT", remaining_nodes)
        self.assertNotIn("GAT", remaining_nodes)
        self.assertNotIn("ATT", remaining_nodes)
        # Check that only query-connected edges remain
        self.assertTrue(dbg.has_successor("ATG", "TGC"))
        self.assertTrue(dbg.has_successor("TGC", "GCA"))
        self.assertTrue(dbg.has_successor("GCA", "CAT"))  
        self.assertTrue(dbg.has_successor("CAT", "ATG"))

    def test_add_kmers_from_fasta_with_min_count(self):
        k = 5
        dbg = deBruijnGraph(k)

        # Create a small FASTA file for testing
        test_fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/test_min_count.fa'))

        # Add k-mers with a minimum count of 2
        dbg.add_kmers_from_fasta(test_fasta, min_count=2)

        # Check that only k-mers with sufficient count are added
        self.assertTrue(dbg.has_successor("ATGC", "TGCA"))  # Edge count is 2
        self.assertTrue(dbg.has_successor("TGCA", "GCAT"))  # Edge count is 2
        self.assertFalse(dbg.has_successor("GCAT", "CATG"))  # Edge count is 1, not added

    def test_add_degrees_attributes_with_isolated_node(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Add an isolated node
        dbg.add_node("ISOL")

        # Add degrees to nodes
        dbg.add_degrees_attributes()

        # Check degrees for the isolated node
        isolated_node_attrs = dbg.get_node_attributes("ISOL")
        self.assertEqual(isolated_node_attrs['out_degree'], 0)
        self.assertEqual(isolated_node_attrs['in_degree'], 0)

    def test_identify_query_subgraph(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT", "CATG", "GATT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Query k-mers
        query_kmers = ["ATG"]

        # Identify query subgraph
        dbg.identify_query_subgraph(query_kmers)

        # Check that only query-connected nodes remain
        remaining_nodes = list(dbg.get_nodes())
        self.assertIn("ATG", remaining_nodes)
        self.assertIn("TGC", remaining_nodes)
        self.assertIn("GCA", remaining_nodes)
        self.assertIn("CAT", remaining_nodes)
        self.assertNotIn("GAT", remaining_nodes)
        self.assertNotIn("ATT", remaining_nodes)

        # Check that only query-connected edges remain
        self.assertTrue(dbg.has_successor("ATG", "TGC"))
        self.assertTrue(dbg.has_successor("TGC", "GCA"))
        self.assertTrue(dbg.has_successor("GCA", "CAT"))
        self.assertTrue(dbg.has_successor("CAT", "ATG"))

        # Check distances in the subgraph
        node_attrs = dbg.get_node_attributes("ATG")
        self.assertEqual(node_attrs['query_distance'], 0)
        node_attrs = dbg.get_node_attributes("TGC")
        self.assertEqual(node_attrs['query_distance'], 1)
        node_attrs = dbg.get_node_attributes("GCA")
        self.assertEqual(node_attrs['query_distance'], 2)
        node_attrs = dbg.get_node_attributes("CAT")
        self.assertEqual(node_attrs['query_distance'], 1)

    def test_identify_query_subgraph_no_query_kmers(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT", "CATG"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Query k-mers not in the graph
        query_kmers = ["AAA"]

        # Attempt to identify query subgraph
        with self.assertRaises(ValueError) as context:
            dbg.identify_query_subgraph(query_kmers)
        self.assertEqual(str(context.exception), "No query k-mers found in the graph.")

    def test_update_all_node_attributes(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Update all node attributes
        dbg.update_all_node_attributes(test_attr=True, value=42)

        # Check that all nodes have the updated attributes
        for kmer in ["ATG", "TGC", "GCA", "CAT"]:
            node_attrs = dbg.get_node_attributes(kmer)
            self.assertTrue(node_attrs['test_attr'])
            self.assertEqual(node_attrs['value'], 42)
            
    def test_identify_query_boundaries(self):
        k = 5
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGCA", "TGCAG", "GCAGT", "CAGTT", "CTGAT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Create a small FASTA file for testing
        test_fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/test_query_boundaries.fa'))

        # Identify query boundaries
        dbg.identify_query_boundaries(test_fasta)

        # Check graph attributes for query start and end
        self.assertEqual(dbg.graph['query_start'], "ATGC")
        self.assertEqual(dbg.graph['query_end'], "AGTT")

        # Check edge attributes for query edges
        self.assertTrue(dbg.get_edge_attributes("ATGC", "TGCA")['query_edge'])
        self.assertTrue(dbg.get_edge_attributes("TGCA", "GCAG")['query_edge'])
        self.assertTrue(dbg.get_edge_attributes("GCAG", "CAGT")['query_edge'])
        self.assertTrue(dbg.get_edge_attributes("CAGT", "AGTT")['query_edge'])

        # Check edge locations
        self.assertEqual(dbg.get_edge_attributes("ATGC", "TGCA")['query_edge_location'], [0])
        self.assertEqual(dbg.get_edge_attributes("TGCA", "GCAG")['query_edge_location'], [1])
        self.assertEqual(dbg.get_edge_attributes("GCAG", "CAGT")['query_edge_location'], [2])
        self.assertEqual(dbg.get_edge_attributes("CAGT", "AGTT")['query_edge_location'], [3])

        # Check visited nodes
        for node in ["ATGC", "TGCA", "GCAG", "CAGT", "AGTT"]:
            self.assertTrue(dbg.get_node_attributes(node)['visited'])
        # Check non visited nodes
        for node in ["CTGA", "TGAT"]:
            self.assertFalse(dbg.get_node_attributes(node)['visited'])
            
    def test_identify_query_boundaries_missing_edge(self):
        k = 5
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGCA", "TGCAG", "GCAGT"]  # Missing "CAGT" -> "AGTT"
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Import a small FASTA file for testing
        test_fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/test_query_boundaries.fa'))

        # Attempt to identify query boundaries
        with self.assertRaises(ValueError) as context:
            dbg.identify_query_boundaries(test_fasta)

        self.assertEqual(
            str(context.exception),
            "Edge CAGT->AGTT does not exist in the graph. Query sequence may not be fully represented."
            )
    
    def test_get_next_node_successor(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)
        
        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Mark nodes as visited
        dbg.set_node_attributes("ATG", visited=False)
        dbg.set_node_attributes("TGC", visited=False)
        dbg.set_node_attributes("GCA", visited=False)

        # Test getting the next successor node
        next_node = dbg.get_next_node("ATG", dir=1)
        self.assertEqual(next_node, "TGC")

        # Mark the successor as visited and test again
        dbg.set_node_attributes("TGC", visited=True)
        next_node = dbg.get_next_node("ATG", dir=1)
        self.assertEqual(next_node, "TGC")
        
        # Mark successor of successor as visited and test again
        dbg.set_node_attributes("GCA", visited=True)
        next_node = dbg.get_next_node("ATG", dir=1)
        self.assertIsNone(next_node)


    def test_get_next_node_predecessor(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)
        
        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Mark nodes as visited
        dbg.set_node_attributes("TGC", visited=False)
        dbg.set_node_attributes("ATG", visited=False)

        # Test getting the next predecessor node
        next_node = dbg.get_next_node("TGC", dir=0)
        self.assertEqual(next_node, "ATG")

        # Mark the predecessor as visited and test again
        dbg.set_node_attributes("ATG", visited=True)
        next_node = dbg.get_next_node("TGC", dir=0)
        self.assertIsNone(next_node)

    def test_get_next_node_no_successors(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add a single k-mer to the graph
        dbg.add_kmer("ATGC")
        
        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Test getting the next successor node when none exist
        next_node = dbg.get_next_node("TGC", dir=1)
        self.assertIsNone(next_node)

    def test_get_next_node_no_predecessors(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add a single k-mer to the graph
        dbg.add_kmer("ATGC")
        
        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Test getting the next predecessor node when none exist
        next_node = dbg.get_next_node("ATG", dir=0)
        self.assertIsNone(next_node)

    def test_get_next_node_multiple_successors(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "TGCT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)
            
        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Mark nodes as visited
        dbg.set_node_attributes("TGC", visited=False)
        dbg.set_node_attributes("GCA", visited=False)
        dbg.set_node_attributes("GCT", visited=False)

        # Test getting the next successor node
        next_node = dbg.get_next_node("TGC", dir=1)
        self.assertIn(next_node, ["GCA", "GCT"])

        # Mark one successor as visited and test again
        dbg.set_node_attributes(next_node, visited=True)
        next_node = dbg.get_next_node("TGC", dir=1)
        self.assertIn(next_node, ["GCA", "GCT"])
        self.assertNotEqual(next_node, dbg.get_node_attributes(next_node)["visited"])

    def test_get_next_node_multiple_predecessors(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "CTGC", "GTGC"]
        for kmer in kmers:
            dbg.add_kmer(kmer)
        
        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Mark nodes as visited
        dbg.set_node_attributes("TGC", visited=False)
        dbg.set_node_attributes("ATG", visited=False)
        dbg.set_node_attributes("CTG", visited=False)
        dbg.set_node_attributes("GTG", visited=False)

        # Test getting the next predecessor node
        next_node = dbg.get_next_node("TGC", dir=0)
        self.assertIn(next_node, ["ATG", "CTG", "GTG"])

        # Mark one predecessor as visited and test again
        dbg.set_node_attributes(next_node, visited=True)
        next_node = dbg.get_next_node("TGC", dir=0)
        self.assertIn(next_node, ["ATG", "CTG", "GTG"])
        self.assertNotEqual(next_node, dbg.get_node_attributes(next_node)["visited"])
        
    def test_walk_forward(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Mark nodes as visited
        dbg.set_node_attributes("ATG", visited=True)
        dbg.set_node_attributes("TGC", visited=False)
        dbg.set_node_attributes("GCA", visited=False)
        dbg.set_node_attributes("CAT", visited=False)

        # Perform a forward walk starting from "ATG"
        sequence = dbg.walk("ATG", dir=1)
        self.assertEqual(sequence, "CAT")

    def test_walk_backward(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGC", "TGCA", "GCAT"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Mark nodes as visited
        dbg.set_node_attributes("CAT", visited=True)
        dbg.set_node_attributes("GCA", visited=False)
        dbg.set_node_attributes("TGC", visited=False)
        dbg.set_node_attributes("ATG", visited=False)

        # Perform a backward walk starting from "CAT"
        sequence = dbg.walk("CAT", dir=0)
        self.assertEqual(sequence, "ATG")

    def test_walk_no_successors(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add a single k-mer to the graph
        dbg.add_kmer("ATGC")

        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Mark the node as visited
        dbg.set_node_attributes("ATG", visited=True)
        # Mark next node as visited
        dbg.set_node_attributes("TGC", visited=True)

        # Perform a forward walk starting from "ATG"
        sequence = dbg.walk("ATG", dir=1)
        self.assertEqual(sequence, "")

    def test_walk_no_predecessors(self):
        k = 4
        dbg = deBruijnGraph(k)

        # Add a single k-mer to the graph
        dbg.add_kmer("ATGC")

        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Mark the node as visited
        dbg.set_node_attributes("TGC", visited=True)
        dbg.set_node_attributes("ATG", visited=True)

        # Perform a backward walk starting from "TGC"
        sequence = dbg.walk("TGC", dir=0)
        self.assertEqual(sequence, "")

    def test_walk_with_branching(self):
        k = 5
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGAT", "TGATT", "TGATG"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Mark start node as visited
        dbg.set_node_attributes("ATGA", visited=True)

        # Perform a forward walk starting from "ATGA"
        sequence = dbg.walk("ATGA", dir=1)
        self.assertEqual(sequence, "TG")
        
        # Add coverage to TT branch, should pick that now
        dbg.set_edge_attributes("TGAT", "GATT", coverage=2) 
        sequence = dbg.walk("ATGA", dir=1)
        self.assertEqual(sequence, "TT")
        
    def test_branch_with_unequal_coverage(self):
        k = 5
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGAT", "TGATT", "TGATG"]
        for kmer in kmers:
            dbg.add_kmer(kmer)
        
        # Add additional coverage to one branch
        dbg.add_kmer("TGATT")

        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Check edge coverage
        edge_attrs_gatc = dbg.get_edge_attributes("TGAT", "GATG")
        edge_attrs_gatt = dbg.get_edge_attributes("TGAT", "GATT")
        self.assertEqual(edge_attrs_gatc['coverage'], 1)
        self.assertEqual(edge_attrs_gatt['coverage'], 2)
        self.assertEqual(dbg.get_node_attributes("GATG")["out_degree"], dbg.get_node_attributes("GATT")["out_degree"])

        # Perform a forward walk starting from "TGAT"
        next_node = dbg.get_next_node("TGAT", dir=1)
        self.assertEqual(next_node, "GATT")  # Higher coverage branch is preferred

    def test_branch_with_unequal_degrees(self):
        k = 5
        dbg = deBruijnGraph(k)

        # Add k-mers to the graph
        kmers = ["ATGAT", "TGATT", "TGATC"]
        for kmer in kmers:
            dbg.add_kmer(kmer)

        # Add another TGATT, so its coverage is higher
        dbg.add_kmer("TGATT")
        
        # Update degree attributes of graph
        dbg.add_degrees_attributes()

        # Check degrees of successors
        succs = dbg.successors("TGAT")
        self.assertEqual(len(list(succs)), 2)
        self.assertEqual(dbg.get_node_attributes("GATT")['out_degree'], 0)
        self.assertEqual(dbg.get_node_attributes("GATC")['out_degree'], 1) # due to reverse complement

        # Perform a forward walk starting from "TGAT"
        next_node = dbg.get_next_node("TGAT", dir=1)
        self.assertEqual(next_node, "GATC")  # Higher out-degree branch is preferred over higher coverage

    def test_compact_except_query(self):
        k = 5
        dbg = deBruijnGraph(k)
        # Add k-mers to the graph
        kmers = ["ATGAT", "TGATT", "TGATC", "GATTA", "ATTAA", "AATGA", "CAATG", "GCAAT", "TCAAT", "GGCAA", "GTCAA"]
        for kmer in kmers:
            dbg.add_kmer(kmer)            
        
        # visualize graph
        dbg.visualize_graph("test_compact.png")
        
        # Import a small FASTA query file for testing
        test_fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/test_compact.fa'))
        dbg.identify_query_boundaries(test_fasta)
        
        # Test query end and start
        self.assertEqual(dbg.graph['query_start'], "ATGA")
        self.assertEqual(dbg.graph['query_end'], "TGAT")
        
        # Add degrees
        dbg.add_degrees_attributes()
        
        # Test compaction
        dbg.compact_except_query()
        
        # visualize graph
        dbg.visualize_graph("test_compacted.png")
        
        # Assert correct compacted sequences
        compact1=dbg.get_node_attributes("GATT")
        compact2=dbg.get_node_attributes("ATCA")
        compact3=dbg.get_node_attributes("CAAT")
        compact4=dbg.get_node_attributes("GGCA")
        compact5=dbg.get_node_attributes("GTCA")
        compact6=dbg.get_node_attributes("TTGC")
        compact7=dbg.get_node_attributes("TTGA")
        norm2=dbg.get_node_attributes("GATC")
        self.assertEqual(compact1["sequence"], "TAATC")
        self.assertEqual(compact2["sequence"], "ATTG")
        self.assertEqual(compact3["sequence"], "TG")
        self.assertEqual(compact4["sequence"], "GGCAA")
        self.assertEqual(compact5["sequence"], "GTCAA")
        self.assertEqual(compact6["sequence"], "CC")
        self.assertEqual(compact7["sequence"], "AC")
        self.assertTrue(compact1["compacted"])
        self.assertFalse(norm2["compacted"])
        
    def test_compact_walk(self):
        k = 5
        dbg = deBruijnGraph(k)
        # Add k-mers to the graph
        kmers = ["ATGAT", "TGATT", "TGATC", "GATTA", "ATTAA", "AATGA", "CAATG", "GCAAT", "TCAAT", "GGCAA", "GTCAA"]
        for kmer in kmers:
            dbg.add_kmer(kmer)            
        
        # visualize graph
        dbg.visualize_graph("test_compact.png")
        
        # Import a small FASTA query file for testing
        test_fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/test_compact.fa'))
        dbg.identify_query_boundaries(test_fasta)
        
        # Test query end and start
        self.assertEqual(dbg.graph['query_start'], "ATGA")
        self.assertEqual(dbg.graph['query_end'], "TGAT")
        
        # Add degrees
        dbg.add_degrees_attributes()
                
        # Compact
        dbg.compact_except_query() 
        
        # Assemble
        assembly = dbg.assemble_from_query()
        self.assertEqual(assembly, "GGCAATGATTAATCATTGAC")
        # there's a couple of branches but it'll choose the one with the lowest alphabet option
    
    def test_not_compact_walk(self):
        # if this graph isn't compacted first, get a smaller output contig
        k = 5
        dbg = deBruijnGraph(k)
        # Add k-mers to the graph
        kmers = ["ATGAT", "TGATT", "TGATC", "GATTA", "ATTAA", "AATGA", "CAATG", "GCAAT", "TCAAT", "GGCAA", "GTCAA"]
        for kmer in kmers:
            dbg.add_kmer(kmer)            
        
        # visualize graph
        dbg.visualize_graph("test_compact.png")
        
        # Import a small FASTA query file for testing
        test_fasta = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/test_compact.fa'))
        dbg.identify_query_boundaries(test_fasta)
        
        # Test query end and start
        self.assertEqual(dbg.graph['query_start'], "ATGA")
        self.assertEqual(dbg.graph['query_end'], "TGAT")
        
        # Add degrees
        dbg.add_degrees_attributes()
        
        # Assembly
        assembly = dbg.assemble_from_query()
        
        # Will miss the bridge option from "GATT"
        self.assertEqual(assembly, "GGCAATGATCATTGAC")

        
        

        
if __name__ == '__main__':
    unittest.main()