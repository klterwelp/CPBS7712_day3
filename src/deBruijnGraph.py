# import from directedGraph
from directedGraph import diGraph
from collections import defaultdict
from collections import deque

# class deBruijnGraph inherits from diGraph
class deBruijnGraph(diGraph):
    # initialize graph 
    def __init__(self, k):
        # call parent class constructor
        super().__init__(k = k, compacted = False)
        # set attributes 
        self.newid = 0
        
# generate reverse complement of sequence
    def reverse_complement(self, seq):
        """
        Generate the reverse complement of a DNA sequence.

        Args:
            seq (str): A string representing the DNA sequence. 
                       The sequence should only contain the characters 'A', 'T', 'C', and 'G'.

        Returns:
            str: The reverse complement of the input DNA sequence.

        Example:
            >>> reverse_complement("ATCG")
            "CGAT"
        """
        # dictionary for reverse complement
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        # reverse complement sequence
        return ''.join([complement[base] for base in reversed(seq)])

# find or generate new de bruijn node 
    def find_or_add_node(self, seq):
        """
        Finds a node in the graph by its sequence or adds it if it does not exist.

        This method checks if a node with the given sequence exists in the graph.
        If the node exists, it retrieves and returns the node's attributes. If the
        node does not exist, it adds the node to the graph, initializes its attributes,
        and then returns the newly added node's attributes.

        Args:
            seq (str): The sequence representing the node.

        Returns:
            dict: A dictionary containing the attributes of the node, either existing
                  or newly added. The attributes include:
                  - 'coverage' (int): Coverage count for the node.
                  - 'sequence' (str): The sequence of the node.
                  - 'id' (int): A unique identifier for the node.
                  - 'visited' (bool): Whether the node has been visited.
                  - 'compacted' (bool): Whether the node has been compacted.
                  - 'sum' (int): Coverage of compacted node.
                  - 'dist' (int): Merged distance of compacted node.
                  - 'query_subgraph' (bool): Whether node belongs to subgraph.
                  - 'query_distance' (int): Distance of node from query. 
                  - 'query_merged' (int): Number of merged query k-mers.
                  - 'tip' (bool): Whether the node is a tip.
                  - 'isolated' (bool): Whether the node is isolated.
                  - 'low_coverage' (bool): Whether the node has low coverage.
        Raises:
                ValueError: If the node was not successfully added to the graph.
        """
        # Check if the node exists
        node_exists = self.has_node(seq)
        if node_exists:
            # Retrieve the node attributes
            node = self.get_node_attributes(seq)
            return node
        else:
            # Add the node to the graph
            self.add_node(seq)
            # Define node attributes
            node_attributes = {
                'coverage': 1,
                'sequence': seq,
                'id': self.newid,
                'visited': False,
                'compacted': False,
                'sum': 0,
                'dist': 0,
                'query_subgraph': False,
                'query_distance': None,
                'query_merged': 0,
                'tip': False,
                'isolated': False,
                'low_coverage': False
            }
            # Set the node attributes
            self.set_node_attributes(seq, **node_attributes)
            # Increment newid
            self.newid += 1
            # Retrieve and return the newly added node
            try: 
                node = self.get_node_attributes(seq)
                return node
            except KeyError:
                raise ValueError(f"Node {seq} was not added.")
            
            
# add de bruijn edge
    def add_nodes_with_edge(self, seq1, seq2):
        """
        Adds nodes and an edge between them in the de Bruijn graph.
        This method ensures that the graph is not compacted before adding nodes or edges.
        It either finds or creates nodes for the given sequences `seq1` and `seq2`, and
        then adds an edge between them. If the edge already exists, its coverage attribute
        is incremented. If the edge does not exist, it is created with an initial coverage
        of 1. If the edge cannot be verified after creation, it is removed, and an error
        is raised.
        Args:
            seq1 (str): The sequence representing the first node.
            seq2 (str): The sequence representing the second node.
        Raises:
            ValueError: If the graph is compacted and edges cannot be added.
            ValueError: If the edge cannot be verified after being added.
        """
        # check if graph is compacted
        if self.graph['compacted']:
            raise ValueError("Cannot add edge to a compacted graph.")
        
        # check that sequences only contain 'A', 'T', 'C', 'G'
        valid_bases = {'A', 'T', 'C', 'G'}
        assert set(seq1) <= valid_bases, f"Invalid base in sequence 1: {seq1}"
        assert set(seq2) <= valid_bases, f"Invalid base in sequence 2: {seq2}"
        
        # find or add node for seq1
        node1 = self.find_or_add_node(seq1)
        # find or add node for seq2
        node2 = self.find_or_add_node(seq2)
        
        # check if edge already exists
        if self.has_successor(seq1, seq2):
            # increment edge coverage
            old_coverage = self.get_edge_attributes(seq1, seq2)['coverage']
            new_coverage = old_coverage + 1
            self.set_edge_attributes(seq1, seq2, coverage=new_coverage)
        else:
            # add new edge
            self.add_edge(seq1, seq2)
            self.set_edge_attributes(seq1, seq2, coverage=1, query_edge=False, query_edge_location=[])
            
            # check that new edge was added
            try: 
                self.get_edge_attributes(seq1, seq2)
            except KeyError:
                self.remove_edge(seq1, seq2)
                raise ValueError(f"Edge {seq1}->{seq2} was not added.")
            

# add k-mer to de bruijn graph
    def add_kmer(self, kmer):
        """
        Adds a k-mer and its reverse complement to the de Bruijn graph.
        This method takes a k-mer, validates its length against the graph's k value,
        and splits it into left and right (k-1)-mers. It then adds the nodes and 
        edges corresponding to the k-mer and its reverse complement to the graph.
        Args:
            kmer (str): The k-mer to be added to the graph.
        Raises:
            AssertionError: If the length of the k-mer does not match the graph's k value.
            ValueError: If there is an issue adding nodes or edges to the graph.
        """
        
        # check that k-mer length matches graph k value
        assert len(kmer) == self.graph['k'], "k-mer length does not match graph k value"
        
        # split k-mer into two sequences
        left_mer = kmer[:-1]
        right_mer = kmer[1:]
        
        # add nodes and edge
        try: 
            self.add_nodes_with_edge(left_mer, right_mer)
        except ValueError as e:
            print(e)
        
        # add reverse complement nodes and edge
        rc_kmer = self.reverse_complement(kmer)
        rc_left_mer = rc_kmer[:-1]
        rc_right_mer = rc_kmer[1:]
        try: 
            self.add_nodes_with_edge(rc_left_mer, rc_right_mer)
        except ValueError as e:
            print(e)
        
# add k-mers from fasta file to de bruijn graph
    def add_kmers_from_fasta(self, fasta_path, min_count=1):
        """
        Reads a FASTA file and adds k-mers from the sequences to the de Bruijn graph
        only if they appear at least `min_count` times.

        Args:
            fasta_path (str): The file path to the FASTA file containing sequences.
            min_count (int): The minimum count a k-mer must have to be added to the graph.

        The method skips lines starting with '>' (FASTA headers) and processes the 
        remaining lines as sequences. It counts k-mer occurrences and adds only those 
        k-mers that meet the `min_count` threshold. This is an initial filter to prevent 
        erronous k-mers from being added to the graph.
        """

        k = self.graph['k']
        kmer_counts = defaultdict(int)

        # First pass: Count k-mer occurrences
        with open(fasta_path, 'r', encoding='utf-8') as file:
            for line in file:
                if line.startswith('>'):  # skip headers
                    continue
                else:
                    sequence = line.strip()
                    # assert that sequence only contains 'A', 'T', 'C', 'G'
                    valid_bases = {'A', 'T', 'C', 'G'}
                    assert set(sequence) <= valid_bases, f"Invalid base in sequence: {sequence}"
                    for i in range(len(sequence) - k + 1):
                        kmer = sequence[i:i+k]
                        kmer_counts[kmer] += 1

        # Second pass: Add k-mers that meet the min_count threshold
        for kmer, count in kmer_counts.items():
            if count >= min_count:
                self.add_kmer(kmer)

# import query fasta file as search queue 
    def query_kmers_from_fasta(self, fasta_path):
        """
        Reads a FASTA file and returns an iterator of k-mers from the sequences.

        Args:
            fasta_path (str): The file path to the FASTA file containing sequences.

        The method skips lines starting with '>' (FASTA headers) and processes the 
        remaining lines as sequences. It yields left and right k-1mers from the 
        forward and reverse complement of the sequence as it reads them. 
        """
        k = self.graph['k']

        with open(fasta_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    continue
                else:
                    sequence = line.strip()
                    # assert that sequence only contains 'A', 'T', 'C', 'G'
                    valid_bases = {'A', 'T', 'C', 'G'}
                    assert set(sequence) <= valid_bases, f"Invalid base in sequence: {sequence}"
                    for i in range(len(sequence) - k + 1):
                        kmer = sequence[i:i+k]
                        left_mer = kmer[:-1]
                        right_mer = kmer[1:]
                        yield left_mer
                        yield right_mer
                        # add reverse complement k-mers 
                        rc_kmer = self.reverse_complement(kmer)
                        rc_left_mer = rc_kmer[:-1]
                        rc_right_mer = rc_kmer[1:]
                        yield rc_left_mer
                        yield rc_right_mer
                        
# identify query subgraph 
## filters de bruijn graph to only include k-mers connected to query k-mers
    def identify_query_connected_nodes(self, query_kmers):
        """
       Identifies k-mers in the De Bruijn Graph connected to query k-mers using a breadth-first search algorithm.

        This method determines all k-mers in the graph that are connected to the input query k-mers. It updates 
        the graph with attributes indicating whether a node is part of the query subgraph and its distance 
        from the query k-mers.

        Args:
            query_kmers (iterable): An iterable of k-mers to be used as the starting points for the search.

        Behavior:
            - Initializes the query k-mers with a distance of 0 and marks them as part of the query subgraph.
            - Performs a breadth-first search (BFS) to traverse the graph, starting from the query k-mers.
            - Updates the graph nodes with attributes:
            - `visited`: Indicates whether the node has been visited during BFS.
            - `query_subgraph`: Indicates whether the node is part of the query subgraph.
            - `query_distance`: The shortest distance from the node to any query k-mer.
            - Tracks distances of connected k-mers in a subgraph distance dictionary.

        Returns:
            dict: A dictionary where keys are k-mers and values are their shortest distances from the query k-mers.
        """
        # initialize subgraph nodes and distance dictionary 
        subgraph_dist = {}
        
        # check that query k-mers are in the graph
        for kmer in query_kmers:
            if self.has_node(kmer):
                # get node attributes from graph
                node_attr = self.get_node_attributes(kmer)
                # set query distance to 0 for query k-mers
                node_attr['query_distance'] = 0
                # update query attributes in graph
                self.set_node_attributes(kmer, **node_attr)
                # add kmer to subgraph_dist 
                subgraph_dist[kmer] = 0
            else:
                print(f"Query k-mer: {kmer} not found in graph.")
                # remove kmer from query_kmers
                query_kmers.remove(kmer)
        
        # if subgraph_dist is empty, raise an error message
        if not subgraph_dist:
            raise ValueError("No query k-mers found in the graph.")
              
        # initialize queue with query k-mers
        # breadth-first search algorithm
        queue = deque(query_kmers)
        # iterate over queue
        while queue:
            # pop the first k-mer (FIFO: first in, first out)
            kmer = queue.popleft()
            # get k-mer attributes 
            node_attr = self.get_node_attributes(kmer)
            # if k-mer is not visited, continue
            if node_attr['visited']:
                continue
            # set visited to True and query subgraph to True
            node_attr['visited'] = True
            node_attr['query_subgraph'] = True        
            # get distance and add to subgraph distance dictionary
            dist = node_attr['query_distance']
            subgraph_dist[kmer] = dist
            # update node attributes in graph
            self.set_node_attributes(kmer, **node_attr)
            # identify adjacent k-mers 
            successors = self.successors(kmer)
            predecessors = self.predecessors(kmer)
            
            # iterate over successors 
            for successor in successors:
                # add successor to queue
                if not self.get_node_attributes(successor)['visited'] and successor not in queue:
                    queue.append(successor)
                # set distance for successor
                if successor not in subgraph_dist:
                    subgraph_dist[successor] = dist + 1
                else : # if already in subgraph_dist, update distance to minimum dist
                    subgraph_dist[successor] = min(subgraph_dist[successor], dist + 1)
                # add successor dist to graph 
                node_attr = self.get_node_attributes(successor)
                node_attr['query_distance'] = subgraph_dist[successor]
                self.set_node_attributes(successor, **node_attr)
                
            # iterate over predecessors
            for predecessor in predecessors:
                # add predecessor to queue if not already visited 
                if not self.get_node_attributes(predecessor)['visited'] and predecessor not in queue:
                    queue.append(predecessor)
                # set distance for predecessor
                if predecessor not in subgraph_dist:
                    subgraph_dist[predecessor] = dist + 1
                else :  #if already in subgraph_dist, update distance to minimum dist
                    subgraph_dist[predecessor] = min(subgraph_dist[predecessor], dist + 1)
                # add predecessor dist to graph 
                node_attr = self.get_node_attributes(predecessor)
                node_attr['query_distance'] = subgraph_dist[predecessor]
                self.set_node_attributes(predecessor, **node_attr)
        
        return(subgraph_dist)
        
# filter graph based on dictionary of nodes 
    def filter_graph_to_subgraph(self, query_kmers):
                """
                Filters the current graph to include only nodes and edges that are part of the query subgraph.

                This method takes a dictionary of query-related k-mers and constructs a subgraph containing only 
                the nodes and edges that are relevant to the query. The subgraph is created efficiently by 
                leveraging a set for membership checks and copying relevant nodes and edges in a single pass.

                Args:
                    query_kmers (dict): A dictionary where keys are k-mers that define the query subgraph, 
                            and values are distances (ignored in this function). Generated with 
                            identify_query_connected_nodes.

                Behavior:
                    - Nodes in the graph that are not part of the query_kmers are excluded.
                    - Edges between nodes are preserved only if both nodes are part of the query_kmers.
                    - The current graph is replaced with the filtered subgraph.

                Notes:
                    - This method assumes that the graph is represented using internal attributes such as 
                      `_node`, `_pred`, and `_succ` for nodes, predecessors, and successors respectively.
                    - The method also updates the graph's metadata, including node and edge counts.

                Raises:
                    KeyError: If any k-mer in query_kmers is not found in the graph.

                """
                # Create a new graph to store the subgraph
                subgraph = diGraph(k=self.graph['k'], compacted=self.graph['compacted'])

                # Use set for faster membership checks
                query_kmers_set = set(query_kmers.keys())

                # Copy nodes and edges to subgraph in a single pass
                for kmer in query_kmers_set:
                    if kmer in self._node:
                        subgraph.add_node(kmer, **self.get_node_attributes(kmer))
                    # Add edges to predecessors
                    for pred in self._pred.get(kmer, {}):
                        if pred in query_kmers_set:
                            subgraph.add_edge(pred, kmer, **self.get_edge_attributes(pred, kmer))
                    # Add edges to successors
                    for succ in self._succ.get(kmer, {}):
                        if succ in query_kmers_set:
                            subgraph.add_edge(kmer, succ, **self.get_edge_attributes(kmer, succ))

                # Replace the current graph with the filtered subgraph
                self.graph = subgraph.graph
                self._node = subgraph._node
                self._succ = subgraph._succ
                self._pred = subgraph._pred
                self.node_count = subgraph.node_count
                self.edge_count = subgraph.edge_count

# identify query connected nodes and generate subgraph 
## combines the above two methods for quicker subgraph generation
    def identify_query_subgraph(self, query_kmers):
        """
        Identifies and extracts a subgraph from the current graph based on the provided query k-mers.
        This method performs a breadth-first search (BFS) starting from the query k-mers, 
        traversing both successors and predecessors, and constructs a subgraph containing 
        all reachable nodes and edges. The subgraph is then used to replace the current graph.
        Args:
            query_kmers (list): A list of k-mers to use as the starting points for identifying the subgraph.
        Raises:
            ValueError: If none of the query k-mers are found in the graph.
        Behavior:
            - Initializes the subgraph and a distance dictionary for tracking distances from query k-mers.
            - Checks if each query k-mer exists in the graph. If not, it is removed from the query list.
            - Performs a BFS to traverse the graph, marking nodes as visited and adding them to the subgraph.
            - Updates node attributes such as `query_distance`, `visited`, and `query_subgraph`.
            - Adds edges between nodes in the subgraph.
            - Replaces the current graph with the newly constructed subgraph.
        Notes:
            - The method assumes that the graph has methods for accessing and modifying node and edge attributes.
            - The graph is expected to have attributes such as `visited` and `query_distance` for nodes.
        Returns:
            None
        """
        
        # initialize subgraph nodes and distance dictionary 
        subgraph_dist = {}
        # Create a new graph to store the subgraph
        subgraph = diGraph(k=self.graph['k'], compacted=self.graph['compacted'])
        
        # Check that query k-mers are in the graph
        for kmer in query_kmers:
            if self.has_node(kmer):
                # get node attributes from graph
                node_attr = self.get_node_attributes(kmer)
                # set query distance to 0 for query k-mers
                node_attr['query_distance'] = 0
                # update query attributes in graph
                self.set_node_attributes(kmer, **node_attr)
                # add kmer to subgraph_dist 
                subgraph_dist[kmer] = 0
            else:
                print(f"Query k-mer: {kmer} not found in graph.")
                # remove kmer from query_kmers
                query_kmers.remove(kmer)
        
        # if subgraph_dist is empty, raise an error message
        if not subgraph_dist:
            raise ValueError("No query k-mers found in the graph.")
        
        # initialize queue with query k-mers
        # breadth-first search algorithm
        queue = deque(query_kmers)
        # iterate over queue
        while queue:
            # pop
            kmer = queue.popleft()
            # get k-mer attributes
            node_attr = self.get_node_attributes(kmer)
            # if k-mer is not visited, continue
            if node_attr['visited']:
                continue
            # set visited to True and query subgraph to True
            node_attr['visited'] = True
            node_attr['query_subgraph'] = True
            # get distance and add to subgraph distance dictionary
            dist = node_attr['query_distance']
            subgraph_dist[kmer] = dist
            # update node attributes in graph
            self.set_node_attributes(kmer, **node_attr)
            # add node to subgraph
            subgraph.add_node(kmer, **node_attr)
            # identify adjacent k-mers
            successors = self.successors(kmer)
            predecessors = self.predecessors(kmer)
            # iterate over successors
            for successor in successors:
                # add successor to queue
                if not self.get_node_attributes(successor)['visited'] and successor not in queue:
                    queue.append(successor)
                # set distance for successor
                if successor not in subgraph_dist:
                    subgraph_dist[successor] = dist + 1
                else:
                    subgraph_dist[successor] = min(subgraph_dist[successor], dist + 1)
                # add successor dist to graph
                node_attr = self.get_node_attributes(successor)
                node_attr['query_distance'] = subgraph_dist[successor]
                self.set_node_attributes(successor, **node_attr)
                # add node to subgraph
                subgraph.add_node(successor, **node_attr)
                # add edge to subgraph
                subgraph.add_edge(kmer, successor, **self.get_edge_attributes(kmer, successor))
            # iterate over predecessors
            for predecessor in predecessors:
                # add predecessor to queue if not already visited
                if not self.get_node_attributes(predecessor)['visited'] and predecessor not in queue:
                    queue.append(predecessor)
                # set distance for predecessor
                if predecessor not in subgraph_dist:
                    subgraph_dist[predecessor] = dist + 1
                else:
                    subgraph_dist[predecessor] = min(subgraph_dist[predecessor], dist + 1)
                # add predecessor dist to graph
                node_attr = self.get_node_attributes(predecessor)
                node_attr['query_distance'] = subgraph_dist[predecessor]
                self.set_node_attributes(predecessor, **node_attr)
                # add node to subgraph
                subgraph.add_node(predecessor, **node_attr)
                # add edge to subgraph
                subgraph.add_edge(predecessor, kmer, **self.get_edge_attributes(predecessor, kmer))
            
        # replace the current graph with the filtered subgraph
        self.graph = subgraph.graph
        self._node = subgraph._node
        self._succ = subgraph._succ
        self._pred = subgraph._pred
        self.node_count = subgraph.node_count
        self.edge_count = subgraph.edge_count
        
    def add_degrees_attributes(self):
        """
        Adds the out-degree and in-degree of each node as attributes in the graph.
        The out-degree of a node is the number of outgoing edges from the node, while
        the in-degree is the number of incoming edges to the node. The method iterates
        over each node in the graph, calculates the out-degree and in-degree, and sets
        these values as attributes of the node.
        """
        for node in self.get_nodes():
            out_degree = self.out_degree(node)
            in_degree = self.in_degree(node)
            self.set_node_attributes(node, out_degree=out_degree, in_degree=in_degree)     

    def update_all_node_attributes(self, **attr):
        """
        Updates the attributes of all nodes in the graph.
        This method updates the attributes of all nodes in the graph with the given
        key-value pairs. It iterates over each node in the graph and updates the
        node's attributes with the provided values.
        Args:
            **attr: Key-value pairs of attributes to update.
        """
        for node in self.get_nodes():
            self.set_node_attributes(node, **attr)
                
    def identify_query_boundaries(self, fasta_path): 
        """
        Identifies query boundaries in the de Bruijn graph based on a single sequence from a FASTA file.
        This method processes a query sequence from a provided FASTA file, validates it, and updates 
        the graph with relevant attributes such as query edges, query start, and query end. It also 
        marks nodes as visited if they are part of the query sequence.
        Args: 
            fasta_path (str): The file path to the FASTA file containing the query sequence. The file 
                              is expected to contain only one sequence.
        Raises:
            AssertionError: If the FASTA file does not contain a sequence or if the sequence contains 
                            invalid characters (not 'A', 'T', 'C', 'G').
            ValueError: If an edge corresponding to a k-mer in the query sequence does not exist in 
                        the graph, indicating that the query sequence is not fully represented.
        Graph Updates:
            - Adds the query sequence as a graph attribute.
            - Marks edges corresponding to k-mers in the query sequence with the "query_edge" attribute 
              and records their positions in the "query_edge_location" attribute.
            - Identifies and sets the start and end nodes of the query sequence in the graph.
            - Marks nodes as visited if they are part of the query sequence.
        """
        k = self.graph['k']
        
        # import fasta sequence from fasta_path
        # assumes there's only one fasta sequence query
        with open(fasta_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    continue
                else:
                    query_sequence = line.strip()
        # if query_sequence is empty error
        assert query_sequence, f"No query sequence in fasta"
        # add query_sequence to graph attributes
        self.set_graph_attributes(query = query_sequence) 
        # assert that sequence only contains 'A', 'T', 'C', 'G'
        valid_bases = {'A', 'T', 'C', 'G'}
        query_sequence = query_sequence.upper()
        assert set(query_sequence) <= valid_bases, f"Invalid base in sequence: {query_sequence}"
        
        for i in range(len(query_sequence) - k + 1):
            kmer = query_sequence[i:i+k]
            left_mer = kmer[:-1]
            right_mer = kmer[1:]
            # check if left_mer-> right_mer edge exists
            if self.has_successor(left_mer, right_mer):
                # get edge attributes
                edge_attr = self.get_edge_attributes(left_mer,right_mer)
                # update edge attributes 
                edge_attr["query_edge"] = True
                # append query_edge_location as i 
                edge_attr["query_edge_location"].append(i)
                # update edge attributes
                self.set_edge_attributes(left_mer, right_mer, **edge_attr) 
                # determine if query_start
                if i == 0: 
                    # update left_mer as query_start
                    self.set_graph_attributes(query_start = left_mer)
                # determine if query_end
                if i == (len(query_sequence) - k):
                    # update right_mer as query_end
                    self.set_graph_attributes(query_end = right_mer)
                        
                # set left and right mers to visited 
                for mer in (left_mer, right_mer): 
                    # get mer attr
                    mer_attr = self.get_node_attributes(mer)
                    # set mer to visited if not already visited
                    if not mer_attr.get("visited", False):
                        mer_attr["visited"] = True
                        # update in graph
                        self.set_node_attributes(mer, **mer_attr)
            else:
                # left_mer-> right_mer edge doesn't exist in graph
                raise ValueError(
                    f"Edge {left_mer}->{right_mer} does not exist in the graph. "
                    "Query sequence may not be fully represented."
                )
                
    # identify next node to select 
    def get_next_node(self, node, dir):
        """
        Determines the next node to visit in a graph traversal based on the specified direction.
        Args:
            node (str): The current node in the graph.
            dir (int): The direction of traversal. Must be 0 (predecessor) or 1 (successor).
        Returns:
            str or None: The next node to visit, or None if no valid next node exists.
        Raises:
            ValueError: If the direction (`dir`) is not 0 or 1.
        Behavior:
            - If the degree of the current node in the specified direction is less than 1, returns None.
            - If the degree is 1, returns the single neighbor unless it is self-referential or leads to a dead-end.
            - If the degree is greater than 1:
                - Prioritizes unvisited neighbors.
                - Among unvisited neighbors, selects the one with the highest degree and edge coverage, breaking ties alphabetically.
                - If no unvisited neighbors exist, considers neighbors with unvisited neighbors, applying the same prioritization.
            - Returns None if all neighbors and their neighbors have been visited.
            - Due to this prioritization, cannot handle loops if all neighbor neighbors are visited. 
        """
      
        # Check that direction is either 0 or 1
        if dir not in (0, 1):
            raise ValueError(f"Invalid direction: {dir}. Expected 0 (predecessor) or 1 (successor).")
        
        # Check that node is set to visisted
        # Node should be visisted if used in the walk function
        # Mark the current node as visited
        node_attr = self.get_node_attributes(node)
        if not node_attr.get("visited", False):
            node_attr["visited"] = True
            self.set_node_attributes(node, **node_attr)
        
        # Define direction-specific variables
        degree_key = "in_degree" if dir == 0 else "out_degree"
        neighbor_func = self.predecessors if dir == 0 else self.successors
        edge_func = (lambda n1, n2: self.get_edge_attributes(n2, n1)) if dir == 0 else self.get_edge_attributes

        # Get the degree of the current node
        degree = self.get_node_attributes(node)[degree_key]

        if degree < 1:
            # Cannot extend path any more
            return None

        if degree == 1:
            # Get the single neighbor
            next_node = list(neighbor_func(node))[0]
            if node == next_node:
                # Skip self-referential node
                return None
            next_node_visited = self.get_node_attributes(next_node)["visited"]
            if next_node_visited:
                # Check if any neighbors of next_node are unvisited
                if any(
                    not self.get_node_attributes(neighbor)["visited"]
                    for neighbor in neighbor_func(next_node)
                    ):
                    # some neighbors not visited, not a dead-end
                    return next_node
                else:
                    # all neighbors visited, dead-end
                    return None
            else:
                # next_node is neither visited or current node
                return next_node

        if degree > 1:
            # Get neighbors
            neighbors = neighbor_func(node)
            # Filter unvisited neighbors
            unvisited_neighbors = [n for n in neighbors if not self.get_node_attributes(n)["visited"]]

            if unvisited_neighbors:
                # Select the next node based on highest coverage and degrees, breaking ties alphabetically
                next_node = sorted(
                    unvisited_neighbors,
                    key=lambda n: (
                        -self.get_node_attributes(n)[degree_key],
                        -edge_func(node, n)["coverage"],
                        n  # Alphabetical order
                        )
                    )[0]
                if next_node == node:
                    # no self referential nodes
                    # this should be avoided by setting node to visited prior to search
                    return None
                else:
                    return next_node
            else:
            # Check if any neighbors have unvisited neighbors
                neighbors_with_unvisited_neighbors = [
                    n for n in neighbors
                    if any(
                    not self.get_node_attributes(neighbor_neighbor)["visited"]
                    for neighbor_neighbor in neighbor_func(n)
                    )
                ]
                if neighbors_with_unvisited_neighbors:
                   # Select the next node based on highest coverage and degrees, breaking ties alphabetically
                    next_node = sorted(
                        neighbors_with_unvisited_neighbors,
                        key=lambda n: (
                            -self.get_node_attributes(n)[degree_key],
                            -edge_func(node, n)["coverage"],
                            n  # Alphabetical order
                            )
                        )[0]
                    if next_node == node:
                        # no self referential nodes
                        # this should be avoided by setting node to visited prior to search
                        return None
                    else: 
                        return next_node
                else:
                    # All neighbors and their neighbors have been visited
                    return None

    def walk(self, start_node, dir): 
        # initiate sequence string
        sequence = ""    
        # check that start_node is visited, raise error if not
            # all query nodes should be visited
        if not self.get_node_attributes(start_node).get("visited", False):
            raise ValueError(f"Start node {start_node} must be marked as visited before calling the walk method.")
        
        # set node to start_node
        node = start_node
        
        # while get next node doesn't return empty
        while True:
            # get the next node
            next_node = self.get_next_node(node, dir)
            if not next_node:
                break
            # update node
            node = next_node
            # get node sequence 
            full_sequence = self.get_node_attributes(node)["sequence"]
            # set node attribute to visited
            self.set_node_attributes(node, visited=True)
            # get node sequence 
            full_sequence = self.get_node_attributes(node)["sequence"]
            # trim sequence based on dir
            if (dir == 0): 
                # get first character of full_sequence
                node_sequence = full_sequence[0]
                # append characters to the beginning of sequence to construct full sequence
                sequence = node_sequence + sequence
            if (dir == 1): 
                # trim sequence to the last character for successor direction
                node_sequence = full_sequence[-1]
                # Append characters from each node's sequence to construct the full sequence
                sequence += node_sequence
        
        # return sequence 
        return sequence
    
    def assemble_from_query(self):
        # intiate empty assembled sequence
        assembled_sequence = ""
        # check that query is available in graph attr
        if not self.graph.get("query"):
            raise ValueError("Query sequence not found in graph attributes. Please import fasta into graph with identify_query_boundaries.")
        
        # get query start and end
        query_start = self.graph['query_start']
        query_end = self.graph['query_end']
        
        # run walk in left direction, starting from query_start
        left_contig = self.walk(query_start, 0)    
        # run walk in right direction, starting from query_end
        right_contig = self.walk(query_end, 1)
        
        # combine into final contig: 
        assembled_sequence += left_contig
        assembled_sequence += self.graph.get("query")
        assembled_sequence += right_contig
        
        return full_sequence
                   
        
        
        
            
            
        
    
    
                    
                
                    
                    
                
                
        
                        
            
  

            
            
                
            



# remove tips from de bruijn graph

# remove isolated nodes from de bruijn graph

# remove low coverage k-mers from de bruijn graph

# remove low coverage edges from de bruijn graph
# bubble removal from de bruijn graph

# eulerian path from de bruijn graph
