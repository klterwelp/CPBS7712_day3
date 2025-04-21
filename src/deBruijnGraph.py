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
                'query_kmer': False
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
    def add_nodes_with_edge(self, seq1, seq2, coverage=1):
        """
        Adds nodes and an edge between them in the de Bruijn graph.
        This method ensures that the graph is not compacted before adding nodes or edges.
        It either finds or creates nodes for the given sequences `seq1` and `seq2`, and
        then adds an edge between them. If the edge already exists, its coverage attribute
        is incremented. If the edge does not exist, it is created with an initial coverage
        of 1. If the edge cannot be verified after creation, it is removed, and an error
        is raised. Can optionally change coverage if using a k-mer count dict. 
        Args:
            seq1 (str): The sequence representing the first node.
            seq2 (str): The sequence representing the second node.
            coverage (int): The number of times an edge was seen. (default is 1)
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
            # only adds one because it assumes that a k-mer dict uses k-mers as keys
            old_coverage = self.get_edge_attributes(seq1, seq2)['coverage']
            new_coverage = old_coverage + 1
            self.set_edge_attributes(seq1, seq2, coverage=new_coverage)
        else:
            # add new edge
            self.add_edge(seq1, seq2)
            self.set_edge_attributes(seq1, seq2, coverage=coverage, query_edge=False, query_edge_location=[])
            
            # check that new edge was added
            try: 
                self.get_edge_attributes(seq1, seq2)
            except KeyError:
                self.remove_edge(seq1, seq2)
                raise ValueError(f"Edge {seq1}->{seq2} was not added.")
            

# add k-mer to de bruijn graph
    def add_kmer(self, kmer, coverage=1):
        """
        Adds a k-mer and its reverse complement to the de Bruijn graph.
        This method takes a k-mer, validates its length against the graph's k value,
        and splits it into left and right (k-1)-mers. It then adds the nodes and 
        edges corresponding to the k-mer and its reverse complement to the graph. If
        adding k-mers from a k-mer count dict, can add the count as coverage to ensure
        coverage passes over.
        Args:
            kmer (str): The k-mer to be added to the graph.
            coverage (int): The amount of coverage for the k-mer edge. (default 1)
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
            self.add_nodes_with_edge(left_mer, right_mer, coverage)
        except ValueError as e:
            print(e)
        
        # add reverse complement nodes and edge
        rc_kmer = self.reverse_complement(kmer)
        rc_left_mer = rc_kmer[:-1]
        rc_right_mer = rc_kmer[1:]
        try: 
            self.add_nodes_with_edge(rc_left_mer, rc_right_mer, coverage)
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
        erronous k-mers from being added to the graph. Then, adds the count as coverage
        for the k-mer edge. 
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
                self.add_kmer(kmer, coverage=count)
        
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
        
        # start empty query_kmer graph attribute
        query_kmers = set()
        for i in range(len(query_sequence) - k + 1):
            kmer = query_sequence[i:i+k]
            left_mer = kmer[:-1]
            right_mer = kmer[1:]
            # check if left_mer-> right_mer edge exists
            if self.has_successor(left_mer, right_mer):
                # add left_mer and right_mer to query_kmers
                query_kmers.add(left_mer)
                query_kmers.add(right_mer)
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
                        mer_attr["query_kmer"] = True
                        # update in graph
                        self.set_node_attributes(mer, **mer_attr)
            else:
                # left_mer-> right_mer edge doesn't exist in graph
                raise ValueError(
                    f"Edge {left_mer}->{right_mer} does not exist in the graph. "
                    "Query sequence may not be fully represented."
                )
        # set query_kmers as graph attribute
        self.set_graph_attributes(query_kmers = query_kmers)
    
    # identify nodes connected to a node using bfs and direction 
    def find_connected_nodes(self, node, visited_nodes=None, connected_nodes=None, dir=0):
        """
        Finds all nodes connected to a given node in a directed graph.
        This method performs a breadth-first search (BFS) to identify all nodes
        connected to the input node, either in the direction of predecessors or
        successors, depending on the specified direction. If `visited_nodes` or
        `connected_nodes` is `None`, they are initialized as empty sets.
        Args:
            node (hashable): The starting node for the search.
            visited_nodes (set, optional): A set of nodes that have already been visited.
                If None, an empty set will be initialized. Defaults to None.
            connected_nodes (set, optional): A set of nodes that are connected to the
                starting node. If None, an empty set will be initialized. Defaults to None.
            dir (int, optional): The direction of the search. Use 0 for predecessors
                and 1 for successors. Defaults to 0.
        Returns:
            connected_nodes (set): A set of all nodes connected to the input node.
        Raises:
            ValueError: If the input node is not found in the graph.
            ValueError: If the direction is not 0 or 1.
        """
        # check if node is in graph
        if not self.has_node(node):
            return set()  # Return an empty set if the node is not found
        # check that direction is a valid option 
        if dir != 0 and dir != 1:
            raise ValueError("Invalid direction: must be 0 (predecessors) or 1 (successors).")
        # determine neighbor function
        neighbor_func = self.predecessors if dir == 0 else self.successors
        # start empty visited nodes set if not set
        if visited_nodes is None:
            visited_nodes = set()
        # start empty connected nodes set if not set
        if connected_nodes is None:
          connected_nodes = set()
        # initiate queue with input node
        queue = deque([node])
        
        # iterate over queue while queue is not empty 
        while queue:
            # pop the first node 
            current_node = queue.popleft()
            # check if current node is already visited
            if current_node in visited_nodes:
                continue    # skips if already visited
            else: 
                # set current node to visited 
                visited_nodes.add(current_node)
            
            if current_node in connected_nodes:
                continue    # skips if already in connected nodes
            else: 
                # add current node to connected nodes
                connected_nodes.add(current_node)
            
            # identify node neighbors
            neighbors = list(neighbor_func(current_node))
            # check if neighbors is empty
            if not neighbors:
                continue
            # iterate over neighbors
            for neighbor in neighbors:
                # check if neighbor is already in connected nodes, visited, or in in queue
                if neighbor not in connected_nodes:
                    if neighbor not in visited_nodes:
                        if neighbor not in queue:
                            queue.append(neighbor)
                else: 
                    continue
        
        return connected_nodes

    # filter out all nodes that are not connected to query k-mers
    def filter_graph(self): 
        
        # identify query start, end, and k-mers
        query_start = self.graph.get("query_start")
        query_end = self.graph.get("query_end")
        query_kmers = self.graph.get("query_kmers")
        
        # check if query start and end are defined
        if not query_start or not query_end:
            raise ValueError("Query start and end not defined. Please run identify_query_boundaries.")
        # check if query k-mers are defined
        if not query_kmers:
            raise ValueError("Query k-mers not defined. Please run identify_query_boundaries.")
        # check if query start and end are in the graph
        if not self.has_node(query_start) or not self.has_node(query_end):
            raise ValueError("Query start and end not found in graph. Please run identify_query_boundaries.")
        # check if graph is already compacted
        if self.graph["compacted"]:
            raise ValueError("Graph is already compacted. Only compact once.")
        
        # identify all nodes connected to query start from the predecessor direction
        connected_start = self.find_connected_nodes(query_start, dir = 0)
        connected_end = self.find_connected_nodes(query_end, dir = 1)
        # combine connected start and end nodes into connected kmers, avoiding duplicates
        connected_kmers = connected_start.union(connected_end)
        # Check that connected start and end are not empty
        if not connected_kmers:
            raise ValueError("No nodes connected to query start or end.")
        # Add all query k-mers to connected_kmers, avoiding duplicates
        connected_kmers.union(query_kmers)
        
        # Remove all nodes that are not in the connected_kmers set
        for node in list(self.get_nodes()):
            if node not in connected_kmers:
                self.remove_node(node)
        
    # compact graph except for query kmers 
    def compact_except_query(self):
        """
        Compacts all non-query k-mers in the graph while preserving query k-mers.

        This method performs a graph compaction process that merges linear paths of k-mers, 
        reducing the number of nodes while maintaining essential connectivity. The method 
        ensures that query k-mers remain unmodified. The process follows these steps:

        1. **Pre-checks**:
            - Ensures the graph is not already compacted.
            - Adds in-degree and out-degree attributes if missing.
            - Validates that query boundaries are defined.

        2. **Compaction Process**:
            - Iteratively identifies compactable nodes (non-query, out-degree = 1).
            - Finds the left-most compactable node in a linear path.
            - Merges consecutive nodes by appending their sequences while preserving essential 
          properties such as distance and coverage.
            - Removes redundant nodes and updates edges accordingly.
    
        3. **Final Updates**:
            - Resets `visited` attributes for all non-query nodes.
            - Updates coverage information for compacted nodes based on `sum/dist`.

        Raises:
            ValueError: If the graph is already compacted or query boundaries are not defined.

        Notes:
            - The function prints intermediate steps for debugging.
            - Query k-mers remain unmodified throughout the process.
        """
        # check if graph is already compacted, raise errors if it is 
        if self.graph["compacted"]:
            raise ValueError("Graph is already compacted. Only compact once.")
        
        # check if nodes have in-degree and out-degree attributes
        if 'in_degree' not in self.get_node_attributes(next(self.get_nodes())):
            self.add_degrees_attributes() # add in-degree and out-degree attributes
        
        # check that query boundaries were established in the graph 
        if not self.graph.get("query"):
            raise ValueError("Query sequence not found in graph attributes. Please import fasta into graph with identify_query_boundaries.")
        
        # get k for graph
        k = self.graph['k']
        
        # iteration counter
        i = 0
        
        while True:
            # update set of compactable nodes
                # do not compact: query_kmers or already compacted nodes
            compactable_nodes = list(self.get_nodes(compacted=False, out_degree=1, query_kmer=False))            
            # if set of compactable nodes is empty, break
            if not compactable_nodes:
                break
           
            # get first node of compactable set 
            node = compactable_nodes.pop(0)
            i += 1 # update iteration counter
            # track whether node visited for left-most node search
            left_walk_visited=set()
            # add node to left_walk_visited
            left_walk_visited.add(node)
            # determine in degree of node and out degree of predecessor
            node_in_degree = self.get_node_attributes(node)['in_degree']
           
            if (node_in_degree == 1):
                pred = list(self.predecessors(node))[0]# get first predecessor
                pred_attr = self.get_node_attributes(pred)
                pred_out_degree = pred_attr['out_degree']
                pred_query = pred_attr['query_kmer']
                pred_visited = pred in left_walk_visited
            else:
                # no predecessor
                pred = 0
                pred_out_degree = 0
                pred_visited = False
                pred_query = False
        
            # walk to left-most compactable node from initial set node
            while (node_in_degree == 1) and (pred_out_degree == 1) and (not pred_visited) and (not pred_query):
                # update node to predecessor, generate new variables with new node
                node = pred
                left_walk_visited.add(node) # add node to visited
                # determine in degree of node and out degree of predecessor
                node_in_degree = self.get_node_attributes(node)['in_degree']
                if (node_in_degree > 0):
                    pred = list(self.predecessors(node))[0] # get first predecessor
                    pred_attr = self.get_node_attributes(pred) # get pred attributes
                    pred_out_degree = pred_attr['out_degree']
                    pred_query = pred_attr['query_kmer']
                    pred_visited = pred in left_walk_visited
                else:
                    # no predecessor
                    pred = 0
                    pred_out_degree = 0
                    pred_visited = False
                    pred_query = False 
        
            # start compacting from left-most node
            node_attr = self.get_node_attributes(node)
            node_out_degree = (self.out_degree(node) == 1)
            node_compacted = node_attr['compacted']
            node_query = node_attr['query_kmer']
        
            # double-check node is compactable
            assert (self.out_degree(node) == 1) and (not node_compacted) and (not node_query), "Error, we hit a non-compactable node."

            first_sequence = True
            # set node attributes to ensure we don't re-compact it
            self.set_node_attributes(node,
                compacted=True,
                sum=0,
                dist=0
                )
        
             # while loop for compaction
            while True:
                # test the following break conditions:
                    # test the following break conditions:
                successors = list(self.successors(node))
                if not successors:
                    #print(f"No successors found for node: {node}")
                    break
                
                # 1) node must have one out_degree
                if (self.out_degree(node) != 1):
                    # more than one successor, no more compaction!
                    #print(f"finished compaction with node {node}")
                    break
                # 2) successor must have one in_degree
                succ = successors[0]
                succ_attr = self.get_node_attributes(succ)
                if (self.in_degree(succ) != 1):
                    #print(f"finished compaction with node {node}")
                    break
                # 3) successor cannot be the node (self-reference loop)
                if succ == node:
                    #print(f"hit a self-reference with node {node}")
                    break
                # 4) successor cannot already be visited
                if succ_attr['visited']:
                    #print(f"hit a loop: {succ} already visited")
                    break
                # 5) successor cannot be query node
                if succ_attr['query_kmer']:
                    #print(f"hit a query kmer: {succ}")
                    break
                # 6) node cannot be a successor or successor's successor
                if node in self.successors(node) or node in self.successors(succ):
                    #print(f"hit a loop: {node} is a successor or a successor's successor")
                    break 
                # set node to visited 
                self.set_node_attributes(node, visited=True)
                # get update node_attributes
                node_attr=self.get_node_attributes(node)
                node_in_degree=node_attr["in_degree"]
                # only keep the full k-mer sequence for first node in a path
                if node_in_degree > 0 and first_sequence:
                    node_sequence = self.get_node_attributes(node)['sequence']
                    trimmed_sequence = node_sequence[k-2:] # trim sequence
                    self.set_node_attributes(node, sequence = trimmed_sequence) # update node sequence
                    first_sequence = False
                else:
                    first_sequence = False  # keeps in_degree 0 node sequence
            
                # add all edges succ.successors to node
                for succ_succ in self.successors(succ):
                    self.add_edge(node, succ_succ)
                    # copy original edge coverage to new edge
                    old_coverage=self.get_edge_attributes(succ, succ_succ)['coverage']                    
                    self.set_edge_attributes(node, succ_succ, coverage=old_coverage)

                # add successor sequence to node sequence
                node_attr = self.get_node_attributes(node) 
            
                if succ_attr['compacted']:
                    # add succ sequence to node
                    node_sequence = node_attr['sequence']
                    succ_sequence = succ_attr['sequence']
                    node_sequence += succ_sequence
                    # add succ distance to node distance
                    node_dist = node_attr['dist']
                    succ_dist = succ_attr['dist']
                    node_dist += succ_dist
                    # add succ sum to node sum
                    node_sum = node_attr['sum']
                    succ_sum = succ_attr['sum']
                    node_sum += succ_sum
                    # update node attributes
                    self.set_node_attributes(node, sequence=node_sequence, 
                                        dist=node_dist,
                                        sum=node_sum
                                        )
                else: 
                    # successor not already compacted, only add last character
                    node_attr["sequence"] += succ_attr["sequence"][-1:]
                    # add edge coverage as the sum (which is the sum of all edges' coverage compacted)
                    node_attr["sum"] += self.get_edge_attributes(node, succ)["coverage"]
                    # moved one further (node->succ)
                    node_attr["dist"] += 1

                # now delete successor and edge between node->succ
                self.remove_edge(node, succ)
                self.remove_node(succ)
                # update node in_degrees, succ_succs should remain the same because edges preserved
                new_out_degree = self.out_degree(node)
                new_in_degree = self.in_degree(node)
                self.set_node_attributes(node, out_degree=new_out_degree, in_degree=new_in_degree)
                        
        # Update coverage of nodes based on sum/dist 
        # Update all nodes except query as not visited for assembly purposes
        reset_list=list(self.get_nodes(query_kmer=False, visited=True))
        for node in reset_list:
            self.set_node_attributes(node, visited=False)
            
        # Update compacted nodes with coverage info
        compacted_list=list(self.get_nodes(compacted=True))
        for node in compacted_list:
            node_attr=self.get_node_attributes(node)
            if (node_attr["dist"] == 0):
                # not actually compacted, just marked
                node_attr["coverage"]=1
                node_attr["compacted"]=False
                self.set_node_attributes(node, **node_attr)
            else:
                # actually compacted, calc coverage
                node_attr["coverage"]=(node_attr["sum"]/node_attr["dist"])
                self.set_node_attributes(node, **node_attr)         
               
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
                        -self.get_node_attributes(n)["coverage"],
                        -edge_func(node, n)["coverage"],
                        -self.get_node_attributes(n)["dist"],
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
                            -self.get_node_attributes(n)["coverage"],
                            -edge_func(node, n)["coverage"],
                            -self.get_node_attributes(n)["dist"],
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
        """
        Walks through the graph starting from the given node and constructs a sequence based on the direction.

        Parameters: 
        - start_node: The node from which to start the walk. Must be marked as visited.
        - dir: The direction of the walk. 0 for predecessor direction, 1 for successor direction.

        Returns:
        - sequence: The constructed sequence from the walk.

        Raises:
        - ValueError: If the start_node is not marked as visited.

        The method performs the following steps:
            1. Initiates an empty sequence string.
            2. Checks if the start_node is visited; raises an error if not.
            3. Sets the current node to start_node.
            4. Iteratively gets the next node in the specified direction until no more nodes are found.
            5. Updates the current node and retrieves its attributes.
            6. Marks the node as visited.
            7. Adds the node's sequence to the constructed sequence based on the direction and whether the node is compacted.
            8. Returns the final constructed sequence.
        """

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
            node_attr = self.get_node_attributes(node)
            full_sequence = node_attr["sequence"]
            dist = node_attr["dist"]
            # set node attribute to visited
            self.set_node_attributes(node, visited=True)
            # add sequence based on dir and compaction
            if node_attr["compacted"]:
                # add full sequence since it is already trimmed in compaction
                if (dir == 0):
                    sequence = full_sequence + sequence
                if (dir == 1):
                    sequence += full_sequence
            else: 
                # node not compacted, need trimming 
                if (dir == 0): 
                    # get first character of sequence for predeccessor direction
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
        """
        Assembles a sequence from the query sequence in the graph by walking in both directions from the query boundaries.
        Returns:
            - assembled_sequence: The final assembled sequence.

        Raises:
            - ValueError: If the query sequence is not found in the graph attributes.

        The method performs the following steps:
            1. Initiates an empty assembled sequence string.
            2. Checks if the query sequence is available in the graph attributes; raises an error if not.
            3. Retrieves the start and end nodes of the query sequence.
            4. Walks in the left direction starting from the query start node to assemble the left contig.
            5. If the node next to the query start is compacted, trims the left contig to remove the start of the query sequence node.
            6. Walks in the right direction starting from the query end node to assemble the right contig.
            7. Combines the left contig, query sequence, and right contig into the final assembled sequence.
            8. Returns the final assembled sequence.
        """
        # intiate empty assembled sequence
        assembled_sequence = ""
        # check that query is available in graph attr
        if not self.graph.get("query"):
            raise ValueError("Query sequence not found in graph attributes. Please import fasta into graph with identify_query_boundaries.")
        
        # get query start and end
        query_start = self.graph['query_start']
        query_end = self.graph['query_end']
        
        # run walk in left direction, starting from query_start
        # prior to starting walk, get the node next to query_start
        next_to_query_start = self.get_next_node(query_start, 0)
        # set next node back to unvisited
        self.set_node_attributes(next_to_query_start, visited=False)
        # assemble left contig
        left_contig = self.walk(query_start, 0)   
        # left contig will be trimmed if next node from query_start is compacted
        if next_to_query_start: 
            # if node exists next to query start
            if self.get_node_attributes(next_to_query_start)["compacted"]:
                # left contig must be trimmed to remove start of query sequence node (k-1 length) from contig end
                contig_len = len(left_contig)
                print(f"contig left length: {contig_len}")
                k = self.graph['k']
                left_contig = left_contig[0:(contig_len-(k-2))]
            
        # run walk in right direction, starting from query_end
        right_contig = self.walk(query_end, 1)
        
        # combine into final contig: 
        assembled_sequence += left_contig
        assembled_sequence += self.graph.get("query")
        assembled_sequence += right_contig
        
        return assembled_sequence

