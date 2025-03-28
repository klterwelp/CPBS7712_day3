# CPBS7712_day3

This repository holds the code for the CPBS7712 programming assignment. The goal of this programming assignment is to identify the sequence context of a given query sequence through targeted assembly.  

## Installation Instructions

Before you can use this repository, you'll need to have conda installed. Follow [conda installation instructions if you need help setting up conda.](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

1. Clone github repository.

  ```sh
  git clone https://github.com/klterwelp/CPBS7712_day3.git
  ```

2. Set up conda environment

```sh
# move into repo environment folder
cd CPBS7712_day3/envs
# create conda environment
conda env create -f dbg.yaml 
# activate conda environment
conda activate dbg_env
```

## Example Usage

### Assemble longest contig containing your query

1. Place your sequence reads and the query sequence into the `data` folder, ensuring they are both fasta files. You can only search one query sequence at a time. 

```sh
# example of move commands
mv reads.fasta data/
mv query.fasta data/
```

2. Import the diGraph and deBruijnGraphs into your python script.

```python

# import from directedGraph
from directedGraph import diGraph
# import from deBruijnGraph
from directedGraph import deBruijnGraph

```

3. Follow this code, changing the paths to your files. You can optionally change k-mer size (k), but I recommend 31 to start with.

```python
# change these variables:
reads_path="data/reads.fa"
query_path="data/query.fa"
contig_path="data/results.fa" # desired path for contig
k=31 # default k-mer size

# Start your empty deBruijn graph
dbg=deBruijnGraph(k)

# Import your reads into the graph
dbg.add_kmers_from_fasta(reads_path)

# Import your query sequence into the graph
dbg.identify_query_boundaries(query_path)

# Compact graph to remove redundant nodes
dbg.compact_except_query()

# Assemble your query containing contig!
contig = dbg.assemble_from_query()

# Write contig as fasta
with open(contig_path, "w") as fasta_file:
  fasta_file.write(f">query_contig\n")
  fasta_file.write(f"{contig}")

```

In the future, I will write this as a simple function someone can call from the src folder. This is a simple example of how I would pull the class functions together.

## Classes

### diGraph Class

The `diGraph` class is a simplified implementation of a directed graph. It provides a lightweight and customizable representation of directed graphs, inspired by the NetworkX `DiGraph` class but with reduced complexity and no external dependencies.

### Key Features

- Nodes and edges can have associated attributes stored as dictionaries.

- Supports adding and removing nodes and edges.
- Provides methods to query successors, predecessors, in-degree, and out-degree of nodes.
- Tracks the total number of nodes and edges in the graph.
- Does not rely on external libraries like NetworkX.

### Graph Structure

- `_succ`: Dictionary representing the successors of each node (outgoing edges).
- `_pred`: Dictionary representing the predecessors of each node (incoming edges).
- `_node`: Dictionary storing node attributes.
- `graph`: Dictionary storing graph-level attributes.
- `node_count`: Tracks the total number of nodes in the graph.
- `edge_count`: Tracks the total number of edges in the graph.

### Methods

- **Node Operations**
  - `add_node`: Add a single node with optional attributes.
  - `remove_node`: Remove a node and its associated edges.
  - `get_nodes`: Retrieve nodes with optional filtering and attribute inclusion.
  - `get_node_attributes`: Retrieve attributes of a specific node.
  - `set_node_attributes`: Set or update attributes of a specific node.

- **Edge Operations**

  - `add_edge`: Add an edge between two nodes with optional attributes.
  - `remove_edge`: Remove an edge between two nodes.
  - `get_edges`: Retrieve edges with optional filtering and attribute inclusion.
  - `get_edge_attributes`: Retrieve attributes of a specific edge.
  - `set_edge_attributes`: Set or update attributes of a specific edge.

- **Graph-Level Operations**

  - `get_graph_attributes`: Retrieve graph-level attributes.
  - `set_graph_attributes`: Set or update graph-level attributes.
  - `clear_graph`: Remove all nodes and edges from the graph.

- **Query Operations**

  - `has_successor`: Check if a directed edge exists between two nodes.
  - `has_predecessor`: Check if a directed edge exists between two nodes.
  - `successors`: Get a list of successors of a node.
  - `predecessors`: Get a list of predecessors of a node.
  - `in_degree`: Get the in-degree of a node.
  - `out_degree`: Get the out-degree of a node.
  - `number_of_nodes`: Get the total number of nodes in the graph.
  - `number_of_edges`: Get the total number of edges in the graph.

This class is the base for the deBruijnGraph class, which is used to construct the de Bruijn graph for targeted assembly.

### deBruijnGraph Class

The `deBruijnGraph` class is a deBruijn directed graph that is based off of the `DiGraph` class. It contains the functions required to assemble reads into the longest contig containing a given query sequence.

### Key Features

- Nodes are k-mers of size `k`
- When k-mers are added, the reverse-complement is also added.
- The graph can be filtered to the subgraph containing paths to the query k-mers.
- The graph can be compacted and assembled.
- This assembly algorithm cannot handle repeat loops.

### Methods

- `reverse_complement`: generates the reverse complement of a sequence.
- `add_kmer`: adds a k-mer and its reverse complement to the de bruijn graph.
- `add_kmers_from_fasta`: adds all k-mers and their reverse complements to the de bruijn graph.
- `query_kmers_from_fasta`: identifies the k-mers from a given query sequence.
- `identify_query_subgraph`: identifies all k-mers in the graph that have a path to a query k-mer.
- `identify_query_boundaries`: identifies the start and end nodes of a given query sequence and determines if these k-mers all exist in the graph.
- `compact_except_query`: compacts the graph except for query nodes. Converts linear stretches of nodes into one node.
- `assemble_from_query`: assembles the graph starting from the left of the query start and then from the right of the query end. Ensures the query sequence remains identical.