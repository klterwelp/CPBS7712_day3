# CPBS7712_day3

This repository holds the code for the CPBS7712 programming assignment. The goal of this programming assignment is to identify the sequence context of a given query sequence through targeted assembly.  

## diGraph Class

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
