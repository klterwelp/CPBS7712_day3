# Description: This file contains the class definition for a directed graph.

class diGraph: 
    """ A simplified implementation of a directed graph class.
This class provides a lightweight and customizable representation of directed graphs, 
inspired by the NetworkX `DiGraph` class but with reduced complexity and dependencies. 
It uses dictionaries to store graph structure and attributes, allowing for efficient 
manipulation of nodes and edges.
Key Features:
- Nodes and edges can have associated attributes stored as dictionaries.
- Supports adding and removing nodes and edges.
- Provides methods to query successors, predecessors, in-degree, and out-degree of nodes.
- Tracks the total number of nodes and edges in the graph.
- Does not rely on external dependencies like NetworkX.
Attributes:
- `node_dict_factory`: Factory function for creating the node dictionary.
- `node_attr_dict_factory`: Factory function for creating node attribute dictionaries.
- `adjlist_outer_dict_factory`: Factory function for creating the outer adjacency dictionary.
- `adjlist_inner_dict_factory`: Factory function for creating the inner adjacency dictionary.
- `edge_attr_dict_factory`: Factory function for creating edge attribute dictionaries.
- `graph_attr_dict_factory`: Factory function for creating graph attribute dictionaries.
Graph Structure:
- `_succ`: Dictionary representing the successors of each node (outgoing edges).
- `_pred`: Dictionary representing the predecessors of each node (incoming edges).
- `_node`: Dictionary storing node attributes.
- `graph`: Dictionary storing graph-level attributes.
- `node_count`: Tracks the total number of nodes in the graph.
- `edge_count`: Tracks the total number of edges in the graph.
Methods:
- `add_node`: Add a single node to the graph with optional attributes.
- `remove_node`: Remove a node from the graph along with its associated edges.
- `add_edge`: Add an edge between two nodes with optional attributes.
- `remove_edge`: Remove an edge between two nodes.
- `has_successor`: Check if a directed edge exists between two nodes.
- `has_predecessor`: Check if a directed edge exists between two nodes.
- `successors`: Get a list of successors of a node. 
- `predecessors`: Get a list of predecessors of a node.
- `in_degree`: Get the in-degree of a node.
- `out_degree`: Get the out-degree of a node.
- `number_of_nodes`: Get the total number of nodes in the graph.
- `number_of_edges`: Get the total number of edges in the graph.
- `get_nodes`: Retrieve nodes with optional filtering and attribute inclusion.
- `get_edges`: Retrieve edges with optional filtering and attribute inclusion.
- `get_graph_attributes`: Retrieve graph-level attributes.
- `set_graph_attributes`: Set or update graph-level attributes.
- `get_node_attributes`: Retrieve attributes of a specific node.
- `set_node_attributes`: Set or update attributes of a specific node.
- `get_edge_attributes`: Retrieve attributes of a specific edge.
- `set_edge_attributes`: Set or update attributes of a specific edge.
- `clear_graph`: Remove all nodes and edges from the graph.
 """
    # add dictionaries to hold attributes/information of graph
    node_dict_factory = dict
        # key = node_id, value = node attributes
    node_attr_dict_factory = dict
        # key = attribute name, value = attribute values
    adjlist_outer_dict_factory = dict
        # key = node_id, value = adjacency info 
    adjlist_inner_dict_factory = dict
        # key = neighbor, value = edges
    edge_attr_dict_factory = dict
        # key = attribute name, value = attribute values
    graph_attr_dict_factory = dict
        # key = attribute name, value = attribute values
    

# initialize graph 
    def __init__(self, **attr):

        self.graph = self.graph_attr_dict_factory()  # dictionary for graph attributes
        self._node = self.node_dict_factory()  # dictionary for node attributes 
        
        # generate two adjacency lists: successor or predecessor 
        self._succ = self.adjlist_outer_dict_factory()  # empty adjacency dict successor
        self._pred = self.adjlist_outer_dict_factory()  # predecessor
        
        # add class attributes to graph
        self.edge_count = 0 # number of edges in the graph
        self.node_count = 0 # number of nodes in the graph



# function to add node to graph 
    def add_node(self, node_for_adding, **attr):
      """Add a single node and its attributes.  
        Parameters
        ----------
        node_for_adding : node
            A node can be any hashable Python object except None.
        attr : keyword arguments, optional
            Set or change node attributes using key=value.
      """  
      if node_for_adding not in self._succ:  #check node existance
          # add node as key to succ pred adj lists
          self._succ[node_for_adding] = self.adjlist_inner_dict_factory()
          self._pred[node_for_adding] = self.adjlist_inner_dict_factory()
          # add node to node attributes list 
          attr_dict = self._node[node_for_adding] = self.node_attr_dict_factory()
          # add attributes to node
          attr_dict.update(attr)
          # increment node count
          self.node_count += 1
      else:  # update attr even if node already exists
          self._node[node_for_adding].update(attr)
      
# function to remove node from graph
    def remove_node(self, n):
      """Remove node n from graph.
        Parameters
        ----------
         n : node
           A node in the graph
      """
      try: 
        nbrs = self._succ[n]
        del self._node[n]
        self.node_count -= 1  # decrement node count
      except KeyError as err: # if node not in graph, raise error
          raise KeyError(f"The node {n} is not in the directed graph.") from err
      for u in nbrs:
        del self._pred[u][n] # remove all edges to n-u from directed graph
        self.edge_count -= 1 # decrement edge count
      del self._succ[n] # remove node from succ dict list 
      for u in self._pred[n]:
        del self._succ[u][n] # remove all edges from n-u from directed graph
        self.edge_count -= 1 # decrement edge count    
      del self._pred[n] # remove node from pred dict list
      
      
# function to add edge to graph
    def add_edge(self, u_of_edge, v_of_edge, **attr):
        """ Add an edge between u and v.
        The nodes u and v will be automatically added if they are
        not already in the graph.
        Parameters
        ----------
        u_of_edge, v_of_edge : nodes
            Nodes can be any hashable Python object. 
        attr : keyword arguments, optional
            Set or change edge attributes using key=value.
        """
        u, v = u_of_edge, v_of_edge
        # add nodes if they don't already exist
        self.add_node(u)
        self.add_node(v)
        # add the edge
        datadict = self._succ[u].get(v, self.edge_attr_dict_factory())
        datadict.update(attr)
        self._succ[u][v] = datadict
        self._pred[v][u] = datadict
        self.edge_count += 1  # increment edge count

# function to remove edge from graph
    def remove_edge(self, u, v):
        """Remove the edge between nodes u and v.
        Parameters
        ----------
        u, v : nodes
            Nodes can be any hashable Python object.
        """
        try:
            del self._succ[u][v]
            del self._pred[v][u]
            self.edge_count -= 1
        except KeyError as err:
            raise KeyError(f"The edge {u}-{v} is not in the directed graph.") from err
        
                   
# function to check if edge u to v exists
    def has_successor(self, u, v):
        """Return True if node v is a successor of node u, False otherwise.
        Parameters
        ----------
        u, v : nodes
            Nodes can be any hashable Python object.
        """
        return v in self._succ[u]        

# function to check if edge v to u exists 
    def has_predecessor(self, u, v):
        """Return True if node v is a predecessor of node u, False otherwise.
        Parameters
        ----------
        u, v : nodes
            Nodes can be any hashable Python object.
        """
        return v in self._pred[u]
    
# function to get the successors of a node
    def successors(self, n):
        """Return an iterator over the successors of node n.
        Parameters
        ----------
        n : node
            A node in the graph.
        """
        try:
            return iter(self._succ[n])
        except KeyError as err:
            raise KeyError(f"The node {n} is not in the directed graph.") from err

# function to get the predecessors of a node
    def predecessors(self, n):
        """Return an iterator over the predecessors of node n.
        Parameters
        ----------
        n : node
            A node in the graph.
        """
        try:
            return iter(self._pred[n])
        except KeyError as err:
            raise KeyError(f"The node {n} is not in the directed graph.") from err
        
# function to get the in_degree of a node
    def in_degree(self, n, weight=None):
        """Return the in-degree of node n.
        Parameters
        ----------
        n : node
            A node in the graph.
        weight : string or None, optional (default=None)
            The edge attribute that holds the numerical value used as a weight.
            If None, the in-degree is returned as the number of edges.
        """
        try:
            if weight is None:
                return len(self._pred[n])
            else:
                return sum(self._pred[n][nbr].get(weight, 1) for nbr in self._pred[n])
        except KeyError as err:
            raise KeyError(f"The node {n} is not in the directed graph.") from err

# function to get the out_degree of a node
    def out_degree(self, n, weight=None):
        """Return the out-degree of node n.
        Parameters
        ----------
        n : node
            A node in the graph.
        weight : string or None, optional (default=None)
            The edge attribute that holds the numerical value used as a weight.
            If None, the out-degree is returned as the number of edges.
        """
        try:
            if weight is None:
                return len(self._succ[n])
            else:
                return sum(self._succ[n][nbr].get(weight, 1) for nbr in self._succ[n])
        except KeyError as err:
            raise KeyError(f"The node {n} is not in the directed graph.") from err

# function to get the number of nodes in graph
    def number_of_nodes(self):
        """Return the number of nodes in the graph.
        Returns
        -------
        nnodes : int
            The number of nodes in the graph.
        """
        return self.node_count
    
# function to get the number of edges in graph
    def number_of_edges(self):
        """Return the number of edges in the graph.
        Returns
        -------
        nedges : int
            The number of edges in the graph.
        """
        return self.edge_count

# function to get the nodes in graph
    def get_nodes(self, include_data=False, **filter_attributes):
        """Return an iterator over the nodes in the graph.
        Parameters
        ----------
        include_data : bool, optional (default=False)
            If True, return a two-tuple of node and node attribute dictionary.
        filter_attributes : keyword arguments, optional
            Filters nodes by specified attributes. Only nodes with matching
            attribute values will be returned.
        Returns
        -------
        node_iterator : iterator
            An iterator over all nodes in the graph, optionally filtered by attributes.
        """
        if include_data:
            if filter_attributes:
                for node, attributes in self._node.items():
                    if all(attributes.get(key) == value for key, value in filter_attributes.items()):
                        yield (node, attributes)
            else:
                for node, attributes in self._node.items():
                    yield (node, attributes)
        else:
            if filter_attributes:
                for node, attributes in self._node.items():
                    if all(attributes.get(key) == value for key, value in filter_attributes.items()):
                        yield node
            else:
                yield from self._succ

# function to get the edges in graph
    def get_edges(self, include_data=False, **filter_attributes):
        """Return an iterator over the edges in the graph.
        Parameters
        ----------
        include_data : bool, optional (default=False)
            If True, return a three-tuple of (source_node, target_node, edge_attribute_dict).
        filter_attributes : keyword arguments, optional
            Filters edges by specified attributes. Only edges with matching
            attribute values will be returned.
        Returns
        -------
        edge_iterator : iterator
            An iterator over all edges in the graph, optionally filtered by attributes.
        """
        for source_node in self._succ:
            for target_node, edge_attributes in self._succ[source_node].items():
                if filter_attributes:
                    if all(edge_attributes.get(attr_key) == attr_value for attr_key, attr_value in filter_attributes.items()):
                        if include_data:
                            yield (source_node, target_node, edge_attributes)
                        else:
                            yield (source_node, target_node)
                else:
                    if include_data:
                        yield (source_node, target_node, edge_attributes)
                    else:
                        yield (source_node, target_node)
                        
# function to get the graph attributes
    def get_graph_attributes(self):
        """Return a dictionary of graph attributes.
        Returns
        -------
        graph_attributes : dict
            A dictionary containing the graph-level attributes.
        """
        return self.graph

# function to set the graph attributes
    def set_graph_attributes(self, **attr):
        """Set or update graph attributes.
        Parameters
        ----------
        attr : keyword arguments
            Set or change graph attributes using key=value.
        """
        self.graph.update(attr)
        
# function to get the node attributes
    def get_node_attributes(self, n):
        """Return a dictionary of attributes for node n.
        Parameters
        ----------
        n : node
            A node in the graph.
        Returns
        -------
        node_attributes : dict
            A dictionary containing the attributes of the specified node.
        """
        try:
            return self._node[n]
        except KeyError as err:
            raise KeyError(f"The node {n} is not in the directed graph.") from err

# function to set the node attributes
    def set_node_attributes(self, n, **attr):
        """Set or update attributes for node n.
        Parameters
        ----------
        n : node
            A node in the graph.
        attr : keyword arguments
            Set or change node attributes using key=value.
        """
        try:
            self._node[n].update(attr)
        except KeyError as err:
            raise KeyError(f"The node {n} is not in the directed graph.") from err

# function to get the edge attributes
    def get_edge_attributes(self, u, v):
        """Return a dictionary of attributes for the edge between nodes u and v.
        Parameters
        ----------
        u, v : nodes
            Nodes can be any hashable Python object.
        Returns
        -------
        edge_attributes : dict
            A dictionary containing the attributes of the specified edge.
        """
        try:
            return self._succ[u][v]
        except KeyError as err:
            raise KeyError(f"The edge {u}-{v} is not in the directed graph.") from err
        
# function to set the edge attributes
    def set_edge_attributes(self, u, v, **attr): 
        """Set or update attributes for the edge between nodes u and v.
        Parameters
        ----------
        u, v : nodes
            Nodes can be any hashable Python object.
        attr : keyword arguments
            Set or change edge attributes using key=value.
        """
        try:
            self._succ[u][v].update(attr)
        except KeyError as err:
            raise KeyError(f"The edge {u}-{v} is not in the directed graph.") from err

# function to clear graph 
    def clear_graph(self): 
        """Remove all nodes and edges from the graph."""
        self.graph.clear()
        self._node.clear() 
        self._succ.clear()
        self._pred.clear()
        self.edge_count = 0
        self.node_count = 0

