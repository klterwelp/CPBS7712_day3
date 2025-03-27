# Tests the directedGraph class in this file.
import unittest
import sys
import os
import shutil
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
try:
    from directedGraph import diGraph
except ModuleNotFoundError:
    raise ModuleNotFoundError("Ensure 'directedGraph.py' exists in the '../src' directory relative to this test file.")

class TestDiGraph(unittest.TestCase):
    """
    Unit tests for the diGraph class.
    This test suite verifies the functionality of the diGraph class, including
    operations for adding, removing, and querying nodes and edges, as well as
    handling attributes and calculating graph properties.
    Test Cases:
    - `test_add_node`: Ensures nodes can be added to the graph and verifies their presence.
    - `test_node_attributes`: Tests adding and updating attributes for nodes.
    - `test_remove_node`: Verifies that nodes can be removed and checks the graph's state.
    - `test_remove_nonexistent_node`: Ensures an exception is raised when attempting to remove a non-existent node.
    - `test_add_edge`: Tests adding edges between nodes and verifies their presence.
    - `test_remove_edge`: Ensures edges can be removed and verifies their absence.
    - `test_has_successor`: Checks if a node has a specific successor.
    - `test_has_predecessor`: Checks if a node has a specific predecessor.
    - `test_successors`: Verifies the list of successors for a given node.
    - `test_predecessors`: Verifies the list of predecessors for a given node.
    - `test_in_degree`: Tests the in-degree calculation for a node.
    - `test_out_degree`: Tests the out-degree calculation for a node.
    - `test_in_degree_weighted`: Verifies the weighted in-degree calculation for a node.
    - `test_out_degree_weighted`: Verifies the weighted out-degree calculation for a node.
    - `test_add_edge_with_attributes`: Tests adding edges with attributes and verifies their values.
    - `test_update_edge_attributes`: Ensures edge attributes can be updated correctly.
    - `test_number_of_nodes`: Verifies the total number of nodes in the graph after various operations.
    - `test_number_of_nodes_after_removing_edge`: Ensures the number of nodes remains consistent after removing edges.
    - `test_number_of_edges_after_removing_node`: Verifies the total number of edges after removing a node.
    - `test_get_nodes`: Tests retrieving nodes with and without attributes, and with attribute filtering.
    - `test_get_edges`: Tests retrieving edges with and without attributes, and with attribute filtering.
    - `test_get_graph_attributes`: Verifies retrieving graph-level attributes.
    - `test_set_graph_attributes`: Ensures graph-level attributes can be set and updated correctly.
    - `test_get_node_attributes`: Tests retrieving attributes for a specific node.
    - `test_set_node_attributes`: Ensures attributes for a specific node can be set and updated.
    - `test_get_edge_attributes`: Tests retrieving attributes for a specific edge.
    - `test_set_edge_attributes`: Ensures attributes for a specific edge can be set and updated.
    - `test_clear_graph`: Verifies that the graph can be cleared of all nodes, edges, and attributes.
    - `test_add_duplicate_node`: Ensures adding a duplicate node updates attributes without overwriting existing ones.
    - `test_add_duplicate_edge`: Ensures adding a duplicate edge updates attributes without overwriting existing ones.
    - `test_remove_nonexistent_edge`: Ensures an exception is raised when attempting to remove a non-existent edge.
    - `test_successors_nonexistent_node`: Ensures an exception is raised when querying successors of a non-existent node.
    - `test_predecessors_nonexistent_node`: Ensures an exception is raised when querying predecessors of a non-existent node.
    - `test_in_degree_nonexistent_node`: Ensures an exception is raised when querying the in-degree of a non-existent node.
    - `test_out_degree_nonexistent_node`: Ensures an exception is raised when querying the out-degree of a non-existent node.
    - `test_get_node_attributes_nonexistent_node`: Ensures an exception is raised when querying attributes of a non-existent node.
    - `test_set_node_attributes_nonexistent_node`: Ensures an exception is raised when setting attributes for a non-existent node.
    - `test_get_edge_attributes_nonexistent_edge`: Ensures an exception is raised when querying attributes of a non-existent edge.
    - `test_set_edge_attributes_nonexistent_edge`: Ensures an exception is raised when setting attributes for a non-existent edge.
    - `test_clear_graph_with_attributes`: Verifies that clearing the graph removes all nodes, edges, and attributes, including graph-level attributes.
    - `test_add_edge_with_nonexistent_nodes`: Ensures edges can be added even if the nodes do not exist, and verifies their creation.
    - `test_str`: Checks the string representation of the graph
    - `test_visualize_graph`: Ensures the graph visualization function produces the expected output (png and pdf files)
    Each test case uses assertions to validate the expected behavior of the diGraph class.
    """

    def setUp(self):
        self.graph = diGraph()

    def test_add_node(self):
        self.graph.add_node(1)
        self.assertIn(1, self.graph._succ)
        self.assertIn(1, self.graph._pred)
        self.assertIn(1, self.graph._node)
        self.assertEqual(self.graph.number_of_nodes(), 1)  
        
    def test_node_attributes(self):
        self.graph.add_node(1, color='blue', size=10)
        self.assertIn(1, self.graph._node)
        self.assertEqual(self.graph._node[1]['color'], 'blue')
        self.assertEqual(self.graph._node[1]['size'], 10)
        # Update attributes
        self.graph.add_node(1, color='red')
        self.assertEqual(self.graph._node[1]['color'], 'red')
        self.assertEqual(self.graph._node[1]['size'], 10)  # Ensure size remains unchanged
    
    def test_has_node(self):
        self.graph.add_node(1)
        self.assertTrue(self.graph.has_node(1))
        self.assertFalse(self.graph.has_node(2))
        
    def test_remove_node(self):
        self.graph.add_node(1)
        self.graph.remove_node(1)
        self.assertNotIn(1, self.graph._succ)
        self.assertNotIn(1, self.graph._pred)
        self.assertNotIn(1, self.graph._node)
        self.assertEqual(self.graph.number_of_nodes(), 0)  

        
    def test_remove_nonexistent_node(self):
        with self.assertRaises(KeyError) as context:
            self.graph.remove_node(99)  # Attempt to remove a node that doesn't exist
        expected_message = f"The node {99} is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)
        
    def test_add_edge(self):
        self.graph.add_edge(1, 2)
        self.assertIn(2, self.graph._succ[1])
        self.assertIn(1, self.graph._pred[2])

    def test_remove_edge(self):
        self.graph.add_edge(1, 2)
        self.graph.remove_edge(1, 2)
        self.assertNotIn(2, self.graph._succ[1])
        self.assertNotIn(1, self.graph._pred[2])

    def test_has_successor(self):
        self.graph.add_edge(1, 2)
        self.assertTrue(self.graph.has_successor(1, 2))
        self.assertFalse(self.graph.has_successor(2, 1))

    def test_has_predecessor(self):
        self.graph.add_edge(1, 2)
        self.assertTrue(self.graph.has_predecessor(2, 1))
        self.assertFalse(self.graph.has_predecessor(1, 2))

    def test_successors(self):
        self.graph.add_edge(1, 2)
        self.graph.add_edge(1, 3)
        self.assertListEqual(list(self.graph.successors(1)), [2, 3])

    def test_predecessors(self):
        self.graph.add_edge(1, 2)
        self.graph.add_edge(3, 2)
        self.assertListEqual(list(self.graph.predecessors(2)), [1, 3])

    def test_in_degree(self):
        self.graph.add_edge(1, 2)
        self.graph.add_edge(3, 2)
        self.assertEqual(self.graph.in_degree(2), 2)

    def test_out_degree(self):
        self.graph.add_edge(1, 2, weight=2)
        self.graph.add_edge(1, 3)
        self.assertEqual(self.graph.out_degree(1), 2)

    def test_in_degree_weighted(self):
        self.graph.add_edge(1, 2, weight=2)
        self.graph.add_edge(3, 2, weight=3)
        self.assertEqual(self.graph.in_degree(2, weight='weight'), 5)

    def test_out_degree_weighted(self):
        self.graph.add_edge(1, 2, weight=2)
        self.graph.add_edge(1, 3, weight=3)
        self.assertEqual(self.graph.out_degree(1, weight='weight'), 5)

    def test_add_edge_with_attributes(self):
        self.graph.add_edge(1, 2, weight=4, color='red')
        self.assertIn(2, self.graph._succ[1])
        self.assertEqual(self.graph._succ[1][2]['weight'], 4)
        self.assertEqual(self.graph._succ[1][2]['color'], 'red')

    def test_update_edge_attributes(self):
        self.graph.add_edge(1, 2, weight=4)
        self.graph.add_edge(1, 2, weight=10, color='blue')
        self.assertEqual(self.graph._succ[1][2]['weight'], 10)
        self.assertEqual(self.graph._succ[1][2]['color'], 'blue')
        
    def test_number_of_nodes(self):
        self.assertEqual(self.graph.number_of_nodes(), 0)
        self.graph.add_node(1)
        self.assertEqual(self.graph.number_of_nodes(), 1)
        self.graph.add_node(2)
        self.assertEqual(self.graph.number_of_nodes(), 2)
        self.graph.remove_node(1)
        self.assertEqual(self.graph.number_of_nodes(), 1)

    def test_number_of_nodes_after_removing_edge(self):
        self.graph.add_edge(1, 2)
        self.graph.add_edge(2, 3)
        self.assertEqual(self.graph.number_of_nodes(), 3)  
        self.graph.remove_edge(1, 2)
        self.assertEqual(self.graph.number_of_nodes(), 3)  

    def test_number_of_edges_after_removing_node(self):
        self.graph.add_edge(1, 2)
        self.graph.add_edge(2, 3)
        self.graph.add_edge(1, 3)
        self.assertEqual(self.graph.number_of_edges(), 3) 
        self.graph.remove_node(2)
        self.assertEqual(self.graph.number_of_edges(), 1) 
        
    def test_get_nodes(self):
        self.graph.add_node(1, color='blue')
        self.graph.add_node(2, size=10)
        self.graph.add_node(3, color='red', size=15)
        
        # Test without attributes
        self.assertListEqual(list(self.graph.get_nodes()), [1, 2, 3])
        
        # Test with include_data=True
        self.assertListEqual(
            list(self.graph.get_nodes(include_data=True)),
            [(1, {'color': 'blue'}), (2, {'size': 10}), (3, {'color': 'red', 'size': 15})]
        )
        
        # Test with attribute filtering
        self.assertListEqual(
            list(self.graph.get_nodes(color='red')),
            [3]
        )
        self.assertListEqual(
            list(self.graph.get_nodes(size=10)),
            [2]
        )
        self.assertListEqual(
            list(self.graph.get_nodes(color='blue', size=10)),
            []
        )
        # Test retrieving nodes with include_data=True and filtering by a single attribute
        self.assertListEqual(
            list(self.graph.get_nodes(include_data=True, color='blue')),
                [(1, {'color': 'blue'})])
            
            # Test retrieving nodes with include_data=True and filtering by multiple attributes
        self.assertListEqual(
            list(self.graph.get_nodes(include_data=True, color='blue', size=10)),
                []
            )

    def test_get_edges(self):
        self.graph.add_edge(1, 2, weight=4, color='red')
        self.graph.add_edge(2, 3, weight=5)
        self.graph.add_edge(3, 1, color='blue')
        
        # Test without attributes
        self.assertListEqual(
            list(self.graph.get_edges()),
            [(1, 2), (2, 3), (3, 1)]
        )
        
        # Test with include_data=True
        self.assertListEqual(
            list(self.graph.get_edges(include_data=True)),
            [
                (1, 2, {'weight': 4, 'color': 'red'}),
                (2, 3, {'weight': 5}),
                (3, 1, {'color': 'blue'})
            ]
        )
        
        # Test with attribute filtering
        self.assertListEqual(
            list(self.graph.get_edges(color='red')),
            [(1, 2)]
        )
        self.assertListEqual(
            list(self.graph.get_edges(weight=5)),
            [(2, 3)]
        )
        self.assertListEqual(
            list(self.graph.get_edges(color='blue', weight=4)),
            []
        )
    def test_get_graph_attributes(self):
        self.graph.set_graph_attributes(name="TestGraph", directed=True)
        attributes = self.graph.get_graph_attributes()
        self.assertEqual(attributes['name'], "TestGraph")
        self.assertEqual(attributes['directed'], True)

    def test_set_graph_attributes(self):
        self.graph.set_graph_attributes(name="Graph1")
        self.assertEqual(self.graph.graph['name'], "Graph1")
        self.graph.set_graph_attributes(name="Graph2", version=1.0)
        self.assertEqual(self.graph.graph['name'], "Graph2")
        self.assertEqual(self.graph.graph['version'], 1.0)

    def test_get_node_attributes(self):
        self.graph.add_node(1, color='blue', size=10)
        attributes = self.graph.get_node_attributes(1)
        self.assertEqual(attributes['color'], 'blue')
        self.assertEqual(attributes['size'], 10)

    def test_set_node_attributes(self):
        self.graph.add_node(1, color='blue')
        self.graph.set_node_attributes(1, size=10)
        self.assertEqual(self.graph._node[1]['color'], 'blue')
        self.assertEqual(self.graph._node[1]['size'], 10)

    def test_set_multiple_node_attributes(self):
        self.graph.add_node(1, color='blue')
        self.graph.set_node_attributes(1, color='red', size=10)
        self.assertEqual(self.graph._node[1]['color'], 'red')
        self.assertEqual(self.graph._node[1]['size'], 10)
    
    def test_set_multiple_node_attributes_from_variable(self):
        node_attributes = {
            'color': 'green',
            'size': 20,
            'shape': 'circle',
            'label': 'Node1',
            'weight': 5.5,
            'active': True,
            'group': 'A',
            'x': 10.0,
            'y': 15.0
        }
        self.graph.add_node(1)
        self.graph.set_node_attributes(1, **node_attributes)
        for key, value in node_attributes.items():
            self.assertEqual(self.graph._node[1][key], value)
        
    def test_get_edge_attributes(self):
        self.graph.add_edge(1, 2, weight=4, color='red')
        attributes = self.graph.get_edge_attributes(1, 2)
        self.assertEqual(attributes['weight'], 4)
        self.assertEqual(attributes['color'], 'red')

    def test_set_edge_attributes(self):
        self.graph.add_edge(1, 2, weight=4)
        self.graph.set_edge_attributes(1, 2, color='blue', weight=10)
        self.assertEqual(self.graph._succ[1][2]['weight'], 10)
        self.assertEqual(self.graph._succ[1][2]['color'], 'blue')

    def test_clear_graph(self):
        self.graph.add_node(1)
        self.graph.add_node(2)
        self.graph.add_edge(1, 2)
        self.graph.clear_graph()
        self.assertEqual(self.graph.number_of_nodes(), 0)
        self.assertEqual(self.graph.number_of_edges(), 0)
        self.assertEqual(len(self.graph._succ), 0)
        self.assertEqual(len(self.graph._pred), 0)
        self.assertEqual(len(self.graph._node), 0)
        self.assertEqual(len(self.graph.graph), 0)
    
    def test_add_duplicate_node(self):
        self.graph.add_node(1, color='blue')
        self.graph.add_node(1, size=10)  # Adding the same node with different attributes
        self.assertEqual(self.graph._node[1]['color'], 'blue')  # Ensure original attribute remains
        self.assertEqual(self.graph._node[1]['size'], 10)  # Ensure new attribute is added

    def test_add_duplicate_edge(self):
        self.graph.add_edge(1, 2, weight=4)
        self.graph.add_edge(1, 2, color='red')  # Adding the same edge with different attributes
        self.assertEqual(self.graph._succ[1][2]['weight'], 4)  # Ensure original attribute remains
        self.assertEqual(self.graph._succ[1][2]['color'], 'red')  # Ensure new attribute is added

    def test_remove_nonexistent_edge(self):
        with self.assertRaises(KeyError) as context:
            self.graph.remove_edge(1, 2)  # Attempt to remove an edge that doesn't exist
        expected_message = f"The edge 1-2 is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)

    def test_successors_nonexistent_node(self):
        with self.assertRaises(KeyError) as context:
            list(self.graph.successors(99))  # Attempt to get successors of a non-existent node
        expected_message = f"The node 99 is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)

    def test_predecessors_nonexistent_node(self):
        with self.assertRaises(KeyError) as context:
            list(self.graph.predecessors(99))  # Attempt to get predecessors of a non-existent node
        expected_message = f"The node 99 is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)

    def test_in_degree_nonexistent_node(self):
        with self.assertRaises(KeyError) as context:
            self.graph.in_degree(99)  # Attempt to get in-degree of a non-existent node
        expected_message = f"The node 99 is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)

    def test_out_degree_nonexistent_node(self):
        with self.assertRaises(KeyError) as context:
            self.graph.out_degree(99)  # Attempt to get out-degree of a non-existent node
        expected_message = f"The node 99 is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)

    def test_get_node_attributes_nonexistent_node(self):
        with self.assertRaises(KeyError) as context:
            self.graph.get_node_attributes(99)  # Attempt to get attributes of a non-existent node
        expected_message = f"The node 99 is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)

    def test_set_node_attributes_nonexistent_node(self):
        with self.assertRaises(KeyError) as context:
            self.graph.set_node_attributes(99, color='blue')  # Attempt to set attributes for a non-existent node
        expected_message = f"The node 99 is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)

    def test_get_edge_attributes_nonexistent_edge(self):
        with self.assertRaises(KeyError) as context:
            self.graph.get_edge_attributes(1, 2)  # Attempt to get attributes of a non-existent edge
        expected_message = f"The edge 1-2 is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)

    def test_set_edge_attributes_nonexistent_edge(self):
        with self.assertRaises(KeyError) as context:
            self.graph.set_edge_attributes(1, 2, color='blue')  # Attempt to set attributes for a non-existent edge
        expected_message = f"The edge 1-2 is not in the directed graph."
        self.assertEqual(str(context.exception).strip("'"), expected_message)

    def test_clear_graph_with_attributes(self):
        self.graph.add_node(1, color='blue')
        self.graph.add_edge(1, 2, weight=4)
        self.graph.set_graph_attributes(name="TestGraph")
        self.graph.clear_graph()
        self.assertEqual(self.graph.number_of_nodes(), 0)
        self.assertEqual(self.graph.number_of_edges(), 0)
        self.assertEqual(len(self.graph.graph), 0)
    
    def test_add_edge_with_nonexistent_nodes(self):
        # Adding an edge where neither node exists
        self.graph.add_edge(1, 2)
        self.assertIn(1, self.graph._succ)
        self.assertIn(2, self.graph._succ)
        self.assertIn(2, self.graph._succ[1])
        self.assertIn(1, self.graph._pred[2])
        self.assertEqual(self.graph.number_of_nodes(), 2)
        self.assertEqual(self.graph.number_of_edges(), 1)

        # Adding an edge where one node exists
        self.graph.add_edge(2, 3)
        self.assertIn(3, self.graph._succ)
        self.assertIn(3, self.graph._succ[2])
        self.assertIn(2, self.graph._pred[3])
        self.assertEqual(self.graph.number_of_nodes(), 3)
        self.assertEqual(self.graph.number_of_edges(), 2)
    
    def test_str(self):
        self.graph.add_node(1, color='blue')
        self.graph.add_edge(1, 2, weight=4)
        self.assertEqual(str(self.graph), "Directed Graph with 2 nodes and 1 edges.")
        
    def test_visualize_graph(self, png_path="../data/test.png"):
        pdf_path="../data/test.png.pdf"
        self.graph.add_edge("ATG", "TGC", weight=4)
        self.graph.visualize_graph(png_path)
        self.assertTrue(os.path.exists(png_path))
        # clean files
        os.remove(png_path)
        os.remove(pdf_path)
        
if __name__ == '__main__':
    unittest.main()