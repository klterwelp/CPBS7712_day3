# Pseudocode

The following is the pseudocode and planning stage of my algorithms. I also include the references that helped me build the code plan.

## References

### Importing fasta files

- [fastapy](https://github.com/aziele/fastapy/blob/main/fastapy.py)

### Graphs
- [NetworkX](https://networkx.org/documentation/)

### De Bruijn Graphs

- [py-debruijn](https://github.com/mitbal/py-debruijn)
- [Eaton Lab Notebook](https://eaton-lab.org/slides/genomics/answers/nb-10.2-de-Bruijn.html)
- [How to apply de Bruijn graphs](https://www.nature.com/articles/nbt.2023)
- [Python bcalm](https://github.com/rchikhi/python-bcalm)
- [miniDBG](https://github.com/mdondrup/miniDBG)

### Targeted Assembly

- [Metacherchant](https://github.com/ctlab/metacherchant)

### General Programming

- [Python Classes](https://docs.python.org/3/tutorial/classes.html)
- [Python classes reddit](https://www.reddit.com/r/learnpython/comments/14y1c31/im_new_to_python_classes_and_objects/)
- [Python data structures](https://docs.python.org/3/tutorial/datastructures.html)

## Algorithm

Import sequences from fasta file

K-mer counter

De Bruijn Graph Class

Identify edge of a graph

Find or create a new node

Compact linear stretches of graph

```python
def import_fasta_sequences(path):
    """ Import fasta and return a named tuple list of sequences """ 

    # read fasta file by line
    lines = readline(path)
    # initiate empty variables
    name = NULL 
    sequence = NULL 
    sequences = () 
    # loop across fasta file
    for line in fasta file 
        if line start ">" 
            # set the name
            name = char after > 
        else 
            if name != NULL 
                # add sequence and name to sequence list
                sequence = line 
                record = named_tuple(
                    name = name, 
                    sequence)
                append record to sequence list 
                set name = NULL 
                set sequence = NULL 
            else 
                log error: fasta file format issue 

```

```python
def get_kmer_list(sequence, k): 
    """ Takes a sequence and generates a list of k-mers """
    i = 0 
    seq_length = num of sequence char
    kmers = list()
    while i is less than seq_length - k + 1: 
        kmer = sequence[i:i+k-1]
        append kmer to kmers
        i += 1
    return kmers
```

### De Bruijn Graph Class
Inspired by NetworkX directed graph 

- `node_dict_factory`
	- key = node ids
	- value = node attributes
	- `self._node` (`__init__`)
- `node_attr_dict_factory`
	- key = attribute name
	- value = attribute values 
- `adjlist_outer_dict_factory`
	- key = node ids
	- value = adjacency info 
	- `self._adj` (`__init__`)
	- `self._pred` (`__init__`)
		- `self._succ` is equal to `_pred`
- `adjlist_inner_dict_factory`
	- key = neighbor
	- value = edge data
- `edge_attr_dict_factory`
	- key = attribute name
	- value = attribute values
- `graph_attr_dict_factory`
	- key = attribute name
	- value = attribute values 
	- `self.graph` (`__init__`)

Add node
```python
def add_node(self, node_for_adding, **attr):
    if node_for_adding not in self._node:




```

Add nodes

Remove node

Remove nodes

Add edge

Add edges