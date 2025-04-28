from deBruijnGraph import deBruijnGraph
import alignment as al
import time
import argparse

parser = argparse.ArgumentParser(
                    prog='Assemble and Align Reads based on a Query Sequence',
                    description='''
                    This program assembles a contig containing the query sequence and the alignments of reads that make up the contig.''',
                    epilog='''Example usage:
                    python src/assemble_align.py -r data/READS.fasta -q data/QUERY.fasta -o data/contig.fa -O data/alignments.tsv -k 31 -K 20 -i 1 -m 1
                    ''')

parser.add_argument('-r', '--read_path', 
                    type=str, required=True,
                    help="Full path to the reads file in fasta format. The file will be used to assemble a contig containing the query sequence.")

parser.add_argument('-q', '--query_path', 
                    type=str, required=True,
                    help="Full path to the query file in fasta format. This fasta file will be used to seed the assembly.")

parser.add_argument('-o', '--output_assembly_path', 
                    type=str, required=True, 
                    help="Full path of the file to write the final contig assembly in fasta format.")

parser.add_argument('-O', '--output_alignment_path',
                    type=str, required=True,
                    help="Full path of the file to write the final alignment in tsv format.")

parser.add_argument('-k', '--k_assembly', 
                    type=int, required=False, default=31, 
                    help="Length of the k-mer to use for assembly. Default is 31.")

parser.add_argument('-K', '--k_alignment', 
                    type=int, required=False, default=20, 
                    help="Length of the k-mer to use for alignment. Default is 20.")

def interval_type(value):
    if value == "default":
        return value
    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Interval must be an integer or 'default'.")

parser.add_argument('-i', '--interval',
                    type=interval_type, required=False, default="default",
                    help="Length of the interval to use for alignment. Default adjusts the interval based on read length.")

parser.add_argument('-m', '--min_count',
                    type=int, required=False, default=1,
                    help="Minimum count of k-mers that appear in reads to include in the assembly. Default is 1.")

# parse the arguments
args = parser.parse_args()
read_path = args.read_path
query_path = args.query_path
output_assembly_path = args.output_assembly_path
output_alignment_path = args.output_alignment_path
k_assembly = args.k_assembly
k_alignment = args.k_alignment
interval = args.interval
min_count = args.min_count

def print_message(message):
    local_time = time.localtime()
    current_time = time.strftime("%H:%M:%S", local_time)
    print(f"{current_time} : {message}")
    
def print_time(start_time, end_time, message):
    elapsed_time = round(end_time - start_time, 2)
    print(f"{message} took: {elapsed_time} seconds")

# initiate empty de bruijn graph
dbg=deBruijnGraph(k_assembly)

# Import reads into the graph
print_message("Adding k-mers to the graph...")
t1 = time.perf_counter()
dbg.add_kmers_from_fasta(read_path, min_count)
t2 = time.perf_counter()
print_time(t1, t2, "Adding read k-mers to the graph")
print_message(f"Pre-filtered graph structure: {str(dbg)}")

# Identifying query boundaries
print_message("Identifying query in the graph...")
t3 = time.perf_counter()
dbg.identify_query_boundaries(query_path)
t4 = time.perf_counter()
print_time(t3, t4, "Identifying query in the graph")
dbg_graph=dbg.get_graph_attributes()
print_message(f"Query start: {dbg_graph['query_start']}")
print_message(f"Query End: {dbg_graph['query_end']}")

# Filter the graph to only include query connected nodes
print_message("Filtering graph to only include query connected k-mers...")
t5 = time.perf_counter()
dbg.filter_graph()
t6 = time.perf_counter()
print_time(t5, t6, "Filtering graph")
print_message(f"Filtered graph structure: {str(dbg)}")

# Compact the graph to remove redundant k-mer nodes
print_message("Compacting graph to remove redundant k-mer nodes...")
t7 = time.perf_counter()
dbg.compact_except_query()
t8 = time.perf_counter()
print_time(t7, t8, "Compacting graph")
print_message(f"Compacted graph structure: {str(dbg)}")

# Assemble contig starting from query sequence
print_message("Assembling contig...")
t9 = time.perf_counter()
assembly=dbg.assemble_from_query()
t10 = time.perf_counter()
print_time(t9, t10, "Assembling contig")

# Save assembly to a file
print_message("Writing assembly to file...")
dbg.write_assembly_to_fasta(output_assembly_path, assembly, "contig")
print_message(f"Assembly written to {output_assembly_path}")

# Import read sequences
print_message("Reading sequences for alignment...")
t11 = time.perf_counter()
sequences = al.read_fasta_with_sequence_name(read_path)
t12 = time.perf_counter()
print_time(t11, t12, "Reading sequences")

# Build index for alignment
print_message("Building index for alignment...")
t13 = time.perf_counter()
contig_index= al.build_index(assembly, k_alignment)
t14 = time.perf_counter()
print_time(t13, t14, "Building index for alignment")

# Align reads to contig
print_message("Aligning reads to contig...")
t15 = time.perf_counter()
alignments = al.align_reads_to_contig(sequences, assembly, "contig", k=k_alignment, I=interval)
t16 = time.perf_counter()
print_time(t15, t16, "Aligning reads to contig")
print_message(f"{len(alignments)} alignments found")

# Export alignments to output file
print_message("Exporting alignments to file...")
t17 = time.perf_counter()
al.export_alignment_to_tsv(alignments, output_alignment_path)
t18 = time.perf_counter()
print_time(t17, t18, f"Exporting alignments to file {output_alignment_path}")

print_message("All tasks completed successfully.")
print_time(t1, t18, "Running full assembly to alignment pipeline")