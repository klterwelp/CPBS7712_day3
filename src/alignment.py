# Simplified BLAST-like seed-and-extend alignment algorithm 
# STEP ONE: BUILD CONTIG INDEX
# build an index for the given contig sequence 
# index is a dictionary where the keys are k-mers and the values are lists of positions

# STEP TWO: SCAN READ SEQUENCE FOR SEEDS
# rather than all continous k-mers, we will use spaced k-mers of interval I 
# where I is interval length as a function of the length of the read, x
# scan the read sequence for spaced k-mers of interval I that match the contig index
# calculate the interval length I using the formula from Bowtie2: 
# I(x) = max(1, floor(1 + 1.15 * √x))
# given a read sequence and a contig index:
# for i in series of 1 to len read sequence with interval I:
# for each k-mer in the read sequence, if it also exists in the contig index:
    # add to the list of seeds
# return the list of seeds as a list of tuples: (contig_index_position, read_index_position)

# STEP THREE: PRIORITIZE SEEDS
# rather than extend all the seeds, we will prioritize the first seed of the most common diagonal. 
# diagonal is defined by the difference between read position - contig position.
# for example, if you compared "ATCG" to "ATCG", the top 2-mer seeds would be:
# (0, 0), (1, 1), (2, 2), (3, 3) which all have a diagonal of 0.

# STEP FOUR: EXTEND SEED
# extend the seed in both directions until the extension fails
# return the extended seed as a tuple: (contig_start_position, read_start_position, length) 

# STEP FIVE: ALIGN READS TO CONTIG
# using the above steps, we will align the forward and reverse complement of a read to the contig
# if the longest alignment is the reverse complement, update the alignment information
# return the alignment information as a tuple: (contig_start, contig_end, read_start, read_end, length, alignment_strand)

import csv

def read_fasta_with_sequence_name(file_path):
    """
    Read a FASTA file and return each sequence and its name.
    Return a dictionary where the keys are the sequence names and the values are the sequences.
    """
    # Open the file and read the lines
    with open(file_path, "r") as file:
        lines = file.readlines()
    # Initialize variables
    sequences = {}
    current_name = None
    current_sequence = []
    # Iterate through the lines in the file
    for line in lines:
        # If the line starts with '>', set as current sequence name
        if line.startswith(">"):
            # If there is a current sequence, save it to the dictionary
            if current_name is not None:
                sequences[current_name] = "".join(current_sequence)
            # Set the new current name and reset the current sequence
            current_name = line[1:].strip() # remove the '>' character
            current_sequence = []
        else:
            # Append the line to the current sequence
            current_sequence.append(line.strip())
    
    # Save the last sequence to the dictionary 
    if current_name is not None:
        sequences[current_name] = "".join(current_sequence)
    # Return the dictionary of sequences 
    return sequences     

def reverse_complement(seq): 
    """
    Return the reverse complement of a given sequence.
    """
    # Create a dictionary to map nucleotides to their complements
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    # Reverse the sequence and get the complement
    return "".join(complement[base] for base in reversed(seq))
        
def build_index(seq, k=15):
    """
    Build an index for the given sequence.
    The index is a dictionary where the keys are k-mers and the values are lists of positions.
    """
    index = {}
    k = k  # Length of k-mer
    # check if k is larger than the sequence length
    if k > len(seq):
        raise ValueError("k cannot be larger than the sequence length")
    # for every k-mer in the sequence, add it to the index and its position
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer not in index:
            index[kmer] = []
        index[kmer].append(i)    
    return index

def scan_read_for_seeds(read_seq, read_len, contig_index, k = 15, I = "default"):
    """
    Scan the read sequence for spaced k-mers of interval I that match the contig index.
    Calculate the interval length I using the formula from Bowtie2:
    I(x) = max(1, floor(1 + 1.15 * √x))
    where x is the length of the read sequence.
    Else, use the provided interval length I. 
    """
    # Calculate the interval length I using the formula from Bowtie2:
    if I == "default":
        I = max(1, int(1 + 1.15 * (read_len ** 0.5)))
    
    # Start empty seeds list
    seeds = []
    # Scan the read sequence at intervals of length I for matching k-mers 
    for r_pos in range(0, read_len - k + 1, I):
        spaced_kmer = read_seq[r_pos:r_pos + k]
        if spaced_kmer in contig_index:
            for q_pos in contig_index[spaced_kmer]:
                seeds.append((q_pos, r_pos))
    return seeds

def prioritize_seeds_by_diagonal(seeds):
    """
    Filter seeds to only those with the most common diagonal value.
    Diagonal is defined by the difference between read position - contig position. 
    For example, if you compared "ATCG" to "ATCG", the top 2-mer seeds would be: 
    (0, 0), (1, 1), (2, 2), (3, 3) which all have a diagonal of 0. 
    Returns the filtered list of seeds with the most common diagonal.
    """
    # Create a dictionary to store the diagonals and their corresponding counts
    diagonal_counts = {}
    for seed in seeds:
        q_pos, r_pos = seed
        diagonal = r_pos - q_pos
        if diagonal not in diagonal_counts:
            diagonal_counts[diagonal] = 0
        diagonal_counts[diagonal] += 1
    
    # Find the diagonal with the most seeds
    max_count = max(diagonal_counts.values())
    if max_count == 1:
        # this means that all diagonals are unique, 
        # which may mean that the seeds are one-off hits
        return None
    max_diagonal = max(diagonal_counts, key=diagonal_counts.get)
    
    # Filter seeds to only those with the most common diagonal
    filtered_seeds = [seed for seed in seeds if (seed[1] - seed[0]) == max_diagonal]
    
    # return the first seed of the most common diagonal
    return filtered_seeds[0]
   
def extending_seed(query, read, q_pos, r_pos, q_len, r_len, dir, max_mismatch):
        """
        Extend the seed in one direction until the extension fails. 
        If dir is -1, extend left. If dir is 1, extend right.
        Allow up to max_mismatch mismatches. 
        Returns the best index (number of positions extended). 
        This is a helper function for the extend_seed function. 
        """
        # initialize score, errors, best_index, best_score, and i (iterator position)
        score = 0 
        errors = 0 
        best_index = 0
        best_score = 0 
        i = 0
        # starting seed: 
        q = q_pos + (dir * i) + dir
        r = r_pos + (dir * i) + dir
        # while extension within bounds of query and read sequence
        while 0 <= q < q_len and 0 <= r < r_len:
            # if query and read match, increment score
            if query[q] == read[r]:
                i += 1
                score += 1
            else:
                # if query and read do not match, increment errors
                i += 1
                score -= 2
                errors += 1
            # calculate new positions
            q = q_pos + (dir * i) + dir
            r = r_pos + (dir * i) + dir
            # if score is greater than best_score, update best index and score
            if score > best_score and errors < max_mismatch:
                best_index = i
                best_score = score
            if errors > max_mismatch:
                # if errors exceed max_mismatch, break the loop 
                break
        
        return (best_index) 

def extend_seed(seed, query, q_len, read, r_len, max_mismatch=1, length_threshold=15):
    """
    Extend the seed in both directions until the extension fails.
    Allow up to max_mismatch mismatches per left or right extension.
    Default is 1 mismatch. 
    If the extended seed is below the length threshold, return None. 
    Return the extended seed as a tuple: (contig_start_position, read_start_position, length)
    """
    # If length threshold larger than read length, set threshold to read length/2
    if r_len < length_threshold: 
        # set length threshold to the read length/2
        length_threshold = r_len/2 
        
    # Unpack the seed
    q_pos, r_pos = seed
    
    # Extend the seed in both directions
    if q_pos != 0 or r_pos != 0:
        left = extending_seed(query, read, q_pos, r_pos, q_len, r_len, -1, max_mismatch)
    else: 
        left = 0 
    if q_pos != (q_len-1) or r_pos != (r_len-1):
        right = extending_seed(query, read, q_pos, r_pos, q_len, r_len, 1, max_mismatch)
    else: 
        right = 0 
        
    # Calculate length of extension 
    length = left + right + 1  # +1 for the seed itself    
    # If the extended seed is below the length threshold, return None
    if length < length_threshold:
        #print(f"Extended seed length of {read} below threshold length: {length_threshold}")
        return None

    # Return the extended seed as a tuple (q_pos, r_pos, length)
    q_start = q_pos - left
    r_start = r_pos - left
    q_end = q_pos + right
    r_end = r_pos + right
    length = left + right + 1  # +1 for the seed itself

    return (r_start, r_end, q_start, q_end, length)

def align_seq_to_contig(read, read_len, contig, contig_len, contig_index, k=15, I="default", max_mismatch=1): 
    """
    Align the read seq to the contig index using a BLAST-like seed-and-extend algorithm. 
    The algorithm consists of the following steps: 
    1. Build the contig index
    2. Scan the read sequence for spaced k-mers of interval I that match the contig index
    3. Prioritize the first seed of the most common diagonal, if there's at least a double-hit
    4. Extend the seed in both directions until the extension fails
    """ 
    # Identify seed matches
    seeds = scan_read_for_seeds(read, read_len, contig_index, k, I)
    # If there are no seeds, return None
    if not seeds:
        #print(f"No seeds found for {read}")
        return None
    # Prioritize the first seed of the most common diagonal 
    top_seed = prioritize_seeds_by_diagonal(seeds)
    # If there is no top seed return None
    if not top_seed:
        #print(f"No top seed found for {read}")
        return None
    # Extend the seed in both directions
    extended_seed = extend_seed(top_seed, contig, contig_len, read, read_len, 
                                max_mismatch=max_mismatch, length_threshold=(k*2))
    # If the extended seed is None, return None
    if not extended_seed:
        #print(f"No extended seed found for {read}")
        return None
    # Return the extended seed
    return extended_seed

def reverse_complement_seed(alignment, read_len): 
    """
    Reverse the extended seed alignment information so that it matches the forward read direction.
    """
    r_start = alignment[0]
    r_end = alignment[1]
    # only reverse alignment found 
    # reverse the alignment information so that it matches read direction
    read_end = read_len - 1
    new_r_start = read_end - r_start
    new_r_end = read_end - r_end
    # add alignment_strand (-1) to the tuple and new read start/ends
    return(new_r_start, new_r_end, alignment[2], alignment[3], alignment[4], -1)

def align_read_to_contig(read, contig, contig_len, contig_index, k=15, I="default", max_mismatch=1):
    """
    Align the read sequence and its reverse complement to the contig.
    If the longest alignment is the reverse complement, update the alignment information.
    Returns: 
    tuple: (contig_start, contig_end, read_start, read_end, length, alignment_strand)
    """
    # Calculate the read length
    #print("Calculating read length...")
    read_len = len(read)
    #print(f"Read length: {read_len}")
    # Generate the reverse complement of the read
    #print("Generating reverse complement of read...")
    rev_comp_read = reverse_complement(read)
    #print(f"Reverse complement read: {rev_comp_read}")
    # Align the read to the contig
    #print("Aligning read to contig...")
    f_alignment = align_seq_to_contig(read, read_len, contig, contig_len, contig_index, k=k, I=I, max_mismatch=max_mismatch)
    #print(f"Forward alignment: {f_alignment}")
    #print("Aligning reverse complement read to contig...")
    r_alignment = align_seq_to_contig(rev_comp_read, read_len, contig, contig_len, contig_index, k=k, I=I, max_mismatch=max_mismatch)
    #print(f"Reverse alignment: {r_alignment}")
    if f_alignment is None: 
        if r_alignment is None: 
            # no alignments found 
            return(None)
        elif r_alignment is not None:
            return(reverse_complement_seed(r_alignment, read_len))
    else: 
        if r_alignment is None:
            # add alignment_strand (1) to the tuple 
            return(f_alignment[0], f_alignment[1], f_alignment[2], f_alignment[3], f_alignment[4], 1)
        else: 
            # If both alignments are found, compare their lengths and return longest
            if f_alignment[4] >= r_alignment[4]:
                # add alignment_strand (1) to the tuple 
                return(f_alignment[0], f_alignment[1], f_alignment[2], f_alignment[3], f_alignment[4], 1)
            elif f_alignment[4] < r_alignment[4]:
                return(reverse_complement_seed(r_alignment, read_len))
                
def align_reads_to_contig(reads, contig, contig_name, k=15, I="default", max_mismatch=1):
    """
    Aligns a list of reads to a contig.
    Inputs:
    -------
    reads : dict 
        a dictionary of reads (read_name: read_sequence)
    contig : str
        contig sequence
    contig_name : str
        name of the contig
    k : int
        length of k-mer
    I : int
        interval length, default is calculated using read length 
    max_mismatch : int
        maximum number of mismatches allowed, default is 1
    Outputs:
    -------
    alignment_tuple : list of tuples
        a list of tuples containing the alignment information for each read:
        (read_name, contig_name, contig_start, contig_end, read_start, read_end, length, alignment_strand)
    """
    # Initialize the alignment list
    alignment_list = []
    # Build the contig index
    contig_index = build_index(contig, k=k)
    # Iterate through the reads and align each read to the contig
    for read_name, read_seq in reads.items():
        #print(f"Aligning {read_name} to {contig_name}...")
        read_len = len(read_seq)
        alignment = align_read_to_contig(read_seq, contig, len(contig), contig_index, k=k, I=I, max_mismatch=max_mismatch)
        if alignment is not None:
            alignment_list.append((read_name, contig_name) + alignment)
    return alignment_list

def export_alignment_to_tsv(alignment_list, output_file):
    """
    Export the alignment information to a TSV file. 
    Inputs:
    -------
    alignment_list : list of tuples
        a list of tuples containing the alignment information for each read:
        (read_name, contig_name, contig_start, contig_end, read_start, read_end, length, alignment_strand)
    output_file : str
        path to the output TSV file
    """
    # Open output file for writing
    with open(output_file, "w") as file: 
        csv.register_dialect("custom", delimiter=" ", skipinitialspace=True)
        writer = csv.writer(file, dialect="custom")
        writer.writerow(("sseqid", "qseqid", "sstart", "send", "qstart", "qend", "length", "strand"))
        # Add each alignment to the file 
        for alignment in alignment_list: 
            writer.writerow(alignment)
            