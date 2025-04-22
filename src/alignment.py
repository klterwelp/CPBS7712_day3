# BLAST like alignment algorithm 
# based on https://community.gep.wustl.edu/repository/course_materials_WU/annotation/Smith_Waterman_to_BLAST_notes.pdf

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

# STEP THREE: EXTEND SEEDS
# extend the seeds in both directions until the extension fails
# given a seed, the sequences of the read and query, and a maximum number of mismatches allowed: 
# extend the seed in both directions until the extension fails
# return the extended seed as a tuple: (contig_start_position, read_start_position, length) 

# STEP FOUR: PRIORITIZE SEEDS
# filter the seeds that are below a certain length threshold
# given a list of seeds and a length threshold:
# 1. filter the seeds that are below the length threshold
# prioritize the seeds based on the following criteria:
# 1. length of the seed
# if longest seed is the same length as the read, return it as a tuple: (contig_start_position, read_start_position, length)
# else return the sorted list of prioritized seeds as a list of tuples: (contig_start_position, read_start_position, length)


# STEP FIVE: ALIGN SEEDS WITH STRIPED SMITH-WATERMAN
# extend the seeds in both directions using dynamic programming approach
# given a seed, the sequences of the read and query: 
# 1. create a scoring matrix of size (query_length + 1) x (read_length + 1)
# 2. fill the scoring matrix only in the banded areas using the following scoring scheme:
#    - match: +1
#    - mismatch: -1
#    - gap: -2
# 3. backtrack the scoring matrix to find the optimal alignment



def build_index(seq, k=15):
    """
    Build an index for the given sequence.
    The index is a dictionary where the keys are k-mers and the values are lists of positions.
    """
    index = {}
    k = k  # Length of k-mer
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer not in index:
            index[kmer] = []
        index[kmer].append(i)
    return index

def scan_read_for_seeds(read_seq, contig_index, k = 15, I = "default"):
    """
    Scan the read sequence for spaced k-mers of interval I that match the contig index.
    Calculate the interval length I using the formula from Bowtie2:
    I(x) = max(1, floor(1 + 1.15 * √x))
    where x is the length of the read sequence.
    Else, use the provided interval length I. 
    """
    # Calculate the interval length I using the formula from Bowtie2:
    len_seq = len(read_seq)
    if I == "default":
        I = max(1, int(1 + 1.15 * (len_seq ** 0.5)))
    
    # Start empty seeds list
    seeds = []
    # Scan the read sequence at intervals of length I for matching k-mers 
    for r_pos in range(0, len(read_seq) - k + 1, I):
        spaced_kmer = read_seq[r_pos:r_pos + k]
        if spaced_kmer in contig_index:
            for q_pos in contig_index[spaced_kmer]:
                seeds.append((q_pos, r_pos))
    return seeds

def extend_seed(seed, query, read, max_mismatch=1):
    """
    Extend the seed in both directions until the extension fails.
    Allow up to max_mismatch mismatches per left or right extension.
    Default is 1 mismatch, same as Bowtie2. 
    Return the extended seed as a tuple: (contig_start_position, read_start_position, length)
    """
    # Unpack the seed
    q_pos, r_pos = seed
    # Get the lengths of the query and read sequences
    q_len = len(query)
    r_len = len(read)
    # Extend left
    left = 0
    mismatches_left = 0
    while q_pos - left - 1 >= 0 and r_pos - left - 1 >= 0:
        if query[q_pos - left - 1] != read[r_pos - left - 1]:
            mismatches_left += 1
            if mismatches_left > max_mismatch:
                break
        left += 1

    # Extend right
    right = 0
    mismatches_right = 0
    while q_pos + right < q_len and r_pos + right < r_len:
        if query[q_pos + right] != read[r_pos + right]:
            mismatches_right += 1
            if mismatches_right > max_mismatch:
                break
        right += 1
    
    # Return the extended seed as a tuple (q_pos, r_pos, length)
    q_pos = q_pos - left
    r_pos = r_pos - left
    length = left + right
    
    return (q_pos, r_pos, length)

    
def prioritize_seeds_by_length(seeds, length_threshold=30):
    """
    Sort the seeds by their length in descending order.
    Filter the seeds that are below the length threshold.
    If the longest seed is the same length as the read, return it as a tuple: (contig_start_position, read_start_position, length)
    Else return the sorted list of prioritized seeds as a list of tuples: (contig_start_position, read_start_position, length)
    """
    # Filter the seeds that are below the length threshold
    seeds = [seed for seed in seeds if seed[2] >= length_threshold]
    
    # If there's only one seed, return it as a tuple: (contig_start_position, read_start_position, length)
    if len(seeds) == 1:
        return seeds[0]
    
    # If there are no seeds above length threshold, return None
    if not seeds:
        return None

    # Sort seeds by their length in descending order
    sorted_seeds = sorted(seeds, key=lambda x: x[2], reverse=True)
    
    # If the longest seed is the same length as the read, return it as a tuple: (contig_start_position, read_start_position, length)
    if sorted_seeds and sorted_seeds[0][2] == len(read):
        return sorted_seeds[0]
    
    # Else return the sorted list of prioritized seeds as a list of tuples: (contig_start_position, read_start_position, length)
    return sorted(seeds, key=lambda x: x[2], reverse=True)    

def banded_smith_waterman_alignment(query, read, seed, bandwidth = 10):
    