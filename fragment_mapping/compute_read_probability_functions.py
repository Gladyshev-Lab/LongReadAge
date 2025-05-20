#!/usr/bin/env python3


# import libraries
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import pysam
from tqdm import tqdm
import pyranges as pr
import matplotlib.pyplot as plt
import re
import concurrent.futures


#################################
#                               #
#           Functions           #
#                               #
#################################



def process_bam_file(bam_file, detailed_markers_subset, num_cell_types, target_cell_type, debug_mode):
    """
    Processes a BAM file to extract read information and calculate probabilities each read originated from
    cell types enumerated in detailed_markers_subset.
    
    Parameters:
    --------
    bam_file (str): Path to the BAM file.
    detailed_markers_subset (pd.DataFrame): DataFrame containing detailed markers subset.
    num_cell_types (int): Number of cell types.
    target_cell_type (str): cell type of interest (e.g. myeloid), probabilities will be calculated for this cell type 
    # and all other cell types will be considered 'background'
    
    Returns:
    --------
    pd.DataFrame: DataFrame containing read probabilities and other generated metrics.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    reads = list(bam.fetch())

    # testing
    # reads = reads[0:30000]
    
    # extract fields of interest from bam file that have an 'ML' tag (meaning they have methylation probabilities)
    bam_qname = [read.query_name for read in reads if read.has_tag('ML')]
    bam_chr = [read.reference_name for read in reads if read.has_tag('ML')]
    bam_pos = [read.reference_start for read in reads if read.has_tag('ML')]
    bam_mapq = [read.mapping_quality for read in reads if read.has_tag('ML')]
    bam_seq = [read.query_sequence for read in reads if read.has_tag('ML')]
    bam_ml = [read.get_tag('ML') for read in reads if read.has_tag('ML')]
    bam_mm = [read.get_tag('MM') for read in reads if read.has_tag('ML')]
    bam_strand = [read.is_reverse for read in reads if read.has_tag('ML')]
    bam_cigar = [read.cigarstring for read in reads if read.has_tag('ML')]

    # ensure all field lengths are equal length
    assert len(bam_qname) == len(bam_chr) == len(bam_pos) == len(bam_mapq) == len(bam_seq) == len(bam_ml) == len(bam_mm) == len(bam_strand) == len(bam_cigar)

    # initialize read info arrays
    prob_target = np.zeros(len(bam_pos))
    prob_background_max = np.zeros(len(bam_pos))
    prob_target_summed = np.zeros(len(bam_pos))
    prob_background_max_summed = np.zeros(len(bam_pos))
    num_overlap = np.zeros(len(bam_pos))
    num_total_methylation_sites = np.zeros(len(bam_pos))

    # filter methylation reads: only consider reads with at least 3 CpG sites, quality scores of ≥20 (≤1% chance mapping is wrong),
    # and found on autosomes
    methylation_site_filter = [
        i for i in range(len(bam_ml)) 
        if len(bam_ml[i]) >= 3 and bam_mapq[i] >= 20 and bam_chr[i] not in ['chrX', 'chrY']
    ]

    # iterate through each read and compute probability it came from particular cell type
    for p in tqdm(methylation_site_filter, desc="Processing reads"):
        process_filtered_read(
            p, bam_chr, bam_pos, bam_seq, bam_ml, bam_mm, bam_strand, bam_cigar, 
            detailed_markers_subset, prob_target, prob_background_max, 
            prob_target_summed, prob_background_max_summed, num_overlap, 
            num_total_methylation_sites, num_cell_types, target_cell_type, debug_mode=debug_mode
        )
        
    assert(len(prob_target) == len(prob_background_max) == len(prob_target_summed) == len(prob_background_max_summed) == len(num_total_methylation_sites) == len(num_overlap) == len(bam_qname))

    # create bam file info df
    bam_file_info = pd.DataFrame({
        'prob_target': prob_target, # multiplicative probability of target
        'prob_background_max': prob_background_max, # max multiplicative probability of background
        'prob_target_summed': prob_target_summed, # sum probability of target
        'prob_background_max_summed': prob_background_max_summed, # max sum probability of background
        'num_total_methylation_sites': num_total_methylation_sites, # number of total methylation sites
        'num_overlap': num_overlap, # number of methylation sites overlapping wth DMRs
        'qnames': bam_qname # list of qnames 
    })
    
    return bam_file_info


def parse_cigar(cigar):
    """
    Parses a CIGAR string into a list of tuples representing the operations and their lengths.
    
    Parameters:
    -----------
    cigar : str
        The CIGAR string representing the alignment of the read to the reference sequence.
        The string is composed of pairs of numbers and characters indicating the operation
        and the length of the operation (e.g., "10M5I3D").
    
    Returns:
    --------
    list of tuples
        Each tuple contains:
            - int: Length of the operation
            - str: Operation character (one of 'M', 'I', 'D', 'X', 'S', '=')
    
    Example:
    --------
    >>> parse_cigar("10M5I3D")
    [(10, 'M'), (5, 'I'), (3, 'D')]
    """

    # This regex will split the CIGAR string into tuples of (number, operation)
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDXS=])', cigar)]


def adjust_positions(cigar, meth_c_positions, first_pos):
    """
    Adjusts the positions of methylated cytosines based on the CIGAR string of an aligned read so they match the reference genome.
    
    Parameters:
    -----------
    cigar : str
        CIGAR string representing the alignment of the read to the reference.
    meth_c_positions : dict
        Dictionary where keys are positions of methylated cytosines in the read (0-based indexing), 
        and values are the methylation probabilities at those positions.
    first_pos : int
        The starting position of the alignment on the reference sequence (1-based indexing).

    Returns:
    --------
    list of tuples
        Each tuple contains the adjusted position of a cytosine on the reference 
        and its methylation probability.
    """
    
    ref_pos = first_pos
    read_pos = 0
    
    # list of tuples in which each tuple contains (adjusted cytosine position, methylation probability at that cytosine)
    adjusted_positions = []
    
    for length, op in parse_cigar(cigar):
        
        if op == '=' or op == 'X':  # Alignment match or mismatch
            for i in range(length):
                
                #print(read_pos + i)
                if (read_pos + i) in meth_c_positions:
                    adjusted_positions.append((ref_pos + i, meth_c_positions[read_pos + i]))
            ref_pos += length
            read_pos += length
        elif op == 'I':  # Insertion to the reference
            read_pos += length
        elif op == 'D':  # Deletion from the reference
            ref_pos += length
        elif op == 'S':  # Soft clipping (not in the reference)
            read_pos += length
        elif op == 'H':  # Hard clipping (not in the reference)
            pass  # Hard clipping does not consume read bases
        else:
            raise ValueError(f"Unknown CIGAR operation: {op}")

    return adjusted_positions


def process_filtered_read(i, bam_chr, bam_pos, bam_seq, bam_ml, bam_mm, bam_strand, bam_cigar, 
                          detailed_markers_subset, prob_target, prob_background_max, 
                          prob_target_summed, prob_background_max_summed, num_overlap, 
                          num_total_methylation_sites, num_cell_types, target_cell_type, debug_mode):
    """
    Processes a single filtered read to calculate probabilities read originated from each cell type
    listed in DMR table.
    
    Parameters:
    --------
    i (int): Index of the read.
    bam_chr (list): List of chromosome names.
    bam_pos (list): List of reference start positions (0-based).
    bam_seq (list): List of query sequences.
    bam_ml (list): List of methylation levels.
    bam_mm (list): List of mismatch positions.
    bam_strand (list): List of strand orientations.
    bam_cigar (list): List of CIGAR strings.
    detailed_markers_subset (pd.DataFrame): DataFrame containing target-cell-type-specific DMRs.
    prob_target (np.ndarray): Array to store multiplicative target probabilities.
    prob_background_max (np.ndarray): Array to store background max multiplicative probabilities.
    prob_target_summed (np.ndarray): Array to store summed target probabilities.
    prob_background_max_summed (np.ndarray): Array to store summed background max probabilities.
    num_overlap (np.ndarray): Array to store number of overlaps.
    num_total_methylation_sites (np.ndarray): Array to store total number of methylation sites.
    num_cell_types (int): Number of cell types.
    target_cell_type (str): cell type of interest (e.g. myeloid), probabilities will be calculated for this cell type 
    # and all other cell types will be considered 'background'.
    debug_mode (bool) : Whether to print additional output and perform additional assertions.

    """
    starting_chr = bam_chr[i]
    starting_position = bam_pos[i]
    seq = bam_seq[i]
    ml = bam_ml[i]
    prob_meth = (np.array(ml) + 0.5) / 256
    mm = [int(x.replace(";", "")) for x in bam_mm[i].split(",")[1:]]
    strand = bam_strand[i]
    cigar = bam_cigar[i]
    
    assert(len(mm) == len(prob_meth))
    
    if debug_mode:
        print('\nStarting chromosome:', starting_chr)
        print('\nStarting position:', starting_position)
        print('\nMM tags:', mm)
        print('\nMethylation probabilities:', prob_meth)
        print('\nStrand:', strand)
        print('\nCIGAR:', cigar)


    # forward strands
    if not strand:
        idx_C_with_methyl_info = get_cytosine_indices(seq, mm)
    

    # reverse strands
    else:
        
        # generate the reverse complement
        seq_rev = seq[::-1].translate(str.maketrans('ATGC', 'TACG'))
        idx_C_with_methyl_info = get_cytosine_indices(seq_rev, mm)
        
        
        # revert indexes to line up with original sequence
        idx_C_with_methyl_info = len(seq) - np.array(idx_C_with_methyl_info) - 2

    
    if debug_mode:
        # determine what nucleotides are present at index C with methylation info
        seq_indexed_C = [seq[i] for i in range(len(seq)) if i in idx_C_with_methyl_info]

        # assert that all nucleotides at index C are indeed C
        assert(all(n == 'C' for n in seq_indexed_C))

        # determine what nucleotides are present at index C with methylation info + 1
        seq_indexed_C_plus_one = [seq[i + 1] for i in range(len(seq)) if i in idx_C_with_methyl_info]

        # assert that all nucleotides at index C +1 are indeed G
        assert(all(n == 'G' for n in seq_indexed_C_plus_one))

        print("Assertions passed!")
    
    
    # create dictionary of methylated cytosine positions and methylation probabilities
    positions = idx_C_with_methyl_info
    statuses = prob_meth

    meth_c_positions = {pos: status for pos, status in zip(positions, statuses)}

    first_pos = starting_position + 1  # Starting position in the reference genome (1-based)

    # correct cytosine positions based on CIGAR string (get list of tuples [(adjusted position, methylation probability), ...])
    adjusted_positions = (adjust_positions(cigar, meth_c_positions, first_pos))

    if debug_mode:
        print('\nOriginal positions:', meth_c_positions.keys())
        print('\nAdjusted positions:', [pos for pos, meth in adjusted_positions])
        print('\nAdjusted positions (minus starting position):', np.array([pos for pos, meth in adjusted_positions]) - first_pos)
        print('\nNumber of adjusted positions:', len(adjusted_positions))
        

    new_pos = np.array([pos for pos, meth in adjusted_positions])
    new_meth = np.array([meth for pos, meth in adjusted_positions])

    
    # create dataframe with all CpG sites in read (indexing is 1-based, C is first coordinate, G is second)
    read_cpg_sites = pd.DataFrame({
        'Chromosome': [starting_chr] * len(new_pos),
        'Start': new_pos,
        'End': new_pos + 1,
        'read_prob_meth': new_meth
    })
    

    # subset DMRs to include only chromosome matching read, merge with read
    overlapping_cpgs = pd.merge(
        detailed_markers_subset[detailed_markers_subset['Chromosome'] == starting_chr],
        read_cpg_sites, on=['Chromosome', 'Start', 'End'], how='inner'
    )

    
    if len(overlapping_cpgs) >= 3:
        
        prob_target[i], prob_background_max[i], prob_target_summed[i], prob_background_max_summed[i] = generate_probabilities(overlapping_cpgs, num_cell_types, target_cell_type)
        
        num_overlap[i] = len(overlapping_cpgs)
        num_total_methylation_sites[i] = len(ml)
        
    else:
        prob_target[i] = 0
        prob_background_max[i] = 0
        prob_target_summed[i] = 0
        prob_background_max_summed[i] = 0
        num_overlap[i] = 0
        num_total_methylation_sites[i] = 0


        
def get_cytosine_indices(seq, mm):
    """
    Determines the indices of cytosines in a sequence of nucleotides with methylation information based on mm positions.
    
    Parameters:
    --------
    seq (str): Sequence of nucleotides.
    mm (list): List of mm positions ("2, 3" means methylated cytosines are the 3rd cytosine and 7th cytosine).
    
    Returns:
    --------
    list: List of indices of cytosines with methylation information.
    """
    
    seq = list(seq)
    idx_C = np.where(np.array(seq) == 'C')[0]
    cumulative_mm = np.cumsum(mm)
    return idx_C[cumulative_mm + np.arange(len(mm))]


def generate_probabilities(overlapping_cpgs, num_cell_types, target_cell_type):
    """
    Generates probabilities a read originated from a target vs background tissue based on overlapping CpG sites
    between read and reference.
    
    Parameters:
    --------
    overlapping_cpgs (pd.DataFrame): DataFrame of overlapping CpG sites.
    num_cell_types (int): Number of cell types.
    target_cell_type (str): cell type of interest (e.g. myeloid), probabilities will be calculated for this cell type 
                            and all other cell types will be considered 'background'.

    """
    
    # convert read methylation to binary values (single-molecule, so only 1 or 0 are accurate)
    overlapping_cpgs.read_prob_meth = np.where(overlapping_cpgs.read_prob_meth >= 0.5, 1, 0)
    
    # isolate reference cell type methylation values
    cell_type_meth = overlapping_cpgs.iloc[:, 3:(3 + num_cell_types)]
    
    # isolate read methylation values
    read_meth = overlapping_cpgs.loc[:, 'read_prob_meth']

    prob_all = np.zeros(num_cell_types)
    prob_all_summed = np.zeros(num_cell_types)

    # compute multiplicative and summed probability read originated from each cell type in reference
    for y in range(num_cell_types):
        cell_type_meth_i = cell_type_meth.iloc[:, y].values
        prob_all[y] = np.prod(
            np.power(cell_type_meth_i, read_meth) * 
            np.power(1 - cell_type_meth_i, 1 - read_meth)
        )
        prob_all_summed[y] = np.sum(1 - np.abs(read_meth - cell_type_meth_i))

    target_index = np.where(cell_type_meth.columns == target_cell_type)[0][0]
    non_target_index = np.where(cell_type_meth.columns != target_cell_type)[0]

    prob_target = prob_all[target_index]
    prob_background_max = np.max(prob_all[non_target_index])
    prob_target_summed = prob_all_summed[target_index]
    prob_background_max_summed = np.max(prob_all_summed[non_target_index])

    return (prob_target, prob_background_max, prob_target_summed, prob_background_max_summed)

    

    
