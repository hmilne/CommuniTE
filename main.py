##########
# CommuniTE- A program to identify transposable element subfamilies using community detection
##########


###### Lots of imports ######

import sys
import os
import argparse
import itertools
from datetime import *
from math import sqrt, log
import random
from igraph import *
import time
import timeit
from Seq import *
import scipy
import json
from scipy import linalg
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import seaborn as sns
import logging






##### Parsing command line inputs #######

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', dest='outfile', nargs=1, help='input file in FASTA format') # stores one argument in outfile
    parser.add_argument('-i', '--infile', dest='infile', nargs=1, help='output file tag') # stores one argument in infile
    # parser.add_argument('-v', dest='verbose', action='store_true', help='show verbose output')

    parser.add_argument('-a', dest='alignment', action='store', nargs=1, help='type of alignment to perform', default=['clu'])
    parser.add_argument('-fragment', dest='fragment_threshold',action='store',nargs=1, help='sequences shorter than this percent of longest sequence will be removed', default=[0.8])
    parser.add_argument('-edge', dest='edge_threshold',action='store',nargs=1, help='edges less than this will be removed', default=[0.75])
    parser.add_argument('-alpha', dest='alpha',action='store',nargs=1, help='relative weight of gap scoring', default=[0.5])
    # parser.add_argument('-dangling', dest='dangling_threshold',action='store',nargs=1, help='nodes beyond this threshold will be considered uncharacterized', default=30)
    parser.add_argument('-c', dest='cd_algorithm', action='store',nargs=1, help='community detection algorithm', default=['fg'])
    parser.add_argument('-s', dest='scoring_algorithm', action='store',nargs=1, help='sequence scoring algorithm',default=['mcic'])
    parser.add_argument('--readalignment', dest='read_alignment', action='store', nargs=1, help='Use an existing fasta alignment', default=[None])
    # options for -a: 'clu' to run Clustal Omega multiple alignment, or 'affine' to run all pairwise affine gap alignments.
    # options for -c enumerated below
    # options for -s: 'fs'/'fifthstate', 'MCIC', 'md'/'missingdata'

    parser.add_argument('--force', dest='force', action='store_true', help='force overwrite of existing files')
    parser.add_argument('--many', dest='many', action='store_true', help='outputs an additional one-line TSV to be concatenated')
    parser.add_argument('-r','--readexisting', dest='read', action='store_true', help='read existing files of given name rather than run segments again')
    args = parser.parse_args()

##### Some global variable thresholds ######

# No idea what these will be quite yet.
fragment_threshold = float(args.fragment_threshold[0]) #This could be a set length or could be rewritten in terms of a percentage of length of the longest element. Say, 10%.
edge_threshold = float(args.edge_threshold[0])
# dangling_threshold = args.dangling_threshold
alpha = args.alpha[0] # Ranges from 0 to 1. 0: All weight on nucleotide scoring. 1: All weight on gap scoring.
if alpha == 'freq' or alpha == 'None':
    beta = 'freq'
else:
    alpha = float(alpha)
    beta = 1 - alpha


# gapscore = "missing data"
# ntscore = None
G = Graph()


#Debugging log
LOG_FILENAME = args.outfile[0] + "_errors.out"
logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG)


###### Parsing input fasta file ########
class Seq:
    def __init__(self, num, seqID):
        self.num = num
        self.ID = seqID
        # self.sequence = sequence
        self.attributes = {}
    def setAttribute(self, key, value):
        self.attributes[key] = value
    def length(self):
        return len(self.sequence)
    def setSequence(self, sequence):
        self.sequence = sequence


def check_file_path(filename):
    '''
    Checks the existence of the directories in the outfile path.
    If the directories do not exist, they will be created.
    '''
    filepath = os.path.dirname(filename)
    if filepath == '':
        filepath = './'

    if os.path.exists(filepath):
        print("Directory " + filepath + " exists!")
    else:
        print("Creating directory " + filepath + "...")
        os.makedirs(filepath)
        print("Done.")

def read_existing(filename):
    '''
    Returns True if it's OK to read an existing file
    of the given name rather than to recompute it, and
    False otherwise.
    '''
    files = os.listdir(".")
    if args.read == True and filename in files:
        print("Using existing file " + filename + "...")
        return True

    else:
        return False


def parse_fasta(filename, add_name_ext):
    '''
    This function parses the input fasta file and returns a list of Seq datatypes (see Seq.py).
    The add_name_ext input is a boolean that determines whether a unique identifier is added to
    the name of each sequence.
    '''
    print("Reading sequences from", filename)
    sequences = [] # this will be of type seq_ID (string) : {dict of sequence information}
    count = -1
    with open(filename, 'r') as infile:
        for line in infile:
            if line[0] == '>': # the first character of the line
                count += 1
                words_in_line = line.split() #split the line into words
                seqID = words_in_line[0][1:]
                if add_name_ext:
                    seqID += "_" + str(count)
                sequences.append(Seq(count, seqID))

                for j in range(1,len(words_in_line)):
                    components = words_in_line[j].split('=')
                    sequences[count].setAttribute(components[0], components[1])
                    # adds each trait on the line as a key:value pair to the inner dictionary
                sequences[count].sequence = '' # initialize seq key with empty string value
            else:
                #if the line doesn't begin with > then append the contents to the previous seq's sequence

                while len(line) > 0 and line[-1] in [" ","\n"]:
                    line = line[:-1]
                sequences[count].sequence += line[:-1]


    print("Done.")
    print("Read", len(sequences), "sequences.")

    return sequences

def overwrite(filename):
    '''
    Returns True if it's OK to write a file with the given
    name (i.e. if no such file exists, or if --force is set
    to true). Otherwise it returns false.
    '''
    files = os.listdir(".")
    if args.force == True:
        return True
    if filename not in files:
        return True
    else:
        print("Refusing to overwrite existing file.")
        print("Use --force to overwrite existing files with the same name.")
        print("Use --readexisting to read existing files with the same name.")
        sys.exit()
        # return False

def write_fasta(seqlist, filename):
    filename = filename + '.fasta'

    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        with open(filename, 'w') as f:
            for seq in seqlist:
                f.write('>')
                f.write(seq.ID + " ")
                for attribute in seq.attributes:
                    f.write(attribute)
                    f.write("=")
                    f.write(seq.attributes[attribute])
                    f.write(" ")

                f.write('\n')
                f.write(seq.sequence)
                f.write('\n')
        return filename



def remove_fragments(seq_list):
    '''
    Takes a seq_list and a threshold value as input, and outputs a copy of the
    list without the sequences at or below the threshold value. I made
    a copy because dictionaries don't allow changes to length while iterating.
    '''
    global fragment_threshold
    maxlen = 0
    maxseq = seq_list[0]
    for seq in seq_list:
        if len(seq.sequence) > maxlen:
            maxlen = len(seq.sequence)
            maxseq = seq
    with open("maxseq42616.txt","w") as f:
        f.write(str(maxlen))
        f.write("\n")
        f.write(str(maxseq.sequence))
        f.write(str(maxseq.ID))
    print(fragment_threshold)
    fraglen = int( maxlen * fragment_threshold )

    print("Removing sequences less than " + str(fraglen) + " (" + str(fragment_threshold * 100) + " % of longest sequence).")

    new_seq_list = [s for s in seq_list if len(s.sequence) >= fraglen]

    print(str(len(new_seq_list)) + " sequences remaining.")
    if len(new_seq_list) == 1:
        print("Only one sequence remaining. Please try again with more data or a lower fragment threshold.")
        sys.exit()
    return new_seq_list

def fragment(seq, threshold):
    if len(seq) <= threshold:
        return True
    else:
        return False



def lengths(seq_list):
    '''
    Finds the average length of all sequences in the dictionary
    '''
    lens = []
    for i in range(len(seq_list)):
        lens.append(len(seq_list[i].sequence))

    avg = sum(lens)//len(lens)
    print("Average length:", avg)
    return lens

def run_clustalomega(force, fasta, outfile):
    '''
    Runs Clustal Omega from the command line, passing 'force' argument.
    Takes infile fasta as input and writes aligned fasta to outfile.
    '''
    outfilename = outfile + ".afasta"

    # Calls the following command with clustal omega.
    if read_existing(outfilename):
        print("Done.")
        return outfilename

    elif args.read_alignment[0] != None:
        return args.read_alignment[0]

    elif overwrite(outfilename):
        print("Performing MSA with Clustal Omega...")
        command = 'clustalo ' + '-i ' + fasta + ' -o ' + outfilename + ' --force'
        print("Writing MSA to " + outfilename + "...")
        os.system(command)
        print("Done.")
        return outfilename



def print_seqs_from_list(seq_list):
    '''
    Just prints all the sequences in a seq dict
    '''
    for i in range(len(seq_list)):
        print(seq_list[i].sequence)
    # Seems to be working
    # Mostly just so I can see the alignments



##################################
####### Alignment Scoring ########
##################################

# Some general use functions for scoring

def remove_gaps(seq1, seq2):
    '''
    Returns the sequences without gapped locations.
    '''
    new1 = ''
    new2 = ''
    if len(seq1) == 0: # If there's no intersection, return a '-'
        return '', '' #
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            new1 += seq1[i]
            new2 += seq2[i]

    # TODO: what if there are no characters left at the end?
    # Also: should I weight the gap/nt parameter by % of gaps vs. nts in the data?
    return new1, new2

def seqstart(seq):
    '''
    Finds the start of a sequence, not including any 5' gap.
    '''
    start = 0
    while seq[start] == '-':
        start += 1
    return start

def seqend(seq):
    '''
    Finds the end of a sequence, not including any 3' gap.
    '''
    end = len(seq)-1
    while seq[end] == '-':
        end -= 1
    return end

def intersection(seq1,seq2):
    '''
    Finds the intersection between two sequences, i.e. the amount
    and position that two sequences overlap. For example,

    AACT----CTG
    ACCTCAG-CTC

    have an intersection of 11, or [0,10] inclusive,

    ---AACT---
    ACTTAGTGTA

    have an intersection of 4, or [3,6] inclusive, and

    ---ACCTCGG
    GTGACCT---

    have an intersection of 4, or [3,6] inclusive.
    '''
    # Should this return the length of the intersection, the position of the intersection,
    # or the actual slices of the sequences themselves where they overlap?
    seq1start = seqstart(seq1)
    seq2start = seqstart(seq2)
    seq1end = seqend(seq1)
    seq2end = seqend(seq2)

    positions = [max(seq1start,seq2start),min(seq1end,seq2end)]
    return seq1[positions[0]:positions[1]+1], seq2[positions[0]:positions[1]+1]

def printTable(table):
    '''
    Just prints a table for easier reading
    '''
    for i in range(len(table)):
        print(table[i])


##### Jukes-Cantor Model #########

def dummy_score(seq1, seq2):
    '''
    A dummy scoring function to allow me to keep working further
    down the pipeline. It gives the similarity a +1 if they share
    the same character in a location, and a 0 else.
    It takes trimmed (sequence intersections only) sequences and 
    normalizes by the length of that intersection
    '''
    if len(seq1) == 0: # If there's no intersection, return a '-'
        return '-'
    score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score += 1
    score = score/len(seq1) # seq1 and seq2 should be same length
    return score


def p_distance(seq1, seq2):
    '''
    Calculates the p-distance
    p = n/L
    where n is the number of sites that differ and L is the number
    of sites over which they were compared, excluding gap sites.
    '''
    p = 0 # Proportion of sites that differ (num sites that differ)/(num compared sites)
    sites = 0 # How many sites were compared? Excludes gaps

    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            sites += 1 # Comparing sites
            if seq1[i] != seq2[i]:
                # Not the same- dist goes up by 1
                p += 1

    p = p/sites # seq1 and seq2 should be same length
    return p


def jukes_cantor(seq1, seq2):
    '''
    Finds the Jukes Cantor distance between two sequences
    (Jukes & Cantor, 1969), for nt scoring of MSA.
    '''
    if len(seq1) == 0 or len(seq2) == 0:
        return '-'
    # NOTE: Need to find a way to normalize this so that the 
    # lowest distances become the highest edge weights...
    p = p_distance(seq1, seq2)
    D = (-0.75) * ( math.log(1-(4/3)*p) )
    # Formula for Jukes-Cantor distance
    return D


######### GTR Model ##########

def GTR_count_mat(seq1,seq2):
    counts = [[0 for i in range(0,4)] for i in range(0,4)]
    nt = {'A':0, 'G':1, 'C':2, 'T':3}
    # across- string2, down- string1
    # counts[seq1char, seq2char]
    # print(seq1)
    # print(seq2)
    for i in range(len(seq1)):
        nt1 = nt[seq1[i]]
        nt2 = nt[seq2[i]]
        counts[nt1][nt2] += 1

    # print(counts)
    return counts


def GTR_count_mat_5thstate(seq1,seq2):
    counts = [[0 for i in range(0,5)] for i in range(0,5)]
    nt = {'A':0, 'G':1, 'C':2, 'T':3, '-':4}
    # across- string2, down- string1
    # counts[seq1char, seq2char]

    for i in range(len(seq1)):
        nt1 = nt[seq1[i]]
        nt2 = nt[seq2[i]]
        counts[nt1][nt2] += 1

    return counts

def transpose(matrix):
    '''
    Takes a matrix as input and returns the transpose of the input matrix.
    '''
    new = [[matrix[i][j] for i in range(len(matrix))] for j in range(len(matrix[0]))]

    return new

def avg_countmat(matrix):
    '''
    Takes a count matrix as input and averages the values with
    those in the transpose of the count matrix. This will only
    work for matrices with symmetrical dimensions (ie 3x3, 4x4).
    '''

    trans = transpose(matrix)
    for i in range(len(trans)):
        for j in range(len(trans[0])):
            trans[i][j] = (trans[i][j] + matrix[i][j])/2

    return trans

def column_sum(matrix, i):
    total = 0
    for row in matrix:
        total += row[i]
    return total

def remove_col(i, matrix):
    '''
    Removes row and column i from the input matrix. If a 4x4 matrix is given as input,
    a 3x3 matrix will be returned. This is to remove columns that sum to zero when calculating
    Phat below.
    '''

    del matrix[i] # Removes row i
    for row in matrix:
        del row[i]

    return matrix

def calc_Phat(counts):
    '''
    Takes an averaged count matrix as input and divides each value by
    the sum of its column to produce
    the estimated P matrix. Returns the P-hat matrix and a list of the
    frequencies of each nucleotide.
    '''
    total = 0
    # rowsums = [sum(counts[i]) for i in range(len(counts))]
    colsums = [column_sum(counts, i) for i in range(len(counts))]
    # print("Colsums", colsums)

    # Check if any colsum is zero
    i = 0
    while i < len(colsums):
        if colsums[i] == 0:
            del colsums[i] # Delete that value from colsums
            remove_col(i, counts) # Delete that row/col from the matrix
        else:
            i += 1

    for i in range(len(counts)):
        for j in range(len(counts[0])):
            total += counts[i][j]
            # print(i, j, counts[i][j], colsums[j]) # i is row and j is column
            counts[i][j] = counts[i][j]/colsums[j]

    freqs = [(colsums[i]/total) for i in range(len(colsums))]

    # print(counts)
    # print(freqs)
    return counts, freqs

def calc_logPhat(Phat):
    '''
    Uses SciPy package to compute the matrix logarithm of the input matrix.
    TODO: Write this function myself, remove SciPy dependency.
    '''
    return scipy.linalg.logm(Phat)

def distance_formula(num):
    d = math.sqrt( ((num.real)**2) + ((num.imag)**2) )
    return d
    
def calc_t_hat(logPhat, freqs):
    '''
    Calculates an estimate of t from the logPhat matrix and nt frequencies.
    t-hat = -trace(logPhat*D)
    '''

    t = 0
    for i in range(len(logPhat)):
        t += (logPhat[i][i] * freqs[i]) # Only adding diagonal, all the rest are zeros.

    if t.imag != 0:
        # If the time is a complex number, the sequences are quite different.
        # So we'll return 1 to say they are divergent beyond our capability of knowing.
        t = distance_formula(t)
        return t

    elif t == 0:
        # We don't want to return -0 because it's sloppy.
        return t

    return t*-1

def GTR(seq1, seq2):
    '''
    1. Create a count matrix of all sites i,j where i=A,j=A, i=A,j=G, etc. XXX
    2. Average the non-diagonal elements with the transpose of the count matrix. XXX
    3. Freqs of A,G,C,T by sum of column / total count XXX
    4. Find P-hat by elements of the averaged count matrix divided by the column sum. XXX Does row vs. column sum matter?
    5. Take matrix log(P-hat) XXX
    6. Diagonal matrix D-hat: freq(A), etc along diagonal
    7. Calculate t-hat = -trace(P-hat*D-hat) (Which is the sum of the product of each diagonal position)
    '''
    if len(seq1) == 0 or len(seq2) == 0:
        return '-'

    counts = GTR_count_mat(seq1, seq2)

    avgcounts = avg_countmat(counts)
    Phat, pivals = calc_Phat(avgcounts)
    logPhat = calc_logPhat(Phat)

    t = calc_t_hat(logPhat, pivals)

    return t

def fifth_state_GTR(seq1, seq2):
    '''
    Performs the same calculations as GTR, but with 5x5 matrix, 
    treating '-' as a 5th nucleotide. This assumes equal chance
    of insertion or deletion, but also provides a natural weighting
    by the frequency of occurence of gaps.
    '''

    # TEST THIS
    if len(seq1) == 0 or len(seq2) == 0:
        return '-'

    counts = GTR_count_mat_5thstate(seq1, seq2)
    avgcounts = avg_countmat(counts)
    Phat, pivals = calc_Phat(avgcounts)
    logPhat = calc_logPhat(Phat)

    t = calc_t_hat(logPhat, pivals)

    # printTable(counts)
    # printTable(avgcounts)
    # printTable(Phat)
    # printTable(logPhat)
    # print(t)
    return t

#####################################

############## MCIC Gap Scoring ###############

def collapse(binseq1,binseq2):
    '''
    A helper function for MCIC gap scoring, this function takes
    sequences converted into lists of binary, and collapses 
    regions that are the same, returning the collapsed bin list.
    Ex:
    111101111101 ->     110101
    110001111101 ->     100101
    '''
    coll1 = [binseq1[0]] # Start by adding the first element
    coll2 = [binseq2[0]] # which repr. the first collapsed block

    i = 1
    while i < len(binseq1):
        # Loop through the entire sequence

        mergestart = i
        mergeend = i

        while mergeend < len(binseq1) and binseq1[mergeend] == binseq1[mergeend-1] and binseq2[mergeend] == binseq2[mergeend-1]:
            # Find the end of the block
            mergeend += 1 

        if mergeend < len(binseq1):
        # Once one same block has ended, append the first column of the next block,
        # unless we've gone over the edge of the sequence.
            coll1.append(binseq1[mergeend])
            coll2.append(binseq2[mergeend])
        i = mergeend + 1

    return coll1, coll2

def remove_same_columns(binseq1, binseq2):
    new1 = []
    new2 = []
    for i in range(len(binseq1)):
        if binseq1[i] != binseq2[i]:
            new1.append(binseq1[i])
            new2.append(binseq2[i])
    return new1, new2

def MCIC_scoring(seq1, seq2):
    '''
    Method by Muller et al (2006) which counts the minimum number of
    transitions that must occur with regard to gap characters in order
    to convert one sequence to another.
    '''
    # NOTE: This represents a DISTANCE, but what I really want is the
    # opposite. How can I normalize and invert this distance value to
    # produce values from 0-1 where 1 is most similar and 0 is least?

    
    # Convert to binary such that 0 repr. a gap and 1 repr. a nt
    if len(seq1) == 0: # If there's no intersection, return a '-'
        return '-'
    else:
        binseq1 = [0 if seq1[i]=='-' else 1 for i in range(len(seq1)) ]
        binseq2 = [0 if seq2[i]=='-' else 1 for i in range(len(seq2)) ]

        coll1, coll2 = collapse(binseq1, binseq2)

        final1, final2 = remove_same_columns(coll1, coll2)

        steps = len(final1)

        return steps

###################################
###################################


###### Reconcile and Normalize Scores #########

def normalize(matrix):
    xmax = 0
    xmin = 10000
    # Find the min and max
    # print(matrix)
    for row in matrix:
        for item in row:
            if item != '-':
                if item > xmax:
                    xmax = item
                if item < xmin:
                    xmin = item
    # Now perform normalization
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] != '-':
                if matrix[i][j] == 0:
                    matrix[i][j] = (matrix[i][j] - xmin)/(xmax - xmin)
                else:
                    matrix[i][j] = (matrix[i][j] - xmin)/(xmax - xmin)                    
    return matrix

def adjacency_matrix(gap_matrix, nt_matrix):

    '''
    Takes the raw Gap Distance Matrix and the raw Nucleotide
    Distance Matrix, then combines them with the equation

    A = 1- (aG + bN)

    where A is the adjacency matrix returned, a is the gap score weight,
    b = 1 - a is the nt score weight, G is the gap distance matrix, and N
    is the nucleotide distance matrix. Represents a weighted average.
    '''
    # Just to get the dimensions
    D = [[0 for i in range(len(gap_matrix[0]))] for j in range(len(gap_matrix))]

    global alpha # alpha is the gap weight
    global beta

    # At this stage, A is actually still the distance matrix (but G and N have
    # been normalized to be combined.


    for i in range(len(D)):
        for j in range(len(D[i])):
            if gap_matrix[i][j] == '-' or nt_matrix[i][j] == '-':
                D[i][j] = '-'
            else:
                D[i][j] = ((alpha * gap_matrix[i][j]) + (beta * nt_matrix[i][j]))

    write_edge_file_from_matrix(D, args.outfile[0] + "_dist")

    A = D.copy()
    for i in range(len(A)):
        for j in range(len(A[i])):
            if A[i][j] != '-':
                A[i][j] = 1 - A[i][j]

    return A

def adjacency_list(A_matrix):
    '''
    Takes an adjacency matrix A and returns all of the non-null values
    in a list, for plotting in a histogram. None of the information about
    which edges correspond to which data is really retained.
    '''
    A_list = []
    for row in A_matrix:
        for item in row:
            if item != '-':
                A_list.append(item)

    return A_list

def nt_count(seq):
    '''
    Returns the number of nucleotides in a sequence.
    '''
    count = 0
    for char in seq:
        if char in ['A','C','G','T']:
            count += 1
    return count


def nt_and_gap_score(seqlist, ntscore, gapscore):
    '''
    Takes the outsequences (aligned sequences) and two scoring functions,
    performs scoring for each unique pair of sequences, calculates G and N matrices,
    then finds matrix A and returns it.
    '''

    print("Scoring all pairwise comparisons...")
    d = len(seqlist)
    N = [['-' for cols in range(d)] for rows in range(d)]
    G = [['-' for cols in range(d)] for rows in range(d)]

    global beta
    global alpha

    if alpha == 'freq':
        total_int_len = 0 # this will keep track of the total length of all intersections we have checked.
        total_nts = 0

    for i in range(len(seqlist)): # Nested loop iterates over each unique pair of sequences
        for j in range(i+1, len(seqlist)):
            # Find the intersection of the sequences
            intseq1, intseq2 = intersection(seqlist[i].sequence, seqlist[j].sequence)
        
            # NT Scoring
            ungapped1, ungapped2 = remove_gaps(intseq1, intseq2)

            if alpha == 'freq': # in this case, alpha = freq(A,C,G,T)/total length of intersections.
                total_int_len += len(intseq1) + len(intseq2)
                total_nts += nt_count(intseq1) + nt_count(intseq2)

            N[i][j] = ntscore(ungapped1,ungapped2)

            # Gap Scoring
            G[i][j] = gapscore(intseq1, intseq2)

    if alpha == 'freq':
        beta = total_nts / total_int_len
        alpha = 1 - beta
    # Normalize the matrices from 0-1
    G = normalize(G)
    N = normalize(N)
    A = normalize(adjacency_matrix(G, N))


    print("Done.")
    return A


def fifth_state_score(seqlist, score):
    '''
    Takes the outsequences (aligned sequences) and two one scoring function
    that can accept gapped sequences, performs scoring for each unique pair
    of sequences, calculates the adjacency matrix A and returns it.
    '''

    print("Scoring all pairwise comparisons...")
    D = [['-' for cols in range(len(seqlist))] for rows in range(len(seqlist))]
    # Distance matrix D

    for i in range(len(seqlist)): # Nested loop iterates over each unique pair of sequences
        for j in range(i+1, len(seqlist)):
            # Find the intersection of the sequences
            intseq1, intseq2 = intersection(seqlist[i].sequence, seqlist[j].sequence)
            
            D[i][j] = score(intseq1,intseq2)

    # Normalize the matrix from 0-1
    dist_hist = histogram(adjacency_list(D), "distances_recon", args.outfile[0])
    write_edge_file_from_matrix(D, args.outfile[0] + "_dist")
    D = normalize(D)

    # Adjacency Matrix
    A = D.copy()
    for i in range(len(A)):
        for j in range(len(A[i])):
            if A[i][j] != '-':
                A[i][j] = 1 - A[i][j]
    # Since we only have one distance matrix, A = 1 - D

    print("Done.")
    return A



def run_seq_scoring(outsequences):
    '''
    Starts the appropriate scoring algorithm based on
    command line input.
    '''
    algorithm = args.scoring_algorithm[0]
    print(algorithm)

    if algorithm.lower() in ['missingdata','md']: # Missing Data
        global alpha
        alpha = 0
        beta = 1

    elif algorithm.lower() in ['mcic']: # Will run MCIC for gaps and GTR for nts
        out = nt_and_gap_score(outsequences, GTR, MCIC_scoring)

    elif algorithm.lower() in ['fs', 'fifthstate']: # Will run GTR with gaps included
        out = fifth_state_score(outsequences, fifth_state_GTR)

    else:
        out = None
    return out


def edge_weight_stats(A):
    alist = adjacency_list(A)
    average = avg(alist)
    variance = var(alist)
    print("Average edge weight:", average)
    print("Variance:", variance)
    return average, variance


####### Write node and edge output files #########

# Want node file to look like:
# Node SeqID Length etc etc
# 1 hg38_AluJ3 32 AACCTCTAG...
# 2 hg38_AluS3X1 175 AACCTCTAG...
# ...

def write_node_file(seq_list, filename):
    '''
    Writes a file with nodes and their attributes to read into
    the network. This will be in the format shown above.
    '''
    filename = filename + "_nodes.txt"
    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        # Write the file if force is true, or if there is no file with that name.
        print("Writing node information to "+ filename + '...')
        with open(filename, 'w') as nodes:
            header = True
            # i = 0
            for i in range(len(seq_list)):
                if header == True:
                    # Write header for first line
                    header = False
                    nodes.write('Node\t')
                    nodes.write('seqID')
                    for key in seq_list[i].attributes:
                        nodes.write('\t')
                        nodes.write(key)
                    nodes.write('\n')

                # Now write the second line with the first sequence
                nodes.write(str(seq_list[i].num))
                nodes.write('\t')
                nodes.write(str(seq_list[i].ID))
                for key in seq_list[i].attributes:
                    nodes.write('\t')
                    nodes.write(seq_list[i].attributes[key])
                nodes.write('\n')
                # i += 1
        print("Done.")
        return filename



# Want edge file to look like:
# Node Node Edge_weight
# 1 2 .94
# 1 3 .76
# 1 4 .65
# 2 3 ...
# ...



def write_edge_file_from_matrix(matrix, filename):
    '''
    Creates an edge file of the format shown above, but starting
    from a matrix of finalized scores and a list of seqIDs that
    go along the axes.
    '''
    filename = filename + "_edges.txt"
    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        # Write the file if force is true, or if there is no file with that name.
        print("Writing edge information to " + filename + "...")

        with open(filename,'w') as edges:
            # First write a header
            edges.write("tail\thead\tweight\n")

            # Now write a line for each unique pair
            for i in range(len(matrix)): # Nested loop iterates over each unique pair of sequences
                for j in range(i+1, len(matrix[i])):
                    # if prune_edge(matrix[i][j]) == False:
                    if matrix[i][j] != '-':

                        edges.write(str(i))
                        edges.write('\t')
                        edges.write(str(j))
                        edges.write('\t')
                        edges.write(str(matrix[i][j]))
                        edges.write('\n')

        print("Done.")
        return filename



######## Construct Network from Node and Edge Files ########

def prune_edge(weight):
    '''
    Returns True when edge either does not exist (seqs have 
    intersection of 0) or the value is at or below the threshold,
    and False when the edge is above the threshold.
    '''

    if weight == '-':
        return True
    elif float(weight) < edge_threshold:
        return True
    else:
        return False

def parse_nodes(filename):
    '''
    Read information from nodes.txt and use this
    to build a network of unconnected nodes with
    the given attributes.
    '''
    print("Reading nodes.txt to build network...")
    with open(filename, 'r') as nodes:
        header = nodes.readline()
        header = header.strip('\n')
        header = header.split('\t')
        # print(header)
        for line in nodes:
            line = line.strip('\n')
            line = line.split('\t')
            global G
            # G.add_node(line[0])
            index = int(line[0])
            seqID = line[1]
            G.add_vertex(name=seqID, id=str(index), k=index)
            # print(G.vs[index])

            i = 2
            while i < len(line):
                G.vs[index][header[i]]=line[i]
                G.vs[index]['k'] = index

                # Add attributes. Header[i] stores name of attribute,
                # line[i] stores value for that node.
                i += 1

    print("Done.")
    print("Network contains " + str(len(G.vs())) + " nodes." )



def parse_edges(filename):
    print("Reading edges.txt to build network...")
    global G
    with open(filename, 'r') as edges:
        header = edges.readline()
        header = header.strip('\n')
        header = header.split('\t')
        # print(header)
        num = 0
        for line in edges:
            line = line.strip('\n')
            line = line.split('\t')
            tail = int(line[0])
            head = int(line[1])
            wt = line[2]

            if prune_edge(wt) == False:
                G.add_edge(tail,head, weight=float(wt))
                G.es[num]['k'] = num
                num+=1

    print("Done.")
    print("Network contains " + str(len(G.es())) + " edges.")




def parse_nodes_networkx(filename):
    '''
    Read information from nodes.txt and use this
    to build a network of unconnected nodes with
    the given attributes.
    '''
    print("Reading nodes.txt to build network...")
    global F
    F = nx.Graph()
    header = True
    with open(filename, 'r') as nodes:

        for line in nodes:
            # We ignore lines commented out with a # symbol
            if line[0] != '#':

                if header == True:
                    # First line of text is the header, containing attribute names
                    header = False
                    firstline = line.strip('\n')
                    firstline = firstline.split('\t')

                else:    
                    line = line.strip('\n')
                    line = line.split('\t')    
                    F.add_node(line[0])
                    i = 1
                    while i < len(line):
                        if firstline[i] != '':
                            F.node[line[0]][firstline[i]]=line[i]
                        # Add attributes. Header[i] stores name of attribute,
                        # line[i] stores value for that node.
                        i += 1

    print("Done.")



def parse_edges_networkx(filename):
    print("Reading edges.txt to build network...")
    global F
    with open(filename, 'r') as edges:
        header = edges.readline()
        header = header.strip('\n')
        header = header.split('\t')

        for line in edges:
            line = line.strip('\n')
            line = line.split('\t')
            tail = line[0]
            head = line[1]
            wt = line[2]

            if prune_edge(wt) == False:
                F.add_edge(tail,head, weight=float(wt))

    print("Done.")

####### Community Detection ########

def run_community_detection(network):
    '''
    Starts the appropriate community detection algorithm
    based on the command line input.
    '''

    algorithm = args.cd_algorithm[0].lower()
    print("Running community detection algorithm " + algorithm + "...")
    if algorithm in ['girvannewman','gn']: # Newman (2004)
        out = girvannewman(network, algorithm)

    elif algorithm in ['labelprop', 'lp']: # Raghavan et al (2007)
        out = labelpropagation(network, algorithm)

    elif algorithm in ['walktrap', 'wt']: # Latapy & Pons (2005)
        out = walktrap(network, algorithm)

    elif algorithm in ['fastgreedy', 'fg']: # Clauset et al (2006)
        out = fastgreedy(network, algorithm)

    elif algorithm in ['eigenvector', 'le']: # Newman 2006
        out = leadingeigenvector(network, algorithm)

    elif algorithm in ['spinglass', 'sg']: # Reichardt & Bornholdt (2006)
        out = spinglass(network, algorithm)

    else:
        out = None
        print("Algorithm not found.")
        print("Please enter a valid algorithm name.")
        sys.exit()
    print("Done.")
    return out

    # Will run the corresponding function below. Need to convert back to network
    # with nodes labeled with a community attribute.

def girvannewman(network, algorithm):
    '''
    Runs the Girvan-Newman algorithm (Newman, 2004) and converts
    the resulting dendrogram to a cluster object by cutting the dendrogram
    at the position that provides the maximum modularity.
    '''
    dendro = network.community_edge_betweenness(None, False, 'weight')

    # Returns a VertexDendrogram object, cut at the point that
    # gives maximum modularity

    clust = dendro_to_cluster(dendro)
    return clust


def fastgreedy(network, algorithm):
    '''
    Runs the FastGreedy algorithm (Clauset et al. 2006) and converts
    the resulting dendrogram to a cluster object by cutting the dendrogram
    at the position that provides the maximum modularity.
    '''
    dendro = network.community_fastgreedy('weight')
    # Returns a VertexDendrogram object

    clust = dendro_to_cluster(dendro)
    return clust

def leadingeigenvector(network, algorithm):
    '''
    Runs the Leading Eigenvector algorithm (Newman 2006) which returns 
    a cluster object.
    '''
    clust = network.community_leading_eigenvector(None, 'weight')

    return clust

def labelpropagation(network, algorithm):
    '''
    Runs the Label Propagation algorithm (Raghavan et al. 2007) which
    returns a cluster object.
    '''
    clust = network.community_label_propagation('weight')
    return clust

def spinglass(network, algorithm):
    '''
    Runs the Spinglass algorithm (Reichardt & Bornholdt 2006) which
    returns a cluster object.
    '''
    clust = network.community_spinglass('weight')
    return clust

def walktrap(network, algorithm):
    '''
    Runs the Walktrap algorithm (Pons & Latapy 2005) and converts
    the resulting dendrogram to a cluster object by cutting the dendrogram
    at the position that provides the maximum modularity.
    '''
    dendro = network.community_walktrap('weight')
    # Returns a VertexDendrogram object
    clust = dendro_to_cluster(dendro)
    return clust



def dendro_to_cluster(dendrogram):
    '''
    Accepts a dendrogram, then uses the optimal number of clusters (maximizing
    modularity) to create a VertexClustering object which can be used for modularity
    calculations and (hopefully) other things down the line.
    '''

    n = dendrogram.optimal_count # Number of clusters providing max modularity
    vc = dendrogram.as_clustering(n)

    return vc



####### Assessment ########

def modularity(vc):
    '''
    I wouldn't mind writing this myself, mostly just for kicks. Or if I can't
    figure out what I am supposed to do with the membership list.
    '''
    membership = vc.graph.vs['cluster']
    Q = vc.graph.modularity(membership,'weight')
    # Q = vc.graph.modularity(vc,'weight')
    # The membership variable can be either a VertexClustering object, which
    # is output by some community detection algorithms, or a "membership list"
    return Q

def betweenness(nodes, network):
    '''
    Returns the betweenness of a given node or nodes, which represents
    the number of shortest paths that pass through that node.
    '''
    return network.betweenness(nodes,'weight') # I think this will work.

def StrAssoc(node, vc):
    '''
    Accepts a node and a VertexClustering object, and returns the
    calculated Strength of Association of that node with its cluster,
    as defined in Greenbaum et al. (2015).
    '''

    origclust = vc.graph.vs[node]['cluster'] # The original cluster of the node under scrutiny 
    Qorig = modularity(vc)

    Qmax = -1 # The lowest possible modularity
    maxclust = origclust

    for i in range(len(vc)):
        if i != origclust: # if this isn't the original clustering assignment
            vc.graph.vs[node]['cluster'] = i # Reassign node     
            Q = modularity(vc)
            if Q > Qmax:
                Qmax = Q
                maxclust = i

    SA = Qorig - Qmax
    vc.graph.vs[node]['SA'] = SA
    vc.graph.vs[node]['cluster'] = origclust # Reset node to original assignment

    return SA

def combine_nodes(attribute_list):
    n = []
    for thing in attribute_list:
        n.append(thing)
    return n

def avg(lst):
    return sum(lst)/len(lst)

def var(lst):
    sqr_avg = avg(lst)**2
    v = 0
    for x in lst:
        v += (x-sqr_avg)
    v = v / len(lst)
    return v


def collapse_network_by_communities(clust):
    '''
    Reduces the graph to a simplified version in which each node is one 
    community, and the edges between communities are the average of the 
    individual edges between nodes in those communities.
    '''
    clust.graph.contract_vertices(clust.graph.vs['cluster'], combine_nodes)

    # We want to normalize the size of our nodes in [a,b] for visualization
    # based on the number of sequences in that cluster.

    xmin = float("inf")
    xmax = 0
    print(clust.graph.vs)
    # Find the min and max cluster sizes
    for node in clust.graph.vs:
        node['cluster'] = node['cluster'][0] 
        if len(node['name']) > xmax:
            xmax = len(node['name'])
        if len(node['name']) < xmin:
            xmin = len(node['name'])

    a = 15
    b = 50
    clust.graph.simplify(True, True, "sum")

    n = 0
    for node in clust.graph.vs:
        x = len(node['name'])
        node['numseqs'] = len(node['name'])
        if xmax == xmin:
            print("Could not size nodes by SA- all values the same.")
            newx = a
        else:
            newx = a + ((x - xmin)*(b-a))/(xmax-xmin)
        node['size'] = newx
        node['names'] = node['name']
        node['id'] = node['cluster']
        node['color'] = node['color'][0]
        node['SA_avg'] = avg(node['SA'])
        node['SA_var'] = var(node['SA'])
        n += 1
    if len(clust.graph.es()) > 0:
        add_edge_thickness_attribute(clust)
    # return collapsed

def add_color_attribute(clust):
    colors = ['#008B8B','#B22222','#228B22','#FFD700','#800080','#FF8C00','#DEB887','#4B0082','#808080','#FF00FF','#A0522D','#7FFFD4','#FF6347','#9ACD32','#DA70D6', '#FF7F50','#CD853F','#D3D3D3','#FFB6C1']
    while len(clust.graph.vs['cluster']) > len(colors):
        colors = colors + [modify_color(hexstring) for hexstring in colors]

    clust.graph.vs['color'] = [colors[cluster] for cluster in clust.graph.vs['cluster']]
    clust.graph.vs['background_color'] = clust.graph.vs['color']

def SAD(vc):
    '''
    Takes a VertexClustering object and calculates the Strength of
    Association Distribution as defined in Greenbaum et al. (2015).
    '''
    SADdict = {}

    for i in range(len(vc.graph.vs)):
        SA = StrAssoc(i, vc)
        if vc.graph.vs[i]['cluster'] in SADdict:
            SADdict[vc.graph.vs[i]['cluster']].append(SA)
        else:
            SADdict[vc.graph.vs[i]['cluster']] = [SA]

    return SADdict
    # After this I need to graph it- SA on x-axis and Frequency on y-axis

def add_size_by_SA(vc):
    '''
    Adds a size attribute, normalized from 15-50px based on the SA of each
    node.
    '''
    xmax = max(vc.graph.vs['SA'])
    xmin = min(vc.graph.vs['SA'])

    a = 15 # Min pixels
    b = 50 # Max pixels

    for node in range(len(vc.graph.vs)):
        x = vc.graph.vs['SA'][node]
        # print("SA", x)
        if xmax == xmin:
            newx = 15
            print("Could not scale node size by SA- all values identical.")
        else:
            newx = a + ( ((x - xmin)*(b-a)) / (xmax-xmin) )
        vc.graph.vs[node]['size'] = newx
        vc.graph.vs[node]['height'] = newx
        vc.graph.vs[node]['width'] = newx

def add_edge_thickness_attribute(vc):
    '''
    Adds a thickness attribute to each edge, scaled from 1-4px
    by the edge weight.
    '''
    xmax = max(vc.graph.es['weight'])
    xmin = min(vc.graph.es['weight'])
    if xmax == xmin:
        for edge in range(len(vc.graph.es)):
            x = vc.graph.es['weight'][edge]

            newx = 2
            vc.graph.es[edge]['thickness'] = newx
            vc.graph.es[edge]['width'] = newx

    else:
        a = 1 # Min pixels
        b = 5 # Max pixels

        for edge in range(len(vc.graph.es)):
            x = vc.graph.es['weight'][edge]

            newx = a + ( ((x - xmin)*(b-a)) / (xmax-xmin) )
            vc.graph.es[edge]['width'] = newx
            vc.graph.es[edge]['thickness'] = newx



def add_cluster_attribute(vc):
    '''
    Assigns the cluster assignment of the VertexClustering object
    to a 'cluster' attribute of each vertex.
    '''
    vc.graph.vs['cluster'] = vc.membership


####### Visualization #########

def histogram(datalist, label, filename):
    '''
    Plots a generic histogram of the given datalist.
    '''
    filename = filename + label + "_hist.png"

    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        print("Plotting Histogram of " + label + " to " + filename + "...")

        xmin = min(datalist)
        xmax = max(datalist)
        # Find the min and max values

        n = max(10, int(len(datalist)/15)) # number of bins
        # At the least, it will be ten, at most, average 15 data points per bin.
        binwidth = (xmax-xmin)/n
        binlist = [xmin + binwidth*i for i in range(n+1)]

        plot = plt.hist(datalist, bins=binlist, normed=False,alpha=0.5)


        plt.xlabel(label)
        plt.ylabel('Frequency')
        plt.xlim(xmin, xmax)
        plt.grid(False)

        plt.savefig(filename, format="png")
        plt.close()
        print("Done.")
        return filename

def modify_color(hexstring):
    newhex = "#"
    for char in hexstring[1:]:
        if char == "0":
            newhex += "5"
        elif char == "1":
            newhex += "6"
        elif char == "2":
            newhex += "7"
        elif char == "3":
            newhex += "8"
        elif char == "4":
            newhex += "9"
        elif char == "5":
            newhex += "a"
        elif char == "6":
            newhex += "b"
        elif char == "7":
            newhex += "c"
        elif char == "8":
            newhex += "d"
        elif char == "9":
            newhex += "e"
        elif char == "a":
            newhex += "f"
        elif char == "b":
            newhex += "0"
        elif char == "c":
            newhex += "1"
        elif char == "d":
            newhex += "2"
        elif char == "e":
            newhex += "3"
        else:
            newhex += "4"
    return newhex

def SAD_histogram(SADdict, filename, numseqs):
    '''
    Takes as input a dictionary which stores community as a key, and a 
    list of SA values of the nodes in that community as a value.
    Creates a plot of the SAD, where each community is a different color,
    and saves it in the file filename_SAD.png. These colors are (unfortunately)
    not correlated to the colors in the networks, but they should be!
    numseqs is the number of datapoints (for bin size calculation).
    '''
    filename = filename + "_SADplot.png"
    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        print("Plotting Strength of Association Distribution to " + filename + "...")

        SAmin = 1.0
        SAmax = -1.0
        # Find the min and max values
        for cluster in SADdict:
            for value in SADdict[cluster]:
                if value < SAmin:
                    SAmin = value
                if value > SAmax:
                    SAmax = value

        n = max(10, int(numseqs/15)) # number of bins
        # At the least, it will be ten, at most, average 15 data points per bin.
        binwidth = (SAmax-SAmin)/n
        binlist = [SAmin + binwidth*i for i in range(n+1)]
        colors = ['#008B8B','#B22222','#228B22','#FFD700','#800080','#FF8C00','#DEB887','#4B0082','#808080','#FF00FF','#A0522D','#7FFFD4','#FF6347','#9ACD32','#DA70D6', '#FF7F50','#CD853F','#D3D3D3','#FFB6C1']
        while len(SADdict) > len(colors):
            colors = colors + [modify_color(hexstring) for hexstring in colors]
        # same colors used for network graph
        for cluster in SADdict:
            plot = plt.hist(SADdict[cluster], facecolor = colors[cluster], bins=binlist, normed=False,alpha=0.5, label=str(cluster))


        plt.xlabel('Strength of Association')
        plt.ylabel('Frequency')
        plt.legend(loc='upper right')
        plt.xlim([SAmin, SAmax])
        # plt.title('Histogram of IQ')
        plt.grid(False)

        plt.savefig(filename, format="png")
        plt.close()
        print("Done.")
        return filename

def is_plottable(ls):
    '''
    Seaborn attempts to plot a density curve for each of the clusters.
    If a cluster has only one node in it, or if all values in the cluster have
    the same value, that is a problem for the plotting function. This function
    checks those conditions before attempting to plot.
    '''
    if len(ls) <= 1:
        return False
    else:
        i = 0
        # Loop through all pairs in the list
        # If any pair is different, then the dist. can be plotted
        while i < len(ls):
            j = i + 1
            while j < len(ls):
                if ls[i] != ls[j]:
                    # print("Returning True, different values at", i, j)
                    return True
                j += 1
            i += 1
        return False


def SAD_histogram_seaborn(SADdict, filename, numseqs):
    '''
    Takes as input a dictionary which stores community as a key, and a 
    list of SA values of the nodes in that community as a value.
    Creates a plot of the SAD, where each community is a different color,
    and saves it in the file filename_SAD.png. These colors are (unfortunately)
    not correlated to the colors in the networks, but they should be!
    numseqs is the number of datapoints (for bin size calculation).
    '''
    try:
        filename = filename + "_SADplot_seaborn.png"
        if read_existing(filename):
            print("Done.")
            return filename

        elif overwrite(filename):
            print("Plotting Strength of Association Distribution to " + filename + "...")

            SAmin = 1.0
            SAmax = -1.0
            # Find the min and max values
            for cluster in SADdict:
                for value in SADdict[cluster]:
                    if value < SAmin:
                        SAmin = value
                    if value > SAmax:
                        SAmax = value

            n = max(10, int(numseqs/15)) # number of bins
            # At the least, it will be ten, at most, average 15 data points per bin.
            binwidth = (SAmax-SAmin)/n
            binlist = [SAmin + binwidth*i for i in range(n+1)]
            # same colors used for network graph
            colors = ['#008B8B','#B22222','#228B22','#FFD700','#800080','#FF8C00','#DEB887','#4B0082','#808080','#FF00FF','#A0522D','#7FFFD4','#FF6347','#9ACD32','#DA70D6', '#FF7F50','#CD853F','#D3D3D3','#FFB6C1']
            while len(SADdict) > len(colors):
                colors = colors + [modify_color(hexstring) for hexstring in colors]
            unplotted = []
            for cluster in SADdict:
                x = SADdict[cluster]

                if is_plottable(x):
                    plot = sns.distplot(x, norm_hist=False, hist=False, color=colors[cluster], kde_kws={"shade": True})
                else:
                    unplotted.append(cluster)

            if len(unplotted) < len(SADdict):
                sns.axlabel("Strength of Association","Frequency")
                plot.set(yticks=[])
                plot.set(xlim=(SAmin-binwidth, SAmax+binwidth))

                plt.savefig(filename)
                plt.close()
                print("Done.")

            if len(unplotted) > 0:
                print("The following clusters could not be plotted because they contained only one node or all nodes had the same value:")
                for item in unplotted:
                    if unplotted.index(item) == len(unplotted) - 1:
                        print(item)
                    else:
                        print(item, end=', ')

            return filename
    except Exception as e:
        print("Could not plot the SAD histogram with Seaborn.")
        print(e)
        logging.exception('Could not plot SAD histogram with Seaborn.')
        # raise

def community_color_palette(clust):
    '''
    Returns a rainbow palette for each community to be coded
    as a different color.
    '''
    n = 0
    for comm in clust:
        n += 1

    palette = RainbowPalette(n)
    return palette
    

def plot_network(clust, filename):
    '''
    Plots clust with nodes colored by their community, and saves
    to outfilename.png.
    '''
    filename = filename + '_network.png'
    add_color_attribute(clust)

    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        print("Plotting network to " + filename + "...")
        plot(clust, filename, layout='kk', vertex_color = clust.graph.vs['color'], bbox = (1000,1000), margin = 20)

        print("Done.")

def plot_collapsed_network(clust, filename):
    '''
    Plots clust with nodes colored by their community, and saves
    to outfilename.png.
    '''

    filename = filename + '_network.png'
    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        print("Plotting network to " + filename + "...")

        plot(clust, filename, layout='kk', vertex_color = clust.graph.vs['color'], bbox = (1000,1000), margin = 50)

        print("Done.")


#################################

######### Output Files ############

def mostfrequent(string):
    '''
    Returns the most frequent character in a string.
    If there is a tie, it returns an arbitrary choice.
    '''
    chardict = {}
    for char in string:
        if char in chardict:
            chardict[char] += 1
        else:
            chardict[char] = 1

    v = list(chardict.values())
    k = list(chardict.keys())
    vote = k[v.index(max(v))]
    return vote

def consensus(seqlist):
    '''
    Takes a list of sequences in one cluster,
    and returns a consensus by voting. Since the sequences
    have already been aligned, the strings should be equal
    lengths.
    '''
    consensus = ''

    for i in range(len(seqlist[0].sequence)):
        chars = ''
        # For each character in each sequence
        for seq in seqlist:
            chars += seq.sequence[i]
        vote = mostfrequent(chars)
        consensus += vote
    return consensus

def community_consensus(seqlist, clust):
    '''
    Takes a list of all sequences in the cluster, and the clustering
    object itself. Then it calculates the consensus of each community
    and returns a list of consensus sequences, which is also assigned
    as the consensus attribute of each collapsed node.
    '''
    consensus_list = []
    for community in clust.graph.vs:
        commseqs = [seqlist[node] for node in community['k']]
        cons = consensus(commseqs)
        community['consensus'] = cons
        consensus_list.append(cons)

    return consensus_list

def output_consensus_fasta(clust, filename):
    '''
    Writes a fasta of the consensus sequences of each cluster
    to filename_cluster_consensus.fasta. No attributes in the 
    top bar right now, except the cluster number.
    '''
    filename = filename + '_cluster_consensus.fasta'
    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        with open(filename, 'w') as f:
            for node in clust.graph.vs:
                f.write('>')
                f.write('cluster_' + str(node['id']) + " ")

                f.write('\n')
                f.write(node['consensus'])
                f.write('\n')
        return filename

def write_one_line_TSV(filename, numseqs_init, numseqs_posttrim, num_comms, q,edge_avg,edge_var):
    # Writes a line like this:
    # Inputfile NumSeqs FragmentThreshold SeqsAfterTrimming ScoringMethod Alpha Beta EdgeThreshold CDAlgorithm NumComms Modularity
    filename = filename + '.tsv'
    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        print("Writing data to tsv...")
        with open(filename, 'w') as f:
            f.write(str(args.infile[0]))
            f.write('\t')
            f.write(str(numseqs_init))
            f.write('\t')
            f.write(str(args.fragment_threshold[0]))
            f.write('\t')
            f.write(str(numseqs_posttrim))
            f.write('\t')
            f.write(str(args.scoring_algorithm[0]))
            f.write('\t')
            global alpha
            global beta
            f.write(str(alpha))
            f.write('\t')
            f.write(str(beta))
            f.write('\t')
            f.write(str(edge_avg))
            f.write('\t')
            f.write(str(edge_var))
            f.write('\t')
            f.write(str(args.edge_threshold[0]))
            f.write('\t')
            f.write(str(args.cd_algorithm[0]))
            f.write('\t')
            f.write(str(num_comms))
            f.write('\t')
            f.write(str(q))
            # f.write('\n')
    return filename


def output_header(numseqs_init, numseqs_posttrim, q):
    global alpha, edge_threshold, fragment_threshold, beta

    header = '# Date: ' + str(datetime.now()) + '\n'
    header += '# Input File: ' + str(args.infile[0]) + '\n'
    header += '# Sequences: ' + str(numseqs_init) + '\n'
    header += '# Sequences after Trimming: ' + str(numseqs_posttrim) + '\n'

    header += '# Scoring Method: '
    algorithm = args.scoring_algorithm[0]
    if algorithm in ['missingdata','md']:
        header += 'GTR (Missing Data)\n'

    elif algorithm in ['mcic']:
        header += 'MCIC/GTR\n'

    elif algorithm in ['fs', 'fifthstate']:
        header += 'GTR (Fifth State)\n'

    header += '# Community Detection: '
    algorithm = args.cd_algorithm[0]
    if algorithm in ['girvannewman','gn']: # Newman (2004)
        header += 'Girvan-Newman (Newman 2004)\n'

    elif algorithm in ['labelprop', 'lp']: # Raghavan et al (2007)
        header += 'Label Propagation (Raghavan et al. 2007)\n'

    elif algorithm in ['walktrap', 'wt']: # Latapy & Pons (2005)
        header += 'Walktrap (Latapy and Pons, 2005)\n'

    elif algorithm in ['fastgreedy', 'fg']: # Clauset et al (2006)
        header += 'FastGreedy (Clauset et al. 2006)\n'

    elif algorithm in ['eigenvector', 'le']: # Newman 2006
        header += 'Leading Eigenvector (Newman 2006)\n'

    elif algorithm in ['spinglass', 'sg']: # Reichardt & Bornholdt (2006)
        header += 'Spinglass (Reichardt and Bornholdt 2006)\n'


    header += '# Modularity: ' + str(q) + '\n'
    header += '# alpha/beta: ' + str(alpha) + '/' + str(beta) + '\n'
    header += '# Edge Threshold: ' + str(edge_threshold) + '\n'
    header += '# Fragment Threshold: ' + str(fragment_threshold) + '\n'
    return header



def node_output(clust, filename, header):
    '''
    Writes an output file with each line containing a node and
    its attributes. Writes to filename_nodeoutput.txt.
    TODO: Think about appending this data (SA, Community)to the
    original node file instead.
    '''
    filename = filename + '_nodeoutput.txt'
    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        print("Writing node information to " + filename + "...")
       
        with open(filename, 'w') as nodes:
            nodes.write(header)

            header = True
            for node in range(len(clust.graph.vs)):
                if header == True:
                    # Modify header for first line
                    header = False
                    nodes.write('num\t')
                    for attr in clust.graph.vs.attributes():
                        nodes.write(attr)
                        nodes.write('\t')
                    nodes.write('\n')


                # Now write the second line with the first sequence
                nodes.write(str(node))
                nodes.write('\t')
                for attr in clust.graph.vs.attributes():
                    nodes.write(str(clust.graph.vs[attr][node]))
                    nodes.write('\t')


                nodes.write('\n')

        print("Done.")
        return filename


def cluster_output(clust, filename, header):
    '''
    Writes an output file with each line containing a community and
    its attributes. Writes to filename_clusteroutput.txt. This is run
    on the collapsed network, where each community is represented as one
    node.
    '''
    filename = filename + '_clusteroutput.txt'
    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        print("Writing cluster information to " + filename + "...")

        with open(filename, 'w') as f:
            f.write(header)

            header = True
            for cluster in range(len(clust.graph.vs)):
                attrs = ['id','numseqs', 'SA_avg','SA_var', 'names']
                if header == True:
                    # Modify header for first line
                    header = False
                    for attr in attrs:
                        f.write(attr)
                        f.write('\t')
                    f.write('\n')


                # Now write the second line with the first sequence

                for attr in attrs:
                    f.write(str(clust.graph.vs[attr][cluster]))
                    f.write('\t')

                f.write('\n')

        print("Done.")
        return filename


def GML_output(clust, filename):
    filename = filename + '.gml'
    if read_existing(filename):
        print("Done.")
        return filename

    elif overwrite(filename):
        print("Writing network information to " + filename + "...")
        clust.graph.write(filename, format='gml')
        print("Done.")

def test_GTR(seq_list, filename):
    '''
    Same as create_dist_mat and write_distmat_to_file but in 
    a simple file format, as shown above.
    '''

    files = os.listdir(".")
    if args.force == True or (filename not in files):
        # Write the file if force is true, or if there is no file with that name.
        print("Writing edge information to " + filename + "...")
        # d = len(seq_dict)
        with open(filename,'w') as edges:
            # First write a header
            edges.write("tail\thead\tGTR\tJC\n")

            # Now write a line for each unique pair
            for i in range(len(seq_list)): # Nested loop iterates over each unique pair of sequences
                for j in range(len(seq_list)):
                    if i > j:
                        # Run a particular scoring function on those sequences
                        int_seq1, int_seq2 = intersection(seq_list[i].sequence, seq_list[j].sequence)
                        int_seq1, int_seq2 = remove_gaps(int_seq1, int_seq2)
                        # print(int_seq1, int_seq2)

                        GTRscore = GTR(int_seq1,int_seq2)
                        JCscore = jukes_cantor(int_seq1, int_seq2)
                        edges.write(str(seq_list[i].ID))
                        edges.write('\t')
                        edges.write(str(seq_list[j].ID))
                        edges.write('\t')
                        edges.write(str(GTRscore))
                        edges.write('\t')
                        edges.write(str(JCscore))
                        edges.write('\n')
                    j += 1
                i += 1
        print("Done.")
    else:
        print("Refusing to overwrite existing file.")
        print("To overwrite existing files with the same name, use --force.")
        print("To read existing files with the same name, use --readexisting.")
        sys.exit()

def convert_igraph_to_networkx(nodes, edges):

    parse_nodes_networkx(nodes)
    parse_edges_networkx(edges)

def write_JSON(clust, filename):
    '''
    Writes the igraph clustering object to a JSON file,
    in the format accepted by igraph. Writes to filename.json.
    However, not currently working. GraphSpace returns an error
    int + str addition.
    '''


    outfile = filename + ".json"
    if read_existing(outfile):
        print("Done.")
        return outfile

    elif overwrite(outfile):
        print("Writing JSON file to " + outfile + "...")
        with open(outfile, "w") as f:

            f.write("{\n")

            # First we'll write all the nodes
            f.write('\t"graph":{\n')
            f.write('\t\t"nodes":[\n')
            for node in clust.graph.vs:
                # print(node.attributes())
                f.write("\t\t\t")
                f.write( json.dumps( {"data" : node.attributes() } ) )
                if node['id'] == str(len(clust.graph.vs) -1):
                    # The last node, no comma
                    f.write("\n")
                else:
                    f.write(",\n")

            f.write("\t\t],\n")

            # Then we'll write all the edges
            f.write('\t\t"edges":[\n')
            for edge in clust.graph.es:

                f.write('\t\t\t')
                attrs = { "source" : str(edge.source), "target": str(edge.target)}
                for item in edge.attributes():
                    attrs[item] = edge.attributes()[item]

                f.write(json.dumps( { "data" : attrs } ))
                if edge.index == len(clust.graph.es) -1:
                    # The last node, no comma
                    f.write("\n")
                else:
                    f.write(",\n")
            f.write("\t\t]\n")

            f.write("\t},\n")

            metadata = { "title" : filename, "description" : "", "tags" : []}

            f.write('\t"metadata":' + json.dumps(metadata) + "\n")

            '''
            }, "metadata": {"title": "Graph Name", "description": "Description of graph.. can also point to an image hosted elsewhere", "tags": ["tutorial"]}
            '''


            f.write("}")

        print("Done.")
        return outfile

def write_runtime(filename, time):
    filename = filename + "_runtime.txt"

    if overwrite(filename):
        minutes, seconds = divmod(time, 60)
        hours, minutes = divmod(minutes, 60)
        runtime = str(int(hours)) + ":" + str(int(minutes)) + ":" + str(seconds)
        print("Runtime: " + runtime)
        with open(filename, "w") as f:
            # f.write(header)
            f.write("Runtime: ")
            f.write(runtime)

def make_hists():
    #### Parsing and alignment
    insequences = parse_fasta(args.infile[0], True)
    length_list = lengths(insequences)
    length_hist = histogram(length_list, "Sequence Length", args.outfile[0])
    trimmed_insequences = remove_fragments(insequences)
    trimmed_fasta = write_fasta(trimmed_insequences, args.outfile[0] + "_trimmed")
    aligned_fasta = run_clustalomega(args.force, trimmed_fasta, args.outfile[0])
    outsequences = parse_fasta(aligned_fasta, False)

    #### Sequence scoring
    A = run_seq_scoring(outsequences)
    A_list = adjacency_list(A)
    A_hist = histogram(A_list, "Sequence Similarities", args.outfile[0])

def run():
    try:
        #### Parsing and alignment
        check_file_path(args.outfile[0])
        insequences = parse_fasta(args.infile[0], True)
        length_list = lengths(insequences)
        length_hist = histogram(length_list, "SequenceLength", args.outfile[0])

        trimmed_insequences = remove_fragments(insequences)
        trimmed_fasta = write_fasta(trimmed_insequences, args.outfile[0] + "_trimmed")
        aligned_fasta = run_clustalomega(args.force, trimmed_fasta, args.outfile[0])
        outsequences = parse_fasta(aligned_fasta, False)

        #### Sequence scoring
        A = run_seq_scoring(outsequences)
        print(A)
        A_list = adjacency_list(A)
        edge_avg = avg(A_list)
        edge_var = var(A_list)
        print("Edge average",edge_avg)
        print("Edge variance", edge_var)
        A_hist = histogram(A_list, "SequenceSimilarities", args.outfile[0])

        nodes = write_node_file(outsequences, args.outfile[0])
        edges = write_edge_file_from_matrix(A, args.outfile[0])

        #### Create the network
        parse_nodes(nodes)
        parse_edges(edges)

        #### Run community detection
        clust = run_community_detection(G)
        add_cluster_attribute(clust)

        SADlist = SAD(clust)
        SAD_histogram(SADlist, args.outfile[0], len(outsequences))
        SAD_histogram_seaborn(SADlist, args.outfile[0], len(outsequences))

        add_size_by_SA(clust)

        add_edge_thickness_attribute(clust)

        q = clust.modularity



        # Output

        global header
        header = output_header(len(insequences),len(outsequences), q)
        tsv = write_one_line_TSV(args.outfile[0], len(insequences),len(outsequences), len(clust), q, edge_avg,edge_var)
        GML_output(clust, args.outfile[0])
        outnodes = node_output(clust,args.outfile[0], header)
        plot_network(clust, args.outfile[0])
        write_JSON(clust, args.outfile[0])
        collapse_network_by_communities(clust)
        consensus(outsequences)
        community_consensus(outsequences, clust)
        output_consensus_fasta(clust, args.outfile[0])
        if len(clust) > 1:
            plot_collapsed_network(clust, args.outfile[0] + "_collapsed")
        outclusters = cluster_output(clust,args.outfile[0], header)

        print(header)
    except:
        logging.exception('Exception raised in main.')
        print("Exception raised.")
        # raise

def run_network(nodefile, edgefile):
    pass
    

def main():
    start = time.time()
    run()
    end = time.time()
    runtime = end - start

    global header
    write_runtime(args.outfile[0], runtime)

    

main()



