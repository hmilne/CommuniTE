from ete3 import Tree
from ete3 import *
import datetime
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib
import seaborn as sns
import math
from scipy.optimize import curve_fit




'''
This script:
1. Generates simulated data with INDELible V1.03 software.
2. Uses the ete3 package to convert the output Newick (PHYLIP) tree in trees.txt to a Tree object.
3. Finds all pairwise distances
4. Performs analysis of performance of the program. ** Need to talk to Anna about this phase.**
'''

def run_sim_data(outname):
    '''
    Dispatches a command to INDELible to create simulated data with the specified parameters.
    Also writes the control file (control.txt) needed to run the simulator.
    '''


    #### TYPE ####
    TYPE = "NUCLEOTIDE 1"

    TYPE_TEXT = "[TYPE] " + TYPE + "\n"

    #### MODEL ####
    MODEL = "GTR"
    MODELNAME = "mGTR"
    # Rates

    a = 0.4 # What does this value actually mean? I know it's the "rate" but...
    b = 0.4
    c = 0.4
    d = 0.4
    e = 0.4
    f = 0.4

    # Equilibrium frequencies
    pi_T = 0.25
    pi_C = 0.25
    pi_A = 0.25
    pi_G = 0.25


    # Indels
    INDELMODEL = "NB"
    q = 0.4
    r = 1 # simple: 1, complex: 10

    INSERTRATE = 0.04
    DELETERATE = 0.06

    MODEL_TEXT = "[MODEL] " + MODELNAME + "\n"
    MODEL_TEXT += "[submodel] " + MODEL + " " + str(a) + " " + str(b) + " " + str(c) + " " + str(d) + " " + str(e) + " " + str(f) + "\n"
    MODEL_TEXT += "[statefreq] " + str(pi_T) + " " + str(pi_C) + " " + str(pi_A) + " " + str(pi_G) + "\n"
    MODEL_TEXT += "[indelmodel] " + INDELMODEL + " " + str(q) + " " + str(r) + "\n"
    MODEL_TEXT += "[insertrate] " + str(INSERTRATE) + "\n"
    MODEL_TEXT += "[deleterate] " + str(DELETERATE) + "\n"

    #### TREE ####
    TREENAME = "tree1"
    ROOTED = False
    NUM_SEQS = 20
    BIRTHRATE = 2.4
    DEATHRATE = 1.1
    SAMPLE_FRAC = 0.2566
    MUTATION_RATE = 0.1 # simple: 0.1, complex: 0.13

    TREE_TEXT = "[TREE] " + TREENAME
    if ROOTED == False:
        TREE_TEXT += " [unrooted] "
    elif ROOTED == True:
        TREE_TEXT += " [rooted] "

    TREE_TEXT += str(NUM_SEQS) + " " + str(BIRTHRATE) + " " + str(DEATHRATE) + " " + str(SAMPLE_FRAC) + " " + str(MUTATION_RATE) + "\n"

    
    #### PARTITIONS ####
    # To be honest I'm not sure what a partition is in the context of the algorithm, but we only need one of them.
    PARTITION_NAME = "pGTR"
    SEQ_LENGTH = 300 # Length of the randomly generated ancestral sequence

    PARTITION_TEXT = "[PARTITIONS] " + PARTITION_NAME + ' [' + TREENAME + " " + MODELNAME + " " + str(SEQ_LENGTH) + "]" + "\n"

    #### SETTINGS ####

    # PRINTANCESTRAL = True
    OUTPUT = "FASTA"
    FASTAEXT = "fas"

    SETTINGS_TEXT = "[SETTINGS] \n"

    # SETTINGS_TEXT += "[ancestralprint] " + str(PRINTANCESTRAL) + "\n"
    SETTINGS_TEXT += "[output] " + str(OUTPUT) + "\n"
    SETTINGS_TEXT += "[fastaextension] " + str(FASTAEXT) + "\n"


    #### EVOLVE ####
    REPETITIONS = 1 #Number of different trees to create
    OUTPUTNAME = outname

    EVOLVE_TEXT = "[EVOLVE] " + PARTITION_NAME + " " + str(REPETITIONS) + " " + OUTPUTNAME + "\n" 



    with open("control.txt", "w") as f:
        currenttime = datetime.datetime.now()
        currenttime = currenttime.strftime("%Y-%m-%d-%H-%M-%S")
        f.write("// " + str(currenttime) + "\n")
        f.write(TYPE_TEXT)
        f.write(MODEL_TEXT)
        f.write(TREE_TEXT)
        f.write(PARTITION_TEXT)
        f.write(SETTINGS_TEXT)
        f.write(EVOLVE_TEXT)


    INDELIBLE_PATH = "INDELible_1.03_Windows.exe"
    # The actual path is in my PATH variable, so I don't have to write out the whole thing in the command.

    os.system(INDELIBLE_PATH)

    """
    Example control file:

    [TYPE] NUCLEOTIDE 1 

    [MODEL]    mGTR
        [submodel]  GTR 0.4 0.4 0.4 0.4 0.4 0.4   //  GTR: alpha beta gamma delta ...
        [statefreq] 0.25 0.25 0.25 0.25           //  pi_T pi_C pi_A pi_G   
        [indelmodel]   NB  0.4 1  //  Geometric indel length distribution (q=0.4, r=1)
        [insertrate]   0.08       //  insertion rate = 0.08 relative to substitution rate of 1
        [deleterate]   0.12       //  deletion rate = 0.12 relative to substitution rate of 1     
      
    [TREE] tree1 [unrooted] 200 2.4 1.1 0.2566 0.34  // nseqs birth-rate, death-rate, sampling fraction, mutation rate

    [PARTITIONS] pGTR [tree1 mGTR 300] // partitionname [ treename modelname seqlen ]  

    [EVOLVE] pGTR 1 simdata1 // partitionname numreps outputname
                                           
    """



def parse_newick_tree(filename):
    '''
    Parses a trees.txt file from INDELible and converts it to an ete3 Tree datatype for further analysis.
    '''
    print("Parsing " + filename + " trees...")
    with open(filename, "r") as f:
        i = 0
        for line in f:
            # print(i)
            # print(line)

            if i == 4:
                # print(line)
                line = line.strip('\n')
                header = line.split('\t')
                # print(header)

            if i == 6:
                # print(line)
                line = line.strip('\n')
                data = line.split('\t')
                # print(data)
            i+=1

    tree_attrs = {}
    for i in range(len(header)):
        tree_attrs[header[i]] = data[i]
    # print(tree_attrs)
    t = Tree(tree_attrs['TREE STRING'], format=5)
    print(t)
    return t


def printTable(table):
    '''
    Just prints a table for easier reading
    '''
    for i in range(len(table)):
        print(table[i])

def tree_dist_mat(tree):
    '''
    Takes the generated tree from the simulated data and calculates distances between individual nodes,
    then fills a distance matrix.
    At this time I'm not sure how to insure that the values are in the same positions in this table
    and the other one (from the pipeline).
    At this time I am also not sure how to make it so it only runs half the matrix (since it's symmetrical anyway).
    '''
    print("Calculating distance matrix from tree...")
    distmat = [["-" for node in tree.iter_leaves()] for node in tree.iter_leaves()] # Init matrix
    i = 0 # Index for node 1
    nodes = []
    for node1 in tree.iter_leaves():
        nodes.append(node1.name)
        j = 0 # Index for node 2

        for node2 in tree.iter_leaves():
            # print(node1)
            # print(node2)
            # print("Distance:", tree.get_distance(node1, node2))
            if i < j:
                distmat[i][j] = tree.get_distance(node1, node2)
            j += 1
        i+=1 

    print("Done.")
    return distmat, nodes

def normalize(matrix):
    xmax = 0
    xmin = 10000
    # Find the min and max
    # print(matrix)
    for row in matrix:
        # print(row)
        for item in row:
            # print(item)
            if item != '-':
                if item > xmax:
                    xmax = item
                if item < xmin:
                    xmin = item
    # Now perform normalization
    # print("MAX:", xmax)
    # print("MIN:", xmin)

    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] != '-':
                if matrix[i][j] == 0:
                    # print(matrix[i][j])
                    matrix[i][j] = (matrix[i][j] - xmin)/(xmax - xmin)
                    # print(matrix[i][j])
                else:
                    matrix[i][j] = (matrix[i][j] - xmin)/(xmax - xmin)                    
                # print()

    return matrix

def adjacency_mat(D):
    '''
    Takes a raw (non-normalized) distance matrix and converts it to a
    normalized adjacency matrix suitable for edge weights.
    '''
    print("Calculating adjacency matrix...")

    D = normalize(D) # Gets all values in the range (0-1)
    # print("Normalized")
    # printTable(D)
    # print()
    A = D.copy() 
    # Now get adjacency rather than distance.
    for i in range(len(A)):
        for j in range(len(A[i])):
            if A[i][j] != '-':
                A[i][j] = 1 - A[i][j]
    print("Done.")
    return A


def write_edge_file(nodelist, A, filename):
    '''
    Same as create_dist_mat and write_distmat_to_file but in 
    a simple file format, as shown above.
    '''
    filename = filename + "_edges.txt"

    # Write the file if force is true, or if there is no file with that name.
    print("Writing edge information to " + filename + "...")
    # d = len(seq_dict)
    with open(filename,'w') as edges:
        # First write a header
        edges.write("tail\thead\tweight\n")

        # Now write a line for each unique pair
        for i in range(len(nodelist)): # Nested loop iterates over each unique pair of sequences
            for j in range(i+1, len(nodelist)):
                # Run a particular scoring function on those sequences
                score = A[i][j]

                edges.write(str(i))
                edges.write('\t')
                edges.write(str(j))
                edges.write('\t')
                edges.write(str(score))
                edges.write('\n')
        print("Done.")
        return filename

def matrix_to_list(matrix):
    '''
    Takes an adjacency matrix A and returns all of the non-null values
    in a list, for plotting in a histogram. None of the information about
    which edges correspond to which data is really retained.
    '''
    lst = []
    for row in matrix:
        for item in row:
            if item != '-':
                lst.append(item)

    # print(lst)
    return lst

def histogram(datalist, label, filename):
    '''
    Plots a generic histogram of the given datalist.
    '''
    filename = filename + label + "_hist.png"

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
    # plt.legend(loc='upper right')
    plt.xlim(xmin, xmax)
    # plt.xticks(np.arange(binlist[0], binlist[-1], binwidth))
    # plt.title('Histogram of IQ')
    plt.grid(False)

    plt.savefig(filename, format="png")
    plt.close()
    # plt.show()
    print("Done.")
    return filename


def check_file_path(filepath):
    '''
    Checks the existence of the directories in the outfile path.
    If the directories do not exist, they will be created.
    '''
    # filepath = os.path.dirname(filename)
    if filepath == '':
        filepath = './'

    if os.path.exists(filepath):
        # print("Directory " + filepath + " exists!")
        pass
    else:
        # print("Creating directory " + filepath + "...")
        os.makedirs(filepath)
        print("Done.")

def change_dir(path):
    '''
    If the path exists, it will navigate to it. If the path
    doesn't exist, then it creates it then moves to that directory.
    '''
    check_file_path(path)
    os.chdir(path)



def run(OUTPUT_PATH, OUTPUTNAME):
    # OUTPUT_PATH = "../Data_outputs/simdata1/simdata1true"
    # Writes to directory simdata1/simdata1true
    # Then the reconstructed ones will write to simdata1/simdata1recon/other stuff
    # OUTPUTNAME = "simdata1true"

    cwd = os.getcwd()

    change_dir(OUTPUT_PATH)
    run_sim_data(OUTPUTNAME)

    t = parse_newick_tree("trees.txt")
    D, nodelist = tree_dist_mat(t)
    # print(nodelist)
    # print(D)

    distedges = write_edge_file(nodelist, D, OUTPUTNAME + "_dist")
    dist_hist = histogram(matrix_to_list(D), "true_distances", OUTPUTNAME)
    # print()
    # print("Distance matrix")
    # printTable(D)
    # print(nodelist)
    # nodelist is just so we know the order of things in the distance matrix.
    # It's the same order as in the fasta.
    # The index of each node in the list is equivalent to its index in the node output file
    # from the main pipeline.
    A = adjacency_mat(D)
    # print()
    # print("Adjacency matrix")
    # printTable(A)
    edges = write_edge_file(nodelist, A, OUTPUTNAME)
    clusters = k_clusters(t, OUTPUTNAME)
    clusterfile = cluster_output(clusters, OUTPUTNAME)

    # dist_scatter(matrix_to_list(A), "simdata1test_edges.txt", "distscatter.png")

    change_dir(cwd)
    return t


##################################
###### Cutting the Tree ##########
##################################



def find_next_merge(tree, clusters):
    '''
    Takes a list of clusters and a tree and returns the indices of
    the next two clusters that should be merged.
    '''
    min_distance=float('inf')
    min_clusters = ()
    i = 0
    while i < len(clusters):
        j = i + 1
        while j < len(clusters):
            if len(clusters[i]) == 1:
                node1 = clusters[i][0]
            else:
                node1 = tree.get_common_ancestor(clusters[i])

            if len(clusters[j]) == 1:
                node2 = clusters[j][0]
            else:
                node2 = tree.get_common_ancestor(clusters[j])

            distance = tree.get_distance(node1, node2)
            if distance < min_distance:
                min_distance = distance
                min_clusters = (i,j)
            j+=1
        i+=1
    
    # print(min_distance)
    # print(min_clusters)
    return min_clusters

def k_clusters(tree, outfile):
    '''
    Returns the tree divided optimally into k clusters. This is 
    performed by hierarchical clustering using distances between nodes.
    '''
    clusters = [] # Will store a list of lists, where each list is a cluster
    # use t.get_distance(node1, node2) to find the distance between two nodes
    # use t.common_ancestor(node1, node2) to get the least common ancestor of two nodes
    # So to get k clusters, the algorithm performs k-1 merges, where at each step it
    # merges the things that are the closest together. When these are just nodes,
    # t.get_distance will find the distance between them. When clusters have more than one
    # node in them, use get_common_ancestor(list of nodes) to get the common ancestor
    for node in tree:
        clusters.append([node.name])
    outfile = outfile + "_k.tsv"

    print("Opened file " + outfile)

    k = len(tree)
    merge = 0
    SDRscores = []
    AICs = []
    AICcs = []
    ks = []
    merges = []
    clusterings = []
    merge += 1
    k -= 1
    # TODO: Choose the smallest k with an SDR below 0.5?
    # TODO: Choose the value where that minimizes the distance to the origin.
    # So it minimizes sqrt(k^2 + score^2)
    mergestart = time.time()
    while merge < len(tree)-1:
        i,j = find_next_merge(tree,clusters)
        clusters[i] = clusters[i] + clusters[j]
        del clusters[j]

        if True:
            sdrstart = time.time()
            score = SDR(tree, clusters)
            sdrend = time.time()
            sdrtime = sdrend - sdrstart
            AIC_score = AIC(k,score)
            if k < len(tree) - 1:
                AICc_score = AICc(k, len(tree),score)
            else:
                AICc_score = "-"
        # else:
        #     score = "-"
        #     AIC_score = "-"
        #     AICc_score = "-"
        #     sdrtime = "-"

        AICs.append(AIC_score)
        SDRscores.append(score)
        ks.append(k)
        merges.append(merge)
        AICcs.append(AICc_score)
        clusterings.append(clusters.copy())

        print("k",k)
        print("merge",merge)
        print("score",score)
        print("AIC",AIC_score)
        print("AICc",AICc_score)
        merge += 1
        k -= 1
        mergeend = time.time()
        # print("Merge time:", mergeend - mergestart - sdrtime)
        # print("SDR time:", sdrtime)

    maxs = max(SDRscores)
    mins = min(SDRscores)

    maxk = max(ks)
    mink = min(ks)

    dists = [ (((ks[i]-mink)/(maxk-mink))**2 + ((SDRscores[i]-mins)/(maxs-mins))**2)**(1/2) for i in range(len(ks))  ]
    print("scores",SDRscores)
    print("ks",ks)
    print("dists",dists)


    with open(outfile, "w") as f:
        f.write("Merges\tk\tSDR\tAIC\tAICc\tDistance\n")
        best = float('inf')
        bestk = 0
        bestcluster = None
        for i in range(len(ks)):
            line = str(merges[i]) + "\t" + str(ks[i]) + "\t" + str(SDRscores[i]) + "\t" + str(AICs[i]) + "\t" + str(AICcs[i]) + "\t" + str(dists[i]) + "\n"
            f.write(line)
            if dists[i] < best:
                best = dists[i]
                bestk = ks[i]
                bestcluster = clusterings[i]


    print("Best", best)
    print("Best k", bestk)
    print("Best clustering", bestcluster)

    return bestcluster

def neighbors(node1, node2, clusters):
    '''
    Returns true if the nodes are in the same cluster
    and False otherwise.
    '''
    for cluster in clusters:
        if node1 in cluster:
            if node2 in cluster:
                return True
                # Found both, must be true
            else:
                return False
                # Found node1 but not node2, must be false
        elif node2 in cluster:
            # already checked for node1, so if we find node2 now,
            # we know we don't have node1 so must be false
            return False
    return False

def SDR(tree, clusters):
    '''
    Returns something equivalent to the subtype diversity ratio
    used in HIV subtype analysis. The mean within-cluster pairwise
    distance over the mean between-cluster pairwise distance. We 
    want to minimize this value.
    '''
    wi = [0,0] # Stores sum of distances in wi[0] and number of pairs examined in wi[1]
    bw = [0,0] # Same here
    nodes = tree.get_leaves()
    i = 0
    while i < len(nodes):
        j = i + 1
        while j < len(nodes):
            # print(wi)
            # print(bw)
            if neighbors(nodes[i].name,nodes[j].name,clusters):
                # print(nodes[i].name,nodes[j].name)
                wi[0] += tree.get_distance(nodes[i].name,nodes[j].name)
                wi[1] += 1
            else:
                bw[0] += tree.get_distance(nodes[i].name,nodes[j].name)
                bw[1] += 1
            j += 1

        i += 1
    # print(clusters)
    # print(wi,bw)
    if wi[1] != 0:
        wi = wi[0]/wi[1]
    else:
        wi = 0
    bw = bw[0]/bw[1]
    score = wi/bw
    return score

def AIC(k, L):
    return 2*k - 2*math.log1p(L)

def AICc(k,n,L):
    return AIC(k,L) + (2*k*(k+1))/(n-k-1)

def plot_curve(x,y,xlab,ylab,label,outfile):
    '''
    Takes an x and y list, plots a scatter plot with a curve, and
    outputs to outfile.
    '''

    matplotlib.pyplot.scatter(x,y)
    matplotlib.pyplot.show()
    filename = filename + label + "_scatter.png"
    print("Plotting Scatter of " + label + " to " + filename + "...")
    plot = plt.hist(datalist, bins=binlist, normed=False,alpha=0.5)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    # plt.xticks(np.arange(binlist[0], binlist[-1], binwidth))
    # plt.title('Histogram of IQ')
    plt.grid(False)

    plt.savefig(filename, format="png")
    plt.close()
    # plt.show()
    print("Done.")
    return filename

def plot_double_curve(x,y1,y2,xlab,y1lab,y2lab,label, outfile):
    y1plot = matplotlib.pyplot.scatter(x,y1,color="blue",label=y1lab)
    y2plot = matplotlib.pyplot.scatter(x,y2,color='red',label=y2lab)
    plt.grid(False)
    plt.xlabel(xlab)
    # plt.ylabel(y1lab)
    plt.legend()
    matplotlib.pyplot.show()
    filename = filename + label + "_scatter.png"
    print("Plotting Scatter of " + label + " to " + filename + "...")

    # plt.xticks(np.arange(binlist[0], binlist[-1], binwidth))
    # plt.title('Histogram of IQ')


    plt.savefig(filename, format="png")
    plt.close()
    # plt.show()
    print("Done.")
    return filename



# edges = [0.2,0.4,0.6,0.8]
# mods = [0.1,0.3,0.6,0.8]
# fowlkes = [1,0.6,.4,0.2]

# plot_double_curve(edges,mods,fowlkes,"Edge Threshold","Modularity","Fowlkes","test","testscatter")

def cluster_output(clusters, filename):
    '''
    Writes an output file with each line containing a community and
    its attributes. Writes to filename_clusteroutput.txt. This is run
    on the collapsed network, where each community is represented as one
    node.
    '''
    filename = filename + '_clusteroutput.txt'

    print("Writing cluster information to " + filename + "...")
   
    #### TODO
    # node = 0
    with open(filename, 'w') as f:
        header = True
        # i = 0
        for i in range(len(clusters)):
            attrs = ['id','numseqs', 'SA_avg','SA_var', 'names']
            if header == True:
                # Modify header for first line
                header = False
                for attr in attrs:
                    f.write(attr)
                    f.write('\t')
                f.write('\n')


            # Now write the second line with the first sequence
            f.write(str(i))
            f.write("\t")
            f.write(str(len(clusters[i])))
            f.write("\t")
            f.write("-")
            f.write("\t")
            f.write("-")
            f.write("\t")
            f.write(str(clusters[i]))
            f.write('\n')

    print("Done.")
    return filename


def main():
    t = run("../Data_outputs/simdatatest", "simdatatest6")
    # k_clusters(t,5)
    return t

# t = main()

# dist_scatter("simdata1_edges.txt", "simdata1test_edges.txt", "distscatter.png")



'''
Two ways to get the precision/recall:
1. Compare distances in the true tree (normalized) to distances found by the pipeline.
2. Run pipeline all the way to community detection, compare how many nodes were placed in the right community.
'''

