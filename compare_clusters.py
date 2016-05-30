import sys

'''
The code in this file takes two clustering output files from CommuniTE (the program I built for my undergraduate thesis)
and compares them to obtain a similarity score called the Fowlkes index.
'''

#### Calculating the Fowlkes index ####

def sqrt(n):
    '''
    Returns the square root of n
    '''
    return n**(1/2)

def choosetwo(n):
    '''
    Returns the value of n choose 2
    '''
    return ((n*(n-1))/2)

def sumchoosetwo(ls):
    '''
    Returns the sum of n choose 2 for all n in ls
    '''
    total = 0
    for item in ls:
        # print(choosetwo(item))
        total += choosetwo(item)
    return total


def wallace(clustering1, clustering2):
    '''
    A clustering is represented as a list of lists where each list represents
    one cluster in a graph.

    This function performs asymmetrical scoring of two clusters: the number of 
    point pairs in the same cluster in both clusterings (npos) divided by
    the sum of the cluster sizes choose two across clustering1.
    '''
    npos = 0
    for cluster1 in clustering1:
        for node1 in cluster1:
            for cluster2 in clustering2:
                for node2 in cluster2:
                    npos += (node1 < node2 and node2 in cluster1 and node1 in cluster2)
                    # if node1 < node2:
                        # print(node1,node2)
                        # if node2 in cluster1 and node1 in cluster2:
                            # print(cluster1,cluster2)
                            # print(node1,node2)
                            # npos += 1
    clust_lengths = [len(cluster) for cluster in clustering1]
    W = npos / sumchoosetwo(clust_lengths)
    # print(W)
    return W



def fowlkes(clustering1, clustering2):
    '''
    Takes the geometric mean of the asymmetrical wallace index with clusterings swapped
    for a symmetrical score of clustering similarity.
    '''
    return sqrt(wallace(clustering1, clustering2)*wallace(clustering2, clustering1))


#### Parsing the clustering files ####


def truncate_at(string, symbol):
    '''
    Truncates a string at (not including) a given symbol.
    '''

    newstring = ''
    found = False
    i = 0
    while (found==False) and i < len(string):
        if string[i] == symbol:
            found = True
        else:
            newstring += string[i]
            i += 1
    return newstring

def get_clusters(clusteroutput):
    '''
    Takes a cluster output file from main.py and returns a list of lists where each list is
    a cluster containing values representing the ID numbers of sequences in that cluster.
    '''
    clusters = []
    with open(clusteroutput, "r") as f:
        header = True
        for line in f:
            if line[0] != "#":
                line = line.strip("\n")
                line = line.split("\t")
                if header == True:
                    header = False
                else:
                    nodes = line[4].lstrip("[")
                    nodes = nodes.rstrip("]")
                    nodes = nodes.replace("'","")
                    nodes = nodes.split(", ")
                    truncated_nodes = []
                    for node in nodes:
                        truncated_nodes.append(node.split("_")[-1])
                    clusters.append(truncated_nodes)
    # print(clusters)
    return clusters



def compare_clusters(clusterfile1, clusterfile2):
    '''
    Takes two clustering output files and performs a pairwise comparison to get a final score of
    clustering similarity.
    '''
    print("COMPARING CLUSTERS")
    clustering1 = get_clusters(clusterfile1)
    clustering2 = get_clusters(clusterfile2)

    return fowlkes(clustering1, clustering2)



##### Stuff for creating true Alu clusters ##### 


def parse_fasta(filename):
    '''
    This function parses the input fasta file and returns a list of Seq datatypes (see Seq.py).
    The add_name_ext input is a boolean that determines whether a unique identifier is added to
    the name of each sequence.
    '''
    print("Reading sequences from", filename)
    sequences = {} # store number: id combos
    count = 0
    with open(filename, 'r') as infile:
        for line in infile:
            if line[0] == '>': # the first character of the line
                words_in_line = line.split() #split the line into words
                seqID = words_in_line[0][1:]
                sequences[count] = seqID
                count += 1

    print("Done.")
    print("Read", len(sequences), "sequences.")

    return sequences


def read_clusters_from_fasta(seqfile):
    '''
    Puts sequences into clusters based on their seq IDs.
    '''
    seqs = parse_fasta(seqfile)
    clusters = {}
    for seq in seqs:
        if seqs[seq] in clusters:
            clusters[seqs[seq]].append(seq)
        else:
            clusters[seqs[seq]] = [seq]
    print(clusters)
    return clusters

def write_clusters_to_file(dct, filename):
    print("Writing clusters to " + filename + "...")
    with open(filename, "w") as f:
        f.write("id\tnumseqs\tSA_avg\tSA_var\tnames\n")
        for cluster in dct:
            f.write(str(cluster))
            f.write("\t")
            f.write(str(len(dct[cluster])))
            f.write("\t")
            f.write("-")
            f.write("\t")
            f.write("-")
            f.write("\t")
            f.write(str(dct[cluster]))
            f.write("\n")
    print("Done.")
    return filename


def cluster_true(seqfile,outfile):
    dct = read_clusters_from_fasta(seqfile)
    write_clusters_to_file(dct, outfile)
    return outfile

# cluster_true("../Data/Alus_500_chr9.fasta", "./test_Alu500cluster.txt")

def main():
    print(sys.argv)
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    f=compare_clusters(file1, file2)
    print(f)

main()