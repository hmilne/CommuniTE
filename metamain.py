import os
import argparse
import simdata
import matplotlib.pyplot as plt 
import numpy as np
import math
import matplotlib
from scipy import stats
import seaborn as sns
from Seq import *
import compare_clusters

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', dest='outfile', nargs=1, help='input file in FASTA format') # stores one argument in outfile
    parser.add_argument('-i', '--infile', dest='infile', nargs=1, help='output file tag') # stores one argument in infile
    # parser.add_argument('-v', dest='verbose', action='store_true', help='show verbose output')

    parser.add_argument('--partone', dest='partone', action='store_true', help='to run part one trials, and part two defaults', default=False)
    parser.add_argument('--parttwo', dest='parttwo', action='store_true', help='to run part one defaults, and part two trials', default=False)
    parser.add_argument("-sim", dest='sim', action='store', nargs=1, help="to run n simulations, then part one for each dataset.", default=[0])
    parser.add_argument("--truealignment", dest='truealignment', action='store_true', help='if true, will run the pipeline on the true alignment', default=False)

    parser.add_argument('--force', dest='force', action='store_true', help='force overwrite of existing files', default=False)
    args = parser.parse_args()

fragment_thresholds = [0,0.2, 0.5, 0.8]
scoring_algorithms = ['MCIC']
alphas = [0, 1, 'freq']

edge_thresholds = [0.5,0.7,0.9] # If edge.weight < edge_threshold, don't include
cd_algorithms = ['fg']

def check_file_path(filename):
    '''
    Checks the existence of the directories in the outfile path.
    If the directories do not exist, they will be created.
    '''
    print("checking file path",filename)
    filepath = os.path.dirname(filename)
    if filepath == '':
        filepath = './'

    if os.path.exists(filepath):
        # print("Directory " + filepath + " exists!")
        pass
    else:
        # print("Creating directory " + filepath + "...")
        os.makedirs(filepath)
        # print("Done.")

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
                # print(seqID)

                for j in range(1,len(words_in_line)):
                    components = words_in_line[j].split('=')
                    # sequences[seqID][components[0]] = components[1]
                    sequences[count].setAttribute(components[0], components[1])
                    # adds each trait on the line as a key:value pair to the inner dictionary
                sequences[count].sequence = '' # initialize seq key with empty string value
            else:
                #if the line doesn't begin with > then append the contents to the previous seq's sequence
                # print(count)
                sequences[count].sequence += line[:-1]

            # sequences[count].setSequence(sequence)
    print("Done.")
    print("Read", len(sequences), "sequences.")

    return sequences

def write_trimmed_edge_file(fasta, edges, outfile, fragment_threshold):
    '''
    Makes edge files with only the sequences above the fragment threshold.
    It does this by reading through the fasta file, identifying sequences that are shorter than the fragment threshold, and
    excluding those sequences when writing the new edge file.
    '''

    seqs = parse_fasta(fasta, False) # Returns a list of seq datatypes
    lengths = [s.length() for s in seqs]
    print(max(lengths))
    threshold = int(max(lengths)*fragment_threshold)
    print(threshold)
    print()
    fragments = []
    for i in range(len(lengths)):
        if lengths[i] < threshold:
            fragments.append(i) #seq.num is the ID number with which the sequence was read in from the fasta


    with open(edges, "r") as original, open(outfile, "w") as new:
        header = original.readline()
        new.write(header)
        # print(header)
        num = 0
        for line in original:
            line = line.strip('\n')
            info = line.split('\t')
            tail = info[0]
            head = info[1]
            wt = info[2]
            if (int(tail) in fragments) or (int(head) in fragments):
                pass
            else:
                new.write(line)
                new.write('\n')


    return outfile


def run_part_one(infile, outfileprefix, simpath=None, truealignment=None):
    '''
    Runs part one of the pipeline with a range of values defined in the global frame.
    The variables include fragment threshold, scoring algorithm, and alpha.
    simpath is the path to the simulated data, if that is what is being used here.
    '''
    global fragment_thresholds
    global scoring_algorithms
    global alpha



    trial = 1 # add this to the end of the outfile prefix
    for f in fragment_thresholds:
        new_f = True # This helps keep track of when to perform a new alignment

        for s in scoring_algorithms:
            if s in ['fifthstate', 'fs']:
                alphas_modified = [None]
                # We don't need alpha for fifth state, so we can
                # save some unnecessary trials by doing this.
            else:
                alphas_modified = alphas
            for a in alphas_modified:
                print()
                print("prefix",outfileprefix)
                filepath = outfileprefix + "\\frag_" + str(f) + "\\" + str(s) + "\\alpha_" + str(a) + "\\"
                outfile = filepath + "\\" + str(extract_filename_from_path(outfileprefix)) + "_" + str(trial)
                print("filepath",filepath)
                print("filename", extract_filename_from_path(outfileprefix))

                print("outfile", outfile)

                check_file_path(outfile)
                if new_f == True:
                    alignment = outfile + ".afasta" # Assign alignment to output file of this trial
                    if simpath != None:
                        sim_adjacency_file = write_trimmed_edge_file(simpath + ".fas", simpath + "_edges.txt", outfile + "true_edges_trimmed.txt", f) #second input True if similarity, false if distance
                        sim_dist_file = write_trimmed_edge_file(simpath + ".fas", simpath + "_dist_edges.txt", outfile + "true_dist_edges_trimmed.txt",f)

                logfile = outfile + ".log"
                mainpipeline = "c:/Users/Heather/Documents/Fall_2015/Thesis/Code/main.py"
                command = "python " + mainpipeline + " -i " + str(infile) + " -o " + str(outfile)
                command += " -fragment " + str(f)
                command += " -s " + str(s)
                command += " -alpha " + str(a)
                if args.truealignment == True and truealignment != None:
                    command += " --readalignment " + str(os.path.abspath(truealignment))
                if args.force == True:
                    command += " --force "

                if new_f == False:
                    command += " --readalignment " + alignment 
                    # Read existing output file from a previous run with this f value
                command += " > "
                command += logfile
                print()
                print()
                print()
                # defaults are used for the edge threshold (0.85) and community detection (fastgreedy)
                print("Running main.py with the following parameters:")
                # print("Input: " + infile)
                # print("Output: " + outfile)
                print("Fragment Threshold: " + str(f))
                print("Scoring Algorithm: " + str(s))
                print("Alpha: " + str(a))
                print("Command:")
                # print("Redirecting stdout to " + logfile + ".")

                print(command)
                os.system(command)
                new_f = False

                if simpath != None:

                    recon_adjacency_file = outfile + "_edges.txt"
                    recon_dist_file = outfile + "_dist_edges.txt"

                    adj_scatter_file = outfile + "_adj"
                    dist_scatter_file = outfile + "_distfs"
                    r_value_adj = scatter(sim_adjacency_file, recon_adjacency_file, "True Normalized Similarity", "Reconstructed Similarity", adj_scatter_file)
                    r_value_dist = scatter(sim_dist_file, recon_dist_file, "True Distance", "Reconstructed Distance", dist_scatter_file)


                    # os.getcwd()
                    with open(outfile + ".tsv", "a") as out:
                        out.write("\t")
                        out.write(str(r_value_dist))
                        out.write("\t")
                        out.write(str(r_value_adj))
                        out.write("\n")
                        out.close()

                tsvs.append(outfile + ".tsv")
                trial += 1

                print("Done.")

    if simpath != None:
        tsv = write_tsv(tsvs, outfileprefix + "/" + str(extract_filename_from_path(outfileprefix)), True)
    elif truealignment != None:
        tsv = write_tsv(tsvs, outfileprefix + "/" + str(extract_filename_from_path(outfileprefix)), "truealign")
    return tsv
    # Finish this!

def run_one_two(infile, outfileprefix, simpath=None, truealignment=None):
    tsvs = []
    global fragment_thresholds
    global scoring_algorithms
    global alpha
    global edge_thresholds
    global cd_algorithms

    trial = 1 # add this to the end of the outfile prefix
    for f in fragment_thresholds:
        new_f = True # This helps keep track of when to perform a new alignment

        for s in scoring_algorithms:
            if s in ['fifthstate', 'fs']:
                alphas_modified = [None]
                # We don't need alpha for fifth state, so we can
                # save some unnecessary trials by doing this.
            else:
                alphas_modified = alphas
            for a in alphas_modified:
                for e in edge_thresholds:
                    for c in cd_algorithms:

                        print()
                        print("prefix",outfileprefix)
                        filepath = outfileprefix + "/frag_" + str(f) + "/" + str(s) + "/alpha_" + str(a) + "/edge_" + str(e) + "/" + str(c)
                        outfile = filepath + "/" + str(extract_filename_from_path(outfileprefix)) + "_" + str(trial)
                        print("filepath",filepath)
                        print("filename", extract_filename_from_path(outfileprefix))

                        print("outfile", outfile)

                        check_file_path(outfile)
                        if new_f == True:
                            alignment = outfile + ".afasta" # Assign alignment to output file of this trial
                            if simpath != None:
                                sim_adjacency_file = write_trimmed_edge_file(simpath + ".fas", simpath + "_edges.txt", outfile + "true_edges_trimmed.txt", f) #second input True if similarity, false if distance
                                sim_dist_file = write_trimmed_edge_file(simpath + ".fas", simpath + "_dist_edges.txt", outfile + "true_dist_edges_trimmed.txt",f)

                        logfile = outfile + ".log"
                        mainpipeline = "c:/Users/Heather/Documents/Fall_2015/Thesis/Code/main.py"
                        command = "python " + mainpipeline + " -i " + str(infile) + " -o " + str(outfile)
                        command += " -fragment " + str(f)
                        command += " -s " + str(s)
                        command += " -alpha " + str(a)
                        command += " -edge " + str(e)
                        command += " -c " + str(c)

                        if args.truealignment == True and truealignment != None:
                            command += " --readalignment " + str(os.path.abspath(truealignment))
                        if args.force == True:
                            command += " --force "

                        if new_f == False:
                            command += " --readalignment " + alignment 
                            # Read existing output file from a previous run with this f value
                        command += " > "
                        command += logfile
                        print()
                        print()
                        print()
                        # defaults are used for the edge threshold (0.85) and community detection (fastgreedy)
                        print("Running main.py with the following parameters:")
                        # print("Input: " + infile)
                        # print("Output: " + outfile)
                        print("Fragment Threshold: " + str(f))
                        print("Scoring Algorithm: " + str(s))
                        print("Alpha: " + str(a))
                        print("Command:")
                        # print("Redirecting stdout to " + logfile + ".")

                        print(command)
                        os.system(command)
                        new_f = False

                        if simpath != None:
                            recon_clust = outfile + "_clusteroutput.txt"
                            true_clust = simpath + "_clusteroutput.txt"
                        else:
                            print("infile",infile)
                            print("outfile",outfile + "_clusteroutput.txt")
                            true_clust = compare_clusters.cluster_true(infile, outfileprefix + "_trueclust.txt")
                            recon_clust = outfile + "_clusteroutput.txt"
                        fowlkes = compare_clusters.compare_clusters(recon_clust, true_clust)
                        print("Fowlkes",fowlkes)
                        print("Writing information to " + outfile + "_fowlkes.txt...")
                        clustercompfile = outfile + "_fowlkes.txt"
                        with open(clustercompfile, "w") as out:
                            out.write("Reconstructed Clusters: " + recon_clust + "\n")
                            out.write("True Clusters: " + true_clust + "\n")
                            out.write("Fowlkes:" +  str(fowlkes) + "\n")
                        with open(outfile + ".tsv", "a") as out:
                            out.write("\t")
                            out.write(str(fowlkes))
                            out.write("\n")
                            out.close()

                        if simpath != None:

                            recon_adjacency_file = outfile + "_edges.txt"
                            recon_dist_file = outfile + "_dist_edges.txt"

                            adj_scatter_file = outfile + "_adj"
                            dist_scatter_file = outfile + "_distfs"
                            r_value_adj = scatter(sim_adjacency_file, recon_adjacency_file, "True Normalized Similarity", "Reconstructed Similarity", adj_scatter_file)
                            r_value_dist = scatter(sim_dist_file, recon_dist_file, "True Distance", "Reconstructed Distance", dist_scatter_file)


                            # os.getcwd()
                            with open(outfile + ".tsv", "a") as out:
                                out.write("\t")
                                out.write(str(r_value_dist))
                                out.write("\t")
                                out.write(str(r_value_adj))
                                out.write("\n")
                                out.close()

                        tsvs.append(outfile + ".tsv")
                        trial += 1

                        print("Done.")

    if simpath != None:
        tsv = write_tsv(tsvs, outfileprefix + "/" + str(extract_filename_from_path(outfileprefix)), True)
    elif truealignment != None:
        tsv = write_tsv(tsvs, outfileprefix + "/" + str(extract_filename_from_path(outfileprefix)), "truealign")
    else:
        print(tsvs)
        tsv = write_tsv(tsvs, outfileprefix + "/" + str(extract_filename_from_path(outfileprefix)),False,False,True)
    return tsv


def write_tsv(file_list, outfile, rval=False, truealign=False, fowlkes=False):
    '''
    Takes all the one-liners from each run and concatenates them into one file.
    '''

    header = "Inputfile\tNumSeqs\tFragmentThreshold\tSeqsAfterTrimming\tScoringMethod\tAlpha\tBeta\tEdgeAvg\tEdgeVar\tEdgeThreshold\tCDAlgorithm\tNumComms\tModularity"
    if rval == True:
        header += "\tRDistance\tRSimilarity"
    if truealign == True:
        header += "\tRTrue\tRAlign"
    if fowlkes == True:
        header += "\tFowlkes"
    header += "\n"

    headerfile = outfile + "_header.txt"
    with open(headerfile, "w") as out:
        out.write(header)

    outfile = outfile + ".tsv"
    command = 'cat '
    command += headerfile + " "
    for f in file_list:
        if os.path.exists(f):
            command += f + " "
    command += " > " + outfile
    os.system(command)

    return outfile

def cat_tsvs(file_list, outfile):
    '''
    Takes all the one-liners from each run and concatenates them into one file.
    '''

    print(file_list)
    with open(outfile, "w") as out, open(file_list[0], "r") as header:
        # all files should have the same header so we'll just take it from the first file in the list
        head = header.readline()
        out.write(head)

    for f in file_list:
        with open(f, "r") as orig, open(outfile, "a") as out:
            header = orig.readline() # just discard it
            for line in orig:
                out.write(line)

    return outfile



def extract_filename_from_path(filepath):
    if "/" in filepath:
        pathlist = filepath.split('/')
    else:
        pathlist = filepath.split("\\")
    return pathlist[-1]


def run_part_two(infile, outfileprefix, simpath=None, a=None, f=None, s=None):
    # Uses the following defaults (from the main program):
    # alpha: freq
    # fragment threshold: 0
    # scoring algorithm: mcic
    tsvs = []
    global edge_thresholds
    global cd_algorithms

    trial = 1 # add this to the end of the outfile prefix
    for e in edge_thresholds:
        for c in cd_algorithms:
            filepath = outfileprefix + "/edge_" + str(e) + "/" + str(c) + "/"
            outfile = filepath + str(extract_filename_from_path(outfileprefix)) + "_" + str(trial)
            print(outfile)
            check_file_path(outfile)

            logfile = outfile + ".log"

            command = "python main.py -i " + str(infile) + " -o " + str(outfile)
            command += " -edge " + str(e)
            command += " -c " + str(c)
            if a != None:
                command += " -alpha " + str(a)
            if f != None:
                command += " -fragment " + str(f)
            if s != None:
                command += " -s " + str(s)
                
            if args.force == True:
                command += " --force "
            command += " > "
            command += logfile
            print()
            print()
            # defaults are used for the edge threshold (0.85) and community detection (fastgreedy)
            print("Running main.py with the following parameters:")
            print("Input: " + infile)
            print("Output: " + outfile)
            print("Edge Threshold: " + str(e))
            print("Community Detection Algorithm: " + str(c))
            print("Redirecting stdout to " + logfile + ".")

            # print(command)
            print()
            print()
            print()

            os.system(command)
            trial += 1
            tsvs.append(outfile + ".tsv")

            if simpath != None:
                recon_clust = outfile + "_clusteroutput.txt"
                true_clust = simpath + "_clusteroutput.txt"
            fowlkes = compare_clusters.compare_clusters(recon_clust, true_clust)
            print("Fowlkes",fowlkes)
            print("Writing information to " + outfile + "_fowlkes.txt...")
            clustercompfile = outfile + "_fowlkes.txt"
            with open(clustercompfile, "w") as out:
                out.write("Reconstructed Clusters: " + recon_clust + "\n")
                out.write("True Clusters: " + true_clust + "\n")
                out.write("Fowlkes:" +  str(fowlkes) + "\n")
            with open(outfile + ".tsv", "a") as out:
                out.write("\t")
                out.write(str(fowlkes))
                out.write("\n")
                out.close()


            print("Done.")


    outfile = write_tsv(tsvs, outfileprefix + "/" + str(extract_filename_from_path(outfileprefix)),False,False,True)
    return outfile




def run_sim_data(n, path, part):
    '''
    Runs the simdata.py file, which will run the simulation and generate the edges.txt file.
    Then runs part one of metamain and for each output it creates a truesimilarity vs. reconstructed
    similarity graph, and a non-normalized truedist vs. recon dist graph.

    Will create n simulated datasets, labeled simdata1... simdatan, then run all part one test cases on 
    '''

    tsvs = []
    for i in range(1,n+1): # Run n simulations
        SIMDATANAME = path.split('/')[-1]  + str(i) + "true" # The last item of the path is the name of the item
        SIMDATAPATH = path + str(i) + "/" + SIMDATANAME # The name of the item will also be the name of the folder containing the outputs

        print("Running simulation " + str(i) + "...")


        # Create the simulated data
        simdata.run(SIMDATAPATH, SIMDATANAME)

        print("Done.")
        print()
        print()

        SIMFASTA = SIMDATAPATH + "/" + SIMDATANAME + ".fas"
        # print(SIMFASTA)

        print("Running reconstruction pipeline on " + SIMFASTA + "...")

        RECONDATANAME = path.split('/')[-1]  + str(i) + "recon"
        RECONDATAPATH = path + str(i) + "/" + RECONDATANAME

        # Run part one analysis on the new data
        SIMALIGN = SIMDATAPATH + "/" + SIMDATANAME + "_TRUE.fas"
        if part == 1:
            tsv = run_part_one(SIMFASTA, RECONDATAPATH, SIMDATAPATH + "/" + SIMDATANAME, SIMALIGN)
            tsvs.append(tsv)
        if part == 2:
            tsv = run_part_two(SIMFASTA, RECONDATAPATH, SIMDATAPATH + "/" + SIMDATANAME, "freq", 0, "MCIC")
            tsvs.append(tsv)

    tsv = cat_tsvs(tsvs, RECONDATAPATH + "_all.tsv")
    print("TSV",tsv)
        # True indicates it is simulated data, so it will run the scatter plots and stuff.

def scatter(truesimfile, reconsimfile, xlab, ylab, outfilename):

    truesimlist = []
    reconsimlist = []
    with open(reconsimfile, "r") as f:
        header = True
        for line in f:
            if header == True:
                header = False
            else:
                line = line.strip("\n")
                line = line.split("\t")
                reconsimlist.append(float(line[2]))
    # print(len(reconsimlist))


    with open(truesimfile, "r") as f:
        header = True
        for line in f:
            if header == True:
                header = False
            else:
                line = line.strip("\n")
                line = line.split("\t")
                truesimlist.append(float(line[2]))

    slope, intercept, r_value, p_value, std_err = stats.linregress(truesimlist,reconsimlist)
    regfile = outfilename + "_regression.txt"
    with open(regfile, "w") as f:
        f.write("Formula:\t")
        formula = "y= " + str(slope) + "x + " + str(intercept)
        f.write(formula)
        f.write("\n")
        f.write("Slope\t")
        f.write(str(slope))
        f.write("\n")
        f.write("Intercept\t")
        f.write(str(intercept))
        f.write("\n")
        f.write("r\t")
        f.write(str(r_value))
        f.write("\n")
        f.write("rsquared\t")
        f.write(str(r_value**2))
        f.write("\n")
        f.write("p\t")
        f.write(str(p_value))
        f.write("\n")
        f.write("Std. Error\t")
        f.write(str(std_err))
        f.write("\n")



    scatterfile = outfilename + "_scatter.png"
    plot = matplotlib.pyplot.scatter(truesimlist, reconsimlist)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    # plt.title('Histogram of IQ')
    plt.grid(False)

    plt.savefig(scatterfile, format="png")
    plt.close()

    return r_value

def scatter_truealign_reconalign(parentdirectory, outfiles, out_tsv):
    '''
    For a given parent directory (like Simplesimdata100), will run main.py on all of the simulated data in the directory,
    then will output two scatter plots for each: true tree similarity vs. true alignment similarity, and true alignment
    similarity vs. Clustal Omega alignment similarity.
    '''
    # Get all the subdirectories
    files = os.listdir(parentdirectory)
    dirs = []
    for f in files:
        fpath = os.path.join(parentdirectory, f)
        if os.path.isdir(fpath):
            dirs.append(f)

    tsvs = []
    with open(out_tsv, "w") as out:
        out.write("Dataset\tInputFile\tR-trues\tR-alignments\n")

        for d in dirs: # at the level d = Simplesimdata100/simplesimdata1001/
            # will create a new directory called, e.g. simplesimdata1001truealign and run part one there on 
            # the existing simulated dataset

            # I will be running the true alignment with fragment=0, scoring=MCIC, and alpha=freq
            # So I need to find the comparable clustal alignment edge file to compare to (recon)

            # each directory d will have one true file to pair will all recon files in the subdirectories.
            dpath = os.path.abspath(os.path.join(parentdirectory, d + "/" + d + "recon"))
            true_dir = parentdirectory + "/" + d + "/" + d + "true/"
            true_edges = os.path.abspath(true_dir + d +  "true_edges.txt")
            truefasta = os.path.abspath(true_dir + d + "true.fas")
            true_alignedfasta = os.path.abspath(true_dir + d + "true_TRUE.fas")

            # print(dpath)
            reconpath = os.path.join(dpath, "frag_0\\MCIC\\alpha_freq")
            # print(reconpath)
            for item in os.listdir(reconpath):
                if item[-10:] == "_edges.txt" and item[-15:] != "_dist_edges.txt":
                    recon_edges = os.path.abspath(os.path.join(reconpath,item)) # This is the reconstructed edge file
                    print(recon_edges)

            # RUN PART ONE on existing simulated data
            out = os.path.abspath(os.path.join(parentdirectory, d + "\\" + d + "truealign\\"))

            print(out)
            print("Running part one")
            tsv = run_part_one(truefasta, out, None, true_alignedfasta)
            print()
            print("TSV")
            print(tsv)
            tsvs.append(tsv)

            truealign_edges_path = os.path.join(out, "frag_0\\MCIC\\alpha_freq\\")
            for item in os.listdir(truealign_edges_path):
                if item[-10:] == "_edges.txt" and item[-15:] != "_dist_edges.txt":
                    truealign_edges = os.path.abspath(os.path.join(truealign_edges_path,item)) # This is the reconstructed edge file
                    print(truealign_edges)

            outfilename = d + "truealign_adj_truetree_truealign"
            print("truealign_edges_path",truealign_edges_path)
            print(outfilename)
            outfilepath = os.path.join(os.path.dirname(truealign_edges_path), outfilename)
            print(outfilepath)
            r_true_true = scatter(true_edges, truealign_edges, "True (Tree) Similarity", "True (Alignment) Similarity", outfilepath)

            outfilename = d + "truealign_adj_truealign_reconalign"
            print(outfilename)
            outfilepath = os.path.join(os.path.dirname(truealign_edges_path), outfilename)
            print(outfilepath)
            r_align = scatter(truealign_edges, recon_edges, "True Alignment Similarity", "Reconstructed Alignment Similarity", outfilepath)
            with open(tsv, "a") as tmp:
                tmp.write("\t")
                tmp.write(str(r_true_true))
                tmp.write("\t")
                tmp.write(str(r_align))
                tmp.write("\n")

    final_tsv = cat_tsvs(tsvs, os.path.join(parentdirectory, d + "truealign.tsv" ))


def scatter_with_threshold(truesimfile, reconsimfile, xlab, ylab, outfilename, threshold, above):
    truesimlist = []
    reconsimlist = []
    with open(reconsimfile, "r") as f:
        header = True
        for line in f:
            if header == True:
                header = False
            else:
                line = line.strip("\n")
                line = line.split("\t")
                reconsimlist.append(float(line[2]))

    with open(truesimfile, "r") as f:
        header = True
        for line in f:
            if header == True:
                header = False
            else:
                line = line.strip("\n")
                line = line.split("\t")
                truesimlist.append(float(line[2]))
                
    # print(truesimlist)
    print(len(truesimlist))
    print(len(reconsimlist))

    i = 0
    while i < len(reconsimlist):
        if above:
            # If we want points above the threshold
            if reconsimlist[i] < threshold:
                del reconsimlist[i]
                del truesimlist[i]
            else:
                i += 1
        else:
            # If we want points below the threshold
            if reconsimlist[i] >= threshold:
                del reconsimlist[i]
                del truesimlist[i]
            else:
                i += 1



    slope, intercept, r_value, p_value, std_err = stats.linregress(truesimlist,reconsimlist)
    regfile = outfilename + "_regression.txt"
    with open(regfile, "w") as f:
        f.write("Formula:\t")
        formula = "y= " + str(slope) + "x + " + str(intercept)
        f.write(formula)
        f.write("\n")
        f.write("Slope\t")
        f.write(str(slope))
        f.write("\n")
        f.write("Intercept\t")
        f.write(str(intercept))
        f.write("\n")
        f.write("r\t")
        f.write(str(r_value))
        f.write("\n")
        f.write("rsquared\t")
        f.write(str(r_value**2))
        f.write("\n")
        f.write("p\t")
        f.write(str(p_value))
        f.write("\n")
        f.write("Std. Error\t")
        f.write(str(std_err))
        f.write("\n")



    scatterfile = outfilename + "_scatter.png"
    plot = matplotlib.pyplot.scatter(truesimlist, reconsimlist)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    # plt.title('Histogram of IQ')
    plt.grid(False)

    plt.savefig(scatterfile, format="png")
    plt.close()

    return r_value

def traverse_to_bottom(parentdirectory):
    """
    Returns a list of the paths to the bottom-most directories
    """
    # print(parentdirectory)
    files = os.listdir(parentdirectory)
    # print(files)
    dirs = []
    for f in files:
        # print(f)
        fpath = os.path.join(parentdirectory, f)
        # print(fpath)
        # print(os.path.isdir(fpath))
        if os.path.isdir(fpath):
            dirs.append(fpath)
        # elif os.path.isfile(fpath):
        #     print("FILE")
    # print("DIRS", dirs)
    if len(dirs) == 0: # Base case - return path 
        # print("ACK")
        return [os.path.abspath(parentdirectory)]

    outs = []
    for d in dirs:
        outs = outs + traverse_to_bottom(d)
    return outs

def cat_subd_tsvs(parentdirectory, out_tsv):
    tsvs = []
    files = os.listdir(parentdirectory)
    dirs = []
    for f in files:
        fpath = os.path.join(parentdirectory, f)
        if os.path.isdir(fpath):
            dirs.append(f)

    for d in dirs: # at the level d = Simplesimdata100/simplesimdata1001/
        # each directory d will have one true file to pair will all recon files in the subdirectories.
        dpath = os.path.abspath(os.path.join(parentdirectory, d + "/" + d + "recon"))

        for filename in os.listdir(dpath):
            if filename[-4:] == ".tsv":
                filename = os.path.join(dpath, filename)
                print(filename)
                tsvs.append(filename)

    outfile = cat_tsvs(tsvs, out_tsv)
    print(outfile)
    return outfile

def run_threshold_scatter(parentdirectory, threshold, out_tsv): # parentdirectory = Simplesimdata100
    """
    Runs the threshold scatter plot on each trial inside the parent directory.
    """
    files = os.listdir(parentdirectory)
    dirs = []
    for f in files:
        fpath = os.path.join(parentdirectory, f)
        if os.path.isdir(fpath):
            dirs.append(f)

    with open(out_tsv, "w") as out:
        out.write("Dataset\tInputFile\tRegThreshold\tFragmentThreshold\tScoringMethod\tAlpha\tR-above\tR-below\n")
        for d in dirs: # at the level d = Simplesimdata100/simplesimdata1001/
            # each directory d will have one true file to pair will all recon files in the subdirectories.
            dpath = os.path.abspath(os.path.join(parentdirectory, d + "/" + d + "recon"))
            all_subds = traverse_to_bottom(dpath)
            for item in all_subds:
                for filename in os.listdir(item):
                    if filename[-10:] == "_edges.txt" and filename[-15:] != "_dist_edges.txt":
                        recon = os.path.abspath(os.path.join(item,filename)) # This is the reconstructed edge file
                        reconpath = recon.split("\\")
                        datasetname = d
                        print(datasetname)
                        print(reconpath)
                        for dirname in reconpath:
                            if "frag_" in dirname:
                                frag = dirname.lstrip("frag_")
                            if dirname == "MCIC" or dirname=="fifthstate":
                                score_alg = dirname
                            if "alpha_" in dirname:
                                alpha = dirname.lstrip("alpha_")

                        true = os.path.join(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(recon))),"fifthstate"),"alpha_None") # Get the directory where trimmed file is stores
                        for thing in os.listdir(true):
                            if "true_edges_trimmed.txt" in thing:
                                true = os.path.join(true, thing)
                                print(true)
                        outfilename = datasetname + "_adj_above_threshold"
                        outfilepath = os.path.join(os.path.dirname(recon), outfilename)
                        r_above = scatter_with_threshold(true, recon, "True Normalized Similarity", "Reconstructed Similarity", outfilepath, threshold, True)
                        # print("r above", r_above)

                        outfilename = datasetname + "_adj_below_threshold"
                        outfilepath = os.path.join(os.path.dirname(recon), outfilename)
                        r_below = scatter_with_threshold(true, recon, "True Normalized Similarity", "Reconstructed Similarity", outfilepath, threshold, False)
                        # print("r below", r_below)
                        out.write(str(datasetname))
                        print(datasetname)
                        out.write("\t")
                        out.write(str(recon))
                        print(recon)
                        out.write("\t")
                        out.write(str(threshold))
                        print(threshold)
                        out.write("\t")
                        out.write(str(frag))
                        print(frag)
                        out.write("\t")
                        out.write(str(score_alg))
                        print(score_alg)
                        out.write("\t")
                        out.write(str(alpha))
                        print(alpha)
                        out.write("\t")
                        out.write(str(r_above))
                        print(r_above)
                        out.write("\t")
                        out.write(str(r_below))
                        print(r_below)
                        out.write("\n")


def log_val(edgefile):
    '''
    Accepts an edge file example.txt and returns a new edge file example_log.txt with the log of each value.
    '''
    edgefilepath = os.path.dirname(edgefile)
    edgefilename = os.path.basename(edgefile)
    edgefileprefix = edgefilename.split(".")[0]
    edgefilesuffix = edgefilename.split(".")[1]

    logedgefile = os.path.join(edgefilepath, edgefileprefix + "_log." + edgefilesuffix)
    print(logedgefile)


    with open(edgefile, "r") as orig, open(logedgefile, "w") as new:
        header = orig.readline()
        new.write(header)
        for line in orig:
            line = line.strip("\n")
            line = line.split("\t")
            head = line[0]
            tail = line[1]
            val = float(line[2])
            if val != 0:
                logval = math.log(val, 10)
            else:
                logval = math.log(0.01, 10)
            new.write(head)
            new.write("\t")
            new.write(tail)
            new.write("\t")
            new.write(str(logval))
            new.write("\n")
    return logedgefile

def log_scatter(file1, file2, xlab, ylab, outfile):
    '''
    Takes two edge files, generates log edge files for both datasets, then creates a scatter plot of those.
    '''

    logfile1 = log_val(file1)
    logfile2 = log_val(file2)
    rval = scatter(logfile1, logfile2, xlab, ylab, outfile)


def main():

    if args.truealignment == True:
        scatter_truealign_reconalign(args.infile[0], args.outfile[0], "./test_truealign_5716.tsv")

    nsims = int(args.sim[0])
    if nsims > 0:
        if args.partone == True:
            tsv = run_sim_data(nsims, args.outfile[0], 1)
        if args.parttwo == True:
            tsv = run_sim_data(nsims, args.outfile[0], 2)
        else:
            print("Please give an option for part I or part II.")

    else:
        if args.partone:
            if args.parttwo:
                tsv = run_one_two(args.infile[0],args.outfile[0])
            else:
                tsv = run_part_one(args.infile[0],args.outfile[0])
        else:
            tsv = run_part_two(args.infile[0],args.outfile[0])


main()

