from Bio import SeqIO

def parse_rmdata(infile, outfile, TEtype, numseqs):
    # takes a fasta file as input with a TE type as a string
    # and saves the sequences containing that string in the id
    # to a new file.

    TEs = []
    with open(infile, "r") as seqfile:
        n = 0
        hits = 0
        for seq in SeqIO.parse(seqfile, "fasta"):
            if hits >= numseqs:
                break
            else: 
                n += 1
                if TEtype.lower() in seq.id.lower():
                    # if the string appears, save to list
                    TEs.append(seq)
                    hits += 1
        print(hits)
        print(numseqs)

    # get name of outfile: TEtype_inputfile
    
    SeqIO.write(TEs, outfile, "fasta")
    # writes seqs to new file.
    print("Sequences: ", n)
    print("Matches: ", hits)
    print("Saved " + TEtype + " seqs to " + outfile + '.')




parse_rmdata('../Data/HSapiens_chr9repeats.fasta','../Data/Alus_500_chr9.fasta', 'Alu', 501)
