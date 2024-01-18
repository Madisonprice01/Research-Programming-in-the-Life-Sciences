#!/usr/bin/env python3
# Name: Rakshya U. Sharma (rausharm)
# Group Members: Kristel Caspe (kcaspe), Gillian Bowers (gfbowers), Vrishab Nukala

import sys
from sequenceAnalysis import NucParams, FastAreader

# Calculate the GC count of the given nucleotide composition
def calculate_gc_content(nuc_composition):
    """
    Calculates the GC content of a given nucleotide composition.
    A dictionary containing the counts of each nucleotide.
    """
    #count the total number of G's and C's
    gc_count = nuc_composition['G'] + nuc_composition['C']
    #count the total number of nucleotides
    total_count = sum(nuc_composition.values())
    #calculate the GC content percentage
    return (gc_count / total_count) * 100

def analyze_genome(filename = None):
    """
    Analyzes a genome by calculating sequence length, GC content,amino acid composition, codon bias, and codon composition. 
    The genome is read from STDIN. A dictionary containing the amino acid composition and a dictionary containing the codon bias.
    """
     # Create a FastAreader object to read from STDIN (no filename specified)
    myReader = FastAreader(filename)
    # Create a NucParams object to analyze the nucleotide and codon composition
    myNuc = NucParams()
    
    ## Read each sequence from the input FASTA file and add it to the NucParams object
    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)
    # Calculate the sequence length in Mb (millions of bases)- Convert bases to Mb
    sequence_length = myNuc.nucCount() / 1000000 
     # Calculate the GC content percentage
    gc_content = calculate_gc_content(myNuc.nucComposition())
    # Calculate the amino acid composition using NucParams's object aaComposition
    aa_composition = myNuc.aaComposition()
    #convert the amino acid compotiion into a dictionary 
    aaComp = dict(aa_composition)
    #Initia;ize an empty dictionary to store codon bias
    codon_bias = {}
    #get the codon composition of the sequence using NucParam's codon composition method
    cc_composition = myNuc.codonComposition()
    #convert the returned compositin into dictionary
    ccComp = dict(cc_composition)

    # Sort the codons in alphabetical order
    #codon_list = sorted(myNuc.rnaCodonTable.items(), key=lambda x: x[0])
    codon_list = sorted(myNuc.rnaCodonTable.items(), key = lambda item:(item[1], item[0])) #Dr.B's method
    # Iterate through each codon and its corresponding amino acid
    for codon, aa in codon_list:
        aa_count = myNuc.aaComposition()[aa]# Get the count of the current amino acid
        codon_count = ccComp[codon]# Get the count of the current codon
        # Calculate the relative codon frequency
        codon_frequency = (codon_count / aa_count) * 100
        #adding the codon_frequency value to codon_bias dictionary
        codon_bias[codon] = codon_frequency
      
    return sequence_length, gc_content, aaComp, codon_bias, ccComp
    
def main():
    """
    Analyze the given genome file and print the results.
    """
    #calling the analyze_genome function 
    sequence_length, gc_content, aaComp, codon_bias, ccComp  = analyze_genome()
    print(f"sequence length = {sequence_length:.2f} Mb\n")#print the length in mega base with two decimal points
    print(f"GC content = {gc_content:.1f}%\n")#print GC content percentage with one decimal place
    # Iterate through each codon and its corresponding frequency in the codon_bias dictionary
    for codon, codon_frequency in codon_bias.items():
        #using rnaCodonTable attribute from NucParams to get the corresponding amino acid
        aa = NucParams().rnaCodonTable[codon]
        #get the count of current amino acid from aaComp dictionary
        aa_count = aaComp[aa]
        #get the frequency of current codon from codon_bias dictionary
        codon_frequency = codon_bias[codon]
        #get the count of current codon from ccComp Dictionary
        codon_count = ccComp[codon]
        
        
        # Print the codon, amino acid, codon frequency, and codon count
        print(f"{codon} : {aa} {codon_frequency:5.1f} ({codon_count:6d})")

if __name__ == "__main__":
    main()