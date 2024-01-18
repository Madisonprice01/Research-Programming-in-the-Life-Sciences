#!/usr/bin/env python3
# Name: Madison Price (maeprice)
# Group Members: Kiana Imani (kimani), Queenie Li (qli86), and Ellizabeth Beer (ebeer)

""" This program imports sequenceAnalysis module and perform an analysis of nucleotide sequences from an input file in fastA format.

    Input: Nucleotide sequence
    Outputs: Sequence Length (Mb)
             GC content (%)
             Statistics on relative codon usage for each codon.
"""
from sequenceAnalysis import FastAreader, NucParams  # import classes from module


def main(fileName=None):
    """Use the sequenceAnalysis module to read a nucleotide sequence from a FastA file, and print an analysis of its content."""
    myReader = FastAreader(fileName)  # initiate fastA reader
    myNuc = NucParams()  # initiate nucParams object
    for head, seq in myReader.readFasta():  # read from fastA reader
        myNuc.addSequence(seq)  # give data to nucParams object

    codons2aa = sorted(myNuc.rnaCodonTable.items(),
                       key=lambda x: (x[1], x[0]))  # sort codons in alpha order, by Amino Acid, using lambda function
    # ^^ creates a list of tuples
    print('sequence length = {:.2f} Mb'.format(myNuc.nucCount() / 1000000))  # mega = 10^6
    print()
    print('GC content  = {:.1f}%'.format(((myNuc.nucComp['G'] + myNuc.nucComp['C']) / myNuc.nucCount()) * 100))
    print()

    for x in codons2aa:
        codon = x[0]
        aa = x[1]
        if myNuc.aaComp[aa] != 0:
            val = myNuc.codonComp[codon]/myNuc.aaComp[aa]  # calculate relative codon usage for each codon and print
            print('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, val * 100, myNuc.codonComp[codon]))
        else:
            print('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, 0, myNuc.codonComp[codon]))

if __name__ == "__main__":
    main('testGenome.fa')                                                 # default is to use stdin