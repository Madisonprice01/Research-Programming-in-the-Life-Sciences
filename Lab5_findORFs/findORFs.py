#!/usr/bin/env python3
# Name: Madison Price
# Group Members: Kiana Imani (kimani), Queenie Li (qli86), and Ellizabeth Beer (ebeer)


# Design:
#
# OrfFinder Class
#  init() method takes optional arguments that specificy requirements of gen candidates
#  findGenes() method generates a list of gene candiates for top and bototm strand as follows:
#  For top strand:
#   -iterate through top strand frames (1 to 3)
#   -for each frame:
#      -record first start position as at beginning of strand
#      -iterate through codons, checking if each is a start codon or stop codon
#      -if start codon, append position to start codon position list (if start codon found at beginning of sequence, don't append position twice)
#      -if stop codon, we have reached the end of an orf and we iterate through all start coodn positions in start codon position list
#      -for each start codon in list, generate gene candidate and append its information to gene candidate list
#      -if longestGene program option is on, only record largest candidate (candidate corresponding to first start codon)
#      -clear stop and start position lists and continue iterating through codons
#      -repeat above process of generating gene candidates for every orf
#      -if we never encounter a stop codon, we have dangling starts (at least one)
#      -iterate through those dangling starts, taking the end of the sequence as a stop codon position, and generate gene candidates
#      -dangling stops are automatically handled since a start position is always recorded at beginning of sequence anyway
#  For bottom strand:
#   -iterate through bottom strand frames (-1 to -3)
#   -for each frame:
#     -execute top strand algorithm, only difference being how we record stop and start positons
#     -record all bottom strand positions with respect to/in terms of the top strand. i.e. use len(top strand)-(where we are in bottom strand) as position
#  Sort the gene candidate list based on gene length, and subsort based on start positions
#  Return sorted list
#
# Main() Function
#  Instantiate a commandLine class object and pass command line 'options' as its input
#  Instantiate a fastAreader class object with main()'s input file as input
#  Call the readFasta() method on the fastAreader object and iterate through each header and sequence in input file
#  For each sequence:
#    -generate an OrfFinder class object and call the findGenes() method that builds a list of its gene candidates
#    -print sequence header and sorted gene candidates
#
# Design Decisions: use same 'queuing' algorithm for top and bottom strand, and record all positions (top and bottom strand) with respect to the top strand


########################################################################
# CommandLine
########################################################################

class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store', nargs='?', const=True, default=False,
                                 help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(0, 100, 200, 300, 500, 1000), default=100,
                                 action='store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action='append', default=['ATG'], nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-t', '--stop', action='append', default=['TAG', 'TGA', 'TAA'], nargs='?',
                                 help='stop Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)



########################################################################
# Main
# Here is the main program
########################################################################

""" This program imports sequenceAnalysis module and fins gene candidates in DNA sequences from an input file in fastA format.

    Input: DNA sequences from a FastA file
    Outputs: Sorted list of possible gene candidates.
"""

from sequenceAnalysis import FastAreader, OrfFinder


def main(inFile='', options=None):
    '''
    For every DNA sequence in a FASTA file, find all possible genes and return them in a sorted list.

    Input: FastA file, program options
    Output: Formatted list of all gene candidates in every sequence of a FastA file
    '''
    thisCommandLine = CommandLine(options)

    ###### replace the code between comments.
    # thisCommandLine.args.longestGene is True if only the longest Gene is desired
    # thisCommandLine.args.start is a list of start codons
    # thisCommandLine.args.stop is a list of stop codons
    # thisCommandLine.args.minGene is the minimum Gene length to include
    #
    #######

    myFile = FastAreader(inFile)  # instantiate fastAreader object with infile
    for header, sequence in myFile.readFasta():  # call readFasta() method to iterate through sequences in file
        Orfs = OrfFinder(sequence, thisCommandLine.args.longestGene, thisCommandLine.args.start,
                         thisCommandLine.args.stop, thisCommandLine.args.minGene)  # instantiate OrfFinder object
        geneList = Orfs.findGenes()  # call findGenes() method to gnerate each sequence's gene candidates
        print(header)
        for gene in geneList:  # print every gene candidate, formatted nicely
            print('{:s} {:>5d}..{:>5d} {:>5d}'.format(gene[0], gene[1], gene[2], gene[3], ))


if __name__ == "__main__":
    main()

