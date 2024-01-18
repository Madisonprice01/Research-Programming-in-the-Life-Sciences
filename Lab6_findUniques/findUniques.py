#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Name: Madison Price (maeprice)
# Group Members: Kiana Imani (kimani), Queenie Li (qli86), and Ellizabeth Beer (ebeer)

import sys


class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence.'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                if not line:  # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


########################################################################

class TRNA:
    '''Analyze a portion of a fastA file containing a tRNA sequence and generate its essential subsequences.

    Input: tRNA header and sequence
    Output: 'Essential' substrings, formatted to be aligned with where they occur in the sequence.

    '''

    tRNAlist = []

    def __init__(self, head, seq):
        '''Initialize a TRNA object with all attributes given at instantiation time.'''
        self.head = head.replace(' ', '')  # get rid of white space in sequence
        self.seq = seq.replace('.', '').replace('_', '').replace('-', '')  # get rid of undesired characters in sequence
        self.powerset = self.getPowerset()
        self.uniques = set()
        self.essentials = set()
        self.nonessentials = set()

    def getPowerset(self):
        '''Generate the powerset of a given tRNA sequence.'''
        powerset = set()
        for nucleotide in range(
                len(self.seq)):  # for each nucleotide, generate all possible subsequences starting at that nucleotide
            for substringEnd in range(nucleotide + 1, len(self.seq) + 1):
                powerset.add(self.seq[nucleotide:substringEnd])
        return powerset

    def getUniques(self):
        '''Generate the unique subsequences of a tRNA sequence.'''
        others = set()  # initialize accumulating union of 'other' powersets
        for tRNA in TRNA.tRNAlist:
            if tRNA != self:  # if tRNA = 'other'
                others = others.union(tRNA.powerset)  # add this 'other's powerset to accumulating union
        self.uniques = self.powerset - others
        return self.uniques

    def getEssentials(self):
        '''Find the shortest unique subsequences aka 'essetials' in a tRNA sequence. '''
        for unique in self.uniques:
            if unique[
               :-1] in self.uniques:  # test if shortening sequence by one nuc yields a subsequnce that is still unique
                self.nonessentials.add(unique)  # if still unique, sequence is nonessential
            elif unique[1:] in self.uniques:
                self.nonessentials.add(unique)
            else:
                pass
        self.essentials = self.uniques - self.nonessentials  # generate uniques by subtracting the nonessential set from the unique set
        return self.essentials


#######################################################################
# Main
# Here is the main program
#
########################################################################
# from sequenceAnalysis import FastAreader

'''
This program uses the FastAreader class and TRNA class to read tRNA headers and sequences from a FastA file and 
output all of the subsequences that are 'essential', or the smallest unique subsequences. 

Input: FastA reader 
Output: Formatted essential subsequences
'''


def main(inCL=None):
    '''For every tRNA header and sequence in a FastA file, generate the essential substrings and output them in
    the desired format.

    Input: FastA file
    Output: Essential substrings for each tRNA sequence
    '''

    myFile = FastAreader()

    for head, seq in myFile.readFasta():  # iterate through tRNA headers and sequences
        myTRNA = TRNA(head, seq)  # instantiate a TRNA class object for every header, sequence
        TRNA.tRNAlist.append(myTRNA)

    TRNA.tRNAlist.sort(key=lambda x: x.head)  # alphabetize tRNA names

    for tRNA in TRNA.tRNAlist:  # iterate through sorted tRNA list
        print(tRNA.head)
        print(tRNA.seq)
        tRNA.getPowerset()
        tRNA.getUniques()
        essentials = tRNA.getEssentials()

        essentialPositions = []  # instantiate a list to store essentials and their positions
        for essential in essentials:
            position = 0
            while position < len(tRNA.seq):  # iterate through tRNA sequence
                try:
                    position = tRNA.seq.index(essential,
                                              position)  # find each occurence of a essential, changing starting position each time you find one
                    essentialPositions.append((position, essential))
                except ValueError:
                    break
                position += 1
            # position = tRNA.seq.find(essential)
            # essentialPositions.append((position,essential))
        sortedEssentials = sorted(essentialPositions, key=lambda x: x[
            0])  # sort list of tuples based on the first field in each tuple (positions)

        for essential in sortedEssentials:  # print the appropriate number of periods for each essential, followed by the esential itself
            print('.' * essential[0], end='')
            print(essential[1])


if __name__ == "__main__":
    main()
