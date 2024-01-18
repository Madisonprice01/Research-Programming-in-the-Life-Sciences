#!/usr/bin/env python3
# Name: Madison Price (maeprice)
# Group Members: Kiana Imani (kimani), Queenie Li (qli86), and Ellizabeth Beer (ebeer)

""" This module contains three classes: ProteinParam, FastAreader, NucParams, and OrfFinder. ProteinParam returns metrics of an
    amino acid sequence's physical and chemical properties. FastAreader reads files of fastA format. NucParams returns
    data containers that depict the content of a nucleotide sequence. OrfFinder finds all of the gene candidates in a given DNA sequence.

    Input: Amino acid sequence, fastA file, or nucleotide sequences
    Outputs: Measurements of an amino acid sequence's properties, header and content of a fastA file, and/or dictionaries
             containing nucleotides, codon, and amino acid frequencies in nucleotide sequences, and sorted lists of gene candidates.
"""

#-------------------------------------------------------------------------------------------------------------------------------------------------
class ProteinParam:
    """Build a ProteinParam object with all of the attributes given at instantiation time."""

    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein):
        """Initialize a ProteinParam object from a protein sequence and build a dictionary object to hold information about its instance arguments."""
        aaSequence = protein.upper()  # uppercase input
        self.aaTable = {}  # initialize table of aa's
        for aa in 'ACDEFGHILKMNPQRSTVYW':
            self.aaTable[aa] = 0
        for aa in aaSequence:
            if aa in self.aaTable:  # check for unexpected letters in protein
                self.aaTable[aa] += 1  # determine freqeuency of each aa in protein

    def aaCount(self):
        """Return a single integer count of valid amino acid characters found in the given protein."""
        count = 0
        for aa in self.aaTable:
            if self.aaTable[aa] != 0:  # ignore aa that weren't found in protein
                count += self.aaTable[aa]  # add the count of each found aa to accumulating sum
        return count

    # def pI (self):                                             # original pI() method
    #     """ Find the  particular pH that yields a neutral net Charge of the given protein."""
    #     if all(aa == 0 for aa in self.aaTable.values()):       # account for empty input
    #         return 0.00
    #     else:
    #         bestCharge = 100000000                             # initialize charge that will not possibly be teh best
    #         for pH100 in range(1400+1):                        # set range to 14x100 then divide by 100 to get 2 decimal places
    #             pH = pH100/100
    #             thisCharge = abs(self._charge_(pH))            # use charge() method to find charge of protein at a certain pH
    #             if thisCharge < bestCharge:                    # as you iterate through range of pH's, test if charge at each pH is closer to zero than the last
    #                 bestCharge = thisCharge
    #                 bestpH = pH
    #         return bestpH

    def pIBinary(self, precision=2):  # extra credit for pI() method
        """Find the  particular pH that yields a neutral net Charge of the given protein using binary search."""
        if all(aa == 0 for aa in self.aaTable.values()):  # account for empty input
            return 0.00
        else:
            high = 14.00  # set upper bound
            low = 0.00  # set lower bound
            while (high - low > 0.01):  # while we havent met degree of precision
                midPoint = (high + low) / 2  # find midpoint within bounds
                charge = self._charge_(midPoint)  # calculate charge at midpoint
                if charge > 0:  # if charge above midpoint, set new range as midpoint to upper bound
                    low = midPoint
                else:  # if charge below midpoint, set new range as midpoint to lower bound
                    high = midPoint
            return midPoint

    def aaComposition(self):
        """Return the dictionary created at instantiation time, keyed by single letter Amino acid code, and having associated values that are the counts of those amino acids in the sequence."""
        return self.aaTable  # dictionary containing aa frequencies was built in init() method

    def _charge_(self, pH):
        """Calculate the net charge on the protein at a specific pH."""
        group1 = ['R', 'K', 'H']
        charge1 = (10 ** ProteinParam.aaNterm / (10 ** ProteinParam.aaNterm + 10 ** pH))
        for aa in group1:  # calculate first term of netcharge equation
            charge1 += self.aaTable[aa] * (
                        10 ** ProteinParam.aa2chargePos[aa] / (10 ** ProteinParam.aa2chargePos[aa] + 10 ** pH))
        group2 = ['D', 'E', 'C', 'Y']
        charge2 = (10 ** pH / (10 ** ProteinParam.aaCterm + 10 ** pH))
        for aa in group2:  # calculate second term of netcharge equation
            charge2 += self.aaTable[aa] * (10 ** pH / (10 ** ProteinParam.aa2chargeNeg[aa] + 10 ** pH))
        netCharge = charge1 - charge2
        return netCharge

    def molarExtinction(self, ):
        """Calculate the molar extinction coefficient that indicates how much light a protein absorbs at a certain wavelength."""
        molarExtinctionCoeff = (self.aaTable['Y'] * ProteinParam.aa2abs280['Y']) + (
                    self.aaTable['W'] * ProteinParam.aa2abs280['W']) + (self.aaTable['C'] * ProteinParam.aa2abs280[
            'C'])  # sum molar extinction coefficens of T,W, and C in protein
        return molarExtinctionCoeff

    def massExtinction(self):
        """Calculate the Mass extinction coefficient from the Molar Extinction coefficient by dividing by the molecularWeight of the given protein."""
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0  # def of mass extinction

    def molecularWeight(self):
        """Calculate the molecular weight of the protein."""
        if all(aa == 0 for aa in self.aaTable.values()):
            return 0.00
        else:
            mwWithWater = 0
            for aa in self.aaTable:  # find total MW of all aa in protein
                mwWithWater += ProteinParam.aa2mw[aa] * self.aaTable[aa]
            mwWithoutWater = mwWithWater - (
                        ProteinParam.mwH2O * (self.aaCount() - 1))  # subtract water MW lost in peptide bond formation
            return mwWithoutWater

#----------------------------------------------------------------------------------------------------------------------------------------------------------------

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
        if self.fname == '':                                                                 # is or == ?
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
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

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

class NucParams:
    """Build a NucParams object with all of the attributes given at instantiation time."""
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}

    def __init__(self, inString=''):
        """Initialize a NucParam object from a nucleotide sequence and create dictionary objects to hold information about its instance arguments."""
        self.nucComp = {}  # instantiate dictionaries for nuc, aa, and codon counts
        for i in 'ACTGNU':  # initiate keys for valid nucs
            self.nucComp[i] = 0
        self.codonComp = {}
        self.aaComp = {}
        for key in NucParams.rnaCodonTable:
            self.codonComp[key] = 0
            self.aaComp[NucParams.rnaCodonTable[key]] = 0
        self.addSequence(inString)  # pass sequence to addSeq to build dictionaries

    def addSequence(self, inSeq):
        """Accept additional nucleotide sequences and add to the data given to the init method."""
        sequence = inSeq.upper()
        rnaSequence = sequence.replace('T', 'U')
        for nuc in rnaSequence:  # iterate through given sequence and test if each character is valid nuc
            if nuc in 'ACTGNU':  # if so, add one to its the count
                self.nucComp[nuc] += 1
        for i in range(0, len(rnaSequence), 3):  # iterate through sequence three characters at a time
            triplet = rnaSequence[i:i + 3]
            if triplet in NucParams.rnaCodonTable:  # check if triplet is a valid codon
                self.codonComp[triplet] += 1
                self.aaComp[NucParams.rnaCodonTable[triplet]] += 1

    def aaComposition(self):
        """Return a dictionary containing counts of the 20 amino acids in the sequence."""
        return self.aaComp  # return aa dictionary

    def nucComposition(self):
        """Return a dictionary of counts of valid nucleotides found in the sequence."""
        return self.nucComp  # return nuc dictionary

    def codonComposition(self):
        """Return a dictionary of counts of valid codons found in the sequence (in RNA format)."""
        return self.codonComp  # return codon dictionary

    def nucCount(self):
        """Return an integer value equal to the number of valid nucleotides found in the sequence."""
        nucCount = sum(self.nucComp.values())  # sum values of nuc dictionary
        return nucCount  # return sum

#-------------------------------------------------------------------------------------------------------------------------------------------------------------


class OrfFinder():
    '''
    Analyze  a sequence of DNA and find the possible ORFs (start and stop codons).

    Input: DNA sequence and program options
    Output: list of gene candidates that adhere to the program options, sorted based on length and start position

    '''

    def __init__(self, sequence, lG=False, s={'ATG'}, t={'TAG', 'TAA', 'TGA'}, mG=100):
        '''Initialize an OrfFinder object with all attributes given at anstantiation time.'''
        self.topStrand = sequence
        self.startCodons = s
        self.stopCodons = t
        self.startCodonPositions = []
        self.stopCodonPositions = []
        self.geneCandidates = []
        self.bottomStrand = self.topStrand.translate(str.maketrans('ATCGN', 'TAGCN'))[::-1]
        self.longestGene = lG
        self.minGene = mG

    def findGenes(self):
        ''' Iterate through a DNA sequence and build a list of all possible gene candidates. '''

        for frame in range(1, 4):  # iterate through top strand frames
            self.startCodonPositions.append(1)  # record start position at beginning of top strand
            for i in range(frame - 1, len(self.topStrand), 3):  # iterate through codons
                codon = self.topStrand[i:i + 3]
                if codon in self.startCodons:  # check if codon is a start codon
                    if i == 0:  # if start codon is at position 1, don't append twice
                        pass
                    else:
                        self.startCodonPositions.append(i + 1)
                if codon in self.stopCodons:  # check if codon is a stop codon
                    self.stopCodonPositions.append(
                        i + 3)  # record position of stop codon as index of its last nucleotide
                    for startPosition in self.startCodonPositions:  # generate gene candidate for every start codon found
                        geneLength = self.stopCodonPositions[0] - startPosition
                        if geneLength > self.minGene:
                            self.geneCandidates.append(['+' + str(frame), startPosition, self.stopCodonPositions[0],
                                                        geneLength + 1])  # record gene candidate
                            if self.longestGene == True:  # if longest gene option is 'on', stop after recording first gene candidate
                                break
                    self.startCodonPositions = []  # clear start and stop position lists before next frame iteration
                    self.stopCodonPositions = []
            for startPosition in self.startCodonPositions:  # record gene candidates for dangling starts (no stop encountered)
                geneLength = len(self.topStrand) - startPosition
                if geneLength > self.minGene:
                    self.geneCandidates.append(['+' + str(frame), startPosition, len(self.topStrand),
                                                geneLength + 1])  # record as a gene candidate
                    if self.longestGene == True:  # if longest gene option is 'on', stop after recording first candidate
                        break
            self.startCodonPositions = []  # clear start and stop position lists before next frame iteration
            self.stopCodonPositions = []

        for frame in range(-1, -4, -1):  # iterate through bottom strand frames
            self.startCodonPositions.append(
                len(self.bottomStrand))  # record start position at beginning of bottom strand (aka end of top strand)
            for i in range(abs(frame) - 1, len(self.bottomStrand), 3):  # iterate through codons
                codon = self.bottomStrand[i:i + 3]
                if codon in self.startCodons:  # check if codon is a start codon
                    if i == 0:  # if start codon is at position 1, don't append twice
                        pass
                    else:
                        self.startCodonPositions.append((len(self.topStrand) - i))
                if codon in self.stopCodons:  # check if codon is a stop codon
                    self.stopCodonPositions.append((
                                                               len(self.topStrand) - i - 2))  # record position of stop codon as index of its last nucleotide, but wrt top strand
                    for startPosition in self.startCodonPositions:  # generate gene candidate for every start codon found
                        geneLength = abs(self.stopCodonPositions[0] - startPosition)
                        if geneLength > self.minGene:
                            self.geneCandidates.append([str(frame), self.stopCodonPositions[0], startPosition,
                                                        geneLength + 1])  # record as a gene candidate
                            if self.longestGene == True:  # if longest gene option is 'on', stop after recording first gene candidate
                                break
                    self.startCodonPositions = []  # clear start and stop position lists before next frame iteration
                    self.stopCodonPositions = []
            for startPosition in self.startCodonPositions:  # record gene candidates for dangling starts (no stop encountered)
                geneLength = startPosition
                if geneLength > self.minGene:
                    self.geneCandidates.append([str(frame), 1, startPosition, geneLength])  # record as a gene candidate
                    if self.longestGene == True:  # if longest gene option is 'on', stop after recording first candidate
                        break
            self.startCodonPositions = []  # clear start and stop position lists before next frame iteration
            self.stopCodonPositions = []

        sortedGenes = sorted(self.geneCandidates, key=lambda x: (x[-1], -x[1]),
                             reverse=True)  # sort gene candidates first by length and then by start codon position
        return sortedGenes


