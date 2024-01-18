#!/usr/bin/env python3
# Name: Mark Forbush (mforbush)
# Group Members: None
'''
Contains the NucParams class

example: 
np = NucParams("LLVMVVC")

Inputs should be formatted similar as to the above example
'''
'''
Contains the FastAreader class

example: 
reader = FastAreader("file.txt")

Inputs should be formatted similar as to the above example
'''
'''
Take in raw protein sequence and return calculations for composition and properties

example: 
pp = ProteinParam("LLVMVVC")

Inputs should be formatted similar as to the above example
'''

########################################################################
# NucParams
########################################################################
class NucParams :
  '''
  Calculate angles and distances among a triad of points.

  Author: Mark Forbush
  Date: May 8, 2023
  Nucleotide sequences are supplied as strings and then cleaned by the class

  Required Modules: math
  initialized: single nucleotide sequence or nothing
            np = NucParams("LLVMVVC")
  attributes: inString : the raw nucleotide string
              codonComp : codon composition dicionary stored in memory
              aaComp : protein composition diciotnary stored in memory
              nucComp : nucleotide composition dictionary stored in memory
  methods:  addSequence(inSeq) : adds a sequence of nucleotide to NucParams
                                 inSeq : string of nucleotides
            nucComposition() : return a dictionary of nucleotides and their count
            aaComposition() : return a dictionary of amino acids and their count
            codonComposition() : return a dictionary of codons and their count
            nucCount() : return the amount of nucleotides
            aaCodonDecomp(dense=False) : return a dictionary of the protein-codon decomposition
                                         dense : boolean indicates should include zero valued entries
  '''
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
  'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
  }
  dnaCodonTable = { key.replace('U','T'):value for key, value in rnaCodonTable.items() }

  def __init__ (self, inString='') :
    ''' Initializes a new NucParams by creating dictionaries and adding the contents
    of inString into them.
    '''
    self.codonComp = { codon:0 for codon in list(self.rnaCodonTable.keys()) }
    self.aaComp = { aa:0 for aa in list(self.rnaCodonTable.values()) }
    self.nucComp = { nuc:0 for nuc in 'ACGTUN' }

    self.addSequence(inString)
      
  def addSequence (self, inSeq) :
    ''' This method must accept new sequences, from the {ACGTUN} alphabet, 
    each sequence can be presumed to start in frame 1. This data must be 
    added to the data that you were given with the init method (if any).
    '''
    cleanSeq = ''.join([ nuc for nuc in inSeq.upper() if nuc in 'ACGTUN' ]).replace('T','U')
    for nuc in 'ACGTUN' : self.nucComp[nuc] += cleanSeq.count(nuc)
    
    codons = [ cleanSeq[i*3:i*3+3] for i in range(int(len(cleanSeq)/3)) ]
    for codon in self.codonComp : self.codonComp[codon] += codons.count(codon)
    
    aminos = [ (self.rnaCodonTable[c] if c in self.rnaCodonTable else '') for c in codons ]
    for aa in self.aaComp : self.aaComp[aa] += aminos.count(aa)

  def aaComposition(self) :
    ''' returns the composition dictionary for amino acids '''
    return self.aaComp

  def nucComposition(self) :
    ''' returns the composition dictionary for nucleotides '''
    return self.nucComp

  def codonComposition(self) :
    ''' returns the composition dictionary for codons '''
    return self.codonComp

  def nucCount(self) :
    ''' returns the count of nucleotides '''
    return sum(self.nucComp.values())

  def aaCodonDecomp(self, dense=False) :
    ''' returns a dictionary of the amino acid codon breakdown. The data is formatted 
    dict = {P{ACT:2,GGG:5}, Y{GGC:3}}, where every amino acid has a dictionary of 
    codons associated with a value. The parameter dense=False indicates entries not 
    present should be left out of the dictionary onn return. Any other value for dense
    will return possible values not present as 0. Used in the extra credit.
    '''
    codons = { c:val for c,val in self.codonComposition().items() if val > 0 or dense}
    aminos = { aa:val for aa,val in self.aaComposition().items() if val > 0 or dense}

    aminoGraph = { aa:{} for aa in aminos } # create the amino graph and sort by alpha
    for codon,val in codons.items() :
      aa = NucParams.rnaCodonTable[codon] # amino acid come from rna codon table always
      aminoGraph[aa][codon] = val # every amino has codons and every codon has a value
    aminoGraph = dict(sorted({ aa:dict(sorted(codons.items())) for aa,codons in aminoGraph.items() }.items()))

    return aminoGraph







########################################################################
# FastAreader
########################################################################
import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence








########################################################################
# ProteinParam
########################################################################
import math
class ProteinParam :
  '''
  Calculate angles and distances among a triad of points.

  Author: Mark Forbush
  Date: April 24, 2023
  Protein sequences are supplied as strings and then cleaned by the class

  Required Modules: math
  initialized: single protein sequence
            pp = ProteinParam("LLVMVVC")
  attributes: protein is the raw protein string
  methods:  aaCount() : returns the count of amino acids
            pI(p=2) : calculate the theoretical isolelectric point
            aaComposition() : return single letter Amino acid code and associated counts
            _charge_(pH) : calculates the net charge on the protein at pH
            molarExtinction() : return how much light a protein absorbs at a certain wavelength
            massExtinction() : return mass extinction coefficient
            molecularWeight() : calculate the molecular weight
  '''

  # These tables are for calculating:
  #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
  #     absorbance at 280 nm (aa2abs280)
  #     pKa of positively charged Amino Acids (aa2chargePos)
  #     pKa of negatively charged Amino acids (aa2chargeNeg)
  #     and the constants aaNterm and aaCterm for pKa of the respective termini
  # As written, these are accessed as class attributes, for example:
  # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

  aa2mw = { # molecular weights of amino acids
    'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
    'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
    'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
    'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

  mwH2O = 18.015 # molecular weight of water
  aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125} # absorbtion of 280nm light

  aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6} # pKa of positive aminos
  aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10} # pKa of negative aminos
  aaNterm = 9.69 # N-terminus
  aaCterm = 2.34 # C-terminus

  alphabet = "ACDEFGHILKMNPQRSTVYW" # valid amino acids

  def __init__ (self, protein=None): # init function
    ''' Holds and processes a sequences of amino acids
    initialized: (protein) single string of amino acids
    attributes: protein the amino acid sequence    
    We're going to clean the protein sequence before it is stored in the 
    class as a string and dictionary. That way, we avoid having to parse it 
    multiple times. '''

    if protein is None : protein = ''
    protein = ''.join([c for c in protein.upper() if c in self.alphabet]) # clean protein
    self.protein = protein # store the protein

    self.proteinCount = {aa : self.protein.count(aa) for aa in self.alphabet} # create amino dictionary

  def aaCount (self): # amino acid composition
    ''' This method will return a single integer count of valid amino acid 
    characters found. Do not assume that this is the length of the input 
    string, since you might have spaces or invalid characters that are 
    required to be ignored. '''
    
    return len(self.protein)

  def pI (self, p=2): # calculate the theoretical isolelectric point
    ''' Calculate the theoretical isolelectric point; can be estimated by 
    finding the particular pH that yields a neutral net Charge (close to 0). 
    By default accurate to 2 decimal places. '''    

    a, b = 0, 14 # pH level ranges

    # bisection algorithm. Converges in 'iters' number of iterations
    iters = math.ceil(math.log((b - a)/(10**(-p)))/math.log(2)) # min iterations

    for i in range(math.ceil(math.log((b - a)/(10**(-p)))/math.log(2))) : 
      c = (a + b)/2 # mid point

      sc = math.copysign(1, self._charge_(c)) # sign of f(c)
      sa = math.copysign(1, self._charge_(a)) # sign of f(a)

      if sc == sa : a = c 
      else : b = c

    return a if abs(a) < abs(b) else b

  def aaComposition (self) : # return single letter Amino acid code and associated counts
    ''' This method returns a dictionary keyed by single letter Amino 
    acid code, and associated values that are the counts of those 
    amino acids in the sequence. Includes all 20 amino acids. 
    Proper amino acids that are not represented in the sequence will have 
    a value of zero. '''

    return self.proteinCount

  def _charge_ (self, pH): # calculates the net charge on the protein at pH
    ''' Calculate the net charge on the protein at a specific pH. 
    This method is used by the pI method. I have marked it with the 
    single _ notation to advise future users of the class that this is not 
    part of the defined interface and it just might change. '''      

    # define a pKa,pH -> charge function to save space
    def pKaToChargePos(pKa) : return 10**pKa / (10**pKa + 10**pH) 
    def pKaToChargeNeg(pKa) : return 10**pH / (10**pKa + 10**pH) 

    # Summation of postivie charges: Arg, Lys, His, Nterm
    # Summation of negative charges: Asp, Glu, Cys, Tyr, Cterm
    posCharge = pKaToChargePos(self.aaNterm) # start with Nterm
    negCharge = pKaToChargeNeg(self.aaCterm) # start with Cterm
    for aa in self.aaComposition() : # go through count and total pos and neg charges
      if aa in self.aa2chargePos : 
        posCharge += self.aaComposition()[aa]*pKaToChargePos(self.aa2chargePos[aa])
      elif aa in self.aa2chargeNeg :
        negCharge += self.aaComposition()[aa]*pKaToChargeNeg(self.aa2chargeNeg[aa])

    # subtract and return result
    return posCharge - negCharge

  def molarExtinction (self, Cystine=True): # return how much light a protein absorbs at 280nm
    ''' Calculate extinction coefficient, which indicates how much light a protein absorbs 
    at a 280nm. It is useful to have an estimation of this 
    coefficient for measuring a protein with a spectrophotometer at a 
    wavelength of 280nm. It has been shown by Gill and von Hippel that it is 
    possible to estimate the molar extinction coefficient of a protein from 
    knowledge of its amino acid composition alone. From the molar extinction 
    coefficient of tyrosine, tryptophan and cystine at a given wavelength, 
    the extinction coefficient of the native protein in water can be computed 
    using the following equation. 
    Note that we will assume for this exercise that all Cysteine residues are 
    represented as Cystine. Under reducing conditions, Cystine does not form 
    however and Cysteine residues do not contribute to absorbance at 280nm. '''
    
    #compute the extinction using the amino acid absorbtion and count
    extinction = sum([c*self.aa2abs280[aa] for (aa,c) in self.aaComposition().items() if aa in self.aa2abs280])
    if not Cystine : extinction -= self.aaComposition()['C']*self.aa2abs280['C']

    return extinction

  def massExtinction (self, Cystine=True) : # return mass extinction coefficient
    ''' Calculates the Mass extinction coefficient from the Molar 
    Extinction coefficient by dividing by the molecularWeight of the 
    corresponding protein. 
    we are assuming that all Cysteine residues are present as Cystine. 
    Provide an optional parameter to both molarExtinction() and 
    massExtinction() to calculate these values assuming reducing conditions. 
    Use an optional parameter with a default of True (Cystine=True) to 
    allow evaluation of both molar and mass extinction under both oxidizing 
    (default) and reducing conditions. 
    oxidizing -> Cystine=True -> Cysteine forms Cystine
    reducing -> Cystine=False -> Cysteine stays Cysteine '''

    myMW = self.molecularWeight()

    return self.molarExtinction(Cystine=Cystine) / myMW if myMW else 0.0

  def molecularWeight (self): # calculate the molecular weight
    ''' This method calculates the molecular weight (MW) of the protein 
    sequence. If we have the composition of our protein, this is done by 
    summing the weights of the individual Amino acids and excluding the waters 
    that are released with peptide bond formation. '''

    totalAminoWeight = -((self.aaCount() - 1)*self.mwH2O)
    for aa in self.aaComposition() : 
      totalAminoWeight += self.aaComposition()[aa]*self.aa2mw[aa]

    return totalAminoWeight
  



########################################################################
# CommandLine
########################################################################
class CommandLine() :
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
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                            epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                            add_help = True, #default is True 
                                            prefix_chars = '-', 
                                            usage = '%(prog)s [options] -option1[default] <input >output'
                                            )
        self.parser.add_argument('-lG', '--longestGene', type=int, nargs='?', default=None, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, nargs='?', default=None, help='minimum Gene length')
        self.parser.add_argument('-g', '--greedy', action = 'store_true', help='greedy genes include multiple starts')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        self.parser.add_argument('-f', '--file', type=str, required=True)
        self.parser.add_argument('-oF', '--outFile', type=str, default=None)
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)




########################################################################
# OrfFinder
########################################################################
class OrfFinder() :
    ''' 
    OrfFiner takes in an input sequence and store it 

    usage:

    seq = 'AaatXCGTTATGggggXgTAA'
    geneFrames = ORF(seq).parseGenes() 
    '''

    def __init__(self, seq=None, args=None) :   
      ''' take in string and CommandLine arguments and clean and store them '''
      
      if seq is None : self.sequence = ''
      else : self.sequence = ''.join([n for n in seq.upper() if n in 'ATCG'])

      if args is None : self.args = CommandLine().parser.parse_args()
      else : self.args = args
    
    def seqLength(self) :
      ''' return the length of the sequence '''
      return len(self.sequence)

    def parseGenes(self) :
      ''' 
      parse the gene into snips of indices and base pair length
        return dict : {+1:[f1,f2,f3...,fn], +2:[f1,f2,f3...,fm], +3:[f1,f2,f3...,fnn]} 
      '''
      genesFrame = {}
      seq = self.sequence
      for d in [1,-1] : # check forward and reverse complement
          for f in range(3) : # check each frame offset
              codons = [seq[i-3:i] for i in range(f+3,self.seqLength()+1,3)] # break seq into codons by frame, then count the start and stop codons
              starts, stops = sum([codons.count(s) for s in self.args.start]), sum([codons.count(s) for s in self.args.stop])

              if starts and stops : # check if the sequence has starts and stops
                  startIndex, stopIndex = 0, min([codons.index(s) for s in self.args.stop if s in codons])

                  for sn in range(stops) : # check each of the stop codons
                      # skip this stop if there are no start codons in it
                      if sum([codons[startIndex:stopIndex].count(s) for s in self.args.start]) > 0 : 
                          if self.args.greedy : # greedy implies multiple start codons per frame
                              startIndex += min([codons[startIndex:stopIndex].index(s) for s in self.args.start if s in codons[startIndex:stopIndex]])
                          else : # non-greedy only has ony start codon in the frame
                              startIndex += (stopIndex-startIndex-1)-min([codons[startIndex:stopIndex][::-1].index(s) for s in self.args.start if s in codons[startIndex:stopIndex]])

                          frameKey = f"{'-'if d==-1 else'+'}{f+1}"
                          geneLength = stopIndex + 1 - startIndex

                          # check to make sure gene length is between defined values
                          if self.args.minGene == None or 3*geneLength >= self.args.minGene :
                            if self.args.longestGene == None or 3*geneLength <= self.args.longestGene :
                              sl, ss, sx = self.seqLength(), startIndex, stopIndex

                              if frameKey not in genesFrame : genesFrame[frameKey] = [] # add each snip

                              if d > 0 : genesFrame[frameKey].append( (3*ss+f+1, 3*(sx+1)+f, 3*geneLength) )
                              else : genesFrame[frameKey].append( (sl-3*(sx+1)+f+1, sl-3*ss+f, 3*geneLength) )

                      if sn + 1 < stops : # no need to check for stops after stop count
                        startIndex = stopIndex + 1
                        stopIndex += 1 + min([codons[startIndex:].index(s) for s in self.args.stop if s in codons[startIndex:]])

          # seq is reversed and changed to lower case compliment, then raised
          seq = seq[::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper() 

      return genesFrame
  