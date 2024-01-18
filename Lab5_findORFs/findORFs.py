#!/usr/bin/env python3
# Name: Mark Forbush 1990886
# Group Members: None

########################################################################
# Main
# Here is the main program
# take in the name of fastA file and optional parameters and output
# gene stippets to file and console

# example 1 : 'the specified output for lab5'
# python findORFs.py -f "tass2.fa" -g -mG 300

# example 2 : 'non greedy gene finder with no minimum'
# python findORFs.py -f "tass2.fa"
########################################################################
from sequenceAnalysis import FastAreader, OrfFinder, CommandLine

def main(inFile=None, options=None):
    ''' Main function for the program
    prints the gene fragments and their frames and stores it in a file
    '''
    thisCommandLine = CommandLine(options)
    print(f"Input params: {thisCommandLine.args}\nreading...")

    if inFile == None : inFile = thisCommandLine.args.file # get the in file from the command line optionally

    reader = FastAreader(inFile) # read the fastA file
    for head, seq in reader.readFasta() : 
        print(f"head: {head}")
        geneFrames = OrfFinder(seq, thisCommandLine.args).parseGenes()
        
        geneStips, geneBases = [], [] # parse gene frames into array of snips and gene length
        for frame,genes in geneFrames.items() :
            for gene in genes : 
                geneStips.append( f"{frame:>2} {gene[0]:>5}..{gene[1]:<5} {gene[2]:>5}" )
                geneBases.append(gene[2])
        
        readOut = '\n'.join(list(zip(*sorted(zip(geneStips,geneBases), key=lambda x:-x[1])))[0])

        outFile = f"{thisCommandLine.args.file.split('.')[0]}ORFdata-{','.join(thisCommandLine.args.start)}-{thisCommandLine.args.minGene}.txt"
        with open(outFile, 'w') as f : f.write(readOut)

        print(readOut) # print the readout

    
if __name__ == "__main__" : main()

