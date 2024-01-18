#!/usr/bin/env python3
# Name: Madison Price (maeprice)
# Group Members: Kiana Imani (kimani), Queenie Li (qli86), and Ellizabeth Beer (ebeer)

""" Read a protein sequence from user input and calculate its physical and chemical properties,
    similar to what ProtParam software does.

    Input: Protein sequence
    Outputs: Number of amino acids and total molecular weight.
             Molar extinction coefficient and Mass extinction coefficient
             Theoretical isoelectric point (pI)
             Amino acid compositon
"""


class ProteinParam:
    """Build a ProteinParam object with all of the attributes given at instantiation time."""
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

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


# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys


def main():
    inString = input('protein sequence?')
    while inString:
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print("Number of Amino Acids: {aaNum}".format(aaNum=myAAnumber))
        print("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print("Theoretical pI: {:.2f}".format(myParamMaker.pIBinary()))
        print("Amino acid composition:")

        if myAAnumber == 0: myAAnumber = 1  # handles the case where no AA are present

        for aa, n in sorted(myParamMaker.aaComposition().items(),
                            key=lambda item: item[0]):
            print("\t{} = {:.2%}".format(aa, n / myAAnumber))

        inString = input('protein sequence?')


if __name__ == "__main__":
    main()