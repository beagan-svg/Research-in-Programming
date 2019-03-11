# /Users/beagannguy/desktop/BME_160
# Name: Beagan Nguy(bnguy)
# Group Members: none
"""
This program contains three classes NucParams, ProteinParams, and FastAreader. Within
each class contains functions and method to analyze the inputed sequence and create a dictionary
storing the counts of aa, nucleotides, and codons. The dictionary will then be used in genomeAnalyzer.py
to output the GC ccontent, sequence length, and codon usage.
"""
import sys

class NucParams:
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
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):  

        # Initializes Dictionary giving it default values of 0 and keys as specifed in the
        # assignment
        # Dictionary to be used to store the counts of aa, nucleotides, and codons
        self.aaComp = {aa: 0 for aa in NucParams.rnaCodonTable.values()}
        self.nucComp = {'A': 0, 'C': 0, 'G': 0, "T": 0, 'U': 0, 'N': 0}
        self.codonComp = {codon: 0 for codon in NucParams.rnaCodonTable.keys()}

        self.addSequence(inString) 

    # Trascribe the raw seq to rna and translate to protein
    # Storing method of codons and protein is commented below
    # GUG -> V   
    def addSequence (self, inSeq):
        # Upper cases the Nucleotides
        rawSeq = str.upper(inSeq)

        # Accepted nuc: {A,C,G,T,N,U}
        accepted = ['A','C','G','T','N','U']

        # To be used to store the string of filtered rawSeq (no spaces or unwanted characters)
        seq = ''

        # Filters the rawSeq 
        for nuc in rawSeq:   
            if nuc in accepted:
                seq += nuc

        # Transcribe T with U
        rna = seq.replace("T", "U")

        # Removes invalid bases
        for aa in rna:
            if 'N' in aa:
                rna = rna.remove(aa)

        # To be used to store as a protein string
        protein = ''
        
        # List to be used to store all the codons  
        codonList = []

        # Adds codon to a list
        for i in range(0,len(rna)-2,3):
            triplet = rna[i:i+3]
            codonList.append(triplet)
            
            # Translate the codon to a aa and adds it to a protein string
            aa = self.rnaCodonTable[triplet]
            protein += aa

        # aaComposition
        """
        for aa in self.rnaCodonTable.values():
            self.aaComp[aa] = self.aaComp.get(aa) + protein.count(aa)
       
        self.aaComp['stop'] = protein.count('stop')
        """
        for codon in codonList:
            if codon in self.rnaCodonTable.keys():
                self.aaComp[NucParams.rnaCodonTable[codon]] += 1

        # nucComposition
        x = 'ACGTNU'
        for nuc in x:
            self.nucComp[nuc] = self.nucComp.get(nuc) + seq.count(nuc)
       
       # codonComposition
        for codon in self.rnaCodonTable.keys():
            self.codonComp[codon] = self.codonComp.get(codon) + codonList.count(codon)


    # Reads the input sequence and counts each aa 
    # {aa: 5}
    def aaComposition(self):    
        return self.aaComp
    
    # Returns counts of nucleotides {ACGTNU} in a dictionary
    # {N: 2}
    def nucComposition(self): 
        return self.nucComp

    # Return counts of the codon in a dictionary
    # {GUU: 5}
    def codonComposition(self):  
        return self.codonComp

    def nucCount(self):
        return sum(self.nucComp.values())

class ProteinParam :
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
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        # Upper case the protein sequence
        upperProtein = str.upper(protein)
        
        # Creates an empty string
        self.str_Protein_Seq = ""
        
        # Filter letters that are not within the dictionary
        # by iterating through the sequence and if it matches with
        # one of the keys in the dictinary, it will be added to a string
        
        for aa in upperProtein:
            if aa in self.aa2mw.keys():
                self.str_Protein_Seq += aa
                
        # print(self.filter_Protein_Seq)
        #['V', 'L', 'S', 'P', 'A', 'D', 'K', 'T', 'N', 'V', 'K', 'A', 'A', 'W']

        # print(self.str_Protein_Seq)
        #VLSPADKTNVKAAW
        
        # Length of Amino Acids
        self.length = 0
        
        # Dictionary where key = aminoAcid and value = count
        self.aaDict = {}
        
        # Value of pH
        self.pH = 0
        
    # Returns count of Amino Acids
    def aaCount (self):
        self.length = len(self.str_Protein_Seq)
        return self.length

    # Returns the theoretical isolelectric point
    def pI (self):
        # Assign values to float and int
        precise_pH = 0.0
        best_Charge = 1000
        
        # Iterates through all the pH
        while (self.pH <= 14.01):
            temp_Charge = self._charge_()
            #print(tempCharge)
            
            # Finds the theoretical isoelectric point closest to 0
            # abs takes the absolute value
            if abs(temp_Charge) <= abs(best_Charge):
                precise_pH = self.pH
                best_Charge = temp_Charge
            self.pH += 0.01
       
        return precise_pH
    
    # Returns the counts of each Amino Acid in a dictionary
    def aaComposition (self) :
        for aminoAcid in self.aa2mw.keys():
            self.aaDict[aminoAcid] = self.str_Protein_Seq.count(aminoAcid)
        
        return self.aaDict
        """
        for key,val in self.aaDict.items():
            print (key, "=>", val)
        """
    # Returns the net charge of the protein at a particular pH    
    def _charge_ (self):
        # Assign values to float
        posCharge = 0.0
        negCharge = 0.0
        numerator = 0.0
        denominator = 0.0
        
        #VLSPADKTNVKAAW
        # Iterates through amino acid sequence and calculates the net sum of all 
        # those charges
        for aminoAcid in self.str_Protein_Seq:
            
            # Calculates the aa with positive charges
            if aminoAcid in self.aa2chargePos.keys():
                numerator = 10 ** self.aa2chargePos[aminoAcid]
                denominator = 10 ** self.aa2chargePos[aminoAcid] + (10 ** self.pH)
                posCharge += (numerator/denominator)
            
            # Calculates the aa with negative charges
            elif aminoAcid in self.aa2chargeNeg.keys():
                numerator = 10 ** self.pH
                denominator = 10 ** self.aa2chargeNeg[aminoAcid] + (10 ** self.pH)
                negCharge += (numerator/denominator)
        
        # Adds the aa(Nterm) to the positive charge
        posCharge += (10 ** self.aaNterm) / (10 ** self.aaNterm + 10 ** self.pH)
        
        # Adds the aa(Cterm) to the negative charge
        negCharge += (10 ** self.pH) / (10 ** self.aaCterm + 10 ** self.pH)
        
        # Calculates the net charge
        netCharge = posCharge - negCharge
        
        return netCharge
    
    # Return the amount of light a protein absorbs at a specific wavelength
    def molarExtinction (self):
        # Calculates the molar coefficient of typrosine, tryptophans, cysteines
        typrosine = self.str_Protein_Seq.count('Y') * self.aa2abs280.get('Y')
        tryptophans = self.str_Protein_Seq.count('W') * self.aa2abs280.get('W')
        cysteines = self.str_Protein_Seq.count('C') * self.aa2abs280.get('C')
        
        # Calculates the molarExtinction by summing the molar coefficient above
        molarExtinction = typrosine + tryptophans + cysteines
        return molarExtinction
        pass
    
    # Returns the Mass extinction coefficient 
    def massExtinction (self):
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    # Returns the total molecular weight of the Amino Acids excluding water
    def molecularWeight (self):
        # Assign values to int
        total_Weight = 0
        aa_Weight = 0
        waterLoss = 0
        
        # Weight of all AA
        for aminoAcid in self.str_Protein_Seq:
            aa_Weight += self.aa2mw.get(aminoAcid)
        
        # Weight of water loss 
        waterLoss = (self.length - 1) * self.mwH2O
       
        # Return net weight
        total_Weight = aa_Weight - waterLoss
        return total_Weight  

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an  entire FastA record and return the sequence header/sequence'''
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
