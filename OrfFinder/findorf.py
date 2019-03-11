# /Users/beagannguy/desktop/BME_160
# Name: Beagan Nguy(bnguy)
# Group Members: none

import sequenceAnalysis
import argparse
import sys 
########################################################################
# CommandLine
########################################################################
    
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
'''
    Main SeudoCode

    We want to find all the ORF that is at least 100 nucleotide long given a sequence. Three functions will be used. Function OrfFinder
    will be created to find all the ORFs given a sequence (positive or negative strand). In the function
    a nest for loop will be used to iterate through each frame and then the codon. The first if
    statement, will be used to find the start dangling ORF if it exist. If not then move on to find
    all the normal ORFs. To find a normal ORF, the program will scan through the sequence until it reaches a start codon.
    From there, it will continue to scan until it finds a stop codon that is 100 nucleotide away from the the start codon.
    Record the position of the start Index, Stop Index, Length, and Frame by calling the second function saveORF which will record the aspects
    of the ORF into a list. The startIndex list will then be be cleared as it is already recorded. Continue this on
    throughout the whole sequence. The last function, findRevOrfs, will be used used to find all the ORFs on the - or complementary
    strand. This function does this by reversing the sequence, then finding the complementary nucleotide throughout the whole reverse
    sequence. Once finish, it will run the reverse sequence into the findOrfs function to find all the ORFs within
    the reverse strand. Lastly, there is a if statement, indicating the end of Orfs search within that strand. If it find a 
    dangling stop ORF it will record it.

'''
    
class CommandLine():

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
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

class OrfFinder():
    def __init__(self, seq):
        self.seq = seq
        self.orfList = [[], [], []] # List for each frame
        self.startCodon = ('ATG')
        self.stopCodon = ('TAG', 'TAA', 'TGA')
        self.startIndexes = []
        self.stopIndexes = []
        self.firstadd = False 
        self.min = 100

    '''
    This function will be used to find all the ORF in the sequence with an minimum of 100 nucleotides long.
    To do this, a nested forloop was implemented where the first for loop iterates through each frame and
    the inner for loop iterates through every three nuc (codons).
    '''
    def findOrfs(self):
        ''' Iterates Through Each Frame (0, 1, 2) '''
        for frame in range(3):
            self.stopIndexes.clear()

   
            ''' Iterates through every third nucleotide '''
            for pos in range(frame, len(self.seq), 3):
                
                ''' codon = the pos -> pos + 2 '''
                codon = self.seq[pos: pos + 3]

                ''' 
                Determine what will be added first to the orf list.
                If the first stop codon in the sequence is 100 nucleotide 
                minimum from the start of the sequence. 
                This will be the dangling start ORF. 
                *** pos + 3 b/c we want to include stop codon
                '''
                if self.firstadd == False:

                    if codon in self.startCodon:
                        self.startIndexes.append(pos) 
                    if codon in self.stopCodon:
                        self.stopIndexes.append(pos)

                        '''
                        If there is nothing in the list at the frame and there is only
                        on stop codon. Determine if it is at least 100 nuc long and if
                        it is, this will be the start dangling ORF
                        '''

                        if not self.orfList[frame] and len(self.stopIndexes) == 1:
                        
                            length = (pos+3) 
                    
                            ''' Checks if the ORF is at least 100 nuc long '''
                            if length >= self.min:
                                self.saveORF(0, pos + 3, length, frame)
                                self.startIndexes.clear()
                                self.firstadd == True
                            else:
                                pass

                        '''
                        If the first stop codon is before the first 100 nucleotide then no 
                        start dangling ORF exist, proceed to find the normal ORFs
                        '''
                        if not self.orfList[frame] and len(self.stopIndexes) > 1:
                            if self.startIndexes:
                                length = (pos + 3) - self.startIndexes[0]

                                if length >= self.min:
                                    self.saveORF(self.startIndexes[0], pos + 3, length, frame)
                                self.startIndexes.clear()
                                self.firstadd == True

                ''' 
                Searches for a start codon after one ORF is found.
                Adds to the startIndexes list
                '''
                if codon in self.startCodon:
                    if self.orfList[frame]:
                        self.startIndexes.append(pos)  

                '''
                Search through frame 1, 2
                '''
                if codon in self.startCodon:
                    if frame != 0:                   
                        self.startIndexes.append(pos)  

                '''
                If it then encounters a stop codon, while having encountered
                a start codon beforehand. Check to see if the length is at least 100
                and add to ORF list
                '''
                if codon in self.stopCodon: 
                    if len(self.startIndexes) != 0:
                        length = (pos + 3) - self.startIndexes[0]
             
                        ''' Checks if the ORF is at least 100 nuc long '''
                        if length >= self.min:
                            self.saveORF(self.startIndexes[0], pos + 3, length, frame)
                            self.startIndexes.clear()
                    
                        else:
                            self.startIndexes.clear()

                '''
                When we reached the end of the sequence.
                If there is a start codon, but no stop, record dangling ORF
                '''
                if pos == len(self.seq) - 4 and len(self.startIndexes) != 0:
                    length = (len(self.seq) - 1) - self.startIndexes[0]
                        
                    '''Check to see if the ORF is at least 100 nuc long'''
                    if length >= self.min:
                        self.saveORF(self.startIndexes[0], (len(self.seq) - 1), length, frame)
        return self.orfList

    '''
    Finds the complementary sequence and runs it through findORFs to find all the ORFs
    ont the negative or complementary strand
    '''
    def findRevOrfs(self):
        '''
        Clear Everything found on the + strand
        '''
        self.orfList = [[],[],[]]
        self.startIndexes.clear()
        self.stopIndexes.clear()

        revSeq = self.seq[::-1]

        '''
        Finds the reverse complementary of the + strand
        '''
        revSeq = revSeq.replace('A', 'X') # A->T
        revSeq = revSeq.replace('G', 'Y') # G->C
        revSeq = revSeq.replace('C', 'G') # C->G
        revSeq = revSeq.replace('T', 'A') # T->A
        revSeq = revSeq.replace('X', 'T')
        revSeq = revSeq.replace('Y', 'C')

        self.seq = revSeq

        return self.findOrfs()

    '''
    Function used to store the ORFS found into a list
    '''
    def saveORF(self, start, stop, length, frame):
        orf = (frame + 1, start + 1, stop, length) # Frame & start + 1 b/c index starts at 0
        self.orfList[frame].append(orf)

def output(inFile, outFile):

    fastaReader = sequenceAnalysis.FastAreader(inFile)
    orfsRecord = []

    f = open(outFile, 'w')
    sys.stdout = f

    ''' Iterates through the input file using the fastaReader'''

    for header, sequence in fastaReader.readFasta():
        print(header)

        '''Find ORFs in the + and - strands of the sequence'''
        
        posList = OrfFinder(sequence).findOrfs()
        negList = OrfFinder(sequence).findRevOrfs()

        '''Iterates through the Orfs found on the positive strands and then iterates through
        each record in the Orf (frame, length, etc.), coverting the list into a tuple and 
        adding it to an overall ORFs list '''
        
        for list in posList:
            for tuple in list:
                frame = tuple[0] 
                startPos = tuple[1] 
                endPos = tuple[2]
                length = tuple[3]
                orfsRecord.append((frame, startPos, endPos, length))

        '''Similar to function above, except performs the coversion on the ORFs found 
        in the negative strand
        '''
        for list in negList:
            for tuple in list:
                frame = tuple[0] 
                startPos = tuple[1] 
                endPos = tuple[2]
                length = tuple[3]
                orfsRecord.append((-frame, len(sequence)-endPos+1, len(sequence)-startPos+1, length))

        '''Sort by decreasing gene and then start position furthest to the left'''
        orfsRecord.sort(key=lambda tup: (tup[3], tup[1]), reverse=True)

        '''
        Prints out the Orfs (Given within the assignment instructions by iterating
        through the list of Orfs 
        '''
        for orf in orfsRecord:
            print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(orf[0], orf[1], orf[2], orf[3]))

    f.close()
########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   

def main(inCL=None):
    '''
    Find some genes.  
    '''
    
    if inCL is None:
        myCommandLine = CommandLine(['tass2.fa',
                                    'tass2ORFdata-ATG-100.txt',
                                    '--longestGene'])

        output(myCommandLine.args.inFile, myCommandLine.args.outFile)
    else :
        myCommandLine = CommandLine(inCL)
    
###### replace the code between comments.
    #print (myCommandLine.args)
        # myCommandLine.args.inFile has the input file name
        # myCommandLine.args.outFile has the output file name
        # myCommandLine.args.longestGene is True if only the longest Gene is desired
        # myCommandLine.args.start is a list of start codons
        # myCommandLine.args.minGene is the minimum Gene length to include
        
#######
    
if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN

