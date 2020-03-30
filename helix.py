
'''
H E L I X
'''

import time
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
from urllib.request import urlopen
from Bio.Seq import Seq
from Bio import SeqIO

# Functions

# param: a protein UnitProt ID
# return: the protein's DNA sequence
def getUniProtFasta(protein_id):
	url = 'http://www.uniprot.org/uniprot/' + protein_id + '.fasta'
	with urlopen(url) as data:
		seq = str(data.read())
		seq = seq[seq.index('\\'):] #get the protein seq without fasta heading
		seq = seq.replace('\\n','') # remove all \n
		seq = seq.replace("'",'') # remove a " ' " at the end
		return seq


# param: a motif/pattern (motif) & a DNA sequence (seq)
# return: all locations of the motif in the seq
def findMotifLocations(motif,seq):
	motifLen = len(motif)
	motifLocations = []

	for i in range(0,len(seq)):
		subseq = seq[i:(i+motifLen)]
		if subseq == motif:
			motifLocations.append(str(i+1))
	return motifLocations


# param: a dna fragment and a list of sequences
# return: a boolean indicating whether the fragment is contained in all sequences
def isGlobalSubfragment(fragment,seqs):
	for seq in seqs:
		if fragment not in seq:
			return False
	return True


# param: a list of dna sequences
# return: all the longest motifs found in the sequences
def findLongestMotif(dnaSeqs):

	# find and extract the smallest dna seq
	currentMinIndex = 0
	currentMinLen = len(dnaSeqs[0])
	for index in range(len(dnaSeqs)):
		seqLen = len(dnaSeqs[index])
		if seqLen < currentMinLen:
			currentMinLen = seqLen
			currentMinIndex = index
	smallestSeq = dnaSeqs.pop(currentMinIndex)

	allLongest = []
	longest = ''
	totalIterations = len(smallestSeq)

	# Iterate over all posible combinations and find the longest shared motif
	for startIndex in range(len(smallestSeq)+1):
		endIndex = len(smallestSeq)
		currentSubSeq = smallestSeq[startIndex:endIndex]

		print(str(startIndex/totalIterations*100)[:4],'%')

		while len(currentSubSeq) >= len(longest):
			if isGlobalSubfragment(currentSubSeq,dnaSeqs):
				if len(currentSubSeq) == len(longest):
					allLongest.append(currentSubSeq)
				else:
					longest = currentSubSeq
					allLongest = []
					allLongest.append(currentSubSeq)
			endIndex -= 1
			currentSubSeq = smallestSeq[startIndex:endIndex]
	return allLongest

'''
proteinIds = ['R1AB_CVHSA','R1AB_BCHK3','R1AB_CVMJH']
allSeq = []

for id in proteinIds:
	allSeq.append(getUniProtFasta(id))
'''

allSequenceRecords = SeqIO.parse('testSeqs.txt','fasta')
allSeq = []
for seq in allSequenceRecords:
	allSeq.append(str(seq.seq))


# Run proceses
startTime = time.time()
longestMotifsFound = findLongestMotif(allSeq)
deltaTime = time.time() - startTime

# Final results
print('_______________________________\n            Report\n_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _\n')
print('- Search time:',str(deltaTime)[:5],'seconds')
print('- Longest motif(s) found:')
for motif in longestMotifsFound:
	print('\t',motif,'\n')
print('_______________________________')


print(findMotifLocations(longestMotifsFound[0],allSeq[0]))






