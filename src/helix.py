
'''
H E L I X
'''

import time
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
from urllib.request import urlopen
from Bio.Seq import Seq
from Bio import SeqIO,motifs
import tkinter as tk
from tkinter import ttk
from PIL import Image


'''
------------------------------------------
				Functions
------------------------------------------
'''


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

		while len(currentSubSeq) >= len(longest):
			if isGlobalSubfragment(currentSubSeq,dnaSeqs):
				if len(currentSubSeq) == len(longest) and currentSubSeq != longest:
					allLongest.append(currentSubSeq)
				else:
					longest = currentSubSeq
					allLongest = []
					allLongest.append(currentSubSeq)
			endIndex -= 1
			currentSubSeq = smallestSeq[startIndex:endIndex]
	return allLongest


# param: a list of dna aminoacid sequences
# return: a list of dna base pair sequences
def aminoAcidToBase(seqs):
	translateDic = {'F':'TTT','L':'CTT','I':'ATT','V':'GTT','F':'TTC','L':'CTC','I':'ATC','V':'GTC','L':'TTA','L':'CTA','I':'ATA','V':'GTA','L':'TTG','L':'CTG','M':'ATG','V':'GTG','S':'TCT','P':'CCT','T':'ACT','A':'GCT','S':'TCC','P':'CCC','T':'ACC','A':'GCC','S':'TCA','P':'CCA','T':'ACA','A':'GCA','S':'TCG','P':'CCG','T':'ACG','A':'GCG','Y':'TAT','H':'CAT','N':'AAT','D':'GAT','Y':'TAC','H':'CAC','N':'AAC','D':'GAC','Stop':'TAA','Q':'CAA','K':'AAA','E':'GAA','Stop':'TAG','Q':'CAG','K':'AAG','E':'GAG','C':'TGT','R':'CGT','S':'AGT','G':'GGT','C':'TGC','R':'CGC','S':'AGC','G':'GGC','Stop':'TGA','R':'CGA','R':'AGA','G':'GGA','W':'TGG','R':'CGG','R':'AGG','G':'GGG'}
	translated = []

	for seq in seqs:
		translatedSeq = ""
		for aa in seq:
			translatedSeq += translateDic[aa]
		translated.append(translatedSeq)

	return translated


# param: a list of dna sequences
# return: a list of dna sequences all with the same lenght
def normalizeSeqs(seqs):
	minLength = min(len(seq) for seq in seqs)
	normalized = []

	for seq in seqs:
		normalized.append(Seq(seq[:minLength]))

	return normalized


def showConfirmation(text,confirmed):

	# Set up window
	confWindow = tk.Tk()
	confWindow.title("Helix")

	canvas2 = tk.Canvas(confWindow,height=125,width=400,bg=appBgColor)
	canvas2.pack()
	frame2 = tk.Frame(confWindow,bg=appBgColor)
	frame2.place(relx=0,rely=0,relwidth=1,relheight=1)

	spacer1 = tk.Label(frame2,text = ' ',font='Arial 30',bg=appBgColor)
	spacer1.pack(side='top')

	confLabel = tk.Label(frame2,text=text,bg=appBgColor,font='Arial 15 bold')
	if confirmed:
		confLabel['fg'] = '#27ae60'
	else:
		confLabel['fg'] = '#e74c3c'
	confLabel.pack(side='top')

	spacer2 = tk.Label(frame2,text = ' ',font='Arial 15',bg=appBgColor)
	spacer2.pack(side='top')

	closeButton = tk.Button(frame2,text='Ok',fg='#575757',width=10,command=lambda: confWindow.destroy())
	closeButton.pack(side='top')


# Home Screen
def mainApp(reload):
	root.title('Helix')

	global canvas, frame

	if reload:
		canvas.destroy()
		frame.destroy()

	canvas = tk.Canvas(root,height=400,width=700,bg=appBgColor)
	canvas.pack()
	frame = tk.Frame(root,bg=appBgColor)
	frame.place(relx=0,rely=0,relwidth=1,relheight=1)

	spacerTop = tk.Label(frame,text='',font='Arial 40',bg=appBgColor)
	spacerTop.pack(side='top')
	helixLogoLabel = tk.Label(frame,image=helixLogo,pady=0, padx=0, borderwidth=0, highlightthickness=0)
	helixLogoLabel.pack(side='top')

	uniProtLogoLabel = tk.Label(frame,image=unitProtLogo,pady=0, padx=0, borderwidth=0, highlightthickness=0)
	uniProtLogoLabel.place(relx=0.25,rely=0.8)
	ncbiLogoLabel = tk.Label(frame,image=ncbiLogo,pady=0, padx=0, borderwidth=0, highlightthickness=0)
	ncbiLogoLabel.place(relx=0.45,rely=0.8)
	bioPythonLogoLabel = tk.Label(frame,image=bioPythonLogo,pady=0, padx=0, borderwidth=0, highlightthickness=0)
	bioPythonLogoLabel.place(relx=0.65,rely=0.8)

	options = ['Find longest shared motif','Get consensus','Get Weblogo','Find motif locations','Get sequence']
	selectedOption = tk.StringVar(frame)
	selectedOption.set(options[0])
	dropDownMenu = tk.OptionMenu(frame,selectedOption,options[0],options[1],options[2],options[3],options[4])
	dropDownMenu.config(bg=appBgColor)
	dropDownMenu.place(relx=0.1,rely=0.5,relwidth=0.5)

	goButton = tk.Button(frame,text='Go',fg=lightLetterColor,command=lambda: redirect(selectedOption.get()))
	goButton.config(bg=appBgColor)
	goButton.place(relx=0.68,rely=0.5,relwidth=0.15)


def redirect(selectedOption):
	if selectedOption == 'Find longest shared motif':
		getLSMPage()
	if selectedOption == 'Get consensus':
		getConsensusPage()
	if selectedOption == 'Get Weblogo':
		getWebLogoPage()
	if selectedOption == 'Find motif locations':
		getMLocPage()
	if selectedOption == 'Get sequence':
		getSeqPage()


def getLSMPage():
	root.title('Longest Shared Motif')

	global canvas, frame
	canvas.destroy()
	frame.destroy()

	canvas = tk.Canvas(root,height=400,width=700,bg=appBgColor)
	canvas.pack()
	frame = tk.Frame(root,bg=appBgColor)
	frame.place(relx=0,rely=0,relwidth=1,relheight=1)

	spacer1 = tk.Label(frame,text='',font='Arial 40',bg=appBgColor)
	spacer1.pack(side='top')
	instruction = tk.Label(frame,text = 'Enter protein IDs or filename',fg='white',bg=appBgColor)
	instruction.pack(side='top')
	entry = tk.Entry(frame,fg='white',bg='#141414',width=50)
	entry.pack(side='top')

	spacer2 = tk.Label(frame,text='',font='Arial 20',bg=appBgColor)
	spacer2.pack(side='top')
	searchButton = tk.Button(frame,text='  Find  ',fg=lightLetterColor,command=lambda: search(entry.get()))
	searchButton.config(bg=appBgColor)
	searchButton.pack(side='top')

	spacer3 = tk.Label(frame,text='',font='Arial 25',bg=appBgColor)
	spacer3.pack(side='top')

	global resultListbox
	resultListbox = tk.Listbox(frame,background='#141414',fg='white',selectbackground='#2f6492',width=50,height=8)
	resultListbox.pack(side='top')

	returnButton = tk.Button(frame,text='  return  ',fg=lightLetterColor,command=lambda: mainApp(reload=True))
	returnButton.config(bg=appBgColor)
	returnButton.pack(side='bottom')

def search(entry):

	if entry != '':
		try:
			if entry[-4:] == '.txt':
				allSequenceRecords = SeqIO.parse(entry,'fasta')
				allSeq = []
				for seq in allSequenceRecords:
					allSeq.append(str(seq.seq))
			else:
				proteinIds = entry.split(',')
				allSeq = []
				for id in proteinIds:
					allSeq.append(getUniProtFasta(id))
			
			# Run process
			startTime = time.time()
			longestMotifsFound = findLongestMotif(allSeq)
			deltaTime = time.time() - startTime

			resultListbox.insert(0,'')
			for motif in longestMotifsFound:
				resultListbox.insert(0,motif)
			resultListbox.insert(0,'> ' + entry + ' T= ' + str(deltaTime)[:5] + 's')
		except:
			resultListbox.insert(0,'')
			resultListbox.insert(0,'> No result: bad input')


def getConsensusPage():
	root.title('Get Consensus')

	global canvas, frame
	canvas.destroy()
	frame.destroy()

	canvas = tk.Canvas(root,height=400,width=700,bg=appBgColor)
	canvas.pack()
	frame = tk.Frame(root,bg=appBgColor)
	frame.place(relx=0,rely=0,relwidth=1,relheight=1)

	spacer1 = tk.Label(frame,text='',font='Arial 40',bg=appBgColor)
	spacer1.pack(side='top')
	instruction = tk.Label(frame,text = 'Enter protein IDs or filename',fg='white',bg=appBgColor)
	instruction.pack(side='top')
	entry = tk.Entry(frame,fg='white',bg='#141414',width=50)
	entry.pack(side='top')

	spacer2 = tk.Label(frame,text='',font='Arial 20',bg=appBgColor)
	spacer2.pack(side='top')
	searchButton = tk.Button(frame,text='  Get Consensus  ',fg=lightLetterColor,command=lambda: getConsensus(entry.get()))
	searchButton.config(bg=appBgColor)
	searchButton.pack(side='top')

	spacer3 = tk.Label(frame,text='',font='Arial 25',bg=appBgColor)
	spacer3.pack(side='top')

	global resultListbox
	resultListbox = tk.Listbox(frame,background='#141414',fg='white',selectbackground='#2f6492',width=50,height=8)
	resultListbox.pack(side='top')

	returnButton = tk.Button(frame,text='  return  ',fg=lightLetterColor,command=lambda: mainApp(reload=True))
	returnButton.config(bg=appBgColor)
	returnButton.pack(side='bottom')

def getConsensus(entry):
	if entry != '':
		try:
			if entry[-4:] == '.txt':
				allSequenceRecords = SeqIO.parse(entry,'fasta')
				allSeq = []
				for seq in allSequenceRecords:
					allSeq.append(str(seq.seq))

				for aa in ['F', 'L', 'I', 'V', 'M', 'S', 'P', 'Y', 'H', 'N', 'D', 'Q', 'K', 'E', 'R', 'W']:
					if aa in allSeq[0]:
						allSeq = aminoAcidToBase(allSeq)
						break

				allSeq = normalizeSeqs(allSeq)
			else:
				proteinIds = entry.split(',')
				allSeq = []
				for id in proteinIds:
					allSeq.append(getUniProtFasta(id))
				allSeq = aminoAcidToBase(allSeq)
				allSeq = normalizeSeqs(allSeq)

			# Run process
			startTime = time.time()
			
			# Find consensus
			motifData = motifs.create(allSeq)
			consensus =  str(motifData.consensus)
			degenerate_consensus = str(motifData.degenerate_consensus)
			
			deltaTime = time.time() - startTime

			resultListbox.insert(0,'')
			resultListbox.insert(0,'Degenerate Consensus: ' + degenerate_consensus)
			resultListbox.insert(0,'Consensus: ' + consensus)
			resultListbox.insert(0,'> ' + entry + ' T= ' + str(deltaTime)[:5] + 's')
		except:
			resultListbox.insert(0,'')
			resultListbox.insert(0,'> No result: bad input')


def getWebLogoPage():
	root.title('Get Weblogo')

	global canvas, frame
	canvas.destroy()
	frame.destroy()

	canvas = tk.Canvas(root,height=400,width=700,bg=appBgColor)
	canvas.pack()
	frame = tk.Frame(root,bg=appBgColor)
	frame.place(relx=0,rely=0,relwidth=1,relheight=1)

	spacer1 = tk.Label(frame,text='',font='Arial 40',bg=appBgColor)
	spacer1.pack(side='top')
	instruction1 = tk.Label(frame,text = 'Enter protein IDs or filename',fg='white',bg=appBgColor)
	instruction1.pack(side='top')
	entry1 = tk.Entry(frame,fg='white',bg='#141414',width=50)
	entry1.pack(side='top')

	spacer2 = tk.Label(frame,text='',font='Arial 40',bg=appBgColor)
	spacer2.pack(side='top')
	instruction2 = tk.Label(frame,text = 'Save to',fg='white',bg=appBgColor)
	instruction2.pack(side='top')
	entry2 = tk.Entry(frame,fg='white',bg='#141414',width=50)
	entry2.pack(side='top')

	spacer3 = tk.Label(frame,text='',font='Arial 20',bg=appBgColor)
	spacer3.pack(side='top')
	getButton = tk.Button(frame,text='  Generate Weblogo  ',fg=lightLetterColor,command=lambda: getWebLogo(entry1.get(),entry2.get()))
	getButton.config(bg=appBgColor)
	getButton.pack(side='top')

	returnButton = tk.Button(frame,text='  return  ',fg=lightLetterColor,command=lambda: mainApp(reload=True))
	returnButton.config(bg=appBgColor)
	returnButton.pack(side='bottom')

def getWebLogo(entry,saveName):
	if entry != '':

		if saveName != '':
			if saveName[-4:] != '.png':
				saveName += '.png'
		else:
			if entry[-4:] == '.txt':
				saveName = entry[:-4] + '-weblogo.png'
			else:
				saveName = entry + '-weblogo.png'

		try:
			if entry[-4:] == '.txt':
				allSequenceRecords = SeqIO.parse(entry,'fasta')
				allSeq = []
				for seq in allSequenceRecords:
					allSeq.append(str(seq.seq))

				for aa in ['F', 'L', 'I', 'V', 'M', 'S', 'P', 'Y', 'H', 'N', 'D', 'Q', 'K', 'E', 'R', 'W']:
					if aa in allSeq[0]:
						allSeq = aminoAcidToBase(allSeq)
						break
						
				allSeq = normalizeSeqs(allSeq)
			else:
				proteinIds = entry.split(',')
				allSeq = []
				for id in proteinIds:
					allSeq.append(getUniProtFasta(id))
				allSeq = aminoAcidToBase(allSeq)
				allSeq = normalizeSeqs(allSeq)

			# Run process
			startTime = time.time()
			
			# Create motif data
			motifData = motifs.create(allSeq)
			motifData.weblogo('weblogos/'+saveName)

			deltaTime = time.time() - startTime

			showConfirmation('Weblogo saved! \n(' + saveName + ')',True)
			
			im = Image.open('weblogos/'+saveName)
			im.show()

		except:
			showConfirmation('Sorry, error creating Weblogo.',False)


def getMLocPage():
	root.title('Motif Locations')

	global canvas, frame
	canvas.destroy()
	frame.destroy()

	canvas = tk.Canvas(root,height=400,width=700,bg=appBgColor)
	canvas.pack()
	frame = tk.Frame(root,bg=appBgColor)
	frame.place(relx=0,rely=0,relwidth=1,relheight=1)

	spacer1 = tk.Label(frame,text='',font='Arial 30',bg=appBgColor)
	spacer1.pack(side='top')
	instruction1 = tk.Label(frame,text = 'Enter protein ID',fg='white',bg=appBgColor)
	instruction1.pack(side='top')
	entry1 = tk.Entry(frame,fg='white',bg='#141414',width=50)
	entry1.pack(side='top')

	spacer2 = tk.Label(frame,text='',font='Arial 20',bg=appBgColor)
	spacer2.pack(side='top')
	instruction2 = tk.Label(frame,text = 'Enter motif sequence',fg='white',bg=appBgColor)
	instruction2.pack(side='top')
	entry2 = tk.Entry(frame,fg='white',bg='#141414',width=50)
	entry2.pack(side='top')

	spacer3 = tk.Label(frame,text='',font='Arial 15',bg=appBgColor)
	spacer3.pack(side='top')
	getButton = tk.Button(frame,text='  Locate  ',fg=lightLetterColor,command=lambda: getMotifLocations(entry1.get(),entry2.get()))
	getButton.config(bg=appBgColor)
	getButton.pack(side='top')

	global resLabel
	spacer4 = tk.Label(frame,text='',font='Arial 20',bg=appBgColor)
	spacer4.pack(side='top')
	resLabel = tk.Label(frame,text = '',font='Arial 20 bold',fg='#3498db',bg=appBgColor)
	resLabel.pack(side='top')
	
	returnButton = tk.Button(frame,text='  return  ',fg=lightLetterColor,command=lambda: mainApp(reload=True))
	returnButton.config(bg=appBgColor)
	returnButton.pack(side='bottom')

def getMotifLocations(protein_id,motif):
	if protein_id != '' and motif != '':
		try:
			protein_seq = getUniProtFasta(protein_id)
			locs = findMotifLocations(motif,protein_seq)

			text = 'Location(s): '

			if len(locs) == 0:
				text = 'No locations found'
			else:
				for loc in locs:
					text += (str(loc) + '  ')
			resLabel['text'] = text

		except:
			resLabel['text'] = 'Error'


def getSeqPage():
	root.title('Get Sequence')

	global canvas, frame
	canvas.destroy()
	frame.destroy()

	canvas = tk.Canvas(root,height=400,width=700,bg=appBgColor)
	canvas.pack()
	frame = tk.Frame(root,bg=appBgColor)
	frame.place(relx=0,rely=0,relwidth=1,relheight=1)

	spacer1 = tk.Label(frame,text='',font='Arial 100',bg=appBgColor)
	spacer1.pack(side='top')
	instruction1 = tk.Label(frame,text = 'Enter protein ID (or IDs separated by ",")',fg='white',bg=appBgColor)
	instruction1.pack(side='top')
	entry1 = tk.Entry(frame,fg='white',bg='#141414',width=50)
	entry1.pack(side='top')

	spacer2 = tk.Label(frame,text='',font='Arial 50',bg=appBgColor)
	spacer2.pack(side='top')
	getButton = tk.Button(frame,text='  Retreive  ',fg=lightLetterColor,command=lambda: getSequence(entry1.get()))
	getButton.config(bg=appBgColor)
	getButton.pack(side='top')

	returnButton = tk.Button(frame,text='  return  ',fg=lightLetterColor,command=lambda: mainApp(reload=True))
	returnButton.config(bg=appBgColor)
	returnButton.pack(side='bottom')

def getSequence(entry):

	try:
		proteinIds = entry.split(',')
		allSeq = []
		for id in proteinIds:
			allSeq.append(getUniProtFasta(id))

		file = open("sequences.txt","w")

		for i in range(len(proteinIds)):
			file.write('> ' + proteinIds[i] + '\n')

			index = 0
			while index <= len(allSeq[i]):
				file.write(allSeq[i][index:index+60] + '\n')
				index += 60

		file.close()

		showConfirmation('Retreived sequences.\nSaved to "sequences.txt"',True)
	except:
		showConfirmation('Error, could not save file.',False)




'''
------------------------------------------
				Run App
------------------------------------------
'''

root = tk.Tk()
root.configure(background='white')

# preload assets
unitProtLogo = tk.PhotoImage(file='assets/logo-uniprot.png')
ncbiLogo = tk.PhotoImage(file='assets/logo-ncbi.png')
bioPythonLogo = tk.PhotoImage(file='assets/logo-biopython.png')
helixLogo = tk.PhotoImage(file='assets/logo-helix.png')
appBgColor = '#181818'
lightLetterColor = '#3d3d3d'

# run home screen
mainApp(reload=False)
root.mainloop()




