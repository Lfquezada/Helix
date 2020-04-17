
'''
H E L I X
'''

import time
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
from urllib.request import urlopen
from Bio.Seq import Seq
from Bio import SeqIO
import tkinter as tk
from tkinter import ttk


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
	#progressBar.config(mode='determinate',maximum=totalIterations)

	# Iterate over all posible combinations and find the longest shared motif
	for startIndex in range(len(smallestSeq)+1):
		endIndex = len(smallestSeq)
		currentSubSeq = smallestSeq[startIndex:endIndex]

		#print(str(startIndex/totalIterations*100)[:4],'%')
		#progressBar['value'] = startIndex
		#progressBar.update()

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
	uniProtLogoLabel.place(relx=0.2,rely=0.8)
	ncbiLogoLabel = tk.Label(frame,image=ncbiLogo,pady=0, padx=0, borderwidth=0, highlightthickness=0)
	ncbiLogoLabel.place(relx=0.35,rely=0.8)
	bioPythonLogoLabel = tk.Label(frame,image=bioPythonLogo,pady=0, padx=0, borderwidth=0, highlightthickness=0)
	bioPythonLogoLabel.place(relx=0.5,rely=0.8)
	pubmedLogoLabel = tk.Label(frame,image=pubmedLogo,pady=0, padx=0, borderwidth=0, highlightthickness=0)
	pubmedLogoLabel.place(relx=0.65,rely=0.8)

	options = ['Find longest shared motif','Find motif locations','Get sequence']
	selectedOption = tk.StringVar(frame)
	selectedOption.set(options[0])
	dropDownMenu = tk.OptionMenu(frame,selectedOption,options[0],options[1],options[2])
	dropDownMenu.config(bg=appBgColor)
	dropDownMenu.place(relx=0.1,rely=0.5,relwidth=0.5)

	goButton = tk.Button(frame,text='Go',fg=lightLetterColor,command=lambda: redirect(selectedOption.get()))
	goButton.config(bg=appBgColor)
	goButton.place(relx=0.68,rely=0.5,relwidth=0.15)


def redirect(selectedOption):
	if selectedOption == 'Find longest shared motif':
		getLSMPage()
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
	searchButton = tk.Button(frame,text='  Search  ',fg=lightLetterColor,command=lambda: search(entry.get()))
	searchButton.config(bg=appBgColor)
	searchButton.pack(side='top')

	spacer3 = tk.Label(frame,text='',font='Arial 25',bg=appBgColor)
	spacer3.pack(side='top')

	#global progressBar
	#progressBar = ttk.Progressbar(frame, length=450)
	#progressBar.pack(side='top')

	global resultListbox
	resultListbox = tk.Listbox(frame,background='#141414',fg='white',selectbackground='#2f6492',width=50,height=8)
	resultListbox.pack(side='top')

	returnButton = tk.Button(frame,text='  return  ',fg=lightLetterColor,command=lambda: mainApp(reload=True))
	returnButton.config(bg=appBgColor)
	returnButton.pack(side='bottom')


# test IDs: R1AB_CVHSA,R1AB_BCHK3,R1AB_CVMJH
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


def getMLocPage():
	root.title('Motif Locations')

	global canvas, frame
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

	returnButton = tk.Button(frame,text='  return  ',fg=lightLetterColor,command=lambda: mainApp(reload=True))
	returnButton.config(bg=appBgColor)
	returnButton.pack(side='bottom')


def getSeqPage():
	root.title('Get Sequence')

	global canvas, frame
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

	returnButton = tk.Button(frame,text='  return  ',fg=lightLetterColor,command=lambda: mainApp(reload=True))
	returnButton.config(bg=appBgColor)
	returnButton.pack(side='bottom')


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
pubmedLogo = tk.PhotoImage(file='assets/logo-pubmed.png')
helixLogo = tk.PhotoImage(file='assets/logo-helix.png')
appBgColor = '#181818'
lightLetterColor = '#3d3d3d'

# run home screen
mainApp(reload=False)
root.mainloop()




