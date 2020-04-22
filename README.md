# Helix
### DNA Motif Analyser

In molecular biology and DNA sequence analysis it is always very useful to identify patterns (motifs) in order to find shared functions between proteins or to identify genes shared by different species (in genomes). For this, Helix, a bioinformatic approach, presents an interactive user interface to analyze different DNA sequences in FASTA format (from UniProt or NCBI) to find shared motifs and their location.

## App
<img src="https://github.com/Lfquezada/Helix/blob/master/assets/helixSS1.png" width="48">

## Tools
Tool | Description
------------ | -------------
Find longest shared motif | Given `n` sequences (by protein IDs or file [FASTA format]), returns the longest(s) global subsequence(s).
Get consensus | Given `n` sequences (by protein IDs or file FASTA format), returns the consensus.
Get Weblogo | Given `n` sequences (by protein IDs or file FASTA format), save to a `.png` file the [Weblogo](https://weblogo.berkeley.edu)
Find motif locations | Given a subsequence `s` and a protein ID or file (FASTA format), returns all the locations of occurrence.
Get sequence | Retrieve a protein aminoacid sequence from UniProt and saves it into a file (FASTA format).

## Requirements
* python 3.x
* BioPython

## Usage
```
python3 helix.py
```

## Author
* **Luis Quezada**
* **Esteban del Valle**