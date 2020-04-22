# Helix
### DNA Motif Analyser

In molecular biology and DNA sequence analysis it is always very useful to identify patterns (motifs) in order to find shared functions between proteins or to identify genes shared by different species (in genomes). For this, Helix, a bioinformatic approach, presents an interactive user interface to analyze different DNA sequences in FASTA format (from UniProt or NCBI) to find shared motifs and their location.

## App
![SS1](/assets/helixSS1.png)

## Tools
Tool | About
------------ | -------------
Find longest shared motif | Given `<n>` sequences (aa or bp), returns the longest(s) global subsequence(s).
Find motif locations | Given a subsequence `<s>` and a protein ID or file name (FASTA format), returns all the locations of occurrence.
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