# Helix
### DNA Motif Analyser

In molecular biology and DNA sequence analysis it is always very useful to identify patterns (motifs) in order to find shared functions between proteins or to identify genes shared by different species (in genomes). For this, Helix, a bioinformatic approach, presents an interactive user interface to analyze different DNA sequences in FASTA format (from UniProt or NCBI) to find shared motifs and their location.

## App
<p align="center">
  <img src="https://github.com/Lfquezada/Helix/blob/master/src/assets/helixSS1.png" width="500">
</p>

## Tools
Tool | Description
------------ | -------------
Find longest shared motif | Given `n` sequences (by protein IDs or file FASTA format), returns the longest global subsequence(s).
Get consensus | Given `n` sequences (by protein IDs or file FASTA format), returns the consensus.
Get Weblogo | Given `n` sequences (by protein IDs or file FASTA format), saves the [Weblogo](https://weblogo.berkeley.edu) to a `.png` file.
Find motif locations | Given a `subsequence` and a protein ID or file (FASTA format), returns all the locations of occurrence.
Get sequence | Retrieves a 1 or more protein aminoacid sequences from UniProt and saves it into a file (FASTA format).

## Requirements
* python 3.x
* [BioPython](https://biopython.org)

## Usage
```
python3 helix.py
```

## Author
* **Luis Quezada 18028**
* **Esteban del Valle 18221**