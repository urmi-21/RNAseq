#!/bin/bash
#output ORFs with start codon but no stop codon
getorf -sequence Homo_sapiens.GRCh38.dna_sm.toplevel.fa -outseq Homo_sapiens.GRCh38.dna_sm.toplevel.fa.orf -minsize 200 -find 3
