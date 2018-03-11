# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 12:43:34 2018

@author: Александр
"""
# import of libraris

from Bio import SeqIO

# class K_mer announcement

class K_mer:
    
    def __init__(self, kmer_name, place, locus):
        self.counter = 1                            # counter
        self.seq = kmer_name                        # sequence of k-mer
        self.pos=[[locus, place]]                   # position locus , position of nucleotide
        self.summary=[]                             # raw of data

# Save the total data in comfortable format
    def summarize(self):
        self.summary = []
        self.summary.append(self.seq)
        self.summary.append(self.pos)
        print(self.summary)
        
# Recieve data about k-mer
    def recieve(self, place, locus):
        self.counter +=1
        self.pos.append([locus, place])
        
# Open fasta file and parse it
with open("C:/Users/Александр/Desktop/seq_y_pestis_2.fasta","r") as fasta:
    sequences = list(SeqIO.parse(fasta, "fasta"))
    
    K = 14                  # Length of k-mer
    kmer_dict = {}          # k-mer dictionary for saving k-mers

    for i in range(len(sequences)):
        sequence = str(sequences[i].seq)
        seq_len = len(sequence)
        Locus = str(sequences[i].name)
        
# k-mer search and writing to dictionary
        for position in range(seq_len-K+1):
            kmer = sequence[position:position+K]
            if kmer in kmer_dict:
                kmer_dict[kmer].recieve(position, Locus)
        
            else:
                kmer_dict[kmer]=K_mer(kmer, position, Locus)

# search of max counter
max_count = 0
for kmer in kmer_dict:
    m = kmer_dict[kmer]
    if m.counter > max_count:
        max_count = m.counter
# writing of k-mers with maximum counter
for kmer in kmer_dict:
    m = kmer_dict[kmer]
    if m.counter > max_count-1:
        m.summarize()