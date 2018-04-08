#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  8 19:57:47 2018

@author: aleksandr
"""

import re
import argparse
import fileinput
from Bio import SeqIO
import numpy as np

if __name__ == "__main__":
     
    parser = argparse.ArgumentParser(description = 'This script allow to find CG/AT regions of maximum lenght')

    parser.add_argument('-i', '--input_file',help='Open file with data', action='store',required=True)
    parser.add_argument('-o', '--output_file', help='Save result in file', action='store',required=True)
    parser.add_argument('-WINDOW','--slidind_window_params', help='Enter the size of sliding window and mean quality', nargs = 2, default = [0,0])
    parser.add_argument('-HEADCROP','--head', help='Enter the number of nucleotides to be trimmed in the head of reads', action='store', default = 0)
    parser.add_argument('-TAILCROP','--tail', help='Enter the number of nucleotides to be trimmed in the tail of reads', action='store', default = 0)

    args = parser.parse_args()

    file = args.input_file
    file_2 = args.ouput_file
    size_of_window = args.window_size[0]
    quality = args.window_size[1]
    head = args.head
    tail = args.tail
        
    class Trim: 
            def __init__(self, file, file_2):
                self.name = file        # ИМЯ ОБРАБАТЫВАЕМОГО (ВХОДНОГО ИЗНАЧАЛЬНО) ФАЙЛА
                self.out = file_2  # ИМЯ ВЫХОДНОГО ФАЙЛА
                self.trimmed=[]         # ВРЕМЕННЫЙ СПИСОК
            
            def head_cut(self, head): #  ОБРЕЗАЕТ ГОЛОВУ И СОХРАНЯЕТ В ФАЙЛ 
                for record in SeqIO.parse(self.name, "fastq"):
                    if len(str(record.seq))>head:
                        self.trimmed.append(record[head:])
                with open(self.out,'w') as trim:
                    for element in self.trimmed:
                        trim.write(element.format("fastq"))
                self.name = self.out
                self.trimmed.clear()
        

            def tail_cut(self, tail): #  ОБРЕЗАЕТ ХВОСТ И СОХРАНЯЕТ В ФАЙЛ
                for record in SeqIO.parse(self.name, "fastq"):
                    if len(str(record.seq))>tail:
                        self.trimmed.append(record[:tail*(-1)])
                with open(self.out,'w') as trim:
                    for element in self.trimmed:
                        trim.write(element.format("fastq"))
                self.name = self.out
                self.trimmed.clear()
            
            def sliding_window(self, size, qual): # СКОЛЬЗЯЩЕЕ ОКНО ПРОВЕРЯЕТ СРЕДНЕЕ ВНУТРИ И ЗАПИСЫВАЕТ НЕНУЛЕВЫЕ И УДВОЛЕТВОРЯЮЩИЕ ПОСЛЕДОВАТЕЛЬНОСТИ
                self.window_size=size
                self.window_quality=qual       
                
                for record in SeqIO.parse(self.name, "fastq"):
                    n_number=len(record.seq)-size+1
                    position=0
                    for nucleotide in range(len(record.seq)-size+1):
                        #window = str(record.seq)[nucleotide:nucleotide+size]
                        quality_of_window = record.letter_annotations["phred_quality"][nucleotide:nucleotide+size]
                        
                        #print(quality_of_window)
                        
                        if np.mean(quality_of_window)>qual:
                            if position==n_number-1:
                                self.trimmed.append(record)
                        else:
                            if position!=0:
                                self.trimmed.append(record[:position*(-1)])
                            break
                        position+=1
                        
                with open(self.out,'w') as trim:
                    for element in self.trimmed:
                        trim.write(element.format("fastq"))
                self.name = self.out
                self.trimmed.clear()
    
    mygod = Trim(file, file_2)
    if head != 0:
        mygod.head_cut(head)
    if tail != 0:
        mygod.tail_cut(tail)
    if size_of_window != 0:
        mygod.sliding_window(size_of_window, quality)
    
    
    