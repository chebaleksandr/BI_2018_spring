
file = "/home/aleksandr/Documents/t.fastq"

from Bio import SeqIO
import numpy as np


class Trim:
    
    def __init__(self, file):
        self.name = file        # ИМЯ ОБРАБАТЫВАЕМОГО (ВХОДНОГО ИЗНАЧАЛЬНО) ФАЙЛА
        self.out = file + '_out'  # ИМЯ ВЫХОДНОГО ФАЙЛА
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
    
natrim = Trim(file)
natrim.head_cut(20)
natrim.sliding_window(10, 45)
    
 