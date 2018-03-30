
file = "/home/aleksandr/Documents/Phylomini/Influenza/test_kmer.fastq"


from Bio import SeqIO
import matplotlib.pyplot as plt


class Kmer_spectre:
    
    def __init__(self):
        pass
    
    def write_kmer_dictionary(self, file):
    # Данный метод позволяет записать словарь k-меров в текстовый документ
        with open(file, 'w') as file:
            for key, value in self.kmer_dict.items():
                file.write('%s:%s' % (key, value))
                file.write('\n')
    
    def read_kmer_dictionary(self, file):
    # Этот метод позволяет прочитать словать k-меров из текстового документа
        self.kmer_dict = {}
        with open(file, 'r') as txt:
            temp=txt.read().splitlines()
            for line in temp:
                key,value = line.split(':')
                self.kmer_dict[key]=int(value)
            self.kmer_size=len(key)
    
    def analyse(self, file, k,q):
    # Составляет словарь k-меров
        self.kmer_size=k
        self.kmer_dict={}
        for record in SeqIO.parse(file, "fastq"):
            
            for nucleotide in range(len(record.seq)-k+1):
                kmer = str(record.seq)[nucleotide:nucleotide+k]
                kmer_quality = record.letter_annotations["phred_quality"][nucleotide:nucleotide+k]
                
                if any(phred < q for phred in kmer_quality):
                    pass
                else:
                    if kmer in self.kmer_dict:
                        self.kmer_dict[kmer]+=1
                    else:
                        self.kmer_dict[kmer]=1
    
    def build_spectre(self):
    # Строит спектр k-меров
        self.quantities = []
        self.spectre_massive = {}
        for kmer in self.kmer_dict:
            quantity = self.kmer_dict[kmer]
            self.quantities.append(quantity)
            if str(quantity) in self.spectre_massive:
                self.spectre_massive[str(quantity)]+=1
            else:
                self.spectre_massive[str(quantity)]=1
    
    def vizualize(self, x, y):
    # Визуализирует спектр ввиде гистограммы
        plt.hist(self.quantities, max(fik.kmer_dict.values())+1, facecolor='green', alpha=0.75)
        plt.axis([0, x, 0, y])
        
    def approx(self, cut):
    # Приблизительно расчитывает размер генома по спектру k-меров
        self.genome_size=0
        for key,value in self.spectre_massive.items():
            if int(key) > cut:
                self.genome_size += int(key)*value
        print(round(self.genome_size/self.kmer_size), 0)
        
        
        

fik = Kmer_spectre()
fik.read_kmer_dictionary("/home/aleksandr/Documents/dictionary_1.txt")
fik.build_spectre()
fik.approx(20)

fik.analyse(file, 11, 35)
fik.kmer_dict
fik.build_spectre()
fik.spectre_massive
