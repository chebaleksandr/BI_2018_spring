# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 00:21:04 2017

"""
import re
import argparse
import fileinput

if __name__ == "__main__":
     
    parser = argparse.ArgumentParser(description = 'This script allow to find CG/AT regions of maximum lenght')

    parser.add_argument('-of', '--open_file',help='Open file with data', action='store')
    parser.add_argument('-sf', '--save_file', help='Save result in text file', action='store', default = '0')
    parser.add_argument('-t','--type', choices=['GC', 'AT'], help='Choice the nucleotide type',default='GC')

    args = parser.parse_args()

    data=args.save_file
    file=args.open_file

    with open(file,'r') as sourse:
        sequence=sourse.readlines()

    if args.type == 'GC':
        reg=re.compile("[G,C,g,c]{1,}")
    
    if args.type == 'AT':
        reg=re.compile("[A,T,a,t]{1,}")
        
    s=reg.findall(str(sequence))
    
    f=0
    for e in s:
        if len(e)>f:
            f=len(e)
            seq=e  
    
    print(seq,len(seq))
    if data != '0':
        with open(data, "w") as txt:
            txt.write(seq+' '+str(len(seq)))