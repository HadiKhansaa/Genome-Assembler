# -*- coding: utf-8 -*-
# @Time    : 23.10.20
# @Authors  : Chenxi N. / Guillaume T.
# @FileName: input.py

## IMPORT

import os
import time

## FUNCTIONS

def read_genome(file):
    '''
    Reads a Text file (.txt) containing a genome sequence
    :param file: full path of the .fa file
    :return: a string
    '''

    time_start = time.time()
    f = open(file)
    genome = ""
    line = f.readline()
    while (line != ""):
        genome += line[:].strip("\n") #changed
        line = f.readline()
    time_stop = time.time()
    print("Loading time : ", str(time_stop-time_start)[0:5], " s")
    print('```````````````````````````````````````````````````````````````````````')
    # print(genome.upper())
    return genome.upper() #changed




def read_fq(file):
    '''
    Reads a text file containing read sequences
    :param file: full path of the text file
    :return: a list of strings
    '''

    time_start = time.time()
    with open(file, 'r') as f:
        reads = [line.strip().upper() for line in f if line.strip()] #changed
    time_stop = time.time()

    print("Loading time : ", str(time_stop - time_start)[0:5], " s")
    return reads


def read_sam(file, withPQline):
    '''
    Reads a Sam file (.sam)
    :param file: full path of the .fq file, boolean whether file contains PQ line
    :return: list of qname, flag, rname, pos, cigar
    '''

    time_start = time.time()
    f = open(file)
    qnames = []
    flags = []
    rname = []
    cigar = []
    pos = []
    print(f)
    f.readline()
    f.readline()
    if (withPQline):
        f.readline()

    while(True):
        line = f.readline()
        if line == "":
            break
        content = line.split('\t')
        #qnames.append(content[0])
        #flags.append(content[1])
        #rname.append(content[2])
        pos.append((content[3]))
        cigar.append(content[5])
    time_stop = time.time()

    print("Loading time : ", str(time_stop-time_start)[0:5], " s")
    return pos, cigar

## EXAMPLES

#x, y = read_sam("/home/adrian/PycharmProjects/BiomedicalEngineering/team09/data_small/output_tiny_30xCov.sam", True)

# Paths

#fa_path = os.getcwd() + "/data_small/genome.chr22.5K.fa"
#fq_path = os.getcwd() + "/data_small/output_tiny_30xCov1.fq"

# Loading

#genome = read_fa(fa_path)
#reads, qnames, sn = read_fq(fq_path)
#print(sn)
