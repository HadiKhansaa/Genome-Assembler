# Authors: Stefan Jokic, Adrian Taubner

import os
import sys
import gzip
import time

from align import getCigarFromIndex, reverse_complement
from bwt_linear import bwt_linear, index_bwt, search, search_in_bwt

# Updated import to match the modified read_fq function
from input import read_fq, read_genome
import numpy as np

from writeOutput import writeToFile

def getFlags(Read_isFwd, isMapped):
    flag = 1 + 64

    if isMapped:
        flag += 2

    if not Read_isFwd:
        flag += 16

    return flag


def run(path, readFile, genomeFile):
    read_file = path + readFile + ".txt"

    # Get files
    if not os.path.exists(read_file):
        sys.exit("ERROR: Could not find specified file: {}".format(read_file))

    print("Loading {} ...".format(genomeFile + ".txt"))
    genome = read_genome(path + genomeFile + ".txt")
    
    print("Loading {} ...".format(read_file))
    reads = read_fq(read_file)

    # initialize variables
    Read_isFwd = True

    print("")
    bwt_txt = input("Is BWT of given genome already stored in a .txt file? If so specify its file name without the .txt ending.\n Otherwise enter anything to continue. \n")
    sa_txt = input("Is suffix array of given genome already stored in a .txt file? If so specify its file name without the .txt ending.\n Otherwise enter anything to continue. \n")
    print("")

    if (os.path.exists(bwt_txt + ".txt") and os.path.exists(sa_txt + ".txt")):
        print("Loading BWT and suffix array from .txt files ...")
        time_start = time.time()
        bwt_file = open(bwt_txt + ".txt", 'r')
        genome_bwt = bwt_file.read()
        sa_file = open(sa_txt + ".txt", 'r')
        genome_sa = list(map(int, list(sa_file.read()[1:-1].split(", "))))
        time_stop = time.time()
        print("Finished loading BWT and suffix array from .txt files in " + str(time_stop - time_start)[0:5] + " s.\n")
    else:
        print("Constructing BWT and suffix array of genome ...")
        time_start = time.time()
        genome_bwt, genome_sa = bwt_linear(genome)
        time_stop = time.time()
        print("Finished constructing BWT and suffix array of genome in " + str(time_stop - time_start)[0:5] + " s.\n")

    c, Occ = index_bwt(genome_bwt)
    min_len_seed = 5
    outputPath = path + "output_" + readFile + ".txt" #changed
    if os.path.exists(outputPath):
        inp = ''
        while inp != 'y' or inp != 'n':
            inp = input("WARNING: Output .sam file already exists! Are you sure you want to continue (File will be overwritten) ? (y/n)\n")
            if inp == 'y':
                os.remove(outputPath)
                break
            elif inp == 'n':
                sys.exit("Aborted.")
            else:
                print("Please enter 'y' to continue or 'n' to abort.")

    file = open(outputPath, 'a')
    header = True

    all_allignments = [] #added
    for i in range(0, len(reads)):
        if i == 1:
            header = False
        print("Processing read #{} out of {}".format(i + 1, len(reads)))
        read_fwd = reads[i]
        # read_rvcmpl = reverse_complement(read_fwd)

        seedLength = len(read_fwd)

        hits_fwd = search(genome_bwt, read_fwd[0:seedLength], genome_sa, c, Occ)
        # hits_rvcmpl = search(genome_bwt, read_rvcmpl[0:seedLength], genome_sa, c, Occ)

        # sp, ep = search_in_bwt(genome_bwt, read_fwd[0:seedLength], c, Occ)

        # while seedLength > min_len_seed and len(hits_fwd) == 0 and len(hits_rvcmpl) == 0:
        #     seedLength = int(seedLength / 2)

        #     if seedLength < 2 * min_len_seed:
        #         seedLength = min_len_seed
        #     for j in range(int(len(read_fwd) / seedLength)):
        #         hits_fwd = search(genome_bwt, read_fwd[seedLength * j:seedLength * (j + 1)], genome_sa, c, Occ)
        #         hits_rvcmpl = search(genome_bwt, read_rvcmpl[seedLength * j:seedLength * (j + 1)], genome_sa, c, Occ)
        #         if len(hits_fwd) != 0 or len(hits_rvcmpl) != 0:
        #             break

        # minCostRead = sys.maxsize
        # minCigarRead = ''
        # minPosRead = 0

        # isMapped = False

        # if len(hits_fwd) or len(hits_rvcmpl):
        #     isMapped = True

        # # forward and reverse complement alignment
        # for hit in hits_fwd:
        #     leftIndex, rightIndex = 0, len(genome)
        #     leftIndex = max(hit - len(read_fwd), leftIndex)
        #     rightIndex = min(hit + len(read_rvcmpl) + 1, rightIndex)
        #     cigar, cost, index = getCigarFromIndex(read_fwd, genome[leftIndex: rightIndex])
        #     if cost < minCostRead:
        #         minCigarRead = cigar
        #         minPosRead = leftIndex + index + 1  # because MAPQ index starts at 1
        #         minCostRead = cost

        # for hit in hits_rvcmpl:
        #     leftIndex, rightIndex = 0, len(genome)
        #     leftIndex = max(hit - len(read_rvcmpl), leftIndex)
        #     rightIndex = min(hit + len(read_rvcmpl) + 1, rightIndex)
        #     cigar, cost, index = getCigarFromIndex(read_rvcmpl, genome[leftIndex: rightIndex])
        #     if cost < minCostRead:
        #         minCigarRead = cigar
        #         minPosRead = leftIndex + index + 1  # because MAPQ index starts at 1
        #         minCostRead = cost
        #         Read_isFwd = False

        # flag = getFlags(Read_isFwd, isMapped)

        # alignments = [str(i + 1), flag, minPosRead, minCigarRead, read_fwd, Read_isFwd, '']
        for hit in hits_fwd:
            alignments = [str(i + 1), 'flag', hit + 1, 'minCigarRead', read_fwd, Read_isFwd, ''] #changed
            all_allignments.append(alignments)
    writeToFile(file, header, '', len(genome), all_allignments) #changed

    file.close()


# Paths
if len(sys.argv) != 4:
    sys.exit("ERROR: Incorrect number of arguments.\n Please specify the following arugments and in this particular order: Name of directory, name of genome txt file and name of txt file with reads.\n Example: python main.py data_small genome.chr22.5K reads.txt")

if not os.path.exists(sys.argv[1]):
    sys.exit("ERROR: Could not find specified directory: {}".format(sys.argv[1]))

fq_path = sys.argv[3]
genome = sys.argv[2]
folder = os.getcwd() + "/" + sys.argv[1] + "/"

run(folder, fq_path, genome)
