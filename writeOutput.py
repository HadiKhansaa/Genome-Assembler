# Authors: Stefan Jokic, Adrian Taubner

# Authors: Stefan Jokic, Adrian Taubner

def alignMismatchesWithGaps(mismatches, gaps):
    '''
    mismatches = [[pos1, read1], [pos2, read2]] to be sorted
    gaps = [g1, g2, g3, g4...] which are sorted ascendingly 
    '''
    mismatches = sorted(mismatches, key=lambda x: x[0])
    i=-1 
    j=0
    new_contigs = []
    while (i< len(mismatches)-1 and j<len(gaps)):
        i+=1
        start,read_with_mismatches = mismatches[i]
        end = start + len(read_with_mismatches)-1 
        if(gaps[j]>end):
            continue
        while(j<len(gaps) and gaps[j]<start):
            j+=1
        if(j>=len(gaps)):
            break
        overlap_start = gaps[j]
        while j<len(gaps) and start<=gaps[j]<=end :
            j+=1
        overlap_end = gaps[j-1]
        # new contig which is read[overlap_start:overlap_end] 
        new_contigs.append([overlap_start,read_with_mismatches[overlap_start-start:overlap_end-start+1]])

    return new_contigs
    

def alignMatches(ln, reads):
    seen_indexes = set()
    #align the reads in the file
    i = 1
    first_index = 1
    contig = ''
    contigs=[]
    gaps = []
    last_index = -1
    reads = sorted(reads, key=lambda x: x[0])
    for read in reads:
        #only look at unique indexes
        if read[0] in seen_indexes:
            continue
        seen_indexes.add(read[0])
        start = 0
        if i >= read[0]: # overlap
            start = i - read[0]
        else: # gap
            if contig != '':
                contigs.append([first_index,contig])
            contig = ''
            gaps += range(i, read[0])
            i+= read[0] - i
            first_index = i
            start = 0
        # interior
        if start > len(read[1]):
            continue
        contig += read[1][start:]
        i+=len(read[1]) - start
    if contig!='':
        contigs.append([first_index,contig])
    gaps += range(i, ln)
    return contigs, gaps

def writeToFile(f, reads):
    for read in reads:
        f.write(f"{read[0]} \t {read[1]}\n")


# def writeToFile(f, header, sn, ln, reads):
#     seen_indexes = set()
#     #align the reads in the file
#     i = 1
#     reads = sorted(reads, key=lambda x: x[0])
#     for read in reads:
#         #only look at unique indexes
#         if read[0] in seen_indexes:
#             continue
#         seen_indexes.add(read[0])

#         # print(read[0])
#         start = 0
#         if i >= read[0]: # overlap
#             start = i - read[0]
#         else: # gap
#             i+= read[0] - i
#             start = 0
#             # f.write('\n')
#         # interior
#         if start > len(read[1]):
#             continue
#         # f.write(read[1][start:])
#         i+=len(read[1]) - start
    