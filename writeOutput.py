# Authors: Stefan Jokic, Adrian Taubner

from align import reverse_complement

def writeToFile(f, header, sn, ln, reads):
    # if(header):
    #     versionNumber = 1.4
    #     s0 = 'unsorted'
    #     # Header
    #     f.write('@HD' + '\t')
    #     f.write('VN:' + str(versionNumber) + '\t')
    #     f.write('SO:' + s0 + '\n')
    #     f.write('@SQ' + '\t')
    #     f.write('SN:' + sn + '\t')
    #     f.write('LN:' + str(ln) + '\n')

    # mapq = 0
    # rnext = '*'
    # pnext = 0
    # tlen = 0
    # seq = '*'
    # qual = '*'

    # # alignment section
    # if reads[5]:
    #     seq = reads[4]
    #     qual = reads[6]
    # else:
    #     seq = reverse_complement(reads[4])
    #     qual = reads[6][::-1]

    # f.write(reads[0] + '\t' + str(reads[1]) + '\t' + sn + '\t' + str(reads[2]) + '\t' + str(mapq) + '\t' + reads[3] + '\t'
    #     + rnext + '\t' + str(pnext) + '\t' + str(tlen) + '\t' + seq + '\t' + qual + '\n')
    
    #open extra_output.txt
    # with open('data_small/extra_output.txt', 'w') as f2:
    #     for i2, read in enumerate(reads):
    #         f2.write("read " + str(i2+1) + " " + str(read) + '\n')

    seen_indexes = set()
    #align the reads in the file
    i = 1
    reads = sorted(reads, key=lambda x: x[2])
    for read in reads:
        #only look at unique indexes
        if read[2] in seen_indexes:
            continue
        seen_indexes.add(read[2])

        print(read[2])
        start = 0
        if i >= read[2]: # overlap
            start = i - read[2]
        else: # gap
            i+= read[2] - i
            start = 0
            f.write('\n')
        # interior
        if start > len(read[4]):
            continue
        f.write(read[4][start:])
        i+=len(read[4]) - start
    