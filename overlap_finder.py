import logging
from itertools import permutations 
import time as t
import itertools

logger = logging.getLogger()

def minEditDistance(ref, samp):
    patternLength = len(samp)
    numAlignments = len(ref) - patternLength + 1
    minDistance = len(samp)
    for start in range(numAlignments):
        editDist = editDistance(ref[start: start+patternLength], samp)
        minDistance = min(minDistance, editDist)
    
    return minDistance

def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    logging.debug("D = %s" %(D))
    return D[-1][-1]

def approximateMatch(ref, samp):
    # Create distance matrix
    D = []
    x = samp
    y = ref 
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = 0
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    logging.debug("D = %s" %(D))
    minEditDistance = min(D[-1])
    logging.debug("Min edit distance:  %d" %(min(D[-1])))
    return minEditDistance





def addKmersForRead(globalKmers, sample, k):
    numAlignments = len(sample) - k + 1
    for start in range(numAlignments):
        currentKmer = sample[start: start+k]
        if not (currentKmer in globalKmers):
            globalKmers[currentKmer] = set()
        globalKmers[currentKmer].add(sample)
    
    return

def findSourceNodes(overlaps):
    sources = set()
    for pair in overlaps:
        sources.add(pair[0])
    return sources

def findOverlapsForReads(reads, k):
    '''
    reads : Array of reads
    k     : Minimum overlap required
    return: Pairs of reads that overlap
    '''
    logging.info("Start: find overlaps")    
    start = t.perf_counter()
    overlaps = []
    globalKmers = {}
    candidatesForRead = {}
    readLengths = []
    for r in reads:
        addKmersForRead(globalKmers, r, k)
        candidatesForRead[r] = set()
    duration = t.perf_counter() - start
    logging.info("Done forming gobal dict of KMERS. Duration: %s secs" %duration)    
    logging.info("KMERS dict size: %d" %(len(globalKmers)))

    start = t.perf_counter()
    for rs in globalKmers.values():
        for r in rs:
            candidatesForRead[r].update(rs)
   
    for r in reads:
        #candidatesForRead[r].remove(r)
        readLengths.append(len(r))

    duration = t.perf_counter() - start
    logging.info("Done forming per read kmer list. Duration: %s secs" %duration)    

    logging.debug("KMERS dict: %s" %(globalKmers))
    
    start = t.perf_counter()
    count = 0
    
    for r in reads:
        if (count % 100 == 0):
            logging.debug("Processing read %d" %count)
        if (count % 1000 == 0):
            logging.info("Processing read %d" %count)
        for cfr  in candidatesForRead[r]:
            #if (overlapMinCheck(r, cfr, k)):
            if (r != cfr):
                ov = overlap_opt(r, cfr, readLengths[count], k)
                if (ov > 0):
                    overlaps.append((r, cfr, ov))
        count += 1
        
    duration = t.perf_counter() - start
    logging.info("Done forming overlap graphs. Duration: %s secs" %duration)    

    return overlaps

def overlap_opt(a, b, readLength, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return readLength-start
        start += 1  # move just past previous match

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def overlapMinCheck(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    start = a.find(b[:min_length])  # look for b's prefix in a
    if start == -1:  # no more occurrences to right
        return 0
    if b.startswith(a[start:]):
        return 1
    return 0

def pickMaximalOverlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    # perms = itertools.permutations(reads, 2)
    perms = findOverlapsForReads(reads, k)
    logging.info("# Permutations: %d" %len(list(perms)))
    for a, b, olen in perms:
        if (olen > best_olen):
            reada, readb = a, b
            best_olen  = olen
    return perms, reada, readb, best_olen

def greedy_scs(reads, k):
    numMerged = 0
    logging.log(21, "Starting with %d reads" %len(reads))
    _, reada, readb, olen = pickMaximalOverlap(reads, k)
    while olen > 0:
        reads.remove(reada)
        reads.remove(readb)
        reads.append(reada + readb[olen:])
        _, reada, readb, olen = pickMaximalOverlap(reads, k) 
        numMerged += 1
        if (numMerged % 10 == 0):
            logging.log(21, "Number of reads merged: %d, num nodes: %d " %(numMerged, len(reads)))

    return ''.join(reads)

def greedy_scs_experminetal(reads, k):
    numMerged = 0
    logging.log(21, "Starting with %d reads" %len(reads))
    perms, reada, readb, best_olen = pickMaximalOverlap(reads, k)


    while best_olen > 0:
        numMerged += 1
        if (numMerged %10 == 0):
            logging.log(21, "Merged %d nodes, graph size: %d edges" %(numMerged, len(perms)))
        merged_node = reada + readb[best_olen:]
        perms_temp = []

        for a, b, olen in perms:
            if (a != reada and b != readb):
                perms_temp.append((a, b, olen))
            if (a == reada):
                continue
            if (a != reada and b == readb):
                perms_temp.append((a, merged_node, olen))
        perms = perms_temp
        best_olen = 0
        for a, b, olen in perms:
            if (olen > best_olen):
                reada, readb = a, b
                best_olen  = olen

    scs = ''
    for a, b, olen in perms:
        scs += a + b[olen]

    return scs, numMerged

def de_bruijn_graph(reads, k):
    edges = []
    nodes = set()
    
    for r in reads:
        for start in range(len(r) -k + 1):
            edges.append( (r[start: start + k -1], r[start + 1: start + k]) )
            nodes.add(r[start: start + k -1]) 
            nodes.add(r[start + 1: start + k])
    return nodes, edges

def de_bruin_graph_scs(reads, k):
    scs = ''
    logging.info("Num reads: %d" %(len(reads)))
    nodes, edges = de_bruijn_graph(reads, k)
    logging.info("DeBruijn graph has %d nodes and %d edges" %(len(nodes), len(edges)))
    #scs = edges[0][0]
    #while (len(edges) > 0):

    wf = open("debruijn.out", "w")
    wf.write("%s" %edges)
    wf.close()
    return scs
    
    
    
def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if len(shortest_sup) == 0:
            shortest_sup.append(sup)  # found shorter superstring
            logging.debug("Len 0: adding %s" %sup )
        elif (len(sup) < len(shortest_sup[0])):
            shortest_sup.clear()
            shortest_sup.append(sup)  # found shorter superstring
            logging.debug("New minimum: clearing and adding %s" %sup )
        elif (len(sup) == len(shortest_sup[0])):
            shortest_sup.append(sup)  # found shorter superstring
            logging.debug("New same minimum: NOT clearing and adding %s" %sup )
            
            
    return shortest_sup  # return shortest
