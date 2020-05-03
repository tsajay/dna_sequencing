import logging
from itertools import permutations 

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
    overlaps = set()
    globalKmers = {}
    kmersForRead = {}
    for r in reads:
        addKmersForRead(globalKmers, r, k)

    for _, rs in globalKmers.items():
        for r in rs:
            if not r in kmersForRead:
                kmersForRead[r] = set()
            kmersForRead[r].update(rs)
   
    logging.debug("KMERS dict: %s" %(globalKmers))
    logging.info("KMERS dict size: %d" %(len(globalKmers)))
    
    count = 0
    for r in reads:
        count += 1
        if (count % 100 == 0):
            logging.debug("Processing read %d" %count)
        for rfk  in kmersForRead[r]:
            if (r == rfk):
                continue
            if (overlapMinCheck(r, rfk, k)):
                overlaps.add((r, rfk))

    return overlaps

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
    start = a.find(b[:min_length], start)  # look for b's prefix in a
    if start == -1:  # no more occurrences to right
        return 0
    if b.startswith(a[start:]):
        return len(a)-start
    return 0
    

