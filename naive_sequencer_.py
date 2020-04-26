import sequence_file_reader as reader
import logging

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def naiveWithNMismatches(p, t, n):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        numMismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                numMismatches += 1
                if numMismatches > n:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def naiveWithReverseComplement(p, t):
    r = reader.reverseComplement(p)
    logging.debug("rc: " + r)
    occurences = naive(p, t)
    if (p != r):
        occurences += naive(r, t)
    occurences.sort()
    return occurences
    

def naiveWithReverseComplement2(p, t):
    r = reader.reverseComplement(p)
    logging.debug("rc: " + r)
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
        match = True
        if not p == r:
            for j in range(len(p)):  # loop over characters
                if t[i+j] != r[j]:  # compare characters
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences
