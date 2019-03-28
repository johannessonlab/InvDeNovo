#!/usr/bin/env python

# InvDeNovo_OneBreakpoint.py
# Version 0.2
# Author: Jesper Svedberg, Department of Organismal Biology, Uppsala University
# E-mail: jesper.svedberg@ebc.uu.se
#
# Companion script to InvDeNovo.py. The script tries to identify inversion in fragmented denovo
# assemblies, but it only looks for single breakpoints instead of both breakpoints of an inversion.
# This makes it possible to identify inversions where one breakpoint coincides with breaks in the
# fragmented assembly (ie, the breakpoint is located between contigs), or when you have complex
# tandem or overlapping inversions, but it also increases the false positive error rate. This
# calls for careful manual curation, and we've found that many candidate inversions are actually
# caused by assembly errors.
# The script will look for single contigs that map at two different locations at the target chromosome
# and output the two alignment locations and the distance between them.
#
# The script takes MUMmer BTAB files as input and outputs to the terminal.
#
# Usage: InvDeNovo_OneBreakpoint.py -c target_chromosome [-s minimum_alignment_block_length -d minimum_alignment_distance -v] input.btab
#
# -s sets minimum alignment block length. Repeats often have very short alignment blocks, and this
# parameter can filter out those.
# -d sets the minimum distance between the two alignment coordinates.
# -v sets verbose output.
#
# Non-verbose output format:
# Inversion candidate: target_chromosome vs. contig	alignment_distance    coord_alignment1    coord_alignment2
#
# ----------------------------
# This software is released under the Creative Commons Attribution-NonCommercial-ShareAlike
# 3.0 Unported (CC BY-NC-SA 3.0) License. See https://creativecommons.org/licenses/by-nc-sa/3.0/
# for details.


import sys
import csv
from operator import itemgetter
import argparse
#import time


SIMPLE_TYPES = "IIIISS"
BTAB_TYPES = "SSISSSIIIIFFIIISISIIIS"
RMOUT_TYPES = "IFFFSIISSSSSSSIS"
MAKER_TYPES = "SSSIISSSS"

QNAME = 0
QLEN = 2
RNAME = 5
QSTART = 6
QEND = 7
RSTART = 8
REND = 9
QDIR = 17
RLEN = 18
ALEN = 12

TRANSLIMIT = 100000
SIZE_LIMIT = 5000


# Converts entries in a table from strings to whatever format that is specified in the type list,
# for instance integers or floats. The type list is now formatted for the MUMmer show-coords
# btab format.
def tableTypeConvert(table, typeList):
    outTable = []
    for line in table:
        newLine = []
        i = 0
        for entry in line:
            if typeList[i] == "S":
                newLine.append(entry)
            elif typeList[i] == "I":
                newLine.append(int(entry))
            elif typeList[i] == "F":
                newLine.append(float(entry))
            i += 1
        
        outTable.append(newLine)
        
    return outTable

def tableTypeConvertIter(table, typeList):
    outTable = []
    for x, line in enumerate(table):
        newLine = []
        i = 0
        for entry in line:
            if typeList[i] == "S":
                newLine.append(entry)
            elif typeList[i] == "I":
                newLine.append(int(entry))
            elif typeList[i] == "F":
                newLine.append(float(entry))
            i += 1
        
        newLine.append(x)
        outTable.append(newLine)
        
    return outTable

# Sorts a table in list or tuple form after several columns in reversed order
def sortTable(table, columns):
    for col in columns:
        table = sorted(table, key=itemgetter(col))
    return table

# Returns column from matrix
def getColumn(matrix,i):
    f = itemgetter(i)
    return map(f,matrix)

# Concatenates a list of string into a string separated by "separator".
def catWords(wordlist, separator=", "):
    return separator.join(wordlist)

# Concatenates a list of arbitrary types into a string, with elements separated by "separator".
def catList(inlist, separator=", "):
    return separator.join([str(x) for x in inlist])

def reversePair(pairString, sep=":"):
    plist = pairString.split(sep)
    
    for p in plist:
        if p[0] == "-":
            p = p[1:]
        else:
            p = "-" + p
    
    return plist[1] + sep + plist[0]

def inversePair(pairString, sep=":"):
    plist = pairString.split(sep)
    olist = []
    
    for p in plist:
        if p[0] == "-":
            p = p[1:]
        else:
            p = "-" + p
        olist.append(p)
    
    return olist[0] + sep + olist[1]


# Creates a dictionary from a BTAB file with with the key query_contig:reference_contig,
# which contains a table of all lines where query_contig aligned with reference_contig.
# Can filter for a specific reference contig, as specified by the second argument.
def createContigPairDict(table, contig="", sizelimit=0):
    pairDict = {}
    for line in table:
        if contig:
            if line[RNAME] != contig:
                continue
            
        if sizelimit:
            if line[ALEN] < sizelimit:
                continue
        
        pairName = line[QNAME] + ":" + line[RNAME]
    
        if pairName in pairDict:
            pairDict[pairName].append(line)
        else:
            pairDict[pairName] = [line]
    
    return pairDict

def createPairDict(pairDict, table, maxLines, pairJump):

    rightmost = 0

    for i in range(maxLines):
        if table[i][17] == "Minus":
            fpname = "-" + table[i][0]
        else:
            fpname = table[i][0]
            
        if table[i+pairJump][17] == "Minus":
            lpname = "-" + table[i+pairJump][0]
        else:
            lpname = table[i+pairJump][0]
        
        if table[i][5] != table[i+pairJump][5]:
            rightmost = 0
            continue
           
        if fpname == lpname:
            continue
        '''
        elif fpname.replace("-","") == lpname.replace("-",""):
            continue
        '''
        isInternal = False
        if fpname.replace("-","") == lpname.replace("-",""):
            if ALLOW_INTERNAL:
                isInternal = True
            else:
                continue
        '''
        if isOverlapped(i, i+pairJump, table):
            continue
        
        overlap = getOverlap(i, i+pairJump, table)

        if overlap > OVERLAP_CUTOFF:
            #print lpname + " " + str(overlap)
            continue
        '''
        if table[i][9] > rightmost:
            rightmost = table[i][9]
        
        if table[i][-1] == "overlapped" or table[i+pairJump][-1] == "overlapped":
            continue

        if table[i+pairJump][9] < rightmost:
            table[i+pairJump].append("overlapped")
            #print lpname + " is overlapping. " + str(i+pairJump+1)
            continue
        else:
            overlap = getOverlap(i, i+pairJump, table)
            if overlap > OVERLAP_CUTOFF:
                #print lpname + " is overlapping with: " + str(overlap) + " " + str(i+pairJump+1)
                table[i+pairJump].append("overlapped")
                continue

        pairName = fpname + ":" + lpname
        pairStats = [i, pairJump-1, table[i][5], table[i][9], scoreEnv(i, i+pairJump, table, isInternal), scoreLength(table[i][12]), scoreLength(table[i+pairJump][12]), scorePercent(table[i]), scorePercent(table[i+pairJump]), scoreGap(pairJump-1)]
        pairStats.append(sum(pairStats[4:]))
        
        if pairName in pairDict:
            pairDict[pairName][1].append(pairStats)
        else:
            pairDict[pairName] = [True, [pairStats]]
    
    return pairDict

# Creates dictionary of all reference contigs with btab line numbers for each alignment
# in a list at each nucleotide position.
def alignCoverage(btab):
    contigDict = {}
    for i, row in enumerate(btab):
        if row[5] in contigDict:
            #print "   " + row[0]
            for j in xrange(row[8],row[9]):
                contigDict[row[5]][j].append(i)
        else:
            print row[5]
            cLen = int(row[18])
            #print cLen
            st=time.time()
            cgList = [[] for _ in xrange(cLen)]
            end=time.time()
            print end - st
            print sys.getsizeof(cgList)
            #print len(cgList)
            for j in xrange(row[8],row[9]):
                cgList[j].append(i)
            contigDict[row[5]] = cgList
    return contigDict

# Calculates average coverage for a portion of a reference contig (between start and
# stop.
def avgCoverage(contiglist, start, stop):
    alist = contiglist[start:stop]
    covSum = 0
    for pos in alist:
        covSum += len(pos)
    return covSum / float(len(alist))

def findMultiCoverage(contiglist, start, stop):
    alist = contiglist[start:stop]
    covSum = 0
    for pos in alist:
        if len(pos) > 1:
            covSum += 1
    return covSum / float(len(alist))

def findOverlap(btab, contig1, contig2, contiglist, start, stop):
    alist = contiglist[start:stop]
    covSum = 0
    for pos in alist:
        if len(pos) < 2:
            continue
        cList = getAlignedContigs(btab, pos)
        if contig1 == contig2:
            if cList.count(contig1) > 1:
                covSum += 1
        else:
            if (contig1 in cList) and (contig2 in cList):
                covSum += 1
    return covSum / float(len(alist))

def getAlignedContigs(btab, row_list):
    outlist = []
    for i in row_list:
        outlist.append(btab[i][0])
    return outlist

def isSameContig(btab, row1, row2):
    if btab[row1][0] == btab[row2][0]:
        return True
    else:
        return False

def printInvList(ll, btab, form="Normal"):
    if form == "Normal":
        out = str(ll[0]) + ", " + str(ll[1]) + ", "
        out += ", ".join([str(x) for x in ll[3:6]]) + ", "+ ", ".join(["%.2f" % x for x in ll[7:]])
        return out

def isOverlapped(row1, row2, btab):
    a1 = btab[row1]
    a2 = btab[row2]
    
    if a1[2] < a2[2]:
        start = a1[8]
        stop = a1[9]
    else:
        start = a2[8]
        stop = a2[9]

    olapRatio = findOverlap(btab, a1[0], a2[0], alignDict[a1[5]], start, stop)

    if olapRatio > 0.1:    
        print a1[0] + ":" + a2[0] + " " + str(olapRatio)

    if olapRatio > OVERLAP_CUTOFF:
        return True
    else:
        return False
    
def getOverlap(row1, row2, btab):
    a1 = btab[row1]
    a2 = btab[row2]
    
    start = a2[8]
    stop = a2[9]

    coverage = 0

    if (a2[8] > a1[8]) and (a2[8] < a1[9]):
        #a1Len = a1[9]-a1[8]
        a2Len = a2[9]-a2[8]

        if (a1[9] > a2[9]):
            coverage = 1
        else:
            coverage = 1-((a2[9]-a1[9])/float(a2Len))

    return coverage


parser = argparse.ArgumentParser(description='Identifies inversions in fragmented denovo assemblies by looking for ONE breakpoints.')
parser.add_argument('filename', help='MUMmer BTAB file')
parser.add_argument('-c','--chromosome',help='Target chromosome.', required=True)
parser.add_argument('-s','--size',help='Minimum alignment size. Default=2000', default=2000)
parser.add_argument('-d','--distance',help='Minimum alignment distance. Default=100000', default=100000)
parser.add_argument('-v','--verbose',help='Verbose output. Default=False', dest='verbose', action='store_true')
args = parser.parse_args()

tabFile = args.filename
size_limit = int(args.size)
isVerbose = args.verbose
target_contig = args.chromosome
TRANSLIMIT = int(args.distance)

'''
isVerbose = False

if len(sys.argv) <= 3:
    print "Usage: findTlInv.py mummer_out.btab target_contig min_alignment_size [verbose]"
    print "Looks for translocations and inversions in a mummer alignment."
    print "'verbose' gives verbose output. Anything else gives simple output."
    print "Output to terminal."
    sys.exit()
elif len(sys.argv) == 4:
    target_contig = sys.argv[2]
    size_limit = sys.argv[3]
elif len(sys.argv) == 5:
    if sys.argv[4] == "verbose":
        isVerbose = True
    target_contig = sys.argv[2]
    size_limit = sys.argv[3]
    
elif len(sys.argv) == 2:
    target_contig = ""
elif len(sys.argv) == 3:
    if sys.argv[2] == "verbose":
        isVerbose = True
        target_contig = ""
    else:
        isVerbose = False
        target_contig = sys.argv[2]


    
tabFile = sys.argv[1]
''' 


table = tableTypeConvertIter([line for line in csv.reader(open(tabFile),delimiter='\t')], BTAB_TYPES)

#alignDict = alignCoverage(table)
#print "alignDict finished"

pairDict = createContigPairDict(table, target_contig, SIZE_LIMIT)

for key, value2 in pairDict.iteritems():
    value = sortTable(value2, [RSTART])
    lenTable = len(value)
    
    if lenTable >= 2:
        gapTable = []
        rlast = 0
        for line in value:
            if rlast == 0:
                rlast = line[REND]
                qlast = line[QEND]
                dirlast = line[QDIR]
                continue
            else:
                if dirlast != line[QDIR]:
                    flip = True
                else:
                    flip = False
                gapTable.append([line[QNAME], line[RNAME], line[QLEN], line[RLEN], line[QSTART]-qlast, line[RSTART]-rlast, flip])
                
                rlast = line[REND]
                qlast = line[QEND]
                dirlast = line[QDIR]

        for i, gap in enumerate(gapTable):
            
            if abs(gap[5]) < TRANSLIMIT:
                #print "test < translimit"
                continue
            else:
                if i+1 < len(gapTable):
                    if abs(value[i+2][RSTART]-value[i][REND]) < TRANSLIMIT:
                        continue
                    else:
                        if gap[-1]:
                            print "Inversion candidate: " + gap[1] + " vs. " + gap[0] + "\t" + str(value[i][REND]) + "\t" + str(value[i+1][RSTART]) + "\t" + str(gap[5])
                        else:
                            print "Transloc. candidate: " + gap[1] + " vs. " + gap[0] + "\t" + str(value[i][REND]) + "\t" + str(value[i+1][RSTART]) + "\t" + str(gap[5])
                        if isVerbose:
                            for line in gapTable:
                                print catList(line, "\t")
                            for line in value:
                                print catList(line, "\t")

print "\n"
