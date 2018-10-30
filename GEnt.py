#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:44:24 2018

@author: skylar
"""

from Bio import AlignIO               # Read fasta files
from Bio import SubsMat
import argparse   
import math
 
PSEUDOCOUNT_MULTIPLIER = 0
def ewhole(alignments,pamMatrix):
    col = len(alignments[0].seq)
    gap = getRaw(alignments[1])
    aaCount = getRaw(alignments[0])
    consensus = getConsensus(alignments)
    aaTotal = getCountPCount(alignments)[0]
    pseudocount = getCountPCount(alignments)[1]
    aaPseudocount = [col][20]
    aaAjusted = [col][20]
    famEntropy = [col]

    
        
    for aa in range(len(aaCount[0])): #i is amino acid               Finds Ajusted values for count and individual pseudocount
        for p in range(col):#j is the positon
            if gap[p]:
                break
            aaPseudocount[p][aa] = pseudocount[p]*pamMatrix[aa]
            aaAjusted[p][aa] = aaTotal[aa]/(aaTotal[aa]+pseudocount[p])*aaCount[p][aa]/aaTotal[p]+pseudocount[p]/(aaTotal[aa]+pseudocount[p])*aaPseudocount/pseudocount[aa]
    for aa in range(col):               #Finds Family Entropy
        for p in range(len(aaCount[0])):
                famEntropy[aa] += math.log((aaCount[aa][p]/aaTotal[aa])/pamMatrix[p],2) 

                
def grpent(inGroup,outGroup,alignments,pamMatrix):
    col = len(alignments[0].seq)
    aaCount = getRaw(alignments[0])
    aaCountIn = getRaw(inGroup[0])
    aaCountOut = getRaw(outGroup[0])
    consensusIn = getConsensus(inGroup)
    consensusOut = getConsensus(outGroup)
    aaTotalIn  = getCountPCount(inGroup)[0]
    aaTotalOut = getCountPCount(outGroup)[0]
    psTotalIn = getCountPCount(inGroup)[1]
    psTotalOut = getCountPCount(outGroup)[1]
    aapsTotalIn = ajustCounts(inGroup,pamMatrix)[0]
    aapsTotalOut = ajustCounts(outGroup,pamMatrix)[0]
    aaAjustedIn = ajustCounts(inGroup,pamMatrix)[1]
    aaAjustedOut = ajustCounts(outGroup,pamMatrix)[1]
    groupEntropy = [col]

    
    for aa in range(col):               #Finds Family Entropy
        for p in range(len(aaCount[0])):
            groupEntropy[aa] += (aaCountIn[aa][p]-aaCountOut[aa]) * math.log((aaCountIn[aa][p]/aaCountOut[aa]))     
def crval(group):
    rows = len(group)
    col = len(group[0].seq)
    gap = getRaw(group)[1]
    aaCount = getRaw(group)[0]
    aapsTotal = getCountPCount(group)[1]
    aaTotal = getCountPCount(group)[0]
    consensus= getConsensus(alignments)
    pseudocount = ajustCounts(group,pamMatrix)[0]
    aaAjusted = ajustCounts(group,pamMatrix)[1]
    
    if(len(group)<3):
        exit
        
    
    
    
def zscore(alignments):
    return ""
def twoent():
    return ""
def main():
    parser = argparse.ArgumentParser(
        description='Calculate group entropy')
    parser.add_argument('settings', type=str,
        help='Settings file')
    f = open(parser.parse_args().settings,'r')
    lines = f.readlines
    f.close()
    pamMatrixName = int(lines[3].remove("Substitution Matrix: pam"))
    pamText = open(pamMatrixName,'r')
    count = 0
    pMat = [20][20]
    for line in pamText:
        pMat[count] = line.split(",")
        count+=1

    alignmentName  = lines[0].remove("Alignment File: ")
    groupsName = lines[1].remove("Groups File: ")
    a = open(alignmentName,'r')
    PSEUDOCOUNT_MULTIPLIER = int(lines[3].remove("Pseudocount Multiplier: "));
    r = AlignIO.read(a,"fasta")
    g = groups(r,open(groupsName,'r'))
    groupNames = g.grpDict.keys
    ewhole(r,pMat)
    for groupN in groupNames:
        grpent(g.getGroup(groupN)[0],g.getGroup(groupN)[1])
        crval(g.getGroup(groupN))
    zscore()
    
    
#Helper Functions
def remNPC(char_AA):
    aaVal = char_AA-65
    if aaVal > 0:
        aaVal-=1
    if aaVal > 7:
        aaVal-=1
    if aaVal > 11:
        aaVal-=1
    if aaVal > 16:
        aaVal-=1
    if aaVal > 18:
        aaVal-=1
    return aaVal

def retAA(intASCII):
    if intASCII > 0:
        intASCII+=1
    if intASCII > 8:
        intASCII+=1
    if intASCII > 13:
        intASCII+=1
    if intASCII > 19:
        intASCII+=1
    if intASCII > 22:
        intASCII+=1
        return chr(intASCII+65)
def getRaw(alignments):
    rows = len(alignments)
    col = len(alignments[0].seq)
    gap = [col]
    aaCount = [col][20]
    for aa in range(col):#i is the length of the sequence
        for p in range (rows):  #j is the number of alignments
            char = (alignments[p].seq)[aa]
            if char == '.':                 #Check for Gaps
                gap[aa] = True
                break
            else:
                gap = False
                aaCount[aa][remNPC(char)]+=1
    return {aaCount,gap}
def getConsensus(alignments):
    col = len(alignments[0])
    aaCount = getRaw(alignments[0])
    gap = getRaw(alignments[1])
    consensus = [col]
    for aa in range(len(aaCount[0])):
        maxV = 0
        temp = ''
        for p in range(col):       #Finds Consensus AA
            if gap[aa]:
                consensus[aa] = '.'
                break
            if maxV < aaCount[p][aa]:
                temp = retAA(p)
                maxV = aaCount[p][aa]
        consensus[aa] = temp
        return consensus
def getCountPCount(alignments):
    col = len(alignments[0])
    aaCount = getRaw(alignments[0])
    gap = getRaw(alignments[1])
    pseudocount = [col]
    aaTotal = [col]

    for aa in range(len(aaCount[0])):        #Finds total AA and pseudocount total
        acids = 0
        for p in range(col):
            if gap[aa]:
                break
            if aaCount[p][aa] > 0:
                acids +=1
        aaTotal[aa] = acids
        pseudocount[aa] = acids*PSEUDOCOUNT_MULTIPLIER
    return {aaTotal,pseudocount}
def ajustCounts(sequences,pamMatrix):
    aaTotal = getCountPCount(sequences)[0]
    aaCount = getRaw(sequences)[0]
    gap = getRaw(sequences)[1]
    aapsTotal = [len(sequences[0])][20]
    aaAjusted = [len(sequences[0])][20]
    psTotal = getCountPCount(sequences)[1]

    for aa in range(len(aaCount[0])): #i is amino acid               Finds Ajusted values for count and individual pseudocount
        for p in range(len(aaCount)):#j is the positon
            if gap[p]:
                break
            aapsTotal[p][aa] = aaCount[p][aa]/aaTotal[p][aa]*pamMatrix[p][aa]/sum(pamMatrix[p])
            aaAjusted[p][aa] = (aaCount[p][aa]+aapsTotal[p][aa])/(aaTotal[p][aa]+psTotal[p][aa])
    return {aapsTotal,aaAjusted}

class groups(object):
    grpDict  = dict()
    def __init__(self,alignments,groupFile):
        temp = ''
        src = ''
        for line in groupFile:
            if line[0:4] == 'Group':
                self.grpDict[line[6:]] = list()
        for line in groupFile:
            if line[0:4] == 'Group':
                temp = line[6:]
            elif line[0:1] != '//' and temp != '':
                src = line
            for i in range(len(alignments)):
                if alignments[i].title == src:
                    self.grpDict[temp].append(alignments[i])
                    break
    def getGroup(self,groupName):
        inGroup = []
        outGroup = []
        for k in self.grpDict.key():
            if k == groupName:
                inGroup == self.grpDict[k]
            else:
                outGroup.append(self.grpDict[k])
        return {inGroup,outGroup}
                
            

    
