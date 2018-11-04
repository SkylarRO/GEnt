#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:44:24 2018

@author: skylar
"""

from Bio import AlignIO               # Read fasta files
import argparse   
import math
import openpyxl                     # Write spreadsheets

global PSEUDOCOUNT_MULTIPLIER 
def ewhole(alignments,pamMatrix,rProb):
    col = len(alignments[0].seq)
    temp =getRaw(alignments) 
    gap = temp[1]
    aaCount = temp[0]
    consensus = getConsensus(alignments)
    temp = getCountPCount(alignments)
    aaTotal = temp[0]
    pseudocount = temp[1]
    aaPseudocount = [[0] * 20 for i in range (col)]
    aaAjusted = [[0] * 20 for i in range (col)]
    famEntropy = [0]*col

        
    for aa in range(col): #i               Finds Ajusted values for count and individual pseudocount
        for p in range(20):
            if gap[aa]:
                break
            aaPseudocount[aa][p] = pseudocount[aa]*float(rProb[p])
            aaAjusted[aa][p] = aaTotal[aa]/(aaTotal[aa]+pseudocount[aa])*aaCount[aa][p]/aaTotal[aa]+pseudocount[aa]/(aaTotal[aa]+pseudocount[aa])*aaPseudocount[aa][p]/pseudocount[aa]
    for aa in range(col):               #Finds Family Entropy
        for p in range(20):
            if gap[aa]:
                break
            odds = aaAjusted[aa][p] / sum(aaAjusted[aa])
            famEntropy[aa] += odds*math.log(odds/float(rProb[p]),2)
                
    return [gap,consensus,famEntropy]

def grpent(inGroup,outGroup,alignments,pamMatrix):
    col = len(alignments[0].seq)
    aaCount = getRaw(alignments)[0]
    aaCountIn = getRaw(inGroup)[0]
    gap = getRaw(inGroup)[1]
    aaCountOut = getRaw(outGroup)[0]
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
    ingroupEntropy = [col]
    outgroupEntropy = [col]
    totalgroupEntropy = [col]

    
    for aa in range(col):               #Finds Family Entropy
        for p in range(len(aaCount[0])):
            if gap[aa]:
                break
            ingroupEntropy[aa] = (aaCountIn[aa][p]-aaCountOut[aa]) * math.log((aaCountIn[aa][p]/aaCountOut[aa]))  
            outgroupEntropy[aa] = (aaCountIn[aa][p]-aaCountOut[aa]) * math.log((aaCountIn[aa][p]/aaCountOut[aa]))  
            totalgroupEntropy[aa] = ingroupEntropy[aa] + -outgroupEntropy[aa]
    return {outgroupEntropy,ingroupEntropy,totalgroupEntropy,consensusIn}
def crval(group):
    rows = len(group)
    col = len(group[0].SeqRecord)
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
    
#Helper Functions
def remNPC(char_AA):
    aaVal = ord(char_AA)-65
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
    gap = [False]*col
    temp = ""
    aaCount = [[0]*20 for i in range(col)]
    for aa in range(col):#i is the length of the sequence
        for p in range (rows):  #j is the number of alignments
            char = (alignments[p].seq)[aa]
            if char == '.':                 #Check for Gaps
                gap[aa] = True
                
    for aa in range(col):#i is the length of the sequence
        for p in range (rows):  #j is the number of alignments
            char = (alignments[p].seq)[aa]

            if(gap[aa] != True):
                temp +=char
                aaCount[aa][remNPC(char)]+=1
                
    return [aaCount,gap]
def getConsensus(alignments):
    col = len(alignments[0].seq)
    aaCount = getRaw(alignments)[0]
    gap = getRaw(alignments)[1]
    consensus = ['']*col
    for aa in range(col):
        maxV = 0
        temp = ''
        for p in range(20):       #Finds Consensus AA
            if gap[aa]:
                consensus[aa] = '.'
                break
            elif maxV < aaCount[aa][p]:
                temp = retAA(p)
                maxV = aaCount[aa][p]
                consensus[aa] = temp        
    return consensus

def getCountPCount(alignments):
    global PSEUDOCOUNT_MULTIPLIER
    col = len(alignments[0].seq)
    aaCount = getRaw(alignments)[0]
    gap = getRaw(alignments)[1]
    pseudocount = [0]*col
    aaTotal = [0]*col

    for aa in range(col):        #Finds total AA and pseudocount total
        acids = 0
        for p in range(20):
            if gap[aa]:
                break
            if aaCount[aa][p] > 0:
                acids +=1
        aaTotal[aa] = acids
        pseudocount[aa] = acids* PSEUDOCOUNT_MULTIPLIER
    return [aaTotal,pseudocount]
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
            aapsTotal[aa][p] = aaCount[aa][p]/aaTotal[aa][p]*pamMatrix[aa][p]/sum(pamMatrix[p])
            aaAjusted[aa][p] = (aaCount[aa][p]+aapsTotal[aa][p])/(aaTotal[aa][p]+psTotal[aa][p])
    return {aapsTotal,aaAjusted}

class groups(object):
    grpDict  = dict()
    def __init__(self,alignments,groupFile):
        temp = ''
        src = ''
        for line in groupFile:
            comp = line[0:5]
            if comp == "Group":
                self.grpDict[line[6:len(line)-1]] = list()
        groupFile.seek(0)
        for line in groupFile:
            if line[0:5] == "Group":
                temp = line[6:len(line)-1]
            elif line[0:1] != '//' and temp != '':
                src = line[:len(line)-1]
                for i in range(len(alignments)):
                    t = alignments[i].name
                    if t == src:
                        self.grpDict[temp].append(alignments[i])
                        break
    def getGroup(self,groupName):
        inGroup = []
        outGroup = []
        
        for k in self.grpDict.keys():
            if k == groupName:
                inGroup.extend(self.grpDict[k])
            else:
                outGroup.extend(self.grpDict[k])
        return [inGroup,outGroup]
    def getAllGrouped(self):
        outGroup = []
        for k in list(self.grpDict.keys()):
            outGroup += self.grpDict[k]
        return outGroup
                
def matMake(matrix, input):
    num = int(input[len(input)-3:])
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j] = pow(float(matrix[i][j]),num)
    return matrix
def main():
    parser = argparse.ArgumentParser(
        description='Calculate group entropy')
    parser.add_argument('settings', type=str,
        help='Settings file')
    f = open(parser.parse_args().settings,'r')
    lines = f.readlines()
    global PSEUDOCOUNT_MULTIPLIER
    PSEUDOCOUNT_MULTIPLIER = float(lines[2][lines[2].rfind(":")+2:len(lines[2])-1])

    f.close()
    alignmentName  = (lines[0][lines[0].rfind(":")+2:len(lines[0])-1])
    a = open(alignmentName,'r')
    r = AlignIO.read(a,"fasta")
    MatrixName = (lines[3][lines[3].rfind(":")+2:len(lines[3])-1])
    matrix = open(MatrixName[:len(MatrixName)-3]+".dat").read()
    matrix = [item.split() for item in matrix.split('\n')[:-1]]
    freq_ran = matrix[20]
    matrix = matrix[0:19]
    matrix = matMake(matrix,MatrixName)
    groupsName = (lines[1][lines[1].rfind(":")+2:len(lines[1])-1])

    g = groups(r,open(groupsName,'r'))
    groupNames = list(g.grpDict.keys())
    fe = ewhole(g.getAllGrouped(),matrix,freq_ran)
    for groupN in groupNames:
        ge = grpent(g.getGroup(groupN)[0],g.getGroup(groupN)[1],r,matrix)
        out = open(groupN + ".csv",'w')
        out.write("Position,Family Entropy,Group Entropy,Partial Group Entropy,Partial Out Group Entropy,Highest Group AA,Highest Family AA")
        for i in range(len(r)):
            if fe[0][i] == False:
                out.write(i+','+fe[1][i] + ','+ge[2][i] + ','+ge[1][i] + ','+ge[0][i] + ','+ge[3][i] +','+fe[2][i])   
                
if __name__ == '__main__':
    main()
