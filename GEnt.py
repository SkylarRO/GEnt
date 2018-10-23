#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:44:24 2018

@author: skylar
"""

from Bio import AlignIO               # Read fasta files
import argparse   
import math
 
PSEUDOCOUNT_MULTIPLIER = 0
def ewhole(alignments,probability): #toDo FIND A PROBABILITY TABLE
    rows = len(alignments)
    col = len(alignments[0].seq)
    gap = [col]
    aaCount = [col][20]
    aaTotal = [col]
    consensus= [col]
    pseudocount = [col]
    aaPseudocount = [col][20]
    aaAjusted = [col][20]
    famEntropy = [col]

    for aa in range(col):#i is the length of the sequence
        for p in range (rows):  #j is the number of alignments
            char = (alignments[p].seq)[aa]
            if char == '.':                 #Check for Gaps
                gap[aa] = True
                break
            else:
                gap = False
                aaCount[aa][remNPC(char)]+=1     #Takes ascii value of char to avoid a much longer if statment complex to figur out what to
                                                #increment If you didn't do stuff with ascii, don't mess with this.
    for aa in range(len(aaCount[0])):
        maxV = 0
        temp = ''
        for p in range(col):       #Finds Consensus AA
            if gap[aa]:
                consensus[aa] = '.'
                break
            if maxV < aaCount[p][aa]:
                temp = retAA(j)
                maxV = aaCount[p][aa]
        consensus[aa] = temp
    for aa in range(len(aaCount[0])):        #Finds total AA and pseudocount total
        acids = 0
        for p in range(col):
            if gap[aa]:
                break
            if aaCount[p][aa] > 0:
                acids +=1
        aaTotal[aa] = acids
        pseudocount[aa] = acids*PSEUDOCOUNT_MULTIPLIER
    for aa in range(len(aaCount[0])): #i is amino acid               Finds Ajusted values for count and individual pseudocount
        for p in range(col):#j is the positon
            if gap[p]:
                break
            aaPseudocount[p][aa] = pseudocount[p]*probability[aa]
            aaAjusted[p][aa] = aaTotal[aa]/(aaTotal[aa]+pseudocount[p])*aaCount[p][aa]/aaTotal[p]+pseudocount[p]/(aaTotal[aa]+pseudocount[p])*aaPseudocount/pseudocount[aa]
    for aa in range(col):               #Finds Family Entropy
        for p in range(len(aaCount[0])):
                famEntropy[aa] += math.log((aaCount[aa][p]/aaTotal[aa])/probability[p],2)                                
def grpent(alignments,group,pamMatrix):
    rows = len(alignments)
    col = len(alignments[0].seq)
    aaCount = [col][20]
    aaCountIn = [col][20]
    aaCountOut = [col][20]
    consensusIn = [col]
    consensusOut = [col]
    aaTotalIn  = [col]
    aaTotalOut = [col]
    psTotalIn = [col]
    psTotalOut = [col]
    aapsTotalIn = [col][20]    
    aapsTotalOut = [col][20]
    gap = [col]
    aaAjustedIn = [col][20]
    aaAjustedOut = [col][20]
    groupEntropy = [col]
    for p in range(col):#i is the length of the sequence
        for a in range (rows):  #j is the number of alignments
            char = (alignments[a].seq)[p]
            if char == '.':                 #Check for Gaps
                gap[p] = True
                break
            else:
                gap = False
                aaCount[p][remNPC(char)]+=1     #Takes ascii value of char to avoid a much longer if statment complex to figur out what to
                                                #increment If you didn't do stuff with ascii, don't mess with this.
    for aa in range(len(aaCountIn[0])):
        maxV = 0
        temp = ''
        for p in range(len(aaCountIn)):
            if gap[aa]:
                consensusIn[aa] = '.'
                break
            if maxV < aaCountIn[p][aa]:
                temp = retAA(j)
                maxV = aaCountIn[p][aa]
        consensusIn[aa] = temp
        for aa in range(len(aaCountOut[0])):
            maxV = 0
            temp = ''
            for p in range(len(aaCountOut)):
                if gap[aa]:
                    consensusOut[aa] = '.'
                    break
                if maxV < aaCountOut[p][aa]:
                    temp = retAA(p)
                    maxV = aaCountOut[p][aa]
            consensusOut[aa] = temp        
    for aa in range(len(aaCountIn[0])):        #Finds total AA and pseudocount total
        acids = 0
        for p in range(len(aaCountIn)):
            if gap[aa]:
                break
            if aaCountIn[p][aa] > 0:
                acids +=1
                
        aaTotalIn[aa] = acids
        psTotalIn[aa] = acids*PSEUDOCOUNT_MULTIPLIER
        
    for aa in range(len(aaCountOut[0])):        #Finds total AA and pseudocount total
        acids = 0
        for p in range(len(aaCountOut)):
            if gap[aa]:
                break
            if aaCountOut[p][aa] > 0:
                acids +=1
        aaTotalOut[aa] = acids
        psTotalOut[aa] = acids*PSEUDOCOUNT_MULTIPLIER
        
    for aa in range(len(aaCountOut[0])): #i is amino acid               Finds Ajusted values for count and individual pseudocount
        for p in range(len(aaCountOut)):#j is the positon
            if gap[p]:
                break
            aapsTotalIn[p][aa] = aaCountIn[p][aa]/aaTotalIn[p][aa]*pamMatrix[p][aa]/sum(pamMatrix[p])
            aapsTotalOut[p][aa] = aaCountOut[p][aa]/aaTotalOut[p][aa]*pamMatrix[p][aa]/sum(pamMatrix[p])
            aaAjustedIn[p][aa] = (aaCountIn[p][aa]+aapsTotalIn[p][aa])/(aaTotalIn[p][aa]+psTotalIn[p][aa])
            aaAjustedOut[p][aa] = (aaCountOut[p][aa]+aapsTotalOut[p][aa])/(aaTotalOut[p][aa]+psTotalOut[p][aa])
    for aa in range(col):               #Finds Family Entropy
        for p in range(len(aaCount[0])):
            groupEntropy[aa] += (aaCountIn[aa][p]-aaCountOut[aa]) * math.log((aaCountIn[aa][p]/aaCountOut[aa]))     
def crval(group):
    if(len(group)<3):
        exit
    
    
    
def zscore():
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
    
    alignmentName  = lines[0].remove("Alignment File: ")
    groupsName = lines[1].remove("Groups File: ")
    a = open(alignmentName,'r')
    PSEUDOCOUNT_MULTIPLIER = int(lines[3].remove("Pseudocount Multiplier: "));
    r = AlignIO.read(alignmentName,"fasta")
    print(r)
    ewhole(r,)
    
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