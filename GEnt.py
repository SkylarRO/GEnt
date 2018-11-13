#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:44:24 2018

@author: skylar
"""

from Bio import AlignIO               # Read fasta files
import argparse   
import math
import time
import sys


global PSEUDOCOUNT_MULTIPLIER 
global ERRORS
def ewhole(alignments,pamMatrix,rProb):             #Calculates the family entropy
    col = len(alignments[0].seq)
    temp =getRaw(alignments) 
    gap = temp[1]
    aaCount = temp[0]
    consensus = getConsensus(alignments)
    temp = getCountPCount(alignments)
    pseudocount = temp[1]
    aaPseudocount = [[0] * 20 for i in range (col)] #Note the proper two dimentional list creation syntax
    aaAjusted = [[0] * 20 for i in range (col)]     #If you no not do it this way, you will have a bunch 
    famEntropy = [0]*col                            
    #of linked columns or nothing at all

    for aa in range(col):                           #Finds Ajusted values for count and individual pseudocount
        for p in range(20):
            if gap[aa]:                             #This prevents the program from calculating results for gapped positions
                break                               #Removing this will cause a division by zero error
            aaPseudocount[aa][p] = pseudocount[aa]*float(rProb[p])
            aaAjusted[aa][p] = len(alignments)/(len(alignments)+pseudocount[aa])*aaCount[aa][p]/len(alignments)+pseudocount[aa]/(len(alignments)+pseudocount[aa])*aaPseudocount[aa][p]/pseudocount[aa]
    for aa in range(col):               #Finds Family Entropy
        for p in range(20):
            if gap[aa]:
                break
            famEntropy[aa] += aaAjusted[aa][p]*math.log(aaAjusted[aa][p]/float(rProb[p]),2)
            
    return [gap,consensus,famEntropy]

def grpent(inGroup,outGroup,alignments,rProb):
    col = len(alignments[0].seq)
    gap = getRaw(alignments)[1]
    consensusIn = getConsensus(inGroup)
    consensusOut = getConsensus(outGroup)
    aaAjustedIn = ajustCounts(inGroup,rProb)[1]
    aaAjustedOut = ajustCounts(outGroup,rProb)[1]
    ingroupEntropy = [0]*col
    outgroupEntropy = [0]*col
    totalgroupEntropy = [0]*col

    
    for aa in range(col):               #Finds Family Entropy
        for p in range(20):
            if gap[aa]:
                break
            ingroupEntropy[aa] += aaAjustedIn[aa][p]*math.log(aaAjustedIn[aa][p]/aaAjustedOut[aa][p],2)     #These calculatipms are from the GENT paper
            outgroupEntropy[aa] += aaAjustedOut[aa][p]*math.log(aaAjustedOut[aa][p]/aaAjustedIn[aa][p],2)   #They were one equation, but the sample results
        totalgroupEntropy[aa] = ingroupEntropy[aa] + outgroupEntropy[aa]                                    #Split between ingroup and outgroup entropy

    return [outgroupEntropy,ingroupEntropy,totalgroupEntropy,consensusIn,consensusOut]  
    
#Helper Functions
def remNPC(char_AA):
    if(char_AA == 'X'):
        char_AA = 'A'
    temp = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    str(char_AA).capitalize()
    for i in range(20):
        if char_AA == temp[i]:
            aaVal = i
    return aaVal

def retAA(intASCII):
    temp = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    return temp[intASCII]
def getRaw(alignments):
    global ERRORS
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
        temp = ""
        for p in range (rows):  #j is the number of alignments
            char = (alignments[p].seq)[aa]

            if(gap[aa] != True):
                temp +=char
                if(temp == 'X'):
                    ERRORS += ''.join["ERROR: At position ", str(aa+1), " in sequence ", str(p+1), " an X was replaced with an A.\n"]
                    
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
            elif maxV <= aaCount[aa][p]:
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
def ajustCounts(sequences,odds):
    col = len(sequences[0].seq)
    aaPCount = getCountPCount(sequences)[1]
    aaCount = getRaw(sequences)[0]
    gap = getRaw(sequences)[1]
    aapsTotal = [[0] * 20 for i in range (col)]
    aaAjusted = [[0] * 20 for i in range (col)]

    for aa in range(col): 
        for p in range(20):#j is the positon
            if gap[aa]:
                break
            aapsTotal[aa][p] = aaPCount[aa]*float(odds[p]) #there is a much better forumla one might be able to use here.
            aaAjusted[aa][p] = len(sequences)/(len(sequences)+aaPCount[aa])*aaCount[aa][p]/len(sequences)+aaPCount[aa]/(len(sequences)+aaPCount[aa])*aapsTotal[aa][p]/aaPCount[aa]
    return [aapsTotal,aaAjusted]

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
    global ERRORS
    ERRORS = ""
    start_time = time.time()
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
    uUnGrouped = (lines[4][lines[4].rfind(":")+2:len(lines[4])-1])    
    freq_ran = matrix[20]
    matrix = matrix[0:19]
    matrix = matMake(matrix,MatrixName)
    groupsName = (lines[1][lines[1].rfind(":")+2:len(lines[1])-1])

    g = groups(r,open(groupsName,'r'))
    groupNames = list(g.grpDict.keys())
    if(uUnGrouped == "false"):
        fe = ewhole(g.getAllGrouped(),matrix,freq_ran)
    else:
        fe = ewhole(r,matrix,freq_ran)
    for groupN in groupNames:
        ge = grpent(g.getGroup(groupN)[0],g.getGroup(groupN)[1],g.getAllGrouped(),freq_ran)
        out = open(groupN + ".csv",'w')
        out.write("Position,Family Entropy,Group Entropy,Partial Group Entropy,Partial Out Group Entropy,Highest Group AA,Highest Family AA\n")
        for i in range(len(fe[0])):
            if fe[0][i] == False:
                out.write(str(i+1)+','+str(fe[2][i]) + ','+str(ge[2][i]) + ','+str(ge[1][i]) + ','+str(ge[0][i]) + ','+str(ge[3][i]) +','+str(fe[1][i])+"\n")   
    sys.stdout.write("GEnt Complete in " + str(time.time()-start_time) + " seconds\n")
    out = open("GEnt.out",'w')
    out.write("GEnt.py\n")
    out.write("Written by Skylar Olson\n")
    out.write("Based on An algorithm for identification and ranking of family-specific residues, applied to the ALDH3 family\n")
    out.write("Written by John Hempel, John Perozich, Troy Wymore, and Hugh B. Nicholas\n")
    out.writelines(ERRORS)
if __name__ == '__main__':
    main()
