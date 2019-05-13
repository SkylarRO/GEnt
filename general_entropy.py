#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:44:24 2018
@author: skylar
Usage: python GEnt.py name.settings
REQUIRES: lg.dat 
"""

import argparse
import math
import time
import os
import sys
from Bio import AlignIO               # Read fasta files


global PC_MUL
global REFRENCE
global GAPS
def ewhole(grp, r_prob, aln):
    """Calculates family entropy."""
    global GAPS
    col = len(aln[0].seq)
    row = len(aln)
    aa_p = [[0] * 20 for i in range(col)]
    aa_ajus = [[0] * 20 for i in range(col)]
    fam_ent = [0]*col

    for i in range(col):
        for j in range(20):
            if grp["gap"][i] > GAPS:
                break
            aa_p[i][j] = grp["ps"][i]*float(r_prob[j])
            aa_ajus[i][j] = row/(row+grp["ps"][i])*grp["aa_ct"][i][j]/row+grp["ps"][i]/(row+grp["ps"][i])*aa_p[i][j]/grp["ps"][i]
    for i in range(col):
        for j in range(20):
            if grp["gap"][i] > GAPS:
                break
            fam_ent[i] += aa_ajus[i][j]*math.log(aa_ajus[i][j]/float(r_prob[j]), 2)
    grp["aa_p_fam"] = fam_ent

def grpent(in_grp, out_grp):
    """Calculates Group Entropy"""
    global REFRENCE
    col = in_grp["col"]
    in_ent = [0]*col
    out_ent = [0]*col
    tot_ent = [0]*col

    for i in range(col):               #Finds Family Entropy
        for j in range(20):
            if in_grp["gap"][i] > GAPS or out_grp["gap"][i] > GAPS:
                break
            in_ent[i] += in_grp["aa_ajus"][i][j]*math.log(in_grp["aa_ajus"][i][j]/out_grp["aa_ajus"][i][j], 2)
            out_ent[i] += out_grp["aa_ajus"][i][j]*math.log(out_grp["aa_ajus"][i][j]/in_grp["aa_ajus"][i][j], 2)
        tot_ent[i] = in_ent[i] + out_ent[i]

    in_grp["in_ent"] =  in_ent
    in_grp["out_ent"] =  out_ent
    in_grp["tot_ent"] =  tot_ent

def ret_num(char_aa):
    """Converts char to a numerical representation of the Amino acid"""
    if char_aa == 'X':
        char_aa = 'A'
    temp = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    str(char_aa).capitalize()
    for i in range(20):
        if char_aa == temp[i]:
            aa_val = i
    return aa_val

def ret_aa(int_ascii):
    """Converts numerical representation back to a char."""
    temp = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    return temp[int_ascii]
def get_raw(aln, ret):
    """Passes the gap percentages and Amino Acid counts back by reference."""
    global GAPS
    rows = len(aln)
    col = ret["col"]
    gap = [0]*col
    temp = ""
    num_gaps = 0
    aa_ct = [[0]*20 for i in range(col)]
    for i in range(col):
        for j in range(rows):
            char = (aln[j].seq)[i]
            if char in ('.', '-'):
                num_gaps += 1
        gap[i] = num_gaps/rows
        num_gaps = 0

    for i in range(col):
        temp = ""
        for j in range(rows):
            char = (aln[j].seq)[i]

            if gap[i] <= GAPS:
                temp += char
                aa_ct[i][ret_num(char)] += 1
    ret["aa_ct"] = aa_ct
    ret["gap"] = gap
def get_consensus(alignments, ret):
    """Passes back consensus Amino Acid by reference."""
    global GAPS
    col = ret["col"]
    aa_ct = ret["aa_ct"]
    gap = ret["gap"]
    con = ['']*col
    for i in range(col):
        max_v = 0
        temp = ''
        for j in range(20):       #Finds Consensus AA
            if gap[i] > GAPS:
                con[i] = '.'
                break
            elif max_v <= aa_ct[i][j]:
                temp = ret_aa(j)
                max_v = aa_ct[i][j]
                con[i] = temp
    ret["con"] = con

def get_counts(alignments, ret):
    """Passes Pseudocount and amino acid totals by reference."""
    global PC_MUL
    global GAPS
    col = ret["col"]
    aa_ct = ret["aa_ct"]
    gap = ret["gap"]
    pseudocount = [0]*col
    aa_tot = [0]*col

    for i in range(col):
        acids = 0
        for j in range(20):
            if gap[i] > GAPS:
                break
            if aa_ct[i][j] > 0:
                acids += 1
        aa_tot[i] = acids
        pseudocount[i] = acids* PC_MUL
    ret["aa_tot"] = aa_tot
    ret["ps"] = pseudocount
def ajust_counts(sequences, odds, ret):
    """Passes back both Pseudocount totals and the Ajusted total"""
    global GAPS
    col = ret["col"]
    aa_pc = ret["ps"]
    aa_ct = ret["aa_ct"]
    gap = ret["gap"]
    aa_ps_tot = [[0] * 20 for i in range(col)]
    aa_ajus = [[0] * 20 for i in range(col)]

    for i in range(col):
        for j in range(20):#j is the positon
            if gap[i] > GAPS:
                break
            aa_ps_tot[i][j] = aa_pc[i]*float(odds[j]) #there is a much better forumla one might be able to use here.
            aa_ajus[i][j] = len(sequences)/(len(sequences)+aa_pc[i])*aa_ct[i][j]/len(sequences)+aa_pc[i]/(len(sequences)+aa_pc[i])*aa_ps_tot[i][j]/aa_pc[i]
    ret["aa_ps_tot"] = aa_ps_tot
    ret["aa_ajus"] = aa_ajus

def ref_position():
    """Finds position within reference sequence."""
    global REFRENCE
    col = len(REFRENCE.seq)
    out = [0]*col
    count = 0
    for i in range(col):
        if REFRENCE.seq[i] != '.' and REFRENCE.seq[i] != '-':
            count += 1
            out[i] = count
    return out

class Groups():
    """Handles data storage."""
    grp_dict = dict()
    dat_dict = dict()
    grp_data_list = dict()
    aln = None
    odds = None
    def __init__(self, alignments, delimiter, groups, odds):
        temp_dict = dict()
        for i in alignments:
            temp_dict[i.name] = i
        temp_int = 0
        global REFRENCE
        for line in groups:
            if line == delimiter + "\n":
                temp = "Group" + str(temp_int)
                temp_int += 1
                self.grp_dict[temp] = list()
            elif line[0] == ">":
                self.grp_dict[temp].append(temp_dict[line[1:len(line)-1]])
        REFRENCE = temp_dict[REFRENCE]
        self.aln = alignments
        self.odds = odds
    def proc_data(self):
        """Calls all calculation functions."""
        self.dat_dict["col"] = len(self.aln[0])
        get_raw(self.aln, self.dat_dict)
        get_consensus(self.aln, self.dat_dict)
        get_counts(self.aln, self.dat_dict)
        ewhole(self.dat_dict, self.odds, self.aln)
        
        for grp in list(self.grp_dict.keys()):
            in_grp = self.get_group(grp)[0]
            out_grp = self.get_group(grp)[1]
            self.grp_data_list[grp + "_in"] = dict()
            self.grp_data_list[grp + "_out"] = dict()
            self.grp_data_list[grp + "_in"]["col"] = len(in_grp[0].seq)
            self.grp_data_list[grp + "_out"]["col"] = len(out_grp[0].seq)
            get_raw(in_grp, self.grp_data_list[grp + "_in"])
            get_raw(out_grp, self.grp_data_list[grp + "_out"])
            get_counts(in_grp, self.grp_data_list[grp + "_in"])
            get_counts(out_grp, self.grp_data_list[grp + "_out"])
            get_consensus(in_grp, self.grp_data_list[grp + "_in"])
            ajust_counts(in_grp, self.odds, self.grp_data_list[grp + "_in"])
            ajust_counts(out_grp, self.odds, self.grp_data_list[grp + "_out"])
            grpent(self.grp_data_list[grp + "_in"], self.grp_data_list[grp + "_out"])

    def get_group(self, group_name):
        """Returns in_group and out_group."""
        in_group = self.grp_dict[group_name]
        out_group = []

        for k in self.grp_dict.keys():
            if k != group_name:
                out_group.extend(self.grp_dict[k])
        return [in_group, out_group]
    def get_all_grouped(self):
        """Returns all group sequences."""
        out_group = []
        for k in list(self.grp_dict.keys()):
            out_group += self.grp_dict[k]
        return out_group

def get_top(grp_ent, ref):
    """Returns top residues."""
    i = sorted(range(len(grp_ent)), key=lambda i: grp_ent[i])[0:]
    ret = []
    i = i[::-1]
    for j in range(10):
        ret.append(ref[i[j]])
    return ret

def main():
    global REFRENCE
    global GAPS
    start_time = time.time()
    parser = argparse.ArgumentParser(
        description='Calculate group entropy')
    parser.add_argument('Alignment', type=str,
                        help='Alignment file in fasta format')
    parser.add_argument('GrDel', type=str,
                        help='What character(s) delimit different groups')
    parser.add_argument('PsMultiplier', type=float,
                        help='Pseudocount Multiplier')
    parser.add_argument('RSequence', type=str,
                        help='Sequence to base positions on')
    parser.add_argument('--gapScore', type=float, nargs='?', const=0,
                        help='Whether or not to count sequences outside of groups')
    global PC_MUL
    REFRENCE = parser.parse_args().RSequence
    PC_MUL = float(parser.parse_args().PsMultiplier)
    matrix = [item.split() for item in open("lg.dat").read().split('\n')[:-1]][20]
    GAPS = parser.parse_args().gapScore
    if GAPS == None:
        GAPS = 0
    de = parser.parse_args().GrDel
    aln = open("aln.temp", "w")
    fasta = open(parser.parse_args().Alignment, 'r')
    file = ""
    for line in fasta:  
        if not de in line:
            file += line
    aln.write(file)
    aln = AlignIO.read("aln.temp", "fasta")
    fasta = open(parser.parse_args().Alignment, 'r')

    g = Groups(aln, de, fasta, matrix)
    g.proc_data()

    for groupN in list(g.grp_dict.keys()):
        out = open(groupN + ".csv", 'w')
        out.write("Position,Family Entropy,Group Entropy,Partial Group Entropy,Partial Out Group Entropy,Highest Group AA,Highest Family AA,Ref. Position\n")
        for i in range(g.dat_dict["col"]):
            if not g.dat_dict["gap"][i] > GAPS:
                out.write(str(i+1)+','+
                          str(g.dat_dict["aa_p_fam"][i]) + ',' +
                          str(g.grp_data_list[groupN + "_in"]["tot_ent"][i]) + ',' +
                          str(g.grp_data_list[groupN + "_in"]["in_ent"][i]) + ',' +
                          str(g.grp_data_list[groupN + "_in"]["out_ent"][i]) + ',' +
                          str(g.grp_data_list[groupN + "_in"]["con"][i]) +',' +
                          str(g.dat_dict["con"][i])+',' +
                          str(ref_position()[i])+"\n")
        
    sys.stdout.write("GEnt Complete in " + str(time.time()-start_time) + " seconds\n")
    out = open("GEnt.out", 'w')
    out.write("GEnt.py\n")
    out.write("Written by Skylar Olson\n")
    out.write("Based on An algorithm for identification and ranking of family-specific residues, applied to the ALDH3 family\n")
    out.write("Written by John Hempel, John Perozich, Troy Wymore, and Hugh B. Nicholas\n")
if __name__ == '__main__':
    main()
