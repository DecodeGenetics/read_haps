#! /usr/bin/python3.4
from sys import stdout, stdin, argv
import string
import os
import re
import statistics
import shutil
import subprocess


for pref in ["raw","res_1","res_2","res_3","res_4","res_5","res_10","res_20"]:
    for suff in ["_read_haps","_verifyBamID.selfSM"]:
        col = 3
        thresh = 0.002
        if suff == "_verifyBamID.selfSM":
            col = 6
            thresh = 0.02
        iter = 100
        if pref == "raw":
            iter = 700
        v = []
        nFail = 0
        for i in range( iter ):
            c = float( open("%1s_%1s%s"%(pref,i+1,suff)).readlines()[1].split()[col])
            v.append(c)
            if c > thresh:
                nFail = nFail+1
        stdout.write("%0.5f & %0.5f & %0.5f & %2.0f"%(statistics.mean(v),statistics.stdev(v),statistics.median(v),float(nFail)/float(iter)*100))
        if suff == "_read_haps":
            stdout.write( " & ")
        else:
            stdout.write( " \\\\ \hline \n")
        
            
    
    
