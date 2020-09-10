#! /usr/bin/python3.4
from sys import stdout, stdin, argv
import string
import os
import re
import tempfile
import shutil
import subprocess


PATH = argv[1]
BAM1 = argv[2]
VCF = argv[3]
RSEED = int(argv[4])
OUT = argv[5]
TARGET_COV = 30.0

try:
    with tempfile.TemporaryDirectory() as tpath:
        print( tpath )
        print( "%1s/read_haps %1s/merged.bam %1s/high_quality_markers_deCODE_2015.txt %1s > %1s"%(PATH,tpath,PATH,VCF,OUT))
        c1 = float( BAM1.split(".")[2].split("x")[0] )
        os.popen( "samtools view -b -s %1s %1s > %1s/f1.bam"%(RSEED+TARGET_COV/c1,BAM1,tpath )).read()
        os.popen( "samtools index %1s/f1.bam"%(tpath )).read()
        os.popen( "ln -s %1s/genome.fa %1s/genome.fa"%(PATH,tpath)).read()
        os.popen( "ln -s %1s/genome.fa.fai %1s/genome.fa.fai"%(PATH,tpath)).read()
        os.popen( "%1s/read_haps %1s/f1.bam %1s/high_quality_markers_deCODE_2015.txt %1s > %1s_read_haps"%(PATH,tpath,PATH,VCF,OUT)).read()
        os.popen( "%1s/verifyBamID --ignoreRG --chip-none --free-full --maxDepth 100 --precise --vcf %1s/1000G_EUR.vcf --bam %1s/f1.bam --out %1s_verifyBamID"%(PATH,PATH,tpath,OUT)).read()

except:
    shutil.rmtree(tpath, ignore_errors=True)
