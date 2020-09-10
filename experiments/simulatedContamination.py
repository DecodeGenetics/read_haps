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
BAM2 = argv[3]
VCF = argv[4]
PERC = float( argv[5] )
RSEED = int(argv[6])
OUT = argv[7]
TARGET_COV = 30.0

try:
    with tempfile.TemporaryDirectory() as tpath:
        print( tpath )
        print( "%1s/read_haps %1s/merged.bam %1s/high_quality_markers_deCODE_2015.txt %1s > %1s"%(PATH,tpath,PATH,VCF,OUT))
        c1 = float( BAM1.split(".")[2].split("x")[0] )
        os.popen( "samtools view -b -s %1s %1s > %1s/f1.bam"%(RSEED+(1.0-PERC)*TARGET_COV/c1,BAM1,tpath )).read()
        c2 = float( BAM2.split(".")[2].split("x")[0] )
        os.popen( "samtools view -b -s %1s %1s > %1s/f2.bam"%(RSEED+PERC*TARGET_COV/c2,BAM2,tpath )).read()
        os.popen( "samtools merge -O bam %1s/merged.bam %1s/f1.bam %1s/f2.bam "%(tpath,tpath,tpath)).read()
        os.popen( "samtools index %1s/merged.bam"%(tpath)).read()
        os.popen( "ln -s %1s/genome.fa %1s/genome.fa"%(PATH,tpath)).read()
        os.popen( "ln -s %1s/genome.fa.fai %1s/genome.fa.fai"%(PATH,tpath)).read()
        os.popen( "%1s/read_haps %1s/merged.bam %1s/high_quality_markers_deCODE_2015.txt %1s > %1s_read_haps"%(PATH,tpath,PATH,VCF,OUT)).read()
        os.popen( "%1s/verifyBamID --ignoreRG --chip-none --free-full --maxDepth 100 --precise --vcf %1s/1000G_EUR.vcf --bam %1s/merged.bam --out %1s_verifyBamID"%(PATH,PATH,tpath,OUT)).read()

except:
    shutil.rmtree(tpath, ignore_errors=True)
