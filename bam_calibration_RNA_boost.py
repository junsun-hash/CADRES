#!/usr/bin/env python

import argparse,os,sys,time,logging,datetime

startTime = time.time()

parser = argparse.ArgumentParser(description='Bam calibration')
parser.add_argument('--bam',help='Input bam file')
parser.add_argument('--output',help='Path and prefix of the output file')
parser.add_argument('--genome',help='genome sequence in fasta format')
parser.add_argument('--KeepTemp', action='store_true', help='Keep tempory files. Disable by default.')

args = parser.parse_args()

bam=args.bam
label=args.output
genome=args.genome
keep=args.KeepTemp

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')
# mark and index bam
com1='picard ReorderSam INPUT='+bam+' OUTPUT='+label+'_reordered.bam S=true R='+genome
com2='picard AddOrReplaceReadGroups INPUT='+label+'_reordered.bam OUTPUT='+label+'_addrg.bam RGID='+label+' RGLB='+label+' RGPL=COMPLETE RGPU=lane1 RGSM='+label
com3='picard MarkDuplicates INPUT='+label+'_addrg.bam OUTPUT='+label+'_dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX=null METRICS_FILE='+label+'_metrics.txt'

# SplitNCigarReads
com4='gatk SplitNCigarReads -R '+genome+' -I '+label+'_dedup.bam -O '+label+'_split.bam'

logging.info('Running command 1: '+com1+'\n')
os.system(com1)

logging.info('Running command 2: '+com2+'\n')
os.system(com2)
if (not keep):
    os.system('rm -f '+label+'_reordered.bam')

logging.info('Running command 3: '+com3+'\n')
os.system(com3)
if (not keep):
    os.system('rm -f '+label+'_addrg.bam')

logging.info('Running command 4: '+com4+'\n')
os.system(com4)
if (not keep):
    os.system('rm -f '+label+'_dedup.bam')

logging.info("Program ended")
currentTime = time.time()
runningTime = currentTime-startTime
logging.info("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))

