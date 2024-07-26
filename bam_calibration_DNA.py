#!/usr/bin/env python
'''
2023.12.7 Grammar modified according to gatk4.
Removed the IndelRealigner step(com5&com6) because gatk4 removed this feature.
'''
import argparse,os,sys,time,logging,datetime

startTime = time.time()

parser = argparse.ArgumentParser(description='Bam calibration for rMATS-DVR')
parser.add_argument('--bam',help='Input bam file')
parser.add_argument('--output',help='Path and prefix of the output file')
parser.add_argument('--genome',help='genome sequence in fasta format')
parser.add_argument('--known', default='NA', help='Known SNVs in vcf format')
parser.add_argument('--KeepTemp', action='store_true', help='Keep tempory files. Disable by default.')

args = parser.parse_args()

bam=args.bam
label=args.output
genome=args.genome
known=args.known
keep=args.KeepTemp

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=label+'.bamCalibration.log'+ str(datetime.datetime.now())+'.txt' ,
                    filemode='w')

# mark and index bam
com1='picard ReorderSam INPUT='+bam+' OUTPUT='+label+'_reordered.bam S=true R='+genome
com2='picard AddOrReplaceReadGroups INPUT='+label+'_reordered.bam OUTPUT='+label+'_addrg.bam RGID='+label+' RGLB='+label+' RGPL=COMPLETE RGPU=lane1 RGSM='+label
com3='picard MarkDuplicates INPUT='+label+'_addrg.bam OUTPUT='+label+'_dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX=null METRICS_FILE='+label+'_metrics.txt'

# SplitNCigarReads
#com5='java -Xmx4g -jar '+directory+'GenomeAnalysisTK.jar -T RealignerTargetCreator -R '+genome+' -I '+label+'_dedup.bam -o '+label+'_realigner.intervals -known /data/scratch/DCT/UCCT/STMOLPHA/czhang/genomes/indel/sorted_final.vcf'
#com6='java -Xmx4g -jar '+directory+'GenomeAnalysisTK.jar -T IndelRealigner -R '+genome+' -I '+label+'_dedup.bam -o '+label+'_split.bam -targetIntervals '+label+'_realigner.intervals -known /data/scratch/DCT/UCCT/STMOLPHA/czhang/genomes/indel/sorted_final.vcf'

# Base Quality Score Recalibration
if (known!='NA'):
    com7='gatk BaseRecalibrator -I '+label+'_dedup.bam -R '+genome+' -O '+label+'_recalibration_report.grp -known-sites '+known
else:
    #com7='java -Xmx4g -jar '+directory+'GenomeAnalysisTK.jar -T BaseRecalibrator -I '+label+'_split.bam -R '+genome+' -o '+label+'_recalibration_report.grp --run_without_dbsnp_potentially_ruining_quality'
    print("NEED KNOWN SITES")
    sys.exit()

com8='gatk ApplyBQSR -R '+genome+' -I '+label+'_dedup.bam -bqsr-recal-file '+label+'_recalibration_report.grp -O '+label+'_recalibration.bam'

logging.debug('Running command 1: '+com1+'\n')
os.system(com1)
logging.debug('Running command 2: '+com2+'\n')
os.system(com2)
if (not keep):
    os.system('rm -f '+label+'_reordered.bam')
logging.debug('Running command 3: '+com3+'\n')
os.system(com3)
if (not keep):
    os.system('rm -f '+label+'_addrg.bam')
#logging.debug('Running command 5: '+com5+'\n')
#os.system(com5)
#logging.debug('Running command 6: '+com5+'\n')
#os.system(com6)


logging.debug('Running command 7: '+com7+'\n')
os.system(com7)
logging.debug('Running command 8: '+com8+'\n')
os.system(com8)
if (not keep):
    #os.system('rm -f '+label+'_split.bam')
    os.system('rm -f '+label+'_split.bai')
    os.system('rm -f '+label+'_dedup.bam')    
logging.debug("Program ended")
currentTime = time.time()
runningTime = currentTime-startTime
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))

