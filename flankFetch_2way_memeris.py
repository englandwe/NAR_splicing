#!/usr/bin/env python

import sys
import os
import subprocess
import argparse
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy import stats as ss

class GtfRec(object):
    def __init__(self,reclist):
         self.seqname=reclist[0]
         self.source=reclist[1]
         self.feature=reclist[2]
         self.start=int(reclist[3])
         self.end=int(reclist[4])
         self.score=reclist[5]
         self.strand=reclist[6]
         self.frame=reclist[7]
         self.attdict={}
         for self.item in reclist[8].strip(';').split('; '):
             self.splitline=self.item.replace('\"','').split(' ')
             self.attdict[self.splitline[0]]=self.splitline[1]

rangefile=sys.argv[1]
fastafile=sys.argv[2]
gtffile=sys.argv[3]
instr=sys.argv[4]
flanklen=int(sys.argv[5])
direction=sys.argv[6]

sys.stderr.write('importing fasta @ %s' % (str(datetime.now())) + '\n')

fasta_dict=SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))

sys.stderr.write('importing gtf @ %s' % (str(datetime.now())) + '\n')

gtfdict={}
with open(gtffile) as infile:
    for line in infile:
        if not line.startswith('#'):
            rectmp=GtfRec(line.strip().split('\t'))
            if rectmp.feature == 'transcript':
                transcript=rectmp.attdict['transcript_id']
                gtfdict[transcript]=rectmp

reclist=[]
cntr=0
with open(rangefile) as infile:
    for line in infile:
        tmpline=line.strip().split()
        if tmpline[4] == 'Pos70':
            cntr+=1
            intron_start=int(tmpline[3])
            tx=tmpline[0]
            chrom=tmpline[2]
            #for 5', intron_start is last base of the upstream exon
            #for 3', it's the first base of the downstream exon 
            if direction == '5':
                endpoints70=[intron_start-69,intron_start+1]
                flanked_endpoints=[endpoints70[0]-flanklen,endpoints70[1]+flanklen]
            elif direction == '3':
                endpoints70=[intron_start,intron_start+70]
                flanked_endpoints=[endpoints70[0]-flanklen,endpoints70[1]+flanklen]
            else:
                sys.stderr.write('invalid direction: %s.  Direction must be "5" or "3".\n' % direction)
            upflank=flanklen
            downflank=flanklen
            #adjust if we're gonna run off the end of the transcript
            txentry=gtfdict[tx]
            if flanked_endpoints[0] < txentry.start:
                flanked_endpoints[0]=txentry.start
                upflank=endpoints70[0]-txentry.start
            #gtfs are right inclusive, so compensate
            elif flanked_endpoints[1] > txentry.end+1:
                flanked_endpoints[1]=txentry.end+1
                downflank=endpoints70[1]-(txentry.end+1)
            range_seq=fasta_dict[chrom].seq[flanked_endpoints[0]-1:flanked_endpoints[1]-1]
            rec=SeqRecord(range_seq,id=instr+'_'+str(cntr),description='_'.join([str(x) for x in [chrom,flanked_endpoints[0],flanked_endpoints[1],upflank,downflank]]))
            reclist.append(rec)
            if cntr % 10 == 0:
                sys.stderr.write('%s records processed...\n' % cntr)

#save a fasta for readin by memeris
#each needs to be in its own subdir, to avoid future tempdir name collisions
drname=instr+'.'+str(flanklen)
os.mkdir(drname)
with open(drname+'/'+drname+'.fa', 'w') as ofile:
    for record in reclist:
        sequence = str(record.seq)
        ofile.write('>'+record.id+' '+record.description+'\n'+sequence+'\n')

#end of step 1
