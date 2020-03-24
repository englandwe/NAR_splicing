#!/usr/bin/env python

import sys
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
                downflank=(txentry.end+1)-endpoints70[1]
            range_seq=fasta_dict[chrom].seq[flanked_endpoints[0]-1:flanked_endpoints[1]-1]
            rec=SeqRecord(range_seq,id=instr+'_'+str(cntr),description='_'.join([str(x) for x in [chrom,flanked_endpoints[0],flanked_endpoints[1],upflank,downflank]]))
            reclist.append(rec)
            if cntr % 10 == 0:
                sys.stderr.write('%s records processed...\n' % cntr)

#save a fasta for future use
with open(instr+'.fa', 'w') as ofile:
    for record in reclist:
        sequence = str(record.seq)
        ofile.write('>'+record.id+' '+record.description+'\n'+sequence+'\n')

big_outs=[]
bypos_dict={}
for record in reclist:
    sequence = str(record.seq)
    plstring='>'+record.id+' '+record.description+'\n'+sequence+'\n'
    #localfold won't do stdin
    #plin=plstring.encode('utf-8')
    tmpname=record.id+'_localfoldtmpfile.fa'
    with open(tmpname,'w') as tmpfile:
        tmpfile.write(plstring)

    lunp_in=record.id+'_lunp'
    #psfile=record.id+'_dp.ps'
    psfile=record.id+'_localfoldtmpfile_W120_L70_skip10.ps'
    flanks=[int(x) for x in record.description.split('_')[3:]]
    #subprocess.run(['RNAplfold', '-u6'],input=plin,check=True)
    subprocess.run(['perl', 'LocalFold-1.0/localfold_fix.pl', '-seqfile='+tmpname, '-u', '6', '-L', '70', '-W', '120'],check=True)
    #retrieve data from lunp file
    #RNAplfold counts from the right side of the hexamer
    #with open(lunp_in) as infile:
    #    next(infile)
    #    next(infile)
    #no lunp for localfold; use .acc instead
    #coding_mouse_invivo_3prime_1_localfoldtmpfile_W100_L70_skip10.acc
    #with open(lunp_in) as infile:
    acc_in=record.id+'_localfoldtmpfile_W120_L70_skip10.acc'
    with open(acc_in) as infile:
        for line in infile:
            tmpline=line.strip().split('\t')
            if int(tmpline[0]) > flanks[0] + 5 and int(tmpline[0]) <= flanks[0] + 70:
                #previous output started from 0, so let's line up with that
                hexpos=str(int(tmpline[0])-flanks[0]-6)
                big_outs.append(record.id+'\t'+hexpos+'\t'+tmpline[6]+'\n')
                #also arrange them by sample and position so we can grab mean and sem
                #shave off the id number
                sampbase='_'.join(record.id.split('_')[0:-1])
                try:
                    bypos_dict[sampbase][hexpos].append(float(tmpline[6]))
                except KeyError:
                    try:
                        bypos_dict[sampbase][hexpos]=[]
                        bypos_dict[sampbase][hexpos].append(float(tmpline[6]))
                    except KeyError:
                        bypos_dict[sampbase]={}
                        bypos_dict[sampbase][hexpos]=[]
                        bypos_dict[sampbase][hexpos].append(float(tmpline[6]))
                     
    #clean up
    subprocess.run(['rm', psfile, acc_in, tmpname],check=True)
    
#and we're done
with open(instr+'_plout.tsv','w') as outfile:
    for entry in big_outs:
        outfile.write(entry)            

#well almost            
with open(instr+'_meanout.tsv','w') as outfile:
    for sample,positions in bypos_dict.items():
        for pos in positions:
            posmean=sum(bypos_dict[sample][pos])/len(bypos_dict[sample][pos])
            possem=ss.sem(bypos_dict[sample][pos])
            outfile.write('\t'.join([str(x) for x in [sample,pos,posmean,possem]])+'\n')
