#!/usr/bin/python

#Usage: intronHunterPC.py icSHAPE_file GTF_file threads output_file

import sys
import os
from multiprocessing import Pool

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

def chunkify(fname,size):
    fileEnd = os.path.getsize(fname)
    with open(fname,'r') as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size,1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break

def child_initialize(_inData):
     global ggtfdict
     ggtfdict = _inData

def genToTx(txname,gtfdict):
    exonlist=[[x.seqname,x.start,x.end,x.strand,int(x.attdict['exon_number']),x.attdict['transcript_biotype']] for x in gtfdict[txname] if x.feature == 'exon']
    exonlist.sort(key=lambda y : y[4])
    rangelist=[]
    strandlist=[]
    chrlist=[]
    tx_start=1
    for exon in exonlist:
        tx_stop=tx_start+(exon[2]-exon[1])
        tmprange=[exon[1],exon[2],tx_start,tx_stop]
        rangelist.append(tmprange)
        strandlist.append(exon[3])
        chrlist.append(exon[0])
        tx_start=tx_stop+1
    return rangelist,strandlist,chrlist

def shapeByRegion(shapefile,chunkStart,chunkSize):
    shapeout=[]
    with open(shapefile) as infile:
        infile.seek(chunkStart)
        lines = infile.read(chunkSize).splitlines()
        for line in lines:
             tmpshape=line.strip().split('\t')
             txname=tmpshape[0]
             bt=[[x.attdict['transcript_biotype'],x.attdict['gene_id']] for x in ggtfdict[txname] if x.feature == 'transcript'][0]
             if bt[0] == 'retained_intron':
                 ret_ranges,ret_strand,ret_chr=genToTx(txname,ggtfdict)
                 if ret_strand[0] == '+':
                     coding_list=[tx for tx in ggtfdict.keys() for entry in ggtfdict[tx] if entry.feature == 'transcript' and entry.attdict['transcript_biotype'] == 'protein_coding' and entry.attdict['gene_id'] == bt[1]]
                     for coder in coding_list:
                         exon_ranges,exon_strands,exon_chr=genToTx(coder,ggtfdict)
                         for possible_ret in ret_ranges:
                             startex=[]
                             stopex=[]
                             for exrange in exon_ranges:
                                 if int(possible_ret[0]) in range(int(exrange[0]-5),int(exrange[1])+6):
                                     startex=exrange
                                 if int(possible_ret[1]) in range(int(exrange[0]-5),int(exrange[1])+6):
                                     stopex=exrange
                             if startex != stopex and len(startex) > 0 and len(stopex) > 0:
                                 #it spanned two exons
                                 #grab genomic endpoint of startex, find that in the retained intron
                                 startex_end=startex[1]
                                 #find distance between startex and genomic start of suspicious exon
                                 ret_dist=int(startex_end)-int(possible_ret[0])
                                 #now add that to the suspicious exon's transcript start position
                                 tx_junct=int(possible_ret[2])+ret_dist
                                 #are there at least 70bp before the previous junction?
                                 if ret_dist > 70:
                                     #then let's do this
                                     subrange=range(tx_junct-69,tx_junct+1)
                                     shapesub=[]
                                     pcsub=[]
				     for pos in subrange: 
                                         shapeval=tmpshape[3:][pos-1]
                                         shapesub.append(shapeval)
                                     shapeout.append([txname,coder,ret_chr[0],str(startex_end)]+shapesub)
    return shapeout


def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final


gtfdict={}
with open(sys.argv[2]) as infile:
    for line in infile:
        if not line.startswith('#'):
            rectmp=GtfRec(line.strip().split('\t'))
            if 'transcript_id' in rectmp.attdict:
                transcript=rectmp.attdict['transcript_id']
                try:
                    gtfdict[transcript].append(rectmp)
                except KeyError:
                    gtfdict[transcript]=[rectmp]

shapefile=sys.argv[1]
cores=int(sys.argv[3])
outfile=sys.argv[4]

shapepool = Pool(cores, initializer = child_initialize, initargs = (gtfdict,))
shapejobs = []
shapejobs=[shapepool.apply_async(shapeByRegion, (shapefile,chunkStart,chunkSize)) for chunkStart,chunkSize in chunkify(shapefile,10*1024*1024)]
shapepool.close()
shapepool.join()
shapeoutput=[shapejob.get() for shapejob in shapejobs]
shape_by_reg=[item for sublist in shapeoutput for item in sublist]
regout=open(outfile,'w')
regout.write(flattenList(shape_by_reg))
regout.close()
