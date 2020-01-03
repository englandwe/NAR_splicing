#!/usr/bin/python

#Usage: sjshapePC_3prime.py icSHAPE_file GTF_file threads output_file

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
    typelist=[]
    chrlist=[]
    tx_start=1
    for exon in exonlist:
        tx_stop=tx_start+(exon[2]-exon[1])
        tmprange=[exon[1],exon[2],tx_start,tx_stop]
        rangelist.append(tmprange)
        chrlist.append(exon[0])
        typelist.append(exon[5])
        strandlist.append(exon[3])
        tx_start=tx_stop+1
    return rangelist,typelist,chrlist,strandlist

def shapeByRegion(shapefile,chunkStart,chunkSize):
    shapeout=[]
    with open(shapefile) as infile:
        infile.seek(chunkStart)
        lines = infile.read(chunkSize).splitlines()
        for line in lines:
             tmpshape=line.strip().split('\t')
             txname=tmpshape[0]
             #translate the genome positions into tx positions using the magic of exons
             ranges,types,chrs,strands=genToTx(txname,ggtfdict)
             #these are the exon ranges, in transcript positions, and their strandedness
             if types[0] == 'protein_coding':
                 for i in range(0,len(ranges)-1):
                     shapesub=[]
                     if ranges[i][3]-ranges[i][2] > 70 and strands[i] == '+':
                         subrange=range(ranges[i][2],ranges[i][2]+70)
                         for pos in subrange: 
                             shapeval=tmpshape[3:][pos-1]
                             shapesub.append(shapeval)
                         shapeout.append([txname,i,chrs[0],ranges[i][1]]+shapesub)
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
