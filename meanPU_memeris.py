#!/usr/bin/env python


import sys
import glob
from scipy import stats as ss

filelist=glob.glob('*/*.pu')

#intron_PC_mouse_invivo_18 11_4177548_4177632_3_11;0.5950;0.5950;0.5311;0.5059;0.5059;0.4897;0.0999;0.0396;0.1253;0.0129;0.0066;0.0069;0.0064;0.0040;0.0038;0.0020;0.0019;0.0179;0.0185;0.0024;0.0002;0.0002;0.0004;0.0004;0.0071;0.0042;0.0009;0.0052;0.0604;0.0150;0.0123;0.0179;0.0152;0.0003;0.0000;0.0001;0.0001;0.0001;0.0001;0.0034;0.2731;0.2687;0.4665;0.6558;0.6453;0.6558;0.6349;0.6666;0.7347;0.7347;0.3966;0.3839;0.2821;0.2438;0.1820;0.1820;0.1820;0.1791;0.0027;0.0017;0.0016;0.0018;0.0085;0.0086;0.0117;0.0097;0.0114;0.0002;0.0003;0.0021;0.0021;0.0000;0.0002;0.0002;0.0028;0.0028;0.0028;0.0030;0.0045;

pudict={}
for file in filelist:
    with open(file) as infile:
        for line in infile:
            tmpline=line.strip('[\n;]').split(';')
            seqid=tmpline[0].split()[0]
            #seqnum=seqid.split('_')[-1]
            #sourcefile='_'.join(seqid.split('_')[:-1])
            flankleft,flankright=[int(x) for x in tmpline[0].split()[1].split('_')[-2:]]
            #we only want the 65 hexamers covering out 70bp of interest, so adjust for flanking regions
            #first hexamer will start at flankleft
            hexcount=0
            for idx,val in enumerate(tmpline[1:]):
                if idx >= flankleft:
                    try:
                        pudict[seqid][hexcount].append(float(val))
                    except KeyError:
                        pudict[seqid]={}
                        for i in range(0,65):
                            pudict[seqid][i]=[]
                        pudict[seqid][hexcount].append(float(val))
                    except ValueError:
                        #catch any remaining '-' out there
                        sys.stderr.write('theres a dash here: %s\n' % tmpline)
                    hexcount+=1
                    if hexcount >= 65:
                        break
metadict={}
with open('meanpu.byseq.out','w') as outfile:
    for id,posdict in pudict.items():
        sourcefile='_'.join(id.split('_')[:-1])
        for pos,puvals in posdict.items():
            try:
                meanpu=sum(puvals)/float(len(puvals))
                sempu=ss.sem(puvals)
            except ZeroDivisionError:
                sys.stderr.write('%s,%s,%s\n' % (id,pos,puvals))
                continue 
            try:
                metadict[sourcefile][pos].append(meanpu)
            except KeyError:
                metadict[sourcefile]={}
                for i in range(0,65):
                     metadict[sourcefile][i]=[]
                     metadict[sourcefile][pos].append(meanpu)
            outfile.write('%s\t%s\t%s\t%s\n' % (id,pos,meanpu,sempu))

with open('meanpu.bypos.out','w') as outfile:    
    for sourcefile,poslist in metadict.items():
        for pos,vals in poslist.items():
            meanpos=sum(vals)/float(len(vals))
            sempos=ss.sem(vals)
            outfile.write('%s\t%s\t%s\t%s\n' % (sourcefile,pos,meanpos,sempos))
                
    
