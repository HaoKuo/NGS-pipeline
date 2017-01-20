#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

import os, sys, time, re, glob, itertools, json
import subprocess, multiprocessing, shlex

work_dir = "/work-z/user/guoh/research/2016-12-28-PedigreeCancerWGS/"
data_dir ="/work-z/user/guoh/research/2016-12-28-PedigreeCancerWGS/RectumCancerPedigreeWGS/NHT160089/Result0513/QC/"
data = data_dir+"*/*.gz"
fq_file_dic = {}
for d in glob.glob(data): 
    fq_name = os.path.basename(d)
    tmp = fq_name.split(".")[0].split("_")
    lane = tmp[4]
    sample = tmp[0]+"-"+tmp[1]
    key = lane+"_"+sample
    fq_index = "R"+str(tmp[-1])
    fq_file_dic.setdefault(key,{}).setdefault(fq_index,d)

f = open("sample_fastqinfo.tsv","wb")
for k in fq_file_dic:
    lane = k.split("_")[0]
    sample = k.split("_")[1]
    tmp_dic = fq_file_dic[k]
    f.write('\t'.join([lane,sample,tmp_dic['R1'],tmp_dic['R2']])+'\n')
f.close()

GATK_json = work_dir + "GATK_workflow_pedigree.json"
jf = file(GATK_json)
data = json.load(jf)
jf.close()

#second_bams = glob.glob(data_dir+"*.bam")
#tmpl = []
#for f in second_bams:
#    name = os.path.basename(f).split(".")[0]
#    bamPath = os.path.abspath(f)
#    baiPath = bamPath+".bai"
#    tmpl.append([name,bamPath,baiPath])

data["GATK_test.sampleMetainfo"]= work_dir + "sample_fastqinfo.tsv"
tmp_dumps = json.dumps(data)
f = open("GATK_workflow_pedigree.json","wb")
f.write(tmp_dumps)
f.close()







