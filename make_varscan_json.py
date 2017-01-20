#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

import os, sys, time, re, glob, itertools, json
import subprocess, multiprocessing, shlex

work_dir ="/work-a/user/zhanggh/Work/06_China1/GATK/"
data_dir = "/work-a/user/zhanggh/Work/06_China1/GATK/data/"

varscan_json = work_dir + "varscan3.json"
jf = file(varscan_json)
data = json.load(jf)
jf.close()
#print data["GATK_test.FinalBams"]

second_bams = glob.glob(data_dir+"*.bam")
#print second_bams
tmpl = []
for f in second_bams:
    name = os.path.basename(f).split(".")[0]
    bamPath = os.path.abspath(f)
    baiPath = bamPath+".bai"
    tmpl.append([name,bamPath])

data["mpileup_varscan.RawBams"]=tmpl
tmp_dumps = json.dumps(data)
f = open("varscan_modified.json","wb")
f.write(tmp_dumps)
f.close()







