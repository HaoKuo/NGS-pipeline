# test VarScan workflow
#
# TASK DEFINITIONS
#
## Align FASTQ reads with BWA-MEM
task BwaMem {
    File BWA
    File input_FASTQ1
    File input_FASTQ2
    String file_basename
    File ref_fasta
    File ref_dict
    File ref_index
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    String THREAD_NUM
    File SAMTOOLS
    String LANE

    command {
        ${BWA} mem -t ${THREAD_NUM} -M ${ref_fasta} ${input_FASTQ1} ${input_FASTQ2} -R '@RG\tID:${file_basename}_${LANE}-\tLB:${file_basename}_${LANE}\tSM:${file_basename}\tPL:ILLUMINA\tPU:flowcell-barcode.lane' | ${SAMTOOLS} view -Sb - > ${LANE}_${file_basename}.bam
        }
    output {
        File aligned_bam = "${LANE}_${file_basename}.bam"
        }
    }

task SortBam {
    File bamFile
    File PICARD
    File ref_fasta
    File ref_dict
    File ref_index
    String sampleName
    command {
        java -Xmx8g -XX:ParallelGCThreads=10 -Djava.io.tmpdir=/work-z/user/guoh/tmp -jar ${PICARD} SortSam \
            INPUT=${bamFile} \
            OUTPUT=${sampleName}.sorted.bam \
            CREATE_INDEX=true \
            SORT_ORDER=coordinate
    }
    output {
        File sorted_Bam = "${sampleName}.sorted.bam"
        File sorted_Bam_index = "${sampleName}.sorted.bai"
        }
    }

task BamFilter {
    String BamFile
    File SAMTOOLS
    File ref_fasta
    File ref_dict
    File ref_index
    String sampleName
    command {
        ${SAMTOOLS} view -F 0x900 ${BamFile} -b -o ${sampleName}.filtered.sorted.bam
                }
    output {
        File bam_filtered = "${sampleName}.filtered.sorted.bam"
        }
    }


task BamFlagstat {
    String BamFile
    File SAMTOOLS
    File ref_fasta
    File ref_dict
    File ref_index
    command {
        ${SAMTOOLS} flagstat ${BamFile} > ${BamFile}.flagstat.txt
                }
    output {
        File bam_index = "${BamFile}.flagstat.txt"
        }
    }


task mpileup_all {
    String BamFile
    File SAMTOOLS
    File ref_fasta
    File ref_dict
    File ref_index
    String sampleName
    command {
        ${SAMTOOLS} mpileup -B -f ${ref_fasta} -q 1 -d 10000000 ${BamFile} | gzip -9 > ${sampleName}.mpileup.gz
                }
    output {
        File mpileup_all_gz = "${sampleName}.mpileup.gz"
        }
    }

task mpileup_target {
    String BamFile
    File SAMTOOLS
    File ref_fasta
    File ref_dict
    File ref_index
    String sampleName
    String panel_bed_path
    command {
        ${SAMTOOLS} mpileup \
        -B -f ${ref_fasta} -q 1 \
        -l ${panel_bed_path} \
        -d 10000000 ${BamFile} | gzip -9 > ${sampleName}.mpileup.gz | parallel
    }
    output {
        File mpileup_target_gz = "${sampleName}.mpileup.gz"

    }

}

task Merge_mpileup_info {
    Array[File] files
    String panel_bed

    command <<<
	python <<CODE
	import os
	mpileup_files="${write_lines(files)}"
	sampleMpileupMap = {}
	with open(mpileup_files, 'r') as vf:
    	    for line in vf:
                tmpFile = line.strip()
                tmpFileDic =  os.path.split(tmpFile)[0]
                tmpFileName = os.path.split(tmpFile)[1]
                tmpName = tmpFileName.split('.')[0]
                sample = tmpName
                sampleMpileupMap.setdefault(sample,[]).append("${panel_bed}")
                sampleMpileupMap.setdefault(sample,[]).append(tmpFile)
	for key in sampleMpileupMap:
    	    fileListStr = '\t'.join(sampleMpileupMap[key])
    	    print (key+'\t'+fileListStr)
	CODE
	>>>
    output {
        File mpileupfileStr = stdout()
        }

    }

task testPrint {
    Array[String] s

    command <<<
        python <<CODE
        tmp = "${write_lines(s)}"
        f= open(tmp,"r")
        for line in f:
            print line
        f.close()
        CODE
        >>>
    output{
        String printScreen = read_string(stdout())
        }
    }


task CovDep{
    File mpileup_info
    String CovDep
    String panel_dic
    command {
        ${CovDep} <(cat ${mpileup_info}) ${panel_dic} > Allsamples.covdep.summary.txt
                }
    output {
        File covdep_summary = "Allsamples.covdep.summary.txt"
        }
    }

task VarscanMpileup2snp {
    File VarScan
    File mpileup
    String snp_output
    command {
        java -Xmx8g -XX:ParallelGCThreads=10 -jar ${VarScan} mpileup2snp \
        <(zcat ${mpileup}) \
        --min-coverage 4 \
        --min-var-freq 0.01 \
        --p-value 0.01 \
        --strand-filter 1 \
        --output-vcf \
        --variants > ${snp_output}.snp.vcf 
    }
    output {
        File snp_vcf = "${snp_output}.snp.vcf"
    }
 }

task VarscanMpileup2indel {
    File VarScan
    File mpileup
    String indel_output
    command {
        java -Xmx8g -XX:ParallelGCThreads=10 -jar ${VarScan} mpileup2indel \
        <(zcat ${mpileup}) \
        --min-coverage 4 \
        --min-var-freq 0.01 \
        --p-value 0.01 \
        --strand-filter 1 \
        --output-vcf \
        --variants > ${indel_output}.indel.vcf
    }
output {
        File indel_vcf = "${indel_output}.indel.vcf"
    }
}

task VarscanSomatic {
    File VarScan
    File normal_mpileup
    File tumor_mpileup
    String somatic_output
    command {
        java -Xmx8g -XX:ParallelGCThreads=10 -jar ${VarScan} somatic \
        <(zcat ${normal_mpileup}) \
        <(zcat ${tumor_mpileup}) \
        ${somatic_output} \
        --min-coverage 10 \
        --strand-filter 1 \
        --output-vcf 1 \
        --somatic-p-value 0.05
    }
    output {
        String somatic_calling = "${somatic_output}"
    }
}



task  MarkDuplicates{
    File bamFile
    File PICARD
    String sampleName
    command {
        java -Xmx8g -XX:ParallelGCThreads=10 -Djava.io.tmpdir=/work-z/user/guoh/tmp -jar ${PICARD} MarkDuplicates \
            INPUT=${bamFile} \
            OUTPUT=${sampleName}.dedup.sorted.bam \
            METRICS_FILE=${sampleName}.dedup.sorted.bam.metrics.txt
    }
    output {
        File deduped_Bam = "${sampleName}.dedup.sorted.bam"
        File dedup_metrics = "${sampleName}.dedup.sorted.bam.metrics.txt"
        }
    }

task BuildBamIndex {
    String bamFile
    File SAMTOOLS
    File ref_fasta
    File ref_dict
    File ref_index
    command {
        ${SAMTOOLS} index ${bamFile}
                }
    output {
        File bam_index = "${bamFile}.bai"
        }
    }


workflow mpileup_varscan {
    File sampleMetainfo
    Array[Array[String]] input_FASTQ = read_tsv(sampleMetainfo)
    File BWA0712
    File PICARD
    File SAMTOOLS131
    File GATK36
    File VarScan

    File ref_fasta
    File ref_dict
    File ref_index
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    String THREAD_NUM

    String panel_bed
    String Panel_bed_path = "/work-a/user/yanll/Tools/genecast-workflow/data/regions/medexome.bed.gz"
		Array[Array[File]] RawBams
    scatter (pair in RawBams){

	#call SortBam as firstSortBam {
    #  input:
    #  bamFile = pair[1],
    #  PICARD = PICARD,
    #  ref_fasta = ref_fasta,
    #  ref_dict = ref_dict,
    #  ref_index = ref_index,
    #  sampleName = pair[0]
#
	 #   }

#	call BamFilter {
#	    input:
#	    BamFile = firstSortBam.sorted_Bam,
#	    SAMTOOLS = SAMTOOLS131,
#      ref_fasta = ref_fasta,
#      ref_dict = ref_dict,
#      ref_index = ref_index,
#	    sampleName = pair[0]
#
#	    }

	call SortBam as secondSortBam {
        input:
        bamFile = pair[1],
        PICARD = PICARD,
        ref_fasta = ref_fasta,
        ref_dict = ref_dict,
        ref_index = ref_index,
        sampleName = pair[0]
        
	    }

	call BamFlagstat {
	    input:
      BamFile = secondSortBam.sorted_Bam,
	    SAMTOOLS = SAMTOOLS131,
      ref_fasta = ref_fasta,
      ref_dict = ref_dict,
      ref_index = ref_index

	    }

    call MarkDuplicates{
        input:
        bamFile = secondSortBam.sorted_Bam,
        PICARD = PICARD,
        sampleName = pair[0]

        }

    call mpileup_all {
	    input:
	    BamFile = MarkDuplicates.deduped_Bam,
        SAMTOOLS = SAMTOOLS131,
        ref_fasta = ref_fasta,
        ref_dict = ref_dict,
        ref_index = ref_index,
	    sampleName = pair[0]
	    }

	call mpileup_target {
        input:
        SAMTOOLS = SAMTOOLS131,
        BamFile = MarkDuplicates.deduped_Bam,
        panel_bed_path = Panel_bed_path,
        ref_fasta = ref_fasta,
        ref_dict = ref_dict,
        ref_index = ref_index,
        sampleName = pair[0]
        }

    call VarscanMpileup2snp {
        input:
        VarScan = VarScan,
        mpileup = mpileup_target.mpileup_target_gz,
        snp_output = pair[0]
        }

    call VarscanMpileup2indel {
        input:
        VarScan = VarScan,
        mpileup = mpileup_target.mpileup_target_gz,
        indel_output = pair[0]
        }


 }
#####
#     set BamFilePairs to call Varscan somatic command
#####
#
#Array[Array[String]] SomaticBams = [["normal1.mpileup","tumor1.mpileup"],["normal2.mpielup","tumor2.mpileup"]]
#scatter (pair in SomaticBams){
#    call VarscanSomatic{
#        input:
#        VarScan = VarScan,
#        normal_mpileup = pair[0],
#        tumor_mpileup = pair[1],
#        somatic_output = "I1_I2"
#        }
#      }


    call Merge_mpileup_info {
    input:
    files = mpileup_all.mpileup_all_gz,
    panel_bed = panel_bed

    }

    call CovDep {
    input:
    CovDep="/work-a/user/yanll/Tools/genecast-workflow/tools/covdepstat",
    mpileup_info = Merge_mpileup_info.mpileupfileStr,
    panel_dic = "/work-a/user/yanll/Tools/genecast-workflow/data/regions"

    }

}
