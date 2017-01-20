# test GATK workflow
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
        ${BWA} mem -t ${THREAD_NUM} -M ${ref_fasta} ${input_FASTQ1} ${input_FASTQ2} -R '@RG\tID:${file_basename}_${LANE}-\tLB:${file_basename}_${LANE}\tSM:${file_basename}\tPL:ILLUMINA\tPU:flowcell-barcode.lane' | ${SAMTOOLS} view -Sb - > ${file_basename}_${LANE}.bam
        }
    output {
        File aligned_bam = "${file_basename}_${LANE}.bam"
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
        java -jar ${PICARD} SortSam \
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

task  MarkDuplicates{
    File bamFile
    File PICARD
    String sampleName
    command {
        java -jar ${PICARD} MarkDuplicates \
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

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File GATK
  File input_bam
  File input_bam_index
  String recal_table_name
  
  File dbSNP_vcf
  File dbSNP_vcf_index
  File gold_indels
  File gold_indels_index

  File ref_dict
  File ref_fasta
  File ref_index
  Array[File]? BQSR_table
  String nct
  Array[String]? BQSR_cmd

  command {
    java -Xmx48g \
      -jar ${GATK} \
      -T BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -o ${recal_table_name} \
      -knownSites ${dbSNP_vcf} \
      -knownSites ${gold_indels} \
      -nct ${nct} ${BQSR_cmd[0]}${BQSR_table[0]}
  }
  output {
    Array[File] recal_table_file = ["${recal_table_name}"]
  }
}

task AnalyzeCovariates {
        File GATK
        File ref_fasta
        File ref_index
        File ref_dict
        String plots_name
        Array[File] before_recal_table
        Array[File] post_recal_table
        command{
        java -Xmx12g \
           -jar ${GATK} \
           -T AnalyzeCovariates \
           -R ${ref_fasta} \
           -before ${before_recal_table[0]} \
           -after ${post_recal_table[0]} \
           -plots ${plots_name}
           }
         output{
           File plots_file = "${plots_name}"
           }
}

task PrintReads {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    String sampleName
    File bamfile
    File bamfile_index
    Array[String] recal_table_file
    String nct
    command{
        java -Xmx48g \
           -jar ${GATK} \
           -T PrintReads \
           -R ${ref_fasta} \
           -I ${bamfile} \
           -BQSR ${recal_table_file[0]} \
      	   -nct ${nct} \
           -o ${sampleName}.recal_reads.bam
    }
    output{
        File recal_bam= "${sampleName}.recal_reads.bam"
        File recal_bam_index= "${sampleName}.recal_reads.bai"
    }
}

task HaplotypeCaller {
    File GATK
    String sampleName
    File bamFile
    File bamFile_index
    File ref_fasta
    File ref_index
    File ref_dict
    File interval_list
    String nct
    command {
        java -Xmx48g \
	-jar ${GATK} \
        -T HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${bamFile} \
        -o ${sampleName}.raw.vcf \
        -L ${interval_list} \
        -nct ${nct} \
        --max_alternate_alleles 5 \
	-stand_emit_conf 10 \
        --read_filter OverclippedRead 
        }
    output {
        File sample_gvcfs = "${sampleName}.raw.vcf"
        File sample_gvcf_indices = "${sampleName}.raw.vcf.idx"
        }
}

task select {
  File GATK
  File ref_fasta
  File ref_dict
  File ref_index
  String sampleName
  String type
  File rawVCF

  command {
    java -Xmx8g -jar ${GATK} \
      -T SelectVariants \
      -R ${ref_fasta} \
      -V ${rawVCF} \
      -selectType ${type} \
      -o ${sampleName}_raw.${type}.vcf
  }
  output {
    File rawSubset = "${sampleName}_raw.${type}.vcf"
  }
}

task hardFilterSNP {
  File GATK
  File ref_fasta
  File ref_dict
  File ref_index
  String sampleName
  File rawSNPs

  command {
    java -Xmx8g -jar ${GATK} \
      -T VariantFiltration \
      -R ${ref_fasta} \
      -V ${rawSNPs} \
      --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
      --filterName "snp_filter" \
      -o ${sampleName}.filtered.snps.vcf
  }
  output {
    File filteredSNPs = "${sampleName}.filtered.snps.vcf"
  }
}

task hardFilterIndel {
  File GATK
  File ref_fasta
  File ref_dict
  File ref_index
  String sampleName
  File rawIndels

  command {
    java -Xmx8g -jar ${GATK} \
      -T VariantFiltration \
      -R ${ref_fasta} \
      -V ${rawIndels} \
      --filterExpression "D < 2.0 || FS	> 200.0 || SOR > 10.0 || InbreedingCoeff < -0.8 || ReadPosRankSum < -20.0" \
      --filterName "indel_filter" \
      -o ${sampleName}.filtered.indels.vcf
  }
  output {
    File filteredIndels = "${sampleName}.filtered.indels.vcf"
  }
}

task combine {
  File GATK
  File ref_fasta
  File ref_dict
  File ref_index
  String sampleName
  File filteredSNPs
  File filteredIndels

  command {
    java -Xmx12g -jar ${GATK} \
      -T CombineVariants \
      -R ${ref_fasta} \
      -V ${filteredSNPs} \
      -V ${filteredIndels} \
      --genotypemergeoption UNSORTED \
      -o ${sampleName}.filtered.snps.indels.vcf
  }
  output {
    File filteredVCF = "${sampleName}.filtered.snps.indels.vcf"
  }
}

task MergeBamStr {
    
    Array[File] files
    
    command <<<
	python <<CODE
	import os
	uMbamlist="${write_lines(files)}"
	sampleBamMap = {}
	with open(uMbamlist, 'r') as vf:
    	    for line in vf:
                tmpBam = line.strip()
                tmpBamDic =  os.path.split(tmpBam)[0]
                tmpBamName = os.path.split(tmpBam)[1]
                tmpName = tmpBamName.split('.')[0]
                sample = tmpName.split('_')[0]
                #lane = tmpName.split('_')[0]
                sampleBamMap.setdefault(sample,[]).append(tmpBam)
	for key in sampleBamMap:
    	    fileListStr = '\t'.join(sampleBamMap[key])
    	    print (key+'\t'+fileListStr)
	CODE        
	>>>
    output {
        Array[Array[String]] bamfileStr = read_tsv(stdout())
        }
    
    }

task MergeBams {
    Array[String]+ bamfileList
    String sampleName
    File PICARD
    command {
        java -Xmx8g -jar ${PICARD} MergeSamFiles \
            I= ${sep=" INPUT= " bamfileList} \
            O= ${sampleName}.merged.bam
        }
    output {
        File mergedBam = "${sampleName}.merged.bam"
        }  
        
    }

task readSampleName {
    Array[String] ll
    command { echo '${ll[0]}'}
    output {String sample = read_string(stdout())}    
    }

task readMergeBamList {
    Array[String] ll
    command <<<
        python <<CODE
        f= open("${write_lines(ll)}","r")
        lines = [i.strip() for i in f.readlines()]
        f.close()
        for line in lines[1:]:
            print line
        CODE
    >>>
    output { 
        Array[String] list = read_lines(stdout())
        }
    }

# WORKFLOW DEFINITIONS

workflow GATK_test {
    File sampleMetainfo
    Array[Array[String]] input_FASTQ = read_tsv(sampleMetainfo)
    File BWA0712
    File PICARD
    File SAMTOOLS131
    File GATK36    

    File ref_fasta
    File ref_dict
    File ref_index
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    String THREAD_NUM

    File dbSNP_vcf
    File dbSNP_vcf_index
    File gold_indels
    File gold_indels_index

    String panel_bed = "/work-a/user/guoh/Data/MedExome_design_files/MedExome_hg19_capture_targets.bed"
#  using filtered Bam files to call GATK workflow
#  Filtered Bams could be generated from VarScan workflow first.
    Array[Array[File]] FinalBams


    scatter (BamPair in FinalBams) {
        call MarkDuplicates {
           input:
           bamFile = BamPair[1],
           PICARD = PICARD,
           sampleName = BamPair[0]
           }
        call BuildBamIndex {
           input:
           bamFile = MarkDuplicates.deduped_Bam,
           SAMTOOLS = SAMTOOLS131,
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_index = ref_index
        
           }
        call BaseRecalibrator as firstBR {
           input:
           GATK = GATK36,
           input_bam= MarkDuplicates.deduped_Bam,
           input_bam_index = BuildBamIndex.bam_index,
           recal_table_name="recal_data.table",
           dbSNP_vcf = dbSNP_vcf ,
           dbSNP_vcf_index = dbSNP_vcf_index,
           gold_indels = gold_indels,
           gold_indels_index = gold_indels_index,
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_index = ref_index,
           nct = "6"
          }

        #call BaseRecalibrator as secondBR {
        #    input:
        #    GATK = GATK36,
        #    input_bam= MarkDuplicates.deduped_Bam,
        #        input_bam_index = BuildBamIndex.bam_index,
        #    recal_table_name="post_recal_data.table",
        #    dbSNP_vcf = dbSNP_vcf,
        #    dbSNP_vcf_index = dbSNP_vcf_index,
        #    gold_indels = gold_indels,
        #    gold_indels_index = gold_indels_index,
        #    ref_fasta = ref_fasta,
        #    ref_dict = ref_dict,
        #    ref_index = ref_index,
        #    nct = "6",
        #    BQSR_cmd = ["-BQSR "],
        #    BQSR_table= firstBR.recal_table_file
        #    }
        #call AnalyzeCovariates {
        #    input:
        #    GATK = GATK36,
        #    ref_fasta = ref_fasta,
        #    ref_dict = ref_dict,
        #    ref_index = ref_index,
        #    before_recal_table = firstBR.recal_table_file,
        #    post_recal_table = secondBR.recal_table_file,
        #    plots_name = "recalibration_plots.pdf"
        #}

        call PrintReads {
           input:
           GATK = GATK36,
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_index = ref_index,
           sampleName = BamPair[0],
           bamfile= MarkDuplicates.deduped_Bam,
           bamfile_index = BuildBamIndex.bam_index,
           recal_table_file = firstBR.recal_table_file,
           nct = "6"
        }
        call HaplotypeCaller {
           input:
           GATK = GATK36,
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_index = ref_index,
           sampleName = BamPair[0],
           nct = "6",
           bamFile = PrintReads.recal_bam,
           bamFile_index = PrintReads.recal_bam_index,
           interval_list = panel_bed
        }
        call select as selectSNPs {
           input:
           sampleName = BamPair[0],
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_index = ref_index,
           GATK=GATK36,
           type="SNP",
           rawVCF=HaplotypeCaller.sample_gvcfs
       }

        call select as selectIndels {
           input:
           sampleName = BamPair[0],
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_index = ref_index,
           GATK=GATK36,
           type="INDEL",
           rawVCF=HaplotypeCaller.sample_gvcfs
       }     

        call hardFilterSNP {
           input:
           sampleName = BamPair[0],
           GATK=GATK36,
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_index = ref_index,
           rawSNPs=selectSNPs.rawSubset
      }
        call hardFilterIndel {
           input:
           sampleName = BamPair[0],
           GATK=GATK36,
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_index = ref_index,
           rawIndels=selectIndels.rawSubset
      }
	    call combine {
           input:
           sampleName = BamPair[0],
           GATK=GATK36,
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_index = ref_index,
           filteredSNPs=hardFilterSNP.filteredSNPs,
           filteredIndels=hardFilterIndel.filteredIndels
      }
 }
}


#####
#Merge bam files of same sample from two lanes
#
#####
#
#scatter (pair in input_FASTQ){
#   call BwaMem {}
#   call SortBam as firstSortBam {}
#}
#
#call MergeBamStr {
#        input: files = firstSortBam.sorted_Bam
#           }
#    scatter (line in MergeBamStr.bamfileStr){
#        call readSampleName{
#           input: ll = line
#           }
#        call readMergeBamList{
#           input: ll = line
#           }
#        call MergeBams{
#           input:
#           bamfileList = readMergeBamList.list,
#           sampleName = readSampleName.sample,
#           PICARD = PICARD
#           }
#	}
