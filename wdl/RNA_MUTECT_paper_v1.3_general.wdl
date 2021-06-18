import 'https://api.firecloud.org/ga4gh/v1/tools/bbecerra:CramToBam.wdl_copy/versions/4/plain-WDL/descriptor' as CramToBam
import 'https://api.firecloud.org/ga4gh/v1/tools/bbecerra:BamRealigner_wdl_v1_0/versions/8/plain-WDL/descriptor' as BamRealigner
import 'https://api.firecloud.org/ga4gh/v1/tools/bbecerra:RemoveAtypicalReads/versions/8/plain-WDL/descriptor' as RemoveAtypicalReads

#The next 3 tasks correspond to the CallSomaticMutations_V131 task under the "Call Somatic Mutations For Capture" workflow.
task CallSomaticMutations_131_Prepare_Task {
		#RUNTIME INPUT PARAMS
		Int preemptible

		#TASK INPUT PARAMS
		File? sampleTargetsIntervalList
		File workspaceTargetsIntervalList
		File refFasta
		File refFastaDict
		File refFastaIdx
		Int nWay
		Float tBamSize
		Float tBaiSize
		Float nBamSize
		Float nBaiSize
		Float dbSNPVCFSize
		Float cosmicVCFSize
		Float normalPanelSize
		Float fastaSize
		Float fastaDictSize
	command <<<
		#increase verbosity
		set -x

		#Calculate disk size for all shards of the mutect run
		SIZE_FILE=split_base_sizes_disk.dat
		awk 'BEGIN { print int((${tBamSize}+${tBaiSize}+${nBamSize}+${nBaiSize}+${dbSNPVCFSize}+${cosmicVCFSize}+${normalPanelSize}+${fastaSize}+${fastaDictSize}+3000000000)/1000000000+1) }' > $SIZE_FILE
		
		#Create list of indices for the scatter job
		seq 0 $((${nWay}-1)) > indices.dat

		#Run the prepare task that splits the .interval_list file into subfiles
		MUTECT_INTERVALS="${if defined(sampleTargetsIntervalList) then sampleTargetsIntervalList else workspaceTargetsIntervalList}"
		java -Xmx2g -jar /usr/local/bin/GatkScatterGatherPrepare.jar . ${nWay} --intervals $MUTECT_INTERVALS --reference_sequence ${refFasta} 

		>>>

	output  {
		Array[File] interval_files=glob("gatk-scatter.*")
		Int mutectDiskGB=read_int("split_base_sizes_disk.dat")
		Array[Int] scatterIndices=read_lines("indices.dat")
		}

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 10 HDD"
		}
	}
task Mutect1_Task {
		#RUNTIME INPUT PARAMS
		Int preemptible
		Int diskGB
		#TASK INPUT PARAMS
		File mutectIntervals 
		File tumorBam
		File normalBam
		File tumorBamIdx
		File normalBamIdx
		File refFasta
		File refFastaDict
		File refFastaIdx
		Float fracContam
		Float contamFloor
		File dbSNPVCF
        File dbSNPVcfIndex
		File cosmicVCF
		Int downsampleToCoverage
		File? readGroupBlackList
		File normalPanel
		String pairName
		String caseName
		String ctrlName

		String? optional1
		String? optional2
		String? optional3
		String? optional4
		String? optional5
	command <<<
	#increase verbosity
	set -x	

	#compute/apply contamination floor for effective contamination
	EFFECTIVE_CONTAMINATION=`/usr/local/bin/python -c 'import sys;print sys.argv[1] if(float(sys.argv[1])>=float(sys.argv[2])) else sys.argv[2]'  ${fracContam} ${contamFloor}` ;
	echo "EFFECTIVE_CONTAMINATION computed to be $EFFECTIVE_CONTAMINATION"

	#mutect 1
	java -jar -Xmx9g /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
	-L ${mutectIntervals}  --normal_sample_name ${ctrlName} -I:normal  ${normalBam}  \
	--tumor_sample_name ${caseName} -I:tumor ${tumorBam}  \
	--reference_sequence ${refFasta} \
	--fraction_contamination $EFFECTIVE_CONTAMINATION --dbsnp ${dbSNPVCF} \
	--cosmic ${cosmicVCF}  ${"--read_group_black_list" + readGroupBlackList}\
	--out ${pairName}.call_stats.txt --coverage_file ${pairName}.coverage.wig.txt \
	--power_file ${pairName}.power.wig.txt  --downsample_to_coverage ${downsampleToCoverage} \
	--normal_panel ${normalPanel} --enable_extended_output ${optional1} ${optional2} ${optional3} ${optional4} ${optional5}

	>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		memory: "12 GB"
		disks: "local-disk ${diskGB} HDD"
		}


	output {
		File mutect1_cs="${pairName}.call_stats.txt"
		File mutect1_pw="${pairName}.power.wig.txt"
		File mutect1_cw="${pairName}.coverage.wig.txt"
		}

	}
task Gather_Task {
		#RUNTIME INPUT PARAMS
		Int preemptible

		#TASK INPUT PARAMS
		Array[File] mutect1_cs
		Array[File] mutect1_pw
		Array[File] mutect1_cw
		String pairName
	command <<<
		#increase verbosity
		set -x

		#mutect1 call_stats merging
		MUTECT1_CS="${pairName}.call_stats.txt"
		head --lines=2 ${mutect1_cs[0]} > $MUTECT1_CS
		cat ${sep =' ' mutect1_cs} | grep -Pv '#'|grep -Pv '^contig' >> $MUTECT1_CS

		MUTECT1_PW="${pairName}.power.wig.txt"
		MUTECT1_PW_GZ="${pairName}.power.wig.txt.gz"
		head --lines=1 ${mutect1_pw[0]} > $MUTECT1_PW
		cat ${sep =' ' mutect1_pw} | grep -Pv '^track' >> $MUTECT1_PW
		tar -zcvf $MUTECT1_PW_GZ $MUTECT1_PW

		MUTECT1_CW="${pairName}.coverage.wig.txt"
		MUTECT1_CW_GZ="${pairName}.coverage.wig.txt.gz"
		head --lines=1 ${mutect1_cw[0]} > $MUTECT1_CW
		cat ${sep =' ' mutect1_cw} | grep -Pv '^track' >> $MUTECT1_CW
		tar -zcvf $MUTECT1_CW_GZ $MUTECT1_CW
		>>>


	output {
		File Mutect_call_stats="${pairName}.call_stats.txt"
		File Mutect_power_wig="${pairName}.power.wig.txt.gz"
		File Mutect_coverage_wig="${pairName}.coverage.wig.txt.gz"
		}


	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 10 HDD"
		}
	}

#This task corresponds to CallstatsToMaflite_V14 task under the "Call Somatic Mutations For Capture" workflow.
task CallstatsToMaflite_14 {
		#RUNTIME INPUT PARAMS
		Int preemptible

		File callstatsFile
		Int genomeBuild
		String mode
		String pairName
		String extraColumns
		String outputFile = "${pairName}.maf"
	command <<<
		#increase verbosity
		set -x
		echo "Calling call_stats_to_maflite.pl on input file ${callstatsFile}. Dumping to ${outputFile}"

		perl /usr/local/bin/call_stats_to_maflite.pl ${callstatsFile} ${genomeBuild} ${mode} ${outputFile} "${extraColumns}"

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 10 HDD"
		}

	output {
		File maflite_capture="${outputFile}"
	}
}

task SummarizeWigFile_1 {
		#RUNTIME INPUT PARAMS
		Int preemptible

		File coverageFile
		String pairName
		String outputFile = "${pairName}.somatic_coverage_summary.txt"
	command <<<
		#increase verbosity
		set -x
		echo "Calling summarizeWigFile.pl on input file ${coverageFile}. Dumping to ${outputFile}"

		
		tar -zxvf ${coverageFile}
		perl /usr/local/bin/summarizeWigFile.pl "${pairName}.coverage.wig.txt" ${outputFile}

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 10 HDD"
		}

	output {
		File somatic_mutations_covered_bases_capture="${outputFile}"
	}
}

task ApplyMafValidation_23 {
		#RUNTIME INPUT PARAMS
		Int preemptible

		File mafliteFile
		String sampleName
		String matchMode
		File valdb
		String valDir = sub(valdb, "\\.tar.gz$", "") 

		String outputFile="${sampleName}.maf.annotated"
	command <<<
		#increase verbosity
		set -x
		echo "Calling ApplyMafValidation.java on input file ${mafliteFile}. Dumping to ${outputFile}"

		VAL_BASENAME=`basename ${valdb}`
		VAL_DIR=`echo $VAL_BASENAME | cut -d'.' -f1`
		tar -zxvf ${valdb}
		
		java -Xmx1g -jar /usr/local/bin/ApplyMAFValidation.jar M=${mafliteFile} OUTPUT_MAF=${outputFile} MATCH_MODE=${matchMode} V=$VAL_DIR

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 10 HDD"
		}

	output {
		File snp_validated_maflite_capture="${outputFile}"
	}
}

task Oncotate_74 {
		#RUNTIME INPUT PARAMS
		Int preemptible
		#TASK INPUT PARAMS
		File snpValidatedMafliteCaptureFile
		String genomeBuild
		File oncoDBTarBall
		String pairName
		String outputFilename = "${pairName}.snp.capture.maf.annotated"
		File defaultConfig
		String inputFormat
		String outputFormat
		String txMode
		String additionalOptions
		File? canonicalTXFile
		String? logFilename
		Int diskGB = ceil(size(snpValidatedMafliteCaptureFile,"G")+size(oncoDBTarBall,"G")*5+size(defaultConfig,"G")+3)

	command <<<
		#increase verbosity
		set -x
		#find TARBALL type
		TYPE=`echo 'if("${oncoDBTarBall}"=~m/z$/) { print "GZ" ; } else { print "TAR" ; } '| perl` ; 

		#obtain the name of the directory for oncodb
		#and unpack based on the TYPE
		if [ "$TYPE" == "GZ" ] ; then
			ONCO_DB_DIR_NAME=`gunzip -c ${oncoDBTarBall} |tar -tf /dev/stdin|head -1` ; 
			tar -xzf ${oncoDBTarBall}
		else
			ONCO_DB_DIR_NAME=`tar -tf ${oncoDBTarBall} |head -1` ;
			tar -xf ${oncoDBTarBall} ;
		fi ;
		echo "ONCO_DB_DIR_NAME is $ONCO_DB_DIR_NAME" ; 

		#Run the merged filtered VCF (from both mutects through Oncotator) 
		/usr/local/lib/python2.7/site-packages/Oncotator-1.8.0.0-py2.7.egg/oncotator/Oncotator.py\
			-v --input_format=${inputFormat} --output_format=${outputFormat} --db-dir=$ONCO_DB_DIR_NAME --no-multicore  --default_config ${defaultConfig}\
			${snpValidatedMafliteCaptureFile} ${outputFilename} ${genomeBuild} --tx-mode=${txMode} ${additionalOptions} ${"--log_name "+logFilename} ${"-c "+canonicalTXFile}


		#Run the adapter from Oncotator_1.8.0 to Oncotator_1.6.1. It changes the names of the output columns to fit with code that expects the 1.6.1 version names.
		python /usr/local/bin/oncotator_backwards_adapter_180_to_161.py ${outputFilename}


		/usr/local/bin/run_checkEmptyMAF.sh /usr/local/MATLAB/MATLAB_Runtime/v901/ ${outputFilename}
		>>>


	output {
		File maf_file_capture = "${outputFilename}"
		Boolean should_keep_running = read_boolean("should_keep_running.txt")
		}


	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		memory: "10 GB"
		disks: "local-disk ${diskGB} HDD"
		}
	}

task MasterFilter_4 {
    #RUNTIME INPUT PARAMS
	Int preemptible

	File mafFileCapture
	String pairName
	File ponFile
	String outputFile = "${pairName}.pon_filtered.txt"
	String outputFile2 = "${pairName}.pon_annotated.txt"
	Int diskGB = ceil(size(mafFileCapture, "G")+size(ponFile, "G")+3)
	command <<<
		#increase verbosity
		set -x
		echo "Calling Matlab script master_filter_wrapper on input file ${mafFileCapture}."

		/usr/local/bin/master_filter/run_master_filter_wrapper.sh /usr/local/MATLAB/MATLAB_Compiler_Runtime/v717/ ${mafFileCapture} ${ponFile} ${pairName}
		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk ${diskGB} HDD"
		}

	output {
		File maf_file_capture_master_filter_removed="${outputFile}"
		File maf_file_capture_master_filter_annotated="${outputFile2}"
	}
}


task CreateOxoGIntervals_22 {
		#RUNTIME INPUT PARAMS
		Int preemptible

		File ponFilteredFile
		String pairName
		String outputFile = "${pairName}.MasterPoN.interval_list"
	command <<<
		#increase verbosity
		set -x
		echo "Calling create createOxoGIntervals.py on input file ${ponFilteredFile}. Dumping to ${outputFile}"

		python /usr/local/bin/createOxoGIntervals.py ${ponFilteredFile} ${outputFile}

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 10 HDD"
		}

	output {
		File maf_PoN_interval_list="${outputFile}"
	}
}

task HiSat2RealignPreprocess {
		#RUNTIME INPUT PARAMS
		Int preemptible

		File ponFilteredFile
		File sampleBAM
		File sampleBAI
		File refFasta
		File refFastaDict
		File refFastaIdx
		String sampleName
		
		String tmp_sequence_1 = "tmp_sequence_1.fastq"
		String tmp_sequence_2 = "tmp_sequence_2.fastq"
		
        Int bufferGB = 10
		Int diskGB = ceil(size(ponFilteredFile, "G")+size(sampleBAM, "G")+size(sampleBAI, "G")+size(refFasta, "G")+size(refFastaDict, "G")+size(refFastaIdx, "G")+bufferGB)
	command <<<
		#increase verbosity
		set -x

		echo "running HiSat_realign_preprocess.sh"

		cp -R /usr/local/bin/HiSat_realign_preprocess/* .
		sh HiSat_realign_preprocess.sh . ${ponFilteredFile} ${sampleBAM} ${refFasta} ${sampleName}

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		memory: "5 GB"
		disks: "local-disk ${diskGB} HDD"
    }

	output {
		File rna_reads_fastq_1 = "${tmp_sequence_1}"
		File rna_reads_fastq_2 = "${tmp_sequence_2}"

	}
}

task HiSat2Realign_3 {
		#RUNTIME INPUT PARAMS
		Int preemptible

		String sampleName
		
		File HISAT_index
		File fastq1
		File fastq2
		String output_prefix = "${sampleName}"
		Int? min_intronlen
		Int? max_intronlen
		File? known_splicesite_infile 
		File? novel_splicesite_infile 
		String? rna_strandness 
		Int? k 
		Int? minins 
		Int? maxins 
		String? mate_orientation 
		String? rg_id 
		String? rg 
		Int? num_threads

		String outputBAM = "${sampleName}.aligned.sorted_by_coord.hisat2.bam"
		String outputBAI = "${sampleName}.aligned.sorted_by_coord.hisat2.bam.bai"
		String outputJunctions = "${sampleName}.HJ.out.tab"

		Int diskGB = ceil(size(HISAT_index, "G")*5+size(fastq1, "G")+size(fastq2, "G")+5)
	command <<<
		echo "running HISAT2 realignment process"
		
		echo -e "${fastq1}\n${fastq2}"   >> ${sampleName}.rna_reads_fastq_list.list

		HISAT_INDEX_DIRECTORY=`gunzip -c ${HISAT_index} |tar -tf /dev/stdin|head -1` ;
		HISAT_INDEX_DIRECTORY=`echo $HISAT_INDEX_DIRECTORY broad_hg19_gencode19 | tr -d " "` ;
		tar -xzvf ${HISAT_index}

		# In preparation for the run we need to make some env changes
		PATH=/usr/local/bin/HISAT2/hisat2-2.0.3-beta:$PATH
		ln -sf /usr/local/bin/samtools-1.3 /usr/local/bin/samtools

		python3 /usr/local/bin/HISAT2/run_HISAT2.py $HISAT_INDEX_DIRECTORY  ${sampleName}.rna_reads_fastq_list.list ${output_prefix} \
		${"--min_intronlen " + min_intronlen} ${"--max_intronlen " + max_intronlen} \
		${"--known_splicesite_infile " + known_splicesite_infile} ${"--novel_splicesite_infile " + novel_splicesite_infile} ${"--rna_strandness " + rna_strandness} \
		${"-k " + k} ${"-I " + minins} ${"-X " + maxins} ${"--mate_orientation " + mate_orientation} ${"--rg_id " + rg_id} ${"--rg " + rg} ${"-p " + num_threads}
		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		memory: "16 GB"
		disks: "local-disk ${diskGB} HDD"
		}

	output {
		File hisat2_bam_file = "${outputBAM}"
		File hisat2_bai_file = "${outputBAI}"
		File hisat2_junctions = "${outputJunctions}"
	}
}

task FilterRNAMutations {
		#RUNTIME INPUT PARAMS
		Int preemptible

		File ponFilteredFile
		File callStatsFile
		String pairName
		File DarnedMat
		File ExacMat
		File RadarMat
		File ponGTEx
		Int minAlt
		Int PoNThr
        File cytoBand

		 Int diskGB = ceil(size(ponFilteredFile, "G")+size(callStatsFile, "G")+size(DarnedMat, "G")+size(ExacMat, "G")+size(RadarMat, "G")+(size(ponGTEx, "G")*25) + 10)
	command <<<
		#increase verbosity
		set -x

        tar -xzf ${ponGTEx}
        ponGTExUnzipped=$(basename ${ponGTEx} .tar.gz)
        echo $ponGTExUnzipped

		/usr/local/bin/run_FilterRNAMutations.sh /usr/local/MATLAB/MATLAB_Runtime/v901/ ${pairName} ${ponFilteredFile} ${callStatsFile} ${minAlt} ${PoNThr} ${DarnedMat} \
		${RadarMat} ${ExacMat} $ponGTExUnzipped ${cytoBand}
		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		memory: "20 GB"
		disks: "local-disk ${diskGB} HDD"
		}

	output {
		File maf_file_rna_post_realign_filter = "${pairName}.intersect.txt"
		File maf_file_rna_post_filtering = "${pairName}.post_filtering.txt"
		File maf_file_rna_prefiltering_with_info = "${pairName}.pre_filtering_plus_info.txt"
	}
}


task FilterRNAMutationsBasedOnDuplicateReads {
		#RUNTIME INPUT PARAMS
		Int preemptible

		String pairName
		File starBAM
		File starBAI
		File postFilteringFile
		Int minAlt

		Int diskGB = ceil((2 * size(starBAM, "G")) + (2 * size(starBAI, "G")) +size(postFilteringFile, "G")+3)
        
        Int memoryGB = 20
        
	command <<<
		#increase verbosity
		set -x
		
		NEWBAM="/usr/input.bam"
		NEWBAI="/usr/input.bam.bai"
		
		cp ${starBAM} $NEWBAM
		cp ${starBAI} $NEWBAI

		/usr/local/bin/run_FilterRNAMutationsBasedOnDuplicateReadsV2.sh /usr/local/MATLAB/MATLAB_Runtime/v901/ ${pairName} $NEWBAM ${postFilteringFile} ${minAlt}

		OUTPUT=${pairName}.post_filtering_remove_duplicates.txt
		if [ -f "$OUTPUT" ]; then
			echo "true" >> is_output.txt
		else 
			# Likely not output due to error or due to no variants left:
			echo "false" >> is_output.txt
			cp ${postFilteringFile} ${pairName}.post_filtering_remove_duplicates.txt
		fi
		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		memory: "${memoryGB} GB"
		disks: "local-disk ${diskGB} HDD"
		}

	output {
		File? maf_file_rna_post_filtering_dup_GTEx = "${pairName}.post_filtering_remove_duplicates.txt"
		Boolean maf_file_rna_post_filtering_dup_GTEx_outputted = read_boolean('is_output.txt')
	}
}

workflow CallingGroup_Workflow {
		#RUNTIME INPUT PARAMS
		Int preemptible

		#WORKFLOW INPUT PARMS
		#General inputs for the whole workflow
		String pairName
		String tumorName
		File tumorBamRNA
		File tumorBamRNAIdx
		String normalName
		File normalBamDNA
		File normalBamDNAIdx
		File refFasta
		File refFastaDict
		File refFastaIdx

		#Input for BamRealignment
		File originalRefFasta
		File originalRefFastaIndex

		File refBwt
		File refSa
		File refAmb
		File refAnn
		File refPac
		File dbSNPVcfIndex

		Array[File] knownIndelsSitesVCFs
		Array[File] knownIndelsSitesIndices

		#Inputs for CallSomaticMutations_131_Prepare_Task
		Int scatterNWay
		File? sampleTargetsIntervalList
		File workspaceTargetsIntervalList

        #Inputs for MasterFilter_4
        File snpFilteringPoNFile

		#Inputs for Mutect1_Task	
		Float fracContam
		File dbSNPVCF
		File cosmicVCF
		Int downsampleToCoverage
		File? readGroupBlackList
		File normalPanelVCF
		String? optional1
		String? optional2
		String? optional3
		String? optional4
		String? optional5

		#Inputs for CallstatsToMaflite
		Int genomeBuild_int
		String mode
		String extraColumns

		# Inputs for ApplyMafValidation
		File valdb
		String matchMode

		#Inputs for Oncotator
		String genomeBuild_string
		File oncoDBTarBall
		File defaultConfig
		String inputFormat
		String outputFormat
		String txMode
		String additionalOptions
		File? canonicalTXFile
		String? logFilename


		#Inputs for HiSat2Realign_3
		File HISAT_index
		Int? min_intronlen
		Int? max_intronlen
		File? known_splicesite_infile 
		File? novel_splicesite_infile 
		String? rna_strandness 
		Int? k 
		Int? minins 
		Int? maxins 
		String? mate_orientation 
		String? rg_id 
		String? rg 
		Int? num_threads

		#Inputs for CallSomaticMutations_131_Prepare_Task after hisat realignment
		Int hisat_scatterNWay
		File? hisat_sampleTargetsIntervalList

		#Inputs for Mutect1_Task after hisat realignment
		Int hisat_downsampleToCoverage
		File? hisat_readGroupBlackList
		String? hisat_optional1
		String? hisat_optional2
		String? hisat_optional3
		String? hisat_optional4
		String? hisat_optional5

		#Inputs for FilterRNAMutations
		File DarnedMat
		File ExacMat
		File RadarMat
        File ponGTEx
		Int minAlt_filter_generic
        Int PoNThr
	    File cytoBand
		#Inputs for FilterRNAMutations_duplicateReads
		Int minAlt_filter_duplicates

		String samtools_docker = "jweinstk/samtools:latest"
        
        Int m2BufferGB = 15

	# PREPARE FOR SCATTER
	call CallSomaticMutations_131_Prepare_Task as CallSomaticMutations_131_Prepare_Task_STAR{
		input: 
			refFasta=refFasta,
			refFastaDict=refFastaDict,
			refFastaIdx=refFastaIdx,
			sampleTargetsIntervalList=sampleTargetsIntervalList,
			workspaceTargetsIntervalList=workspaceTargetsIntervalList,
			nWay=scatterNWay,
			preemptible=preemptible,
			tBamSize=size(tumorBamRNA),
			tBaiSize=size(tumorBamRNAIdx),
			nBamSize=size(normalBamDNA),
			nBaiSize=size(normalBamDNAIdx),
			dbSNPVCFSize=size(dbSNPVCF),
			cosmicVCFSize=size(cosmicVCF),
			normalPanelSize=size(normalPanelVCF),
			fastaSize=size(refFasta),
			fastaDictSize=size(refFastaDict)

	}

	#SCATTER AND ANALYZE
	scatter (idx in CallSomaticMutations_131_Prepare_Task_STAR.scatterIndices) {
			
			call Mutect1_Task as Mutect1_Task_STAR {
				input:
					contamFloor=0.01,
					tumorBam=tumorBamRNA,
					normalBam=normalBamDNA,
					tumorBamIdx=tumorBamRNAIdx,
					pairName=pairName,
					caseName=tumorName,
					ctrlName=normalName,
					normalBamIdx=normalBamDNAIdx,
					mutectIntervals=CallSomaticMutations_131_Prepare_Task_STAR.interval_files[idx],
					refFasta=refFasta,
					refFastaDict=refFastaDict,
					refFastaIdx=refFastaIdx,
					fracContam=fracContam,
					dbSNPVCF=dbSNPVCF,
                    dbSNPVcfIndex=dbSNPVcfIndex,
					cosmicVCF=cosmicVCF,
					downsampleToCoverage=downsampleToCoverage,
					readGroupBlackList=readGroupBlackList,
					normalPanel=normalPanelVCF,
					preemptible=preemptible,
					diskGB=CallSomaticMutations_131_Prepare_Task_STAR.mutectDiskGB + m2BufferGB,
					optional1=optional1,
					optional2=optional2,
					optional3=optional3,
					optional4=optional4,
					optional5=optional5
					
				}

			
			}


	call Gather_Task as Gather_Task_STAR {
		input:
			pairName=pairName,
			mutect1_cs=Mutect1_Task_STAR.mutect1_cs,
			mutect1_pw=Mutect1_Task_STAR.mutect1_pw,
			mutect1_cw=Mutect1_Task_STAR.mutect1_cw,
			preemptible=preemptible,

	}


	call CallstatsToMaflite_14 {
		input:
			preemptible=preemptible,
			callstatsFile=Gather_Task_STAR.Mutect_call_stats,
			genomeBuild=genomeBuild_int,
			mode=mode,
			pairName=pairName,
			extraColumns=extraColumns
	}

	call SummarizeWigFile_1 {
		input:
			preemptible=preemptible,
			coverageFile=Gather_Task_STAR.Mutect_coverage_wig,
			pairName=pairName
	}

	call ApplyMafValidation_23 {
		input:
			preemptible=preemptible,
			mafliteFile=CallstatsToMaflite_14.maflite_capture,
			sampleName=tumorName,
			valdb=valdb,
			matchMode=matchMode
	}

	call Oncotate_74 {
		input:
			preemptible=preemptible,
			snpValidatedMafliteCaptureFile=ApplyMafValidation_23.snp_validated_maflite_capture,
			genomeBuild=genomeBuild_string,
			oncoDBTarBall=oncoDBTarBall,
			pairName=pairName,
			defaultConfig=defaultConfig,
			inputFormat=inputFormat,
			outputFormat=outputFormat,
			txMode=txMode,
			additionalOptions=additionalOptions,
            canonicalTXFile=canonicalTXFile,
            logFilename=logFilename,


	}


	if (Oncotate_74.should_keep_running) {
        call MasterFilter_4 {
                input:
                    preemptible=preemptible,
                    mafFileCapture=Oncotate_74.maf_file_capture,
                    ponFile=snpFilteringPoNFile,
                    pairName=pairName
        }

		call CreateOxoGIntervals_22 {
			input:
				preemptible=preemptible,
				ponFilteredFile=MasterFilter_4.maf_file_capture_master_filter_removed,
				pairName=pairName,
		}

		call HiSat2RealignPreprocess as HiSat2RealignPreprocess_tumor {
			input:
				preemptible=preemptible,
				ponFilteredFile=MasterFilter_4.maf_file_capture_master_filter_removed,
				sampleBAM=tumorBamRNA,
				sampleBAI=tumorBamRNAIdx,
				refFasta=refFasta,
				refFastaDict=refFastaDict,
				refFastaIdx=refFastaIdx,
				sampleName=tumorName

		}

		call RemoveAtypicalReads.RemoveAtypicalReads as RemoveAtypicalReads_tumor {
			input:
				fastq1 = HiSat2RealignPreprocess_tumor.rna_reads_fastq_1,
				fastq2 = HiSat2RealignPreprocess_tumor.rna_reads_fastq_2,
				preemptible = preemptible
		}

		call HiSat2Realign_3 as HISAT2_tumor {
			input:
				preemptible=preemptible,
				sampleName=tumorName,
				HISAT_index=HISAT_index,
				fastq1=RemoveAtypicalReads_tumor.fastq1_removed,
				fastq2=RemoveAtypicalReads_tumor.fastq2_removed,
				min_intronlen=min_intronlen,
				max_intronlen=max_intronlen,
				known_splicesite_infile=known_splicesite_infile,
				novel_splicesite_infile=novel_splicesite_infile,
				rna_strandness=rna_strandness,
				k=k,
				minins=minins,
				maxins=maxins,
				mate_orientation=mate_orientation,
				rg_id = rg_id,
				rg = rg,
				num_threads = num_threads 
		
		}

		call HiSat2RealignPreprocess as HiSat2RealignPreprocess_normal {
			input:
				preemptible=preemptible,
				ponFilteredFile=MasterFilter_4.maf_file_capture_master_filter_removed,
				sampleBAM=normalBamDNA,
				sampleBAI=normalBamDNAIdx,
				refFasta=refFasta,
				refFastaDict=refFastaDict,
				refFastaIdx=refFastaIdx,
				sampleName=normalName

		}

		call RemoveAtypicalReads.RemoveAtypicalReads as RemoveAtypicalReads_normal {
			input:
				fastq1 = HiSat2RealignPreprocess_normal.rna_reads_fastq_1,
				fastq2 = HiSat2RealignPreprocess_normal.rna_reads_fastq_2,
				preemptible = preemptible
		}

		call HiSat2Realign_3 as HISAT2_normal {
			input:
				preemptible=preemptible,
				sampleName=normalName,
				HISAT_index=HISAT_index,
				fastq1=RemoveAtypicalReads_normal.fastq1_removed,
				fastq2=RemoveAtypicalReads_normal.fastq2_removed,
				min_intronlen=min_intronlen,
				max_intronlen=max_intronlen,
				known_splicesite_infile=known_splicesite_infile,
				novel_splicesite_infile=novel_splicesite_infile,
				rna_strandness=rna_strandness,
				k=k,
				minins=minins,
				maxins=maxins,
				mate_orientation=mate_orientation,
				rg_id = rg_id,
				rg = rg,
				num_threads = num_threads 
		
		}

		call CallSomaticMutations_131_Prepare_Task as CallSomaticMutations_131_Prepare_Task_HISAT2 {
			input: 
				refFasta=refFasta,
				refFastaDict=refFastaDict,
				refFastaIdx=refFastaIdx,
				sampleTargetsIntervalList=hisat_sampleTargetsIntervalList,
				workspaceTargetsIntervalList=CreateOxoGIntervals_22.maf_PoN_interval_list,
				nWay=hisat_scatterNWay,
				preemptible=preemptible,
				tBamSize=size(HISAT2_tumor.hisat2_bam_file),
				tBaiSize=size(HISAT2_tumor.hisat2_bai_file),
				nBamSize=size(HISAT2_normal.hisat2_bam_file),
				nBaiSize=size(HISAT2_normal.hisat2_bai_file),
				dbSNPVCFSize=size(dbSNPVCF),
				cosmicVCFSize=size(cosmicVCF),
				normalPanelSize=size(normalPanelVCF),
				fastaSize=size(refFasta),
				fastaDictSize=size(refFastaDict)

		}


		#SCATTER AND ANALYZE
		scatter (idx in CallSomaticMutations_131_Prepare_Task_HISAT2.scatterIndices) {
				
				call Mutect1_Task as Mutect1_Task_HISAT2 {
					input:
						contamFloor=0.01,
						tumorBam=HISAT2_tumor.hisat2_bam_file,
						normalBam=HISAT2_normal.hisat2_bam_file,
						tumorBamIdx=HISAT2_tumor.hisat2_bai_file,
						pairName=pairName,
						caseName=tumorName,
						ctrlName=normalName,
						normalBamIdx=HISAT2_normal.hisat2_bai_file,
						mutectIntervals=CallSomaticMutations_131_Prepare_Task_HISAT2.interval_files[idx],
						refFasta=refFasta,
						refFastaDict=refFastaDict,
						refFastaIdx=refFastaIdx,
						fracContam=fracContam,
						dbSNPVCF=dbSNPVCF,
                        dbSNPVcfIndex=dbSNPVcfIndex,
						cosmicVCF=cosmicVCF,
						downsampleToCoverage=hisat_downsampleToCoverage,
						readGroupBlackList=hisat_readGroupBlackList,
						normalPanel=normalPanelVCF,
						preemptible=preemptible,
						diskGB=CallSomaticMutations_131_Prepare_Task_HISAT2.mutectDiskGB,
						optional1=hisat_optional1,
						optional2=hisat_optional2,
						optional3=hisat_optional3,
						optional4=hisat_optional4,
						optional5=hisat_optional5
					}

				
				}

		call Gather_Task as Gather_Task_HISAT2 {
			input:
				pairName=pairName,
				mutect1_cs=Mutect1_Task_HISAT2.mutect1_cs,
				mutect1_pw=Mutect1_Task_HISAT2.mutect1_pw,
				mutect1_cw=Mutect1_Task_HISAT2.mutect1_cw,
				preemptible=preemptible,

		}

		call FilterRNAMutations {
			input:
				preemptible=preemptible,
				ponFilteredFile=MasterFilter_4.maf_file_capture_master_filter_removed,
				callStatsFile=Gather_Task_HISAT2.Mutect_call_stats,
				pairName=pairName,
				DarnedMat=DarnedMat,
				ExacMat=ExacMat,
				RadarMat=RadarMat,
				minAlt=minAlt_filter_generic,
                ponGTEx=ponGTEx,
				PoNThr=PoNThr,
                cytoBand=cytoBand
		}

		call FilterRNAMutationsBasedOnDuplicateReads{
			input:
				preemptible=preemptible,
				pairName=pairName,
				starBAM=tumorBamRNA,
				starBAI=tumorBamRNAIdx,
				postFilteringFile=FilterRNAMutations.maf_file_rna_post_filtering,
				minAlt=minAlt_filter_duplicates
		}

	}
	output {
		File Mutect_call_stats_STAR = Gather_Task_STAR.Mutect_call_stats
		File Mutect_power_wig_STAR = Gather_Task_STAR.Mutect_power_wig
		File Mutect_coverage_wig_STAR = Gather_Task_STAR.Mutect_coverage_wig
		File maflite_capture_STAR = CallstatsToMaflite_14.maflite_capture
		File somatic_mutations_covered_bases_capture_STAR = SummarizeWigFile_1.somatic_mutations_covered_bases_capture
		File snp_validated_maflite_capture_STAR = ApplyMafValidation_23.snp_validated_maflite_capture
		File maf_file_capture_STAR = Oncotate_74.maf_file_capture
        Boolean oncotate_should_keep_running = Oncotate_74.should_keep_running
		File? maf_PoN_interval_list_STAR = CreateOxoGIntervals_22.maf_PoN_interval_list
        File? maf_file_capture_master_filter_removed = MasterFilter_4.maf_file_capture_master_filter_removed
		File? maf_file_capture_master_filter_annotated = MasterFilter_4.maf_file_capture_master_filter_annotated
		File? rna_reads_fastq_1_tumor = HiSat2RealignPreprocess_tumor.rna_reads_fastq_1
		File? rna_reads_fastq_2_tumor = HiSat2RealignPreprocess_tumor.rna_reads_fastq_2
		File? hisat2_bam_file_tumor = HISAT2_tumor.hisat2_bam_file
		File? hisat2_bai_file_tumor = HISAT2_tumor.hisat2_bai_file
		File? hisat2_junctions_tumor = HISAT2_tumor.hisat2_junctions
		File? rna_reads_fastq_1_normal = HiSat2RealignPreprocess_normal.rna_reads_fastq_1
		File? rna_reads_fastq_2_normal = HiSat2RealignPreprocess_normal.rna_reads_fastq_2
		File? hisat2_bam_file_normal = HISAT2_normal.hisat2_bam_file
		File? hisat2_bai_file_normal = HISAT2_normal.hisat2_bai_file
		File? hisat2_junctions_normal = HISAT2_normal.hisat2_junctions
		File? Mutect_call_stats_HISAT2 = Gather_Task_HISAT2.Mutect_call_stats
		File? Mutect_power_wig_HISAT2 = Gather_Task_HISAT2.Mutect_power_wig
		File? Mutect_coverage_wig_HISAT2 = Gather_Task_HISAT2.Mutect_coverage_wig
		File? maf_file_rna_post_realign_filter = FilterRNAMutations.maf_file_rna_post_realign_filter
		File? maf_file_rna_post_filtering = FilterRNAMutations.maf_file_rna_post_filtering
		File? maf_file_rna_prefiltering_with_info = FilterRNAMutations.maf_file_rna_prefiltering_with_info
		File? maf_file_rna_post_filtering_dup_GTEx = FilterRNAMutationsBasedOnDuplicateReads.maf_file_rna_post_filtering_dup_GTEx
		Boolean? maf_file_rna_post_filtering_dup_GTEx_outputted = FilterRNAMutationsBasedOnDuplicateReads.maf_file_rna_post_filtering_dup_GTEx_outputted
	}

}
