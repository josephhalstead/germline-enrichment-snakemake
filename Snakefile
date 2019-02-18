
"""

Germline Enrichment Pipeline Using GATK4

To DO:

Gaps calculations
-variants to text (VariantReporterSpark)
trusight cancer specific


"""

from pathlib import Path 


#-----------------------------------------------------------------------------------------------------------------#
# Configuration variables
#-----------------------------------------------------------------------------------------------------------------#

# Which YAMl config to use

config_location = "testvm.yaml"
configfile: config_location

# How many lanes do we have?
folder, name, lanes = glob_wildcards("{folder}/{sample_number}_{lane}_R1_001.fastq.gz")
lanes =  list((set(lanes)))

# What are our samples and sample numbers?
folder, sample_names, sample_numbers = glob_wildcards("{folder}/{sample_name}_{sample_number}_L001_R1_001.fastq.gz")

# Get other data from config file
chromosomes = config["chromosomes"]
panel = config["panel"]
seq_id = config["seqId"]

#-----------------------------------------------------------------------------------------------------------------#
# Utility Input Functions For Getting Files
#-----------------------------------------------------------------------------------------------------------------#

def get_fastqc(wildcards):
	"""	
	Function to return the fastqc file input into multiqc
	https://groups.google.com/forum/#!topic/snakemake/xxdADOSK7mY

	"""
	file_list = []
	for lane in lanes:

		for sample_name, sample_number in zip(sample_names,sample_numbers ):

			file_list.append("output/qc_reports/fastqc/" + sample_name + "_" + sample_number + "_" + lane + "_R1_001.qfilter_fastqc.html" )
			file_list.append("output/qc_reports/fastqc/" + sample_name + "_" + sample_number + "_" + lane + "_R2_001.qfilter_fastqc.html" )

	return file_list

def get_all_custom_coverage(wildcards):
	"""
	Return a list of strings with all the custom coverage files in.

	"""

	file_list = []
	beds = glob_wildcards(config["hotspot_bed_dir"] + "{bedfile}.bed")

	for sample_name, sample_number in zip(sample_names,sample_numbers ):

		for bedfile in beds[0]:

			my_input = "output/depth/hotspot_coverage/" + sample_name + "_" + sample_number + "_" + bedfile + ".coverage"
			file_list.append(my_input)

	return file_list

#-----------------------------------------------------------------------------------------------------------------#
# Starting Rule
#-----------------------------------------------------------------------------------------------------------------#

# All function pulls all rules together
rule all:
	input:
		expand("output/pipeline_finished/{seq_id}.finished" , seq_id=seq_id),

#-----------------------------------------------------------------------------------------------------------------#
# Initial QC and Preprocessing
#-----------------------------------------------------------------------------------------------------------------#

# Create a PED file for downstream analysis
rule create_ped_file:
	input:
		config_location
	output:
		"output/config/{seq_id}.ped"
	shell:
		"python scripts/make_ped.py --config {input} > {output}"

# Run the fastp program to generate a read quality report and trim read adapters and by quality
rule fastp:
	input:
		fwd = "{sample_name}/{sample_name}_{sample_number}_{lane}_R1_001.fastq.gz",
		rev = "{sample_name}/{sample_name}_{sample_number}_{lane}_R2_001.fastq.gz"
	output:
		html = "output/qc_reports/fastp/{sample_name}_{sample_number}_{lane}_fastp.html",
		json = "output/qc_reports/fastp/{sample_name}_{sample_number}_{lane}_fastp.json",
		fwd = temp("output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R1_001.qfilter.fastq.gz"),
		rev = temp("output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R2_001.qfilter.fastq.gz")
	threads:
		config["fastp_threads"]
	shell:
		"fastp -i {input.fwd} "
		"-I {input.rev} "
		"-o {output.fwd} "
		"-O {output.rev} "
		"-h {output.html} "
		"-j {output.json} "
		"--disable_quality_filtering "
		"--detect_adapter_for_pe "
		"-w {threads}"

# Check for inter-species contamination
rule fastq_screen:
	input:
		fwd = "{sample_name}/{sample_name}_{sample_number}_{lane}_R1_001.fastq.gz",
		rev = "{sample_name}/{sample_name}_{sample_number}_{lane}_R2_001.fastq.gz"
	output:
		"output/qc_reports/fastq_screen/{sample_name}_{sample_number}_{lane}_R1_001_screen.html",
		temp("output/qc_reports/fastq_screen/{sample_name}_{sample_number}_{lane}_R1_001_screen.png"),
		"output/qc_reports/fastq_screen/{sample_name}_{sample_number}_{lane}_R1_001_screen.txt",
		"output/qc_reports/fastq_screen/{sample_name}_{sample_number}_{lane}_R2_001_screen.html",
		temp("output/qc_reports/fastq_screen/{sample_name}_{sample_number}_{lane}_R2_001_screen.png"),
		"output/qc_reports/fastq_screen/{sample_name}_{sample_number}_{lane}_R2_001_screen.txt"	
	threads:
		config["fastq_screen_threads"]
	params:
		fastq_screen_config = config["fastq_screen_config"]
	shell:
		"fastq_screen "
		"--aligner bwa "
		"--threads {threads} "
		"--outdir output/qc_reports/fastq_screen/ "
		"--conf {params.fastq_screen_config} "
		"{input.fwd} "
		"{input.rev}"

# Run fastqc
rule fastqc:
	input:
		fwd = "output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R1_001.qfilter.fastq.gz",
		rev = "output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R2_001.qfilter.fastq.gz"
	output:
		"output/qc_reports/fastqc/{sample_name}_{sample_number}_{lane}_R1_001.qfilter_fastqc.html",
		"output/qc_reports/fastqc/{sample_name}_{sample_number}_{lane}_R1_001.qfilter_fastqc/summary.txt",	
		"output/qc_reports/fastqc/{sample_name}_{sample_number}_{lane}_R2_001.qfilter_fastqc.html",
		"output/qc_reports/fastqc/{sample_name}_{sample_number}_{lane}_R2_001.qfilter_fastqc/summary.txt"
	threads:
		config["fastqc_threads"]
	params:
		temp_dir = config["fastqc_temp_dir"]
	shell:
		"fastqc "
		"--threads {threads} "
		"--dir {params.temp_dir} "
		"--outdir output/qc_reports/fastqc "
		"--extract "
		"{input.fwd} "
		"{input.rev}"

#-----------------------------------------------------------------------------------------------------------------#
# Alignment and Further Preprocessing
#-----------------------------------------------------------------------------------------------------------------#

# Align reads with bwa mem, pipe into samtools and sort
rule bwa_align:
	input:
		fwd = "output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R1_001.qfilter.fastq.gz",
		rev = "output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R2_001.qfilter.fastq.gz"
	output:
		temp("output/alignments/{sample_name}_{sample_number}_{lane}.bam")
	threads:
		config["bwa_threads"]
	params:
		ref = config["bwa_reference"],
		seq_id = config["seqId"],
		centre = config["centre"],
		samtools_temp_dir = config["samtools_temp_dir"]
	shell:
		"bwa mem "
		"-t {threads} "
		"-M "
		"-R '@RG\\tID:{params.seq_id}.{params.seq_id}.{wildcards.lane}\\tCN:{params.centre}\\tSM:{wildcards.sample_name}\\tLB:{params.seq_id}\\tPL:ILLUMINA' "
		"{params.ref} {input.fwd} {input.rev} | "
		"samtools view -Sb - | "
		"samtools sort -T {params.samtools_temp_dir}.temp -O bam > {output}"

# Index the bam file
rule index_original_bam:
	input:
		"output/alignments/{sample_name}_{sample_number}_{lane}.bam"
	output:
		temp("output/alignments/{sample_name}_{sample_number}_{lane}.bam.bai")
	shell:
		"samtools index {input}"

# Merge the bams and mark duplicates
rule merge_and_remove_duplicates:
	input:
		bams = expand("output/alignments/{{sample_name}}_{{sample_number}}_{lane}.bam", lane=lanes),
		bam_indexes = expand("output/alignments/{{sample_name}}_{{sample_number}}_{lane}.bam.bai", lane=lanes),
	output:
		bam = temp("output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam"),
		index = temp("output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai"),
		metrics = "output/qc_reports/mark_duplicates/{sample_name}_{sample_number}_MarkDuplicatesMetrics.txt"
	params:
		temp = config["picard_temp_dir"],
		merge_duplicates_max_records = config["merge_duplicates_max_records"],
		files = lambda wildcards, input: " I=".join(input.bams),
		java_options = config["picard_memory_options"],
		java_home = config["java_home"]

	shell:
		"export JAVA_HOME={params.java_home}; picard {params.java_options} MarkDuplicates I={params.files} "
		"O={output.bam} "
		"METRICS_FILE={output.metrics} "
		"CREATE_INDEX=true "
		"MAX_RECORDS_IN_RAM={params.merge_duplicates_max_records} "
		"VALIDATION_STRINGENCY=SILENT "
		"TMP_DIR={params.temp} "


if config["perform_bqsr"] == True:

	# Create the BQSR report to later apply
	rule create_base_quality_report:
		input:
			bam = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
			index = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai",
		output:
			"output/bqsr_tables/{sample_name}_{sample_number}_recal_data.table"
		params:
			ref = config["reference"],
			known_sites_dbsnp = config["dbsnp_vcf"],
			known_sites_indels = config["indels_1k_vcf"],
			known_sites_gold = config["gold_standard_indels"],
			bed_file = config["capture_bed_file"],
			padding = config["interval_padding_bqsr"],
			java_options = config["gatk_bqsr_java_options"]
		shell:
			"gatk --java-options '{params.java_options}' BaseRecalibrator "
			"-I {input.bam} "
			"-O {output} "
			"-R {params.ref} "
			"--known-sites {params.known_sites_dbsnp} "
			"--known-sites {params.known_sites_indels} "
			"--known-sites {params.known_sites_gold} "
			"-L {params.bed_file} "
			"--interval-padding {params.padding}"


	# Apply BQSR Report
	rule apply_base_quality_report:
		input:
			bam = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
			bam_index = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai",
			bqsr_report = "output/bqsr_tables/{sample_name}_{sample_number}_recal_data.table"
		output:
			bam_file = "output/final_bam/{sample_name}_{sample_number}_final.bam",
			bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai"
		params:
			ref = config["reference"],
			java_options = config["gatk_bqsr_java_options"]
		shell:
			"gatk --java-options '{params.java_options}' "
			" ApplyBQSR -R {params.ref} "
			" -I {input.bam} "
			"-bqsr {input.bqsr_report} "
			"-O {output.bam_file} "

	# Second pass to get table for AnalyseCovariates
	# Do we really need this?
	rule create_base_quality_report_second_pass:
		input:
			bam = "output/final_bam/{sample_name}_{sample_number}_final.bam",
			bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai"
		output:
			"output/bqsr_tables/{sample_name}_{sample_number}_post_recal_data.table"
		params:
			ref = config["reference"],
			known_sites_dbsnp = config["dbsnp_vcf"],
			known_sites_indels = config["indels_1k_vcf"],
			known_sites_gold = config["gold_standard_indels"],
			bed_file = config["capture_bed_file"],
			padding = config["interval_padding_bqsr"],
			java_options = config["gatk_bqsr_java_options"]
		shell:
			"gatk --java-options '{params.java_options}' "
			"BaseRecalibrator "
			"-I {input.bam} "
			"-O {output} "
			"-R {params.ref} "
			"--known-sites {params.known_sites_dbsnp} "
			"--known-sites {params.known_sites_indels} "
			"--known-sites {params.known_sites_gold} "
			"-L {params.bed_file} "
			"--interval-padding {params.padding}"

	# Analyse co-variates before and after BQSR
	rule analyse_covariates:
		input:
			before = "output/bqsr_tables/{sample_name}_{sample_number}_recal_data.table",
			after = "output/bqsr_tables/{sample_name}_{sample_number}_post_recal_data.table"
		output:
			pdf = "output/qc_reports/bqsr/{sample_name}_{sample_number}_bqsr_covariation.pdf",
			csv = "output/qc_reports/bqsr/{sample_name}_{sample_number}_bqsr_covariation.csv",
		params:
			java_options = config["gatk_analyse_covariates_java_options"]
		shell:
			"gatk --java-options '{params.java_options}' "
			"AnalyzeCovariates "
			"-before {input.before} "
			"-after {input.after} "
			"-csv {output.csv} "
			"-plots {output.pdf} "

else:

	rule move_final_bam:
		input:
			bam = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
			bam_index = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai",
		output:
			bam = "output/final_bam/{sample_name}_{sample_number}_final.bam",
			bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai"
		shell:
			"cp {input.bam} {output.bam}; cp {input.bam_index} {output.bam_index} "



#-----------------------------------------------------------------------------------------------------------------#
# Post Alignment QC
#-----------------------------------------------------------------------------------------------------------------#


# Create an interval file from the BED file for use in Picard tools such as CollectHsMetrics
rule create_interval_file:
	input:
		config["capture_bed_file"]
	output:
		temp("output/config/" + Path(config["capture_bed_file"]).name.split(".")[0] + ".interval_list")
	params:
		sequence_dict = config["reference_sequence_dict"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard BedToIntervalList I={input} O={output} SD={params.sequence_dict}" 

# Collect the insert size metrics using picard
rule collect_insert_size_metrics:
	input:
		bam = "output/final_bam/{sample_name}_{sample_number}_final.bam",
		bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai"
	output:
		txt="output/qc_reports/insert_size_metrics/{sample_name}_{sample_number}_InsertSizeMetrics.txt",
		pdf="output/qc_reports/insert_size_metrics/{sample_name}_{sample_number}_InsertSizeMetrics.pdf"
	params:
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard CollectInsertSizeMetrics "
		"I={input.bam} "
		"O={output.txt} "
		"HISTOGRAM_FILE={output.pdf} "


# Collect the HS metrics using picard
rule collect_hs_metrics:
	input:
		bam = "output/final_bam/{sample_name}_{sample_number}_final.bam",
		bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai",
		intervals = "output/config/" + Path(config["capture_bed_file"]).name.split(".")[0] + ".interval_list"
	output:
		"output/qc_reports/hs_metrics/{sample_name}_{sample_number}_HsMetrics.txt"
	params:
		ref = config["reference"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard CollectHsMetrics "
		"I={input.bam} "
		"O={output} "
		"R={params.ref} "
		"BAIT_INTERVALS={input.intervals} "
		"TARGET_INTERVALS={input.intervals}"

# Collect alignment summary metrics using picard
rule collect_alignment_metrics:
	input:
		bam = "output/final_bam/{sample_name}_{sample_number}_final.bam",
		bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai"
	output:
		"output/qc_reports/alignment_metrics/{sample_name}_{sample_number}_AlignmentSummaryMetrics.txt"
	params:
		ref = config["reference"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard CollectAlignmentSummaryMetrics "
		"I={input.bam} "
		"O={output} "
		"R={params.ref}"

# Extract high confidence SNPs for contamination analysis
rule extract_high_confidence_snps:
	input:
		config["1000k_high_confidence_snps"]
	output:
		vcf = temp("output/config/snps/1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf"),
		index = temp("output/config/snps/1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf.idx")
	params:
		ref = config["reference"],
		bed = config["capture_bed_file"],
		java_options = config["gatk_variants_java_options"]
	shell:
		"gatk --java-options '{params.java_options}' "
		"SelectVariants "
		"-R {params.ref} "
		"-V {input} "
		"-L {params.bed} "
		"-O {output.vcf} "
		"-select-type SNP "
		"--restrict-alleles-to BIALLELIC "
		"--exclude-intervals X "
		"--exclude-intervals Y "
		"--exclude-intervals MT "
		"--exclude-non-variants "
		"--exclude-filtered "

# Check for contamination using verify bam id
rule verify_bam_id:
	input:
		bam = "output/final_bam/{sample_name}_{sample_number}_final.bam",
		bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai",
		vcf = "output/config/snps/1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf"
	output:
		depthSM = "output/qc_reports/verify_bam_id/{sample_name}_{sample_number}_verify_bam_id.depthSM",
		log = "output/qc_reports/verify_bam_id/{sample_name}_{sample_number}_verify_bam_id.log",
		selfSM = "output/qc_reports/verify_bam_id/{sample_name}_{sample_number}_verify_bam_id.selfSM"
	shell:
		"verifyBamID "
		"--vcf {input.vcf} "
		"--bam {input.bam} "
		"--out output/qc_reports/verify_bam_id/{wildcards.sample_name}_{wildcards.sample_number}_verify_bam_id "
		"--verbose "
		"--ignoreRG "
		"--chip-none "
		"--minMapQ 20 "
		"--maxDepth 1000 "
		"--precise "

# Calculate per base coverage - replace with sambamba or mosdepth for speed?
rule generate_per_base_coverage:
	input:
		bam = "output/final_bam/{sample_name}_{sample_number}_final.bam",
		bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai"
	output:
		depth = temp("output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage"),
		summary = "output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.sample_summary"
	params:
		ref = config["reference"],
		bed = config["capture_bed_file"],
		min_coverage = config["minimum_coverage"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; gatk3 -T DepthOfCoverage "
		"-R {params.ref} "
		"-L {params.bed} "
		"-I {input.bam} "
		"-o {output.depth} " 
		"--countType COUNT_FRAGMENTS "
		"--minMappingQuality 20 "
		"--minBaseQuality 10 "
		"-ct {params.min_coverage} "
		"--omitLocusTable "
		"-rf MappingQualityUnavailable "
		"-dt NONE"

# Create a compressed and indexed version of the depth file
rule compress_and_index_coverage:
	input:
		"output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage"
	output:
		depth = "output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.gz",
		index = "output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.gz.tbi"
	shell:
		"sed 's/:/\t/g' {input} | grep -v 'Locus' | sort -k1,1 -k2,2n | bgzip > {output.depth} && "
		"tabix -b 2 -e 2 -s 1 {output.depth} "


#-----------------------------------------------------------------------------------------------------------------#
# Find Coverage Gaps
#-----------------------------------------------------------------------------------------------------------------#

rule calculate_coverage_metrics:
	input:
		depth = "output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.gz",
		index = "output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.gz.tbi"	
	output:
		"output/depth/metrics/" + seq_id + "_{sample_name}_{sample_number}_Gaps.bed"
	params:
		roi_bed = config["capture_bed_file"],
		panel = panel,
		hotspots = config["germline_hotspots"],
		min_coverage = config["minimum_coverage"],
		seq_id = seq_id,
		refseq = config["refseq"]
	shell:
		"bash scripts/get_coverage.sh "
		"{params.roi_bed} "
		"{params.panel} "
		"{params.hotspots} "
		"{input.depth} "
		"{params.seq_id} "
		"{wildcards.sample_name}_{wildcards.sample_number} "
		"{params.min_coverage} "
		"{params.refseq} "
		"output/depth/metrics/ "
	
#-----------------------------------------------------------------------------------------------------------------#
# SNP and Small Indel Calling with GATK Haplotype Caller
#-----------------------------------------------------------------------------------------------------------------#

# Sort ROI bed for splitting by bedextract
rule sort_capture_bed:
	input:
		config["capture_bed_file"]
	output:
		temp("output/config/sorted_beds/{{panel}}_sorted.bed".format(panel=panel))
	shell:
		"sort-bed {input} > {output}"

# Split the bed by chromosome for input into create_gvcfs
rule split_bed_by_chromosome:
	input:
		"output/config/sorted_beds/{panel}_sorted.bed".format(panel=panel)
	output:
		temp(expand("output/config/split_capture_bed/{chr}.bed", chr=chromosomes))
	params:
		chromosomes = chromosomes
	shell:
		"for chr in {params.chromosomes}; do bedextract $chr {input} > output/config/split_capture_bed/$chr.bed; done"

# Create GVCF using Haplotype Caller for each sample chromosome combination
rule create_gvcfs:
	input:
		bam_file = "output/final_bam/{sample_name}_{sample_number}_final.bam",
		bam_index= "output/final_bam/{sample_name}_{sample_number}_final.bai",
		bed = "output/config/split_capture_bed/{chr}.bed",
	output:
		gvcf_file = temp("output/gvcfs/{sample_name}_{sample_number}_chr{chr}.g.vcf"),
		index = temp("output/gvcfs/{sample_name}_{sample_number}_chr{chr}.g.vcf.idx")
	params:
		ref = config["reference"],
		padding = config['interval_padding_haplotype_caller'],
		java_options = config['gatk_hc_java_options']
	shell:
		"gatk --java-options '{params.java_options}' HaplotypeCaller -R {params.ref} "
		"-I {input.bam_file} "
		"--emit-ref-confidence GVCF "
		"-O {output.gvcf_file} "
		"-L {input.bed} "
		"--interval-padding {params.padding}"


# Consolidate all samples into a genomics db for joint genotyping
rule create_genomics_db:
	input:
		gvcfs = expand("output/gvcfs/{sample_name}_{sample_number}_chr{{chr}}.g.vcf" , zip, sample_name=sample_names, sample_number=sample_numbers),
		index = expand("output/gvcfs/{sample_name}_{sample_number}_chr{{chr}}.g.vcf.idx", zip, sample_name=sample_names, sample_number=sample_numbers),
		bed = "output/config/split_capture_bed/{chr}.bed"
	output:
		temp(directory("output/genomicdbs/{seq_id}_chr{chr}"))
	params:
		files = lambda wildcards, input: " -V ".join(input.gvcfs),
		java_options = config["gatk_genomics_db_java_options"],
		padding = config['interval_padding_haplotype_caller']
	shell:
		"gatk --java-options '{params.java_options}' "
		" GenomicsDBImport -V {params.files} "
		"--genomicsdb-workspace-path {output} "
		"-L {input.bed} "
		"--interval-padding {params.padding}"


# Genotype the gvcfs and produce a joint vcf
rule genotype_gvcfs:
	input:
		db = "output/genomicdbs/{seq_id}_chr{chr}",
		bed = "output/config/split_capture_bed/{chr}.bed",
	output:
		vcf = temp("output/jointvcf_per_chr/{seq_id}_chr{chr}.vcf"),
		index = temp("output/jointvcf_per_chr/{seq_id}_chr{chr}.vcf.idx")
	params:
		ref = config["reference"],
		java_options = config['gatk_hc_java_options'],
		padding = config['interval_padding_haplotype_caller']
	shell:
		"gatk --java-options '{params.java_options}'  GenotypeGVCFs -R {params.ref} "
		"-V gendb://{input.db} "
		"-G StandardAnnotation "
		"-O {output.vcf} "
		"-L {input.bed} "
		"--interval-padding {params.padding} "

# Combine the chromsome vcfs into one final vcf with all samples and all chromosomes
rule collect_vcfs:
	input:
		vcf = expand("output/jointvcf_per_chr/{{seq_id}}_chr{chr}.vcf", chr= chromosomes),
		index = expand("output/jointvcf_per_chr/{{seq_id}}_chr{chr}.vcf.idx", chr= chromosomes)
	output:
		vcf = temp("output/jointvcf/{seq_id}_all_chr.vcf"),
		index = temp("output/jointvcf/{seq_id}_all_chr.vcf.idx")
	params:
		files = lambda wildcards, input: " I=".join(input.vcf),
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard GatherVcfs "
		"I={params.files} "
		"O={output.vcf}"


#-----------------------------------------------------------------------------------------------------------------#
# Filter Variants on Quality (Hard Filtering)
#-----------------------------------------------------------------------------------------------------------------#

# Select just the SNPs
rule select_snps_for_filtering:
	input:
		vcf = "output/jointvcf/{seq_id}_all_chr.vcf",
		index = "output/jointvcf/{seq_id}_all_chr.vcf.idx"
	output:
		vcf = temp("output/jointvcf_snps/{seq_id}_all_chr_snps.vcf"),
		index = temp("output/jointvcf_snps/{seq_id}_all_chr_snps.vcf.idx")
	params:
		ref = config["reference"],
		padding = config["interval_padding_haplotype_caller"],
		bed = config["capture_bed_file"],
		java_options = config["gatk_variants_java_options"]
	shell:
		"gatk --java-options '{params.java_options}' SelectVariants "
		"-R {params.ref} "
		"-V {input.vcf} "
		"-select-type SNP "
		"-L {params.bed} "
		"--interval-padding {params.padding} "
		"-O {output.vcf} "

# Filter the SNPs on quality
rule filter_snps:
	input:
		vcf = "output/jointvcf_snps/{seq_id}_all_chr_snps.vcf",
		index = "output/jointvcf_snps/{seq_id}_all_chr_snps.vcf.idx"
	output:
		vcf = temp("output/jointvcf_snps_filtered/{seq_id}_all_chr_snps_filtered.vcf"),
		index = temp("output/jointvcf_snps_filtered/{seq_id}_all_chr_snps_filtered.vcf.idx")
	params:
		ref = config["reference"],
		padding = config["interval_padding_haplotype_caller"],
		bed = config["capture_bed_file"],
		min_qual = ["snp_min_qual"],
		min_QD = ["snp_min_QD"],
		max_FS = ["snp_max_FS"],
		min_MQ = ["snp_min_MQ"],
		min_MQRankSum = ["snp_min_MQRankSum"],
		min_ReadPosRankSum = ["snp_min_ReadPosRankSum"],
		java_options = config["gatk_variants_java_options"]
	shell:
		"gatk --java-options '{params.java_options}' "
		"VariantFiltration "
		"-R {params.ref} "
		"-V {input.vcf} "
		"-L {params.bed} "
		"--interval-padding {params.padding} "
		"-O {output.vcf} "
		"--filter-expression 'QUAL < {params.min_qual}' "
		"--filter-name 'LowQual' "
		"--filter-expression 'QD <  {params.min_QD}' "
		"--filter-name 'QD' "
		"--filter-expression 'FS >  {params.max_FS}' "
		"--filter-name 'FS' "
		"--filter-expression 'MQ <  {params.min_MQ}' "
		"--filter-name 'MQ' "
		"--filter-expression 'MQRankSum <  {params.min_MQRankSum}' "
		"--filter-name 'MQRankSum' "
		"--filter-expression 'ReadPosRankSum <  {params.min_ReadPosRankSum}' "
		"--filter-name 'ReadPosRankSum' "

# Select all non SNPs
rule select_non_snps_for_filtering:
	input:
		vcf = "output/jointvcf/{seq_id}_all_chr.vcf",
		index = "output/jointvcf/{seq_id}_all_chr.vcf.idx",
	output:
		vcf = temp("output/jointvcf_indels/{seq_id}_all_chr_indels.vcf"),
		index = temp("output/jointvcf_indels/{seq_id}_all_chr_indels.vcf.idx")
	params:
		ref = config["reference"],
		padding = config["interval_padding_haplotype_caller"],
		bed = config["capture_bed_file"],
		java_options = config["gatk_variants_java_options"]
	shell:
		"gatk --java-options '{params.java_options}' SelectVariants "
		"-R {params.ref} "
		"-V {input.vcf} "
		"--select-type-to-exclude SNP "
		"-L {params.bed} "
		"--interval-padding {params.padding} "
		"-O {output.vcf} "

# Filter the non SNPs
rule filter_non_snps:
	input:
		vcf = "output/jointvcf_indels/{seq_id}_all_chr_indels.vcf",
		index = "output/jointvcf_indels/{seq_id}_all_chr_indels.vcf.idx"
	output:
		vcf = temp("output/jointvcf_indels_filtered/{seq_id}_all_chr_indels_filtered.vcf"),
		index = temp("output/jointvcf_indels_filtered/{seq_id}_all_chr_indels_filtered.vcf.idx")
	params:
		ref = config["reference"],
		padding = config["interval_padding_haplotype_caller"],
		bed = config["capture_bed_file"],
		min_qual = config["indel_min_qual"],
		min_QD = config["indel_min_QD"],
		max_FS = config["indel_max_FS"],
		max_SOR = config["indel_max_SOR"],
		min_ReadPosRankSum = config["indel_min_ReadPosRankSum"],
		min_inbreeding_coeff = config["indel_min_inbreeding_coeff"],
		java_options = config["gatk_variants_java_options"]
	shell:
		"gatk --java-options '{params.java_options}' "
		"VariantFiltration "
		"-R {params.ref} "
		"-V {input.vcf} "
		"-L {params.bed} "
		"--interval-padding {params.padding} "
		"-O {output.vcf} "
		"--filter-expression 'QUAL < {params.min_qual}' "
		"--filter-name 'LowQual' "
		"--filter-expression 'QD <  {params.min_QD}' "
		"--filter-name 'QD' "
		"--filter-expression 'FS >  {params.max_FS}' "
		"--filter-name 'FS' "
		"--filter-expression 'SOR >  {params.max_SOR}' "
		"--filter-name 'MQ' "
		"--filter-expression 'ReadPosRankSum <  {params.min_ReadPosRankSum}' "
		"--filter-name 'ReadPosRankSum' "
		"--filter-expression 'InbreedingCoeff != 'NaN' && InbreedingCoeff < {params.min_inbreeding_coeff}' "
		"--filter-name 'ReadPosRankSum' "

# Combine the filtered SNPs and Indels into a single VCF
rule combine_filtered_snps_and_indels:
	input:
		snps = "output/jointvcf_snps_filtered/{seq_id}_all_chr_snps_filtered.vcf",
		indels = "output/jointvcf_indels_filtered/{seq_id}_all_chr_indels_filtered.vcf"
	output:
		vcf = temp("output/jointvcf_all_variants_filtered/{seq_id}_all_variants_filtered.vcf"),
		index = temp("output/jointvcf_all_variants_filtered/{seq_id}_all_variants_filtered.vcf.idx")
	params:
		java_options = config["gatk_variants_java_options"]
	shell:
		"gatk --java-options '{params.java_options}' "
		"MergeVcfs "
		"-I {input.snps} "
		"-I {input.indels} "
		"-O {output.vcf}"		

# Filter on genotype quality
rule filter_genotypes:
	input:
		vcf = "output/jointvcf_all_variants_filtered/{seq_id}_all_variants_filtered.vcf",
		index = "output/jointvcf_all_variants_filtered/{seq_id}_all_variants_filtered.vcf.idx"
	output:
		vcf = "output/jointvcf_all_variants_filtered_genotype/{seq_id}_all_variants_filtered_genotype.vcf",
		index = temp("output/jointvcf_all_variants_filtered_genotype/{seq_id}_all_variants_filtered_genotype.vcf.idx")
	params:
		ref = config["reference"],
		padding = config["interval_padding_haplotype_caller"],
		bed = config["capture_bed_file"],
		min_genotype_depth = config["min_genotype_depth"],
		java_options = config["gatk_variants_java_options"]
	shell:
		"gatk --java-options '{params.java_options}' "
		"VariantFiltration "
		"-R {params.ref} "
		"-V {input.vcf} "
		"-L {params.bed} "
		"--interval-padding {params.padding} "
		"-O {output.vcf} "
		"--genotype-filter-expression 'DP < {params.min_genotype_depth}' "
		"--genotype-filter-name 'LowDP' "


# Compress and index the filtered vcf
rule compress_and_index_filtered_vcf:
	input:
		"output/jointvcf_all_variants_filtered_genotype/{seq_id}_all_variants_filtered_genotype.vcf"
	output:
		"output/jointvcf_all_variants_filtered_genotype/{seq_id}_all_variants_filtered_genotype.vcf.gz",
		"output/jointvcf_all_variants_filtered_genotype/{seq_id}_all_variants_filtered_genotype.vcf.gz.tbi"
	shell:
		"bgzip {input} && tabix {input}.gz"	


# Filter out variants outside of ROI
rule filter_by_roi:
	input:
		vcf = "output/jointvcf_all_variants_filtered_genotype/{seq_id}_all_variants_filtered_genotype.vcf.gz",
		index = "output/jointvcf_all_variants_filtered_genotype/{seq_id}_all_variants_filtered_genotype.vcf.gz.tbi"
	output:
		temp("output/jointvcf_all_variants_filtered_genotype_roi/{seq_id}_all_variants_filtered_genotype_roi.vcf")
	params:
		bed = config["capture_bed_file"],
		ref = config["reference"]
	shell:
		"bcftools view -R {params.bed} {input.vcf} | "
		"vt normalize -r {params.ref} - > {output} "

# Use vt to split multiallelics and normalise variants
rule decompose_and_normalise:
	input:
		"output/jointvcf_all_variants_filtered_genotype_roi/{seq_id}_all_variants_filtered_genotype_roi.vcf"
	output:
		temp("output/jointvcf_all_variants_filtered_genotype_roi_norm/{seq_id}_all_variants_filtered_genotype_roi_norm.vcf")
	params:
		ref = config["reference"]
	shell:
		"cat {input} | "
		"vt decompose -s - | "
		"vt normalize -r {params.ref} - > {output}"

	
#-----------------------------------------------------------------------------------------------------------------#
# Variant Annotation
#-----------------------------------------------------------------------------------------------------------------#


# Annotate the vcf using VEP
rule annotate_vep:
	input:
		"output/jointvcf_all_variants_filtered_genotype_roi_norm/{seq_id}_all_variants_filtered_genotype_roi_norm.vcf"
	output:
		"output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf"
	params:
		vep_cache = config["vep_cache_location"],
		ref = config["reference"],
		gnomad_genomes = config["gnomad_genomes"],
		gnomad_exomes = config["gnomad_exomes"],
		ccr_bed = config["ccr_bed"],
		spliceai = config["spliceai"]
	threads:
		config["vep_threads"]
	shell:
		"vep --verbose "
		"--format vcf "
		"--everything "
		"--fork {threads} "
		"--species homo_sapiens "
		"--assembly GRCh37  "
		"--input_file {input}  "
		"--output_file {output} "
		"--force_overwrite "
		"--cache "
		"--dir  {params.vep_cache} "
		"--fasta {params.ref} "
		"--offline "
		"--cache_version 94 "
		"--no_escape "
		"--shift_hgvs 1 "
		"--vcf "
		"--refseq "
		"--flag_pick "
		"--pick_order biotype,canonical,appris,tsl,ccds,rank,length "
		"--exclude_predicted "
		"--custom {params.gnomad_genomes},gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_POPMAX "
		"--custom {params.gnomad_exomes},gnomADe,vcf,exact,0,AF_POPMAX "
		"--custom {params.ccr_bed},ccrs,bed,overlap,0 "
		"--custom {params.spliceai},SpliceAI,vcf,exact,0,DS_AG,DS_AL,DS_DG,DS_DL,SYMBOL "

rule compress_and_index_vep:
	input:
		"output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf"
	output:
		"output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz",
		"output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz.tbi"
	shell:
		"bgzip {input} && tabix {input}.gz"	

# Convert vcf to csv
rule convert_to_csv:
	input:
		vcf = "output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz",
		index = "output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz.tbi"
	output:
		temp("output/vcf_csv/{seq_id}_vcf.csv")
	shell:
		"gatk VariantsToTable -V {input.vcf} "
		"-O {output} -F CHROM -F POS -F REF -F ALT -F ID -F QUAL -F FILTER -F CSQ -F AC "
		"-GF GT -GF GQ -GF DP"

# Get the CSQ string which describes the VEP fields	
rule get_csq_string:
	input:
		vcf = "output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz",
		index = "output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz.tbi"
	output:
		temp("output/config/csq/{seq_id}_csq.txt")
	shell:
		"zcat {input.vcf} | grep \"^##INFO=<ID=CSQ\" | awk 'BEGIN {{ FS = \":\" }} ; {{ print $2 }}' | tr -d '>\" ' > {output} "

# Run the germline variant filter program
rule create_variant_reports:
	input:
		csv = "output/vcf_csv/{seq_id}_vcf.csv",
		ped = "output/config/{seq_id}.ped",
		csq = "output/config/csq/{seq_id}_csq.txt"
	output:
		"output/variant_reports/{seq_id}_finished.txt"
	params:
		germline_variant_filter = config["germline_variant_filter"],
		filter_config = config["filter_config"],
		local_panel_app_dump = config["local_panel_app_dump"],
	shell:
		"python {params.germline_variant_filter} "
		"--config {params.filter_config} "
		"--ped {input.ped} "
		"--input {input.csv} "
		"--panelapp "
		"--local-panel-app-dump {params.local_panel_app_dump} "
		"--spliceai "
		"--smart-synonymous "
		"--gnomad-constraint-scores "
		"--worksheet {wildcards.seq_id} "
		"--results-dir output/variant_reports/ "
		"--csq $(cat {input.csq}) "
		"--add-ccrs && touch {output}"

# Merge the metadata from all samples into a single file
rule collect_meta_data:
	input:
		meta = expand("output/vcf_meta/{sample_name}_{sample_number}_meta.txt", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		meta = temp("output/vcf_meta/merged_meta.txt")
	shell:
		"cat {input.meta} > {output.meta}"


# Add metadata to vcf header for variant database import
rule add_meta_to_vcf:
	input:
		vcf = "output/jointvcf_all_variants_filtered_genotype_roi/{seq_id}_all_variants_filtered_genotype_roi.vcf",
		meta = "output/vcf_meta/merged_meta.txt"
	output:
		temp("output/jointvcf_all_variants_filtered_genotype_roi_meta/{seq_id}_all_variants_filtered_genotype_roi_meta.vcf")
	shell:
		"""
		 grep '^##' {input.vcf} > {output}
		 cat {input.meta} >> {output}
		 grep -v '^##' {input.vcf} >> {output}
		"""

# Exclude Mitochrondial variants for inclusion in variant database
rule create_vcf_without_mt:
	input:
		"output/jointvcf_all_variants_filtered_genotype_roi_meta/{seq_id}_all_variants_filtered_genotype_roi_meta.vcf"
	output:
		"output/jointvcf_all_variants_filtered_genotype_roi_meta_nomt/{seq_id}_all_variants_filtered_genotype_roi_meta_nomt.vcf"
	shell:
		"awk '$1 !~ /^MT/ {{ print $0 }}' {input} > {output}"


# Check the final vcf in valid
rule validate_vcf:
	input:
		vcf = "output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz",
		index = "output/jointvcf_all_variants_filtered_genotype_roi_norm_vep/{seq_id}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz.tbi"
	output:
		"output/validated_vcf/{seq_id}.validated"
	params:
		ref = config["reference"]
	shell:
		"gatk ValidateVariants "
		"-V {input.vcf} "
		"-R {params.ref} "
		"&& touch {output}"


#-----------------------------------------------------------------------------------------------------------------#
# Final QC checks
#-----------------------------------------------------------------------------------------------------------------#

# Evaluate variant calling - need sex and rawSequenceQuality inputs
rule variant_evaluation:
	input:
		"output/jointvcf_all_variants_filtered_genotype/{seq_id}_all_variants_filtered_genotype.vcf.gz"
	output:
		"output/qc_reports/variant_calling_metrics/{seq_id}_CollectVariantCallingMetrics.variant_calling_detail_metrics",
		"output/qc_reports/variant_calling_metrics/{seq_id}_CollectVariantCallingMetrics.variant_calling_summary_metrics"
	params:
		dbsnp_vcf = config["dbsnp_vcf_eval"],
		java_home = config["java_home"]
	threads:
		config["variant_evaluation_threads"],
	shell:
		"export JAVA_HOME={params.java_home}; picard CollectVariantCallingMetrics "
		"I={input} "
		"O=output/qc_reports/variant_calling_metrics/{wildcards.seq_id}_CollectVariantCallingMetrics "
		"DBSNP={params.dbsnp_vcf} "
		"THREAD_COUNT={threads}"

# Use the Y chromosome coverage to calculate the sex
rule calculate_sex:
	input:
		depth_summary = "output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.sample_summary",
		bam = "output/final_bam/{sample_name}_{sample_number}_final.bam",
		bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai"
	output:
		sex = "output/qc_reports/sex/{sample_name}_{sample_number}_sex.txt",
		y_bed = temp("output/config/y_bed/{sample_name}_{sample_number}_y.bed"),
		y_cov = temp("output/depth/y_coverage/{sample_name}_{sample_number}_DepthOfCoverage.sample_summary"),
		y_cov_stats = temp("output/depth/y_coverage/{sample_name}_{sample_number}_DepthOfCoverage.sample_statistics")
	params:
		bed = config["capture_bed_file"],
		ref = config["reference"],
		output_name = "output/depth/y_coverage/{sample_name}_{sample_number}_DepthOfCoverage",
		java_home = config["java_home"]
	shell:
		"""
		#sex analysis using Y chrom coverage

		export JAVA_HOME={params.java_home};
		awk '{{if ($1 == "Y") print $0}}' {params.bed} > {output.y_bed}

		meanOnTargetCoverage=$(head -n2 {input.depth_summary} | tail -n1 | cut -s -f3)

		if [ $(wc -l {output.y_bed} | cut -d' ' -f1) -gt 0 ] && [ $(awk -v meanOnTargetCoverage="$meanOnTargetCoverage" 'BEGIN{{printf "%3.0f", meanOnTargetCoverage}}') -gt 10 ]; then

		    #calc Y coverage
		    gatk3 \
		    -T DepthOfCoverage \
		    -R {params.ref} \
		    -o {params.output_name} \
		    --omitDepthOutputAtEachBase \
		    --omitIntervalStatistics \
		    --omitLocusTable \
		    -L {output.y_bed} \
		    -XL Y:10000-2649520 \
		    -XL Y:59034049-59363566 \
		    -I {input.bam} \
		    --countType COUNT_FRAGMENTS \
		    --minMappingQuality 20 \
		    -rf MappingQualityUnavailable \
		    -dt NONE

		    #extract Y mean coverage
		    meanYCov=$(head -n2 {output.y_cov} | tail -n1 | cut -s -f3)
		    calcSex=$(awk -v meanOnTargetCoverage="$meanOnTargetCoverage" -v meanYCov="$meanYCov" 'BEGIN {{if (meanYCov > 10 && (meanYCov / meanOnTargetCoverage) > 0.1){{print "MALE"}} else if (meanYCov < 10 && (meanYCov / meanOnTargetCoverage) < 0.1) {{print "FEMALE" }} else {{print "UNKNOWN"}} }}')

		else
		    calcSex="UNKNOWN"
		    touch {output.y_cov}
		fi


		echo $calcSex > {output.sex}

		"""

# Gather the QC metrics into a single file
rule gather_qc_metrics:
	input:
		bam = "output/final_bam/{sample_name}_{sample_number}_final.bam",
		bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai",
		insert_size_metrics = "output/qc_reports/insert_size_metrics/{sample_name}_{sample_number}_InsertSizeMetrics.txt",
		duplicate_metrics = "output/qc_reports/mark_duplicates/{sample_name}_{sample_number}_MarkDuplicatesMetrics.txt",
		hs_metrics = "output/qc_reports/hs_metrics/{sample_name}_{sample_number}_HsMetrics.txt",
		depth_summary = "output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.sample_summary",
		alignments_summary = "output/qc_reports/alignment_metrics/{sample_name}_{sample_number}_AlignmentSummaryMetrics.txt",
		contamination = "output/qc_reports/verify_bam_id/{sample_name}_{sample_number}_verify_bam_id.selfSM",
		fastqc_fwd = expand("output/qc_reports/fastqc/{{sample_name}}_{{sample_number}}_{lane}_R1_001.qfilter_fastqc/summary.txt", lane =lanes),
		fastqc_rev = expand("output/qc_reports/fastqc/{{sample_name}}_{{sample_number}}_{lane}_R2_001.qfilter_fastqc/summary.txt", lane =lanes),
		sex = "output/qc_reports/sex/{sample_name}_{sample_number}_sex.txt"
	output:
		summary = "{sample_name}/{sample_name}_{sample_number}_QC.txt",
		high_coverage = "output/qc_reports/high_coverage_bams/{sample_name}_{sample_number}_high_coverage.txt",
		low_coverage = "output/qc_reports/low_coverage_bams/{sample_name}_{sample_number}_low_coverage.txt",
		meta = "output/vcf_meta/{sample_name}_{sample_number}_meta.txt"
	params:
		seq_id = seq_id,
		worksheet_id  = lambda w: config["sample_info"][w.sample_name]["worklistId"],
		panel = panel,
		pipeline_name = config["pipelineName"],
		pipeline_version = config["pipelineVersion"],

	shell:
		"""
		countQCFlagFails() {{
		grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
		grep -v ^PASS | \
		grep -v ^WARN | \
		wc -l | \
		sed 's/^[[:space:]]*//g'
		}}	

		rawSequenceQuality=PASS

		for report in  {input.fastqc_fwd} {input.fastqc_rev}; do

			 if [ $(countQCFlagFails $report) -gt 0 ]; then
				rawSequenceQuality=FAIL
			 fi

		done

		meanInsertSize=$(head -n8 {input.insert_size_metrics} | tail -n1 | cut -s -f5) #mean insert size
		sdInsertSize=$(head -n8 {input.insert_size_metrics}  | tail -n1 | cut -s -f6) #insert size standard deviation
		duplicationRate=$(head -n8 {input.duplicate_metrics}  | tail -n1 | cut -s -f9) #The percentage of mapped sequence that is marked as duplicate.
		totalReads=$(head -n8 {input.hs_metrics} | tail -n1 | cut -s -f6) #The total number of reads in the SAM or BAM file examine.
		pctSelectedBases=$(head -n8 {input.hs_metrics} | tail -n1 | cut -s -f19) #On+Near Bait Bases / PF Bases Aligned.
		totalTargetedUsableBases=$(head -n2 {input.depth_summary} | tail -n1 | cut -s -f2) #total number of usable bases. NB BQSR requires >= 100M, ideally >= 1B
		meanOnTargetCoverage=$(head -n2 {input.depth_summary} | tail -n1 | cut -s -f3) #avg usable coverage
		pctTargetBasesCt=$(head -n2 {input.depth_summary} | tail -n1 | cut -s -f7) #percentage panel covered with good enough data for variant detection
		freemix=$(tail -n1 {input.contamination} | cut -s -f7) #percentage DNA contamination. Should be <= 0.02
		pctPfReadsAligned=$(grep ^PAIR {input.alignments_summary} | awk '{{print $7*100}}') #Percentage mapped reads
		atDropout=$(head -n8 {input.hs_metrics} | tail -n1 | cut -s -f51) #A measure of how undercovered <= 50% GC regions are relative to the mean
		gcDropout=$(head -n8 {input.hs_metrics} | tail -n1 | cut -s -f52) #A measure of how undercovered >= 50% GC regions are relative to the mean
		calcSex=$(head -n1 {input.sex})

		echo -e "TotalReads\tRawSequenceQuality\tTotalTargetUsableBases\tDuplicationRate\tPctSelectedBases\tPctTargetBasesCt\tMeanOnTargetCoverage\tSex\tEstimatedContamination\tMeanInsertSize\tSDInsertSize\tPercentMapped\tAtDropout\tGcDropout" > {output.summary}
		echo -e "$totalReads\t$rawSequenceQuality\t$totalTargetedUsableBases\t$duplicationRate\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage\t$calcSex\t$freemix\t$meanInsertSize\t$sdInsertSize\t$pctPfReadsAligned\t$atDropout\t$gcDropout" >> {output.summary}

		if [ $(echo "$meanOnTargetCoverage" | awk '{{if ($1 > 20) print "true"; else print "false"}}') = true ]; then
			echo {input.bam} > {output.high_coverage} 
		else

			touch {output.high_coverage} 

		fi

		if [ $(echo "$meanOnTargetCoverage" | awk '{{if ($1 < 20) print "true"; else print "false"}}') = true ]; then
			echo {input.bam} > {output.low_coverage} 
		else

			touch {output.low_coverage} 

		fi		

		RemoteVcfFilePath=$(dirname $PWD)/output/jointvcf_all_variants_filtered_genotype_vep_gnomad_roi_meta/{params.seq_id}_jointvcf_all_variants_filtered_genotype_vep_gnomad_roi_meta.vcf

		RemoteBamFilePath=$(dirname $PWD)/{input.bam}

		# print vcf meta data needed for db import
		echo \#\#SAMPLE\=\<ID\="{wildcards.sample_name}",Tissue\=Germline,WorklistId\={params.worksheet_id},SeqId\={params.seq_id},,Assay\={params.panel},PipelineName\={params.pipeline_name},PipelineVersion\={params.pipeline_version},RawSequenceQuality\="$rawSequenceQuality",PercentMapped\="$pctPfReadsAligned",ATDropout\="$atDropout",GCDropout\="$gcDropout",MeanInsertSize\="$meanInsertSize",SDInsertSize\="$sdInsertSize",DuplicationRate\="$duplicationRate",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",PctTargetBasesCt\="$pctTargetBasesCt",EstimatedContamination\="$freemix",GenotypicGender\="$calcSex",TotalTargetedUsableBases\="$totalTargetedUsableBases",RemoteVcfFilePath\="$RemoteVcfFilePath",RemoteBamFilePath\="$RemoteBamFilePath"\> > {output.meta}

		"""


# Relatedness Testing
rule relatedness_test:
	input:
		"output/jointvcf_all_variants_filtered_genotype_roi/{seq_id}_all_variants_filtered_genotype_roi.vcf"
	output:
		"output/qc_reports/relatedness/{seq_id}.relatedness2",
	shell:
		"vcftools --relatedness2 "
		"--out output/qc_reports/relatedness/{wildcards.seq_id} "
		"--vcf {input} "

# Merge all the sample QC summaries into a single file
rule create_merged_qc:
	input:
		expand("{sample_name}/{sample_name}_{sample_number}_QC.txt", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		"output/qc_reports/combined_qc/" + seq_id + "_combined_QC.txt"
	shell:
		"python scripts/merge_qc_files.py . ; mv combined_QC.txt {output}"

# Multiqc to compile all qc data into one file
rule multiqc:
	input:
		expand("output/qc_reports/insert_size_metrics/{sample_name}_{sample_number}_InsertSizeMetrics.txt", zip, sample_name=sample_names, sample_number=sample_numbers),
		expand("output/qc_reports/mark_duplicates/{sample_name}_{sample_number}_MarkDuplicatesMetrics.txt",zip, sample_name=sample_names, sample_number=sample_numbers),
		expand("output/qc_reports/hs_metrics/{sample_name}_{sample_number}_HsMetrics.txt",zip, sample_name=sample_names, sample_number=sample_numbers),
		expand("output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.sample_summary",zip, sample_name=sample_names, sample_number=sample_numbers),
		expand("output/qc_reports/alignment_metrics/{sample_name}_{sample_number}_AlignmentSummaryMetrics.txt",zip, sample_name=sample_names, sample_number=sample_numbers),
		expand("output/qc_reports/verify_bam_id/{sample_name}_{sample_number}_verify_bam_id.selfSM",zip, sample_name=sample_names, sample_number=sample_numbers),
		expand("output/qc_reports/variant_calling_metrics/{seq_id}_CollectVariantCallingMetrics.variant_calling_detail_metrics", seq_id = seq_id),
		expand("output/qc_reports/variant_calling_metrics/{seq_id}_CollectVariantCallingMetrics.variant_calling_summary_metrics", seq_id = seq_id),
		expand("output/qc_reports/relatedness/{seq_id}.relatedness2",  seq_id = seq_id),
		fastqc = get_fastqc
	output:
		html = "output/qc_reports/multiqc/" + seq_id + ".html",
		data = directory("output/qc_reports/multiqc/" + seq_id + "_data")
	params:
		seq_id = seq_id
	shell:
		"multiqc --filename {params.seq_id} "
		"--outdir output/qc_reports/multiqc/ "
		"--exclude fastp output/qc_reports"	



#-----------------------------------------------------------------------------------------------------------------#
# Structural Variant and CNV Calling
#-----------------------------------------------------------------------------------------------------------------#

# Create the manta runWorkflow.py script
rule run_manta_config:
	input:
		bam_file = "output/final_bam/{sample_name}_{sample_number}_final.bam"
	output:
		temp("output/manta/{sample_name}_{sample_number}/runWorkflow.py")
	params:
		ref = config["reference"]
	conda:
		"envs/python2.yaml"
	group: "manta"
	shell:
		"configManta.py "
		"--bam {input} "
		"--referenceFasta {params.ref} "
		"--exome "
		"--runDir output/manta/{wildcards.sample_name}_{wildcards.sample_number}"

# Execute the manta script
rule run_manta:
	input:
		"output/manta/{sample_name}_{sample_number}/runWorkflow.py"
	output:
		"output/manta/{sample_name}_{sample_number}/results/variants/diploidSV.vcf.gz",
		"output/manta/{sample_name}_{sample_number}/results/variants/diploidSV.vcf.gz.tbi",
	conda:
		"envs/python2.yaml"
	threads:
		config["manta_threads"]
	group: "manta"
	shell:
		"{input} "
		"--quiet "
		"-m local "
		"-j {threads}"

# Just keep the Manta
rule copy_manta_results:
	input:
		vcf = "output/manta/{sample_name}_{sample_number}/results/variants/diploidSV.vcf.gz",
		index = "output/manta/{sample_name}_{sample_number}/results/variants/diploidSV.vcf.gz.tbi"
	output:
		vcf = "output/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz",
		index = "output/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz.tbi"
	shell:
		"cp {input.vcf} {output.vcf} && cp {input.index} {output.index} && rm -r output/manta/{wildcards.sample_name}_{wildcards.sample_number}/"			

# Collect all the high coverage bams and put them in a single file.
rule high_coverage_bam_list:
	input:
		expand("output/qc_reports/high_coverage_bams/{sample_name}_{sample_number}_high_coverage.txt", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		"output/qc_reports/final_high_coverage_bams/high_coverage_bams.txt"
	shell:
		"cat {input} > {output}"

# Collect all the low coverage bams and put them in a single file.
rule low_coverage_bam_list:
	input:
		expand("output/qc_reports/low_coverage_bams/{sample_name}_{sample_number}_low_coverage.txt", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		"output/qc_reports/final_low_coverage_bams/low_coverage_bams.txt"
	shell:
		"cat {input} > {output}"

# Create the bed file for CNV analysis
rule make_cnv_bed:
	input:
		config["capture_bed_file"]
	output:
		"output/config/cnv_bed/" + panel + "_ROI_b37_CNV.bed"
	params:
		genomic_superdups_bed = config["genomic_superdups_bed"],
		ref_index = config["reference_index"],
	shell:
		"bedtools slop -i {input} "
		"-g {params.ref_index} "
		"-b 250 | "
		"grep -v ^X | "
		"grep -v ^Y | "
		"grep -v ^MT | "
		"bedtools sort "
		"-faidx {params.ref_index} | "
		"bedtools merge | "
		"bedtools subtract "
		"-A -a - -b {params.genomic_superdups_bed} | "
		"awk -F\"\t\" '{{print $1\"\t\"$2\"\t\"$3\"\tr\"NR}}' > {output} "

# Exome depth requires bam indexes as .bam.bai rather than .bai
rule change_bam_indexes:
	input:
		bam_index = "output/final_bam/{sample_name}_{sample_number}_final.bai"
	output:
		bam_index = temp("output/final_bam/{sample_name}_{sample_number}_final.bam.bai")
	shell:
		"cp {input.bam_index} {output.bam_index}"


# Run the R ExomeDepth program - create empty files for samples with low depth
rule call_cnvs:
	input:
		high_bam_list = "output/qc_reports/final_high_coverage_bams/high_coverage_bams.txt",
		low_bam_list = "output/qc_reports/final_low_coverage_bams/low_coverage_bams.txt",
		bed = "output/config/cnv_bed/" + panel + "_ROI_b37_CNV.bed",
		bam_indexes = expand("output/final_bam/{sample_name}_{sample_number}_final.bam.bai", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		expand("output/exome_depth/{sample_name}_{sample_number}_final_cnv_fixed.vcf.gz", zip, sample_name=sample_names, sample_number=sample_numbers),
		expand("output/exome_depth/{sample_name}_{sample_number}_final_cnv_fixed.vcf.gz.tbi", zip, sample_name=sample_names, sample_number=sample_numbers),
		"output/exome_depth/ExomeDepth.log",
		"output/exome_depth/" + seq_id + "_ExomeDepth_Metrics.txt",
	params:
		ref = config["reference"],
		seq_id = seq_id,
		prefix = "output/exome_depth/",
		sequence_dict = config["reference_sequence_dict"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; bash scripts/call_cnvs.sh {input.high_bam_list} {params.ref} {input.bed} {params.seq_id} {params.prefix} {params.sequence_dict} {input.low_bam_list}"

#-----------------------------------------------------------------------------------------------------------------#
# Panel Specific Rules
#-----------------------------------------------------------------------------------------------------------------#
if panel == "IlluminaTruSightCancer":


	# Generate single bedfile from gene beds in order to enable custom reporting of gaps and coverage for TSC panel
	rule create_hotspots_bed_file:
		input:
			config["hotspot_bed_dir"]
		output:
			"output/config/hotspot_bed/IlluminaTruSightCancer_CustomROI_b37.bed"
		shell:
			"cat {input}*.bed | sort -k1,1 -k2,2n > {output}"

	# Merge the Manta and CNV calls into a single csv file along with QC data
	rule create_combined_sv_report:
		input:
			bed = "output/config/hotspot_bed/IlluminaTruSightCancer_CustomROI_b37.bed",
			manta = expand("output/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz", zip, sample_name=sample_names, sample_number=sample_numbers),
			manta_index = expand("output/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz.tbi", zip, sample_name=sample_names, sample_number=sample_numbers),
			exome_depth = expand("output/exome_depth/{sample_name}_{sample_number}_final_cnv_fixed.vcf.gz.tbi", zip, sample_name=sample_names, sample_number=sample_numbers),
			high_coverage_bams = "output/qc_reports/final_high_coverage_bams/high_coverage_bams.txt",
			exome_depth_metrics = "output/exome_depth/" + seq_id + "_ExomeDepth_Metrics.txt",
		output:
			"output/combined_sv_report/" + seq_id + "_cnvReport.csv"
		params:
			depth_folder = "output/depth/depth_of_coverage/",
			exome_depth_folder = "output/exome_depth/",
			manta_folder = "output/manta/"			
		shell:
			"Rscript --vanilla scripts/generateCnvReport.R "
			"{input.bed} "
			"{input.exome_depth_metrics} "
			"{input.high_coverage_bams} "
			"{params.depth_folder} "
			"{params.exome_depth_folder} "
			"{params.manta_folder} "
			"{output} "

	# Create custom coverage data for each sample
	rule get_custom_coverage:
		input:
			depth = "output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.gz",
			index = "output/depth/depth_of_coverage/{sample_name}_{sample_number}_DepthOfCoverage.gz.tbi",
			bed = config["hotspot_bed_dir"] + "{bedfile}.bed",
		output:
			"output/depth/hotspot_coverage/{sample_name}_{sample_number}_{bedfile}.coverage",
			"output/depth/hotspot_coverage/{sample_name}_{sample_number}_{bedfile}.gaps",
			"output/depth/hotspot_coverage/{sample_name}_{sample_number}_{bedfile}.missing",
			"output/depth/hotspot_coverage/{sample_name}_{sample_number}_{bedfile}.totalCoverage"
		params:
			coverage_calculator = config["coverage_calculator"],
			min_depth = config["minimum_coverage"]
		group:
			"custom_coverage"
		shell:
			"python  {params.coverage_calculator} "
			"-B {input.bed} "
			"-D {input.depth} "
			"--depth {params.min_depth} "
			"--padding 0 "
			"--outdir output/depth/hotspot_coverage/ "
			"--outname {wildcards.sample_name}_{wildcards.sample_number}_{wildcards.bedfile} "

	rule collect_custom_coverage:
		input:
			get_all_custom_coverage
		output:
			"output/depth/hotspot_coverage/custom.finished"
		group:
			"custom_coverage"
		shell:
			"touch {output}"

#-----------------------------------------------------------------------------------------------------------------#
# Final Rules
#-----------------------------------------------------------------------------------------------------------------#

# Rules to create outputs. We have a seperate final rule for each assay.

if config["perform_bqsr"] == True:

	rule final:
		input:
			"output/validated_vcf/{seq_id}.validated",
			expand("output/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz", zip, sample_name=sample_names, sample_number=sample_numbers),
			"output/variant_reports/{seq_id}_finished.txt",
			"output/jointvcf_all_variants_filtered_genotype_roi_meta_nomt/{seq_id}_all_variants_filtered_genotype_roi_meta_nomt.vcf",
			expand("output/qc_reports/bqsr/{sample_name}_{sample_number}_bqsr_covariation.csv", zip, sample_name=sample_names, sample_number=sample_numbers),
			"output/qc_reports/multiqc/" + seq_id + ".html",
			"output/qc_reports/combined_qc/" + seq_id + "_combined_QC.txt"
		output:
			"output/pipeline_finished/{seq_id}.finished"
		shell:
			"touch {output}"

else:

	# Collect custom coverage and SV combined report
	if panel == "IlluminaTruSightCancer":

		rule final:
			input:
				"output/validated_vcf/{seq_id}.validated",
				expand("output/manta/{sample_name}_{sample_number}/results/variants/diploidSV.vcf.gz", zip, sample_name=sample_names, sample_number=sample_numbers),
				"output/depth/hotspot_coverage/custom.finished",
				"output/vcf_csv/{seq_id}_vcf.csv",
				"output/variant_reports/{seq_id}_finished.txt",
				"output/combined_sv_report/" + seq_id + "_cnvReport.csv",
				"output/jointvcf_all_variants_filtered_genotype_roi_meta_nomt/{seq_id}_all_variants_filtered_genotype_roi_meta_nomt.vcf",
				"output/qc_reports/multiqc/" + seq_id + ".html",
				"output/qc_reports/combined_qc/" + seq_id + "_combined_QC.txt",
			output:
				"output/pipeline_finished/{seq_id}.finished"
			shell:
				"touch {output}"

	else:


		rule final:
			input:
				"output/validated_vcf/{seq_id}.validated",
				expand("output/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz", zip, sample_name=sample_names, sample_number=sample_numbers),
				"output/vcf_csv/{seq_id}_vcf.csv",
				"output/variant_reports/{seq_id}_finished.txt",
				"output/jointvcf_all_variants_filtered_genotype_roi_meta_nomt/{seq_id}_all_variants_filtered_genotype_roi_meta_nomt.vcf",
				"output/qc_reports/multiqc/" + seq_id + ".html",
				"output/qc_reports/combined_qc/" + seq_id + "_combined_QC.txt"
			output:
				"output/pipeline_finished/{seq_id}.finished"
			shell:
				"touch {output}"




















