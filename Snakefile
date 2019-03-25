
"""
Germline Enrichment Pipeline Using GATK

Run:

snakemake -p

"""

from pathlib import Path 


#-----------------------------------------------------------------------------------------------------------------#
# Configuration variables
#-----------------------------------------------------------------------------------------------------------------#

# Which YAMl config to use
config_location = "config.yaml"
configfile: config_location

# How many lanes do we have?
folder, name, lanes = glob_wildcards("{folder}/{sample_number}_{lane}_R1_001.fastq.gz")
lanes =  list((set(lanes)))

# What are our samples and sample numbers?
folder, sample_names, sample_numbers = glob_wildcards("{folder}/{sample_name}_{sample_number}_L001_R1_001.fastq.gz")

# Get other data from config file
chromosomes = config["chromosomes"]
panel = config["panel"]
seqid = config["seqId"]

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

			file_list.append(sample_name + "/" + sample_name + "_" + sample_number + "_" + lane + "_R1_001_qfilter_fastqc.html" )
			file_list.append(sample_name + "/" + sample_name + "_" + sample_number + "_" + lane + "_R2_001_qfilter_fastqc.html" )

	return file_list

def get_all_custom_coverage(wildcards):
	"""
	Return a list of strings with all the custom coverage files in.

	"""

	file_list = []
	beds = glob_wildcards(config["hotspot_bed_dir"] + "{bedfile}.bed")

	for sample_name, sample_number in zip(sample_names,sample_numbers ):

		for bedfile in beds[0]:

			my_input = sample_name + "/hotspot_coverage/" + sample_name + "_" + sample_number + "_" + bedfile + ".coverage"
			file_list.append(my_input)

	return file_list

#-----------------------------------------------------------------------------------------------------------------#
# Starting Rule
#-----------------------------------------------------------------------------------------------------------------#

# All function pulls all rules together
rule all:
	input:
		expand("{seqid}.finished" , seqid=seqid),

#-----------------------------------------------------------------------------------------------------------------#
# Initial QC and Preprocessing
#-----------------------------------------------------------------------------------------------------------------#

# Create a PED file for downstream analysis
rule create_ped_file:
	input:
		config_location
	output:
		"{seqid}.ped"
	shell:
		"python scripts/make_ped.py --config {input} > {output}"

# Run the fastp program to generate a read quality report and trim read adapters 
rule fastp:
	input:
		fwd = "{sample_name}/{sample_name}_{sample_number}_{lane}_R1_001.fastq.gz",
		rev = "{sample_name}/{sample_name}_{sample_number}_{lane}_R2_001.fastq.gz"
	output:
		html = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_fastp.html",
		json = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_fastp.json",
		fwd = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R1_001_qfilter.fastq.gz"),
		rev = temp("{sample_name}/" + seqid  + "_{sample_name}_{sample_number}_{lane}_R2_001_qfilter.fastq.gz")
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
		"--length_required 35 "
		"-w {threads}"

# Check for inter-species contamination
rule fastq_screen:
	input:
		fwd = "{sample_name}/{sample_name}_{sample_number}_{lane}_R1_001.fastq.gz",
		rev = "{sample_name}/{sample_name}_{sample_number}_{lane}_R2_001.fastq.gz"
	output:
		"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R1_001_screen.html",
		temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R1_001_screen.png"),
		"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R1_001_screen.txt",
		"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R2_001_screen.html",
		temp("{sample_name}" + seqid + "_/{sample_name}_{sample_number}_{lane}_R2_001_screen.png"),
		"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R2_001_screen.txt"	
	threads:
		config["fastq_screen_threads"]
	params:
		fastq_screen_config = config["fastq_screen_config"]
	shell:
		"fastq_screen "
		"--aligner bwa "
		"--threads {threads} "
		"--outdir {wildcards.sample_name} "
		"--conf {params.fastq_screen_config} "
		"{input.fwd} "
		"{input.rev}"

# Run fastqc
rule fastqc:
	input:
		fwd = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R1_001_qfilter.fastq.gz",
		rev = "{sample_name}/" + seqid  + "_{sample_name}_{sample_number}_{lane}_R2_001_qfilter.fastq.gz"
	output:
		html_r1 = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R1.html",
		summary_r1 = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R1.txt",	
		htmlr2 = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R2.html",
		summary_r2 = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R2.txt"
	threads:
		config["fastqc_threads"]
	params:
		temp_dir = config["fastqc_temp_dir"],
		seqid = seqid
	shell:
		"""
		fastqc --threads {threads} \
		--dir {params.temp_dir} \
		--outdir {wildcards.sample_name} \
		--extract \
		{input.fwd} \
		{input.rev} \

		# move summary file
		mv {wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R1_001_qfilter_fastqc/summary.txt \
		{wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R1.txt
		
		# move summary file
		mv {wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R1_001_qfilter_fastqc.html \
		{wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R1.html
		
		# move summary file
		mv {wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R2_001_qfilter_fastqc/summary.txt \
		{wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R2.txt
		
		# move summary file
		mv {wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R2_001_qfilter_fastqc.html \
		{wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R2.html
		
		# remove unneeded files
		rm -r {wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R1_001_qfilter_fastqc
		rm -r {wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R2_001_qfilter_fastqc
		rm {wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R1_001_qfilter_fastqc.zip
		rm {wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_R2_001_qfilter_fastqc.zip

		"""
	

#-----------------------------------------------------------------------------------------------------------------#
# Alignment and Further Preprocessing
#-----------------------------------------------------------------------------------------------------------------#

# Align reads with bwa mem, pipe into samtools and sort
rule bwa_align:
	input:
		fwd = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_{lane}_R1_001_qfilter.fastq.gz",
		rev = "{sample_name}/" + seqid  + "_{sample_name}_{sample_number}_{lane}_R2_001_qfilter.fastq.gz"
	output:
		temp("temp/bam/" + seqid + "_{sample_name}_{sample_number}_{lane}.bam")
	threads:
		config["bwa_threads"]
	params:
		ref = config["bwa_reference"],
		seqid = config["seqId"],
		centre = config["centre"],
		samtools_temp_dir = config["samtools_temp_dir"]
	shell:
		"bwa mem "
		"-t {threads} "
		"-M "
		"-R '@RG\\tID:{params.seqid}.{wildcards.lane}\\tCN:{params.centre}\\tSM:{wildcards.sample_name}\\tLB:{params.seqid}\\tPL:ILLUMINA' "
		"{params.ref} {input.fwd} {input.rev} | "
		"samtools view -Sb - | "
		"samtools sort -T {params.samtools_temp_dir}.temp -O bam > {output}"

# Index the bam file
rule index_original_bam:
	input:
		"temp/bam/" + seqid + "_{sample_name}_{sample_number}_{lane}.bam"
	output:
		temp("temp/bam/" + seqid + "_{sample_name}_{sample_number}_{lane}.bam.bai")
	shell:
		"samtools index {input}"

# Merge the bams and mark duplicates
rule merge_and_remove_duplicates:
	input:
		bams = expand("temp/bam/" + seqid + "_{{sample_name}}_{{sample_number}}_{lane}.bam", lane=lanes),
		bam_indexes = expand("temp/bam/" + seqid + "_{{sample_name}}_{{sample_number}}_{lane}.bam.bai", lane=lanes),
	output:
		bam = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups.bam"),
		index = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups.bai"),
		metrics = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_MarkDuplicatesMetrics.txt"
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

# identify regions requiring indel realignment
rule identify_realignment_regions:
	input:
		bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups.bam",
		index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups.bai"
	output:
		temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_realign.intervals")
	params:
		ref = config["reference"],
		bed = config["capture_bed_file"],
		known_sites_indels = config["indels_1k_vcf"],
		known_sites_gold = config["gold_standard_indels"],
		java_options = config['gatk_hc_java_options'],
		java_home = config["java_home"],
		padding = config["interval_padding_bqsr"]
	shell:
		"export JAVA_HOME={params.java_home}; gatk3 "
		"{params.java_options} "
		"-T RealignerTargetCreator  "
		"-R {params.ref} "
		"-known {params.known_sites_indels} "
		"-known {params.known_sites_gold} "
		"-I {input.bam} "
		"-o {output} "
		"-L {params.bed} "
		"-ip {params.padding} "
		"-dt NONE "

# Realign around indels
rule realign_indels:
	input:
		bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups.bam",
		index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups.bai",
		intervals = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_realign.intervals"
	output:
		bam = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups_realigned.bam"),
		index = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups_realigned.bai")
	params:
		ref = config["reference"],
		known_sites_indels = config["indels_1k_vcf"],
		known_sites_gold = config["gold_standard_indels"],
		java_options = config['gatk_hc_java_options'],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; gatk3 "
		"{params.java_options} "
		"-T IndelRealigner "
		"-R {params.ref} "
		"-known {params.known_sites_indels} "
		"-known {params.known_sites_gold} "
		"-targetIntervals {input.intervals} "
		"-I {input.bam} "
		"-o {output.bam} "
		"-dt NONE "


if config["perform_bqsr"] == True:

	# Create the BQSR report to later apply
	rule create_base_quality_report:
		input:
			bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups_realigned.bam",
			index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups_realigned.bai"
		output:
			"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_recal_data.table"
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
			bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups_realigned.bam",
			index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups_realigned.bai",
			bqsr_report = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_recal_data.table"
		output:
			bam_file = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
			bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai"
		params:
			ref = config["reference"],
			java_options = config["gatk_bqsr_java_options"]
		shell:
			"gatk --java-options '{params.java_options}' "
			"ApplyBQSR -R {params.ref} "
			"-I {input.bam} "
			"-bqsr {input.bqsr_report} "
			"-O {output.bam_file} "

	# Second pass to get table for AnalyseCovariates
	# Do we really need this?
	rule create_base_quality_report_second_pass:
		input:
			bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
			bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai"
		output:
			"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_post_recal_data.table"
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
			before = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_recal_data.table",
			after = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_post_recal_data.table"
		output:
			pdf = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_recalibration_plots.pdf",
			csv = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_recalibration.csv",
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

	# Move bam to final folder if we don't perform bqsr
	rule move_final_bam:
		input:
			bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups_realigned.bam",
			bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_merged_nodups_realigned.bai"
		output:
			bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
			bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai"
		shell:
			"cp {input.bam} {output.bam}; cp {input.bam_index} {output.bam_index} "

#-----------------------------------------------------------------------------------------------------------------#
# Post Alignment QC
#-----------------------------------------------------------------------------------------------------------------#


# Create an interval file from the BED file for use in Picard tools such as CollectHsMetrics
rule create_interval_file:
	input:
		ancient(config["capture_bed_file"])
	output:
		temp(Path(config["capture_bed_file"]).name.split(".")[0] + ".interval_list")
	params:
		sequence_dict = config["reference_sequence_dict"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard BedToIntervalList I={input} O={output} SD={params.sequence_dict}" 

# Collect the insert size metrics using picard
rule collect_insert_size_metrics:
	input:
		bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
		bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai"
	output:
		txt = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_InsertSizeMetrics.txt",
		pdf = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_InsertSizeMetrics.pdf"
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
		bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
		bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai",
		intervals = Path(config["capture_bed_file"]).name.split(".")[0] + ".interval_list"
	output:
		"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_HsMetrics.txt"
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
		bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
		bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai"
	output:
		"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_AlignmentSummaryMetrics.txt"
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
		ancient(config["1000k_high_confidence_snps"])
	output:
		vcf = temp("1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf"),
		index = temp("1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf.idx")
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
		bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
		bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai",
		vcf = "1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf"
	output:
		depthSM = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_Contamination.depthSM",
		log = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_Contamination.log"),
		selfSM = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_Contamination.selfSM"
	params:
		seqid = seqid
	shell:
		"verifyBamID "
		"--vcf {input.vcf} "
		"--bam {input.bam} "
		"--out {wildcards.sample_name}/{params.seqid}_{wildcards.sample_name}_{wildcards.sample_number}_Contamination "
		"--verbose "
		"--ignoreRG "
		"--chip-none "
		"--minMapQ 20 "
		"--maxDepth 1000 "
		"--precise "

# Calculate per base coverage
rule generate_per_base_coverage:
	input:
		bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
		bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai"
	output:
		depth = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage"),
		summary = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage.sample_summary"
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
		"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage"
	output:
		depth = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage.gz",
		index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage.gz.tbi"
	shell:
		"sed 's/:/\t/g' {input} | grep -v 'Locus' | sort -k1,1 -k2,2n | bgzip > {output.depth} && "
		"tabix -b 2 -e 2 -s 1 {output.depth} "


#-----------------------------------------------------------------------------------------------------------------#
# Find Coverage Gaps
#-----------------------------------------------------------------------------------------------------------------#

rule calculate_coverage_metrics:
	input:
		depth = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage.gz",
		index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage.gz.tbi"	
	output:
		"{sample_name}/" + seqid + "_{sample_name}_{sample_number}_Gaps.bed"
	params:
		roi_bed = config["capture_bed_file"],
		panel = panel,
		hotspots = config["germline_hotspots"],
		min_coverage = config["minimum_coverage"],
		seqid = seqid,
		refseq = config["refseq"]
	shell:
		"bash scripts/get_coverage.sh "
		"{params.roi_bed} "
		"{params.panel} "
		"{params.hotspots} "
		"{input.depth} "
		"{params.seqid} "
		"{wildcards.sample_name}_{wildcards.sample_number} "
		"{params.min_coverage} "
		"{params.refseq} "
		"temp/depth_metrics "
		"{wildcards.sample_name}"
	
#-----------------------------------------------------------------------------------------------------------------#
# SNP and Small Indel Calling with GATK Haplotype Caller
#-----------------------------------------------------------------------------------------------------------------#

# Sort ROI bed for splitting by bedextract
rule sort_capture_bed:
	input:
		ancient(config["capture_bed_file"])
	output:
		temp("temp/sorted_beds/{{panel}}_sorted.bed".format(panel=panel))
	shell:
		"sort-bed {input} > {output}"

# Split the bed by chromosome for input into create_gvcfs
rule split_bed_by_chromosome:
	input:
		"temp/sorted_beds/{panel}_sorted.bed".format(panel=panel)
	output:
		temp(expand("temp/split_capture_bed/{chr}.bed", chr=chromosomes))
	params:
		chromosomes = chromosomes
	shell:
		"for chr in {params.chromosomes}; do bedextract $chr {input} > temp/split_capture_bed/$chr.bed; done"


# Create GVCF using Haplotype Caller for each sample chromosome combination
rule create_gvcfs:
	input:
		bam_file = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
		bam_index= "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai",
		bed = "temp/split_capture_bed/{chr}.bed",
		ped = seqid + ".ped"
	output:
		gvcf_file = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_chr{chr}.g.vcf"),
		index = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_chr{chr}.g.vcf.idx")
	params:
		ref = config["reference"],
		padding = config['interval_padding_haplotype_caller'],
		java_options = config['gatk_hc_java_options'],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; gatk3 "
		"{params.java_options} "
		"-T HaplotypeCaller "
		"-R {params.ref} "
		"-I {input.bam_file} "
		"-L {input.bed} "
		"-ip {params.padding} "
		"-o {output.gvcf_file} "
		"-ped {input.ped} "
		"--genotyping_mode DISCOVERY "
		"--emitRefConfidence GVCF "
		"-dt NONE "

# Genotype the gvcfs and produce a joint vcf
rule genotype_gvcfs:
	input:
		gvcfs = expand("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_chr{{chr}}.g.vcf" , zip, sample_name=sample_names, sample_number=sample_numbers),
		index = expand("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_chr{{chr}}.g.vcf.idx", zip, sample_name=sample_names, sample_number=sample_numbers),
		bed = "temp/split_capture_bed/{chr}.bed",
		ped = seqid + ".ped"
	output:
		vcf = temp("temp/jointvcf_per_chr/{seqid}_chr{chr}.vcf"),
		index = temp("temp/jointvcf_per_chr/{seqid}_chr{chr}.vcf.idx")
	params:
		ref = config["reference"],
		java_options = config['gatk_hc_java_options'],
		padding = config['interval_padding_haplotype_caller'],
		files = lambda wildcards, input: " -V ".join(input.gvcfs),
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; gatk3 "
		"{params.java_options} "
		"-T GenotypeGVCFs "
		"-R {params.ref} "
		"-V {params.files} "
		"-L {input.bed} "
		"-ip {params.padding} "
		"-o {output.vcf} "
		"-ped {input.ped} "
		"-dt NONE "

# Combine the chromsome vcfs into one final vcf with all samples and all chromosomes
rule collect_vcfs:
	input:
		vcf = expand("temp/jointvcf_per_chr/{{seqid}}_chr{chr}.vcf", chr= chromosomes),
		index = expand("temp/jointvcf_per_chr/{{seqid}}_chr{chr}.vcf.idx", chr= chromosomes)
	output:
		vcf = temp("temp/jointvcf/{seqid}_all_chr.vcf"),
		index = temp("temp/jointvcf/{seqid}_all_chr.vcf.idx")
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
		vcf = "temp/jointvcf/{seqid}_all_chr.vcf",
		index = "temp/jointvcf/{seqid}_all_chr.vcf.idx"
	output:
		vcf = temp("temp/jointvcf_snps/{seqid}_all_chr_snps.vcf"),
		index = temp("temp/jointvcf_snps/{seqid}_all_chr_snps.vcf.idx")
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
		vcf = "temp/jointvcf_snps/{seqid}_all_chr_snps.vcf",
		index = "temp/jointvcf_snps/{seqid}_all_chr_snps.vcf.idx"
	output:
		vcf = temp("temp/jointvcf_snps_filtered/{seqid}_all_chr_snps_filtered.vcf"),
		index = temp("temp/jointvcf_snps_filtered/{seqid}_all_chr_snps_filtered.vcf.idx")
	params:
		ref = config["reference"],
		padding = config["interval_padding_haplotype_caller"],
		bed = config["capture_bed_file"],
		min_qual = config["snp_min_qual"],
		min_QD = config["snp_min_QD"],
		max_FS = config["snp_max_FS"],
		min_MQ = config["snp_min_MQ"],
		min_MQRankSum = config["snp_min_MQRankSum"],
		min_ReadPosRankSum = config["snp_min_ReadPosRankSum"],
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
		vcf = "temp/jointvcf/{seqid}_all_chr.vcf",
		index = "temp/jointvcf/{seqid}_all_chr.vcf.idx",
	output:
		vcf = temp("temp/jointvcf_indels/{seqid}_all_chr_indels.vcf"),
		index = temp("temp/jointvcf_indels/{seqid}_all_chr_indels.vcf.idx")
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

# Filter the non SNPs e.g. INDEL, MIXED, MNP, SYMBOLIC, NO_VARIATION)
rule filter_non_snps:
	input:
		vcf = "temp/jointvcf_indels/{seqid}_all_chr_indels.vcf",
		index = "temp/jointvcf_indels/{seqid}_all_chr_indels.vcf.idx"
	output:
		vcf = temp("temp/jointvcf_indels_filtered/{seqid}_all_chr_indels_filtered.vcf"),
		index = temp("temp/jointvcf_indels_filtered/{seqid}_all_chr_indels_filtered.vcf.idx")
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
		"--filter-name 'SOR' "
		"--filter-expression 'ReadPosRankSum <  {params.min_ReadPosRankSum}' "
		"--filter-name 'ReadPosRankSum' "
		"--filter-expression \"InbreedingCoeff != 'NaN' && InbreedingCoeff < {params.min_inbreeding_coeff}\" "
		"--filter-name 'InbreedingCoeff' "

# Combine the filtered SNPs and Indels into a single VCF
rule combine_filtered_snps_and_indels:
	input:
		snps = "temp/jointvcf_snps_filtered/{seqid}_all_chr_snps_filtered.vcf",
		indels = "temp/jointvcf_indels_filtered/{seqid}_all_chr_indels_filtered.vcf"
	output:
		vcf = temp("temp/jointvcf_all_variants_filtered/{seqid}_all_variants_filtered.vcf"),
		index = temp("temp/jointvcf_all_variants_filtered/{seqid}_all_variants_filtered.vcf.idx")
	params:
		java_options = config["gatk_variants_java_options"]
	shell:
		"gatk --java-options '{params.java_options}' "
		"MergeVcfs "
		"-I {input.snps} "
		"-I {input.indels} "
		"-O {output.vcf}"		

# Filter on genotype depth - just mark genotypes with low depth
rule filter_genotypes:
	input:
		vcf = "temp/jointvcf_all_variants_filtered/{seqid}_all_variants_filtered.vcf",
		index = "temp/jointvcf_all_variants_filtered/{seqid}_all_variants_filtered.vcf.idx"
	output:
		vcf = "temp/jointvcf_all_variants_filtered_genotype/{seqid}_all_variants_filtered_genotype.vcf",
		index = temp("temp/jointvcf_all_variants_filtered_genotype/{seqid}_all_variants_filtered_genotype.vcf.idx")
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
		"temp/jointvcf_all_variants_filtered_genotype/{seqid}_all_variants_filtered_genotype.vcf"
	output:
		"temp/jointvcf_all_variants_filtered_genotype/{seqid}_all_variants_filtered_genotype.vcf.gz",
		"temp/jointvcf_all_variants_filtered_genotype/{seqid}_all_variants_filtered_genotype.vcf.gz.tbi"
	shell:
		"bgzip {input} && tabix {input}.gz"	


# Filter out variants outside of ROI
rule filter_by_roi:
	input:
		vcf = "temp/jointvcf_all_variants_filtered_genotype/{seqid}_all_variants_filtered_genotype.vcf.gz",
		index = "temp/jointvcf_all_variants_filtered_genotype/{seqid}_all_variants_filtered_genotype.vcf.gz.tbi"
	output:
		temp("temp/jointvcf_all_variants_filtered_genotype_roi/{seqid}_all_variants_filtered_genotype_roi.vcf")
	params:
		bed = config["capture_bed_file"],
		ref = config["reference"]
	shell:
		"bcftools view -R {params.bed} {input.vcf} > {output} "

# Use vt to split multiallelics and normalise variants
rule decompose_and_normalise:
	input:
		"temp/jointvcf_all_variants_filtered_genotype_roi/{seqid}_all_variants_filtered_genotype_roi.vcf"
	output:
		temp("temp/jointvcf_all_variants_filtered_genotype_roi_norm/{seqid}_all_variants_filtered_genotype_roi_norm.vcf")
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
		"temp/jointvcf_all_variants_filtered_genotype_roi_norm/{seqid}_all_variants_filtered_genotype_roi_norm.vcf"
	output:
		vcf = "{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf",
		summary = temp("{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf_summary.html"),
		warnings = temp("{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf_warnings.txt")
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
		"--output_file {output.vcf} "
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

# Compress and index the VEP annotated vcf
rule compress_and_index_vep:
	input:
		"{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf"
	output:
		"{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz",
		"{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz.tbi"
	shell:
		"bgzip {input} && tabix {input}.gz"	

# Convert vcf to csv
rule convert_to_csv:
	input:
		vcf = "{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz",
		index = "{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz.tbi"
	output:
		temp("temp/vcf_csv/{seqid}_vcf.csv")
	shell:
		"gatk VariantsToTable -V {input.vcf} "
		"-O {output} -F CHROM -F POS -F REF -F ALT -F ID -F QUAL -F FILTER -F CSQ -F AC "
		"-GF GT -GF GQ -GF DP"

# Get the CSQ string which describes the VEP fields	
rule get_csq_string:
	input:
		vcf = "{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz",
		index = "{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz.tbi"
	output:
		temp("temp/csq/{seqid}_csq.txt")
	shell:
		"zcat {input.vcf} | grep \"^##INFO=<ID=CSQ\" | awk 'BEGIN {{ FS = \":\" }} ; {{ print $2 }}' | tr -d '>\" ' > {output} "

# Run the germline variant filter program
rule create_variant_reports:
	input:
		csv = "temp/vcf_csv/{seqid}_vcf.csv",
		ped = "{seqid}.ped",
		csq = "temp/csq/{seqid}_csq.txt"
	output:
		"variant_reports/{seqid}_finished.txt"
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
		"--gnomad-constraint-scores "
		"--worksheet {wildcards.seqid} "
		"--results-dir variant_reports/ "
		"--csq $(cat {input.csq}) "
		"--add-ccrs && touch {output}"

# Merge the metadata from all samples into a single file
rule collect_meta_data:
	input:
		meta = expand("temp/vcf_meta/{sample_name}_{sample_number}_meta.txt", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		meta = temp("temp/vcf_meta/merged_meta.txt")
	shell:
		"cat {input.meta} > {output.meta}"


# Add metadata to vcf header for variant database import
rule add_meta_to_vcf:
	input:
		vcf = "temp/jointvcf_all_variants_filtered_genotype_roi/{seqid}_all_variants_filtered_genotype_roi.vcf",
		meta = "temp/vcf_meta/merged_meta.txt"
	output:
		temp("temp/jointvcf_all_variants_filtered_genotype_roi_meta/{seqid}_all_variants_filtered_genotype_roi_meta.vcf")
	shell:
		"""
		 grep '^##' {input.vcf} > {output}
		 cat {input.meta} >> {output}
		 grep -v '^##' {input.vcf} >> {output}
		"""

# Exclude Mitochrondial variants for inclusion in the variant database
rule create_vcf_without_mt:
	input:
		"temp/jointvcf_all_variants_filtered_genotype_roi_meta/{seqid}_all_variants_filtered_genotype_roi_meta.vcf"
	output:
		"{seqid}_all_variants_filtered_genotype_roi_meta_nomt.vcf"
	shell:
		"awk '$1 !~ /^MT/ {{ print $0 }}' {input} > {output}"

# Check the final vcf is valid
rule validate_vcf:
	input:
		vcf = "{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz",
		index = "{seqid}_all_variants_filtered_genotype_roi_norm_vep.vcf.gz.tbi"
	output:
		"temp/validated_vcf/{seqid}.validated"
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

# Evaluate variant calling
rule variant_evaluation:
	input:
		"temp/jointvcf_all_variants_filtered_genotype/{seqid}_all_variants_filtered_genotype.vcf.gz"
	output:
		"{seqid}_CollectVariantCallingMetrics.variant_calling_detail_metrics",
		"{seqid}_CollectVariantCallingMetrics.variant_calling_summary_metrics"
	params:
		dbsnp_vcf = config["dbsnp_vcf_eval"],
		java_home = config["java_home"]
	threads:
		config["variant_evaluation_threads"],
	shell:
		"export JAVA_HOME={params.java_home}; picard CollectVariantCallingMetrics "
		"I={input} "
		"O={wildcards.seqid}_CollectVariantCallingMetrics "
		"DBSNP={params.dbsnp_vcf} "
		"THREAD_COUNT={threads}"

# Use the Y chromosome coverage to calculate the sex
rule calculate_sex:
	input:
		depth_summary = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage.sample_summary",
		bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
		bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai"
	output:
		sex = temp("temp/sex/{sample_name}_{sample_number}_sex.txt"),
		y_bed = temp("temp/y_bed/{sample_name}_{sample_number}_y.bed"),
		y_cov = temp("temp/y_coverage/{sample_name}_{sample_number}_DepthOfCoverage.sample_summary")
	params:
		bed = config["capture_bed_file"],
		ref = config["reference"],
		output_name = "temp/y_coverage/{sample_name}_{sample_number}_DepthOfCoverage",
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
		bam = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
		bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai",
		insert_size_metrics = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_InsertSizeMetrics.txt",
		duplicate_metrics = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_MarkDuplicatesMetrics.txt",
		hs_metrics = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_HsMetrics.txt",
		depth_summary = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage.sample_summary",
		alignments_summary = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_AlignmentSummaryMetrics.txt",
		contamination = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_Contamination.selfSM",
		fastqc_fwd = expand("{{sample_name}}/" + seqid + "_{{sample_name}}_{{sample_number}}_{lane}_R1.txt", lane =lanes),
		fastqc_rev = expand("{{sample_name}}/" + seqid + "_{{sample_name}}_{{sample_number}}_{lane}_R2.txt", lane =lanes),
		sex = "temp/sex/{sample_name}_{sample_number}_sex.txt"
	output:
		summary = "{sample_name}/{sample_name}_{sample_number}_QC.txt",
		high_coverage = "temp/high_coverage_bams/{sample_name}_{sample_number}_high_coverage.txt",
		low_coverage = "temp/low_coverage_bams/{sample_name}_{sample_number}_low_coverage.txt",
		meta = "temp/vcf_meta/{sample_name}_{sample_number}_meta.txt"
	params:
		seqid = seqid,
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

		# Create high coverage bam list
		if [ $(echo "$meanOnTargetCoverage" | awk '{{if ($1 > 20) print "true"; else print "false"}}') = true ]; then
			echo {input.bam} > {output.high_coverage} 
		else

			touch {output.high_coverage} 

		fi
		
		# Create low coverage bam list
		if [ $(echo "$meanOnTargetCoverage" | awk '{{if ($1 < 20) print "true"; else print "false"}}') = true ]; then
			echo {input.bam} > {output.low_coverage} 
		else

			touch {output.low_coverage} 

		fi		

		RemoteVcfFilePath=$(dirname $PWD)/{params.seqid}_all_variants_filtered_genotype_roi_meta_nomt.vcf

		RemoteBamFilePath=$(dirname $PWD)/{input.bam}

		# print vcf meta data needed for db import
		echo \#\#SAMPLE\=\<ID\="{wildcards.sample_name}",Tissue\=Germline,WorklistId\={params.worksheet_id},seqid\={params.seqid},Assay\={params.panel},PipelineName\={params.pipeline_name},PipelineVersion\={params.pipeline_version},RawSequenceQuality\="$rawSequenceQuality",PercentMapped\="$pctPfReadsAligned",ATDropout\="$atDropout",GCDropout\="$gcDropout",MeanInsertSize\="$meanInsertSize",SDInsertSize\="$sdInsertSize",DuplicationRate\="$duplicationRate",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",PctTargetBasesCt\="$pctTargetBasesCt",EstimatedContamination\="$freemix",GenotypicGender\="$calcSex",TotalTargetedUsableBases\="$totalTargetedUsableBases",RemoteVcfFilePath\="$RemoteVcfFilePath",RemoteBamFilePath\="$RemoteBamFilePath"\> > {output.meta}

		"""

# Relatedness Testing
rule relatedness_test:
	input:
		"temp/jointvcf_all_variants_filtered_genotype_roi/{seqid}_all_variants_filtered_genotype_roi.vcf"
	output:
		"{seqid}.relatedness2",
		temp("{seqid}.log")
	shell:
		"vcftools --relatedness2 "
		"--out {wildcards.seqid} "
		"--vcf {input} "

# Merge all the sample QC summaries into a single file
rule create_merged_qc:
	input:
		expand("{sample_name}/{sample_name}_{sample_number}_QC.txt", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		seqid + "_combined_QC.txt"
	shell:
		"python scripts/merge_qc_files.py . ; mv combined_QC.txt {output}"


#-----------------------------------------------------------------------------------------------------------------#
# Structural Variant and CNV Calling
#-----------------------------------------------------------------------------------------------------------------#

# Collect all the high coverage bams and put them in a single file.
rule high_coverage_bam_list:
	input:
		expand("temp/high_coverage_bams/{sample_name}_{sample_number}_high_coverage.txt", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		"temp/final_high_coverage_bams/high_coverage_bams.txt"
	shell:
		"cat {input} > {output}"

# Collect all the low coverage bams and put them in a single file.
rule low_coverage_bam_list:
	input:
		expand("temp/low_coverage_bams/{sample_name}_{sample_number}_low_coverage.txt", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		"temp/final_low_coverage_bams/low_coverage_bams.txt"
	shell:
		"cat {input} > {output}"

# Run manta
rule run_manta:
	input:
		bam_file = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam",
		bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai",
		high_bam_list = "temp/final_high_coverage_bams/high_coverage_bams.txt",
		low_bam_list = "temp/final_low_coverage_bams/low_coverage_bams.txt"
	output:
		vcf = "temp/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz",
		index = "temp/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz.tbi"
	params:
		ref = config["reference"]
	threads:
		config["manta_threads"]
	conda:
		"envs/python2.yaml"
	shell:
		"bash scripts/run_manta.sh "
		"{input.high_bam_list} "
		"{input.bam_file} "
		"{params.ref} "
		"temp/manta/{wildcards.sample_name}_{wildcards.sample_number} "
		"{threads} "
		"{wildcards.sample_name}_{wildcards.sample_number} "
		"temp/manta "
		"{input.low_bam_list}"

# Create the bed file for CNV analysis
rule make_cnv_bed:
	input:
		config["capture_bed_file"]
	output:
		"temp/cnv_bed/" + panel + "_ROI_b37_CNV.bed"
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
		bam_index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bai"
	output:
		bam_index = temp("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam.bai")
	shell:
		"cp {input.bam_index} {output.bam_index}"

# Run the R ExomeDepth program - create empty files for samples with low depth
rule call_cnvs:
	input:
		high_bam_list = "temp/final_high_coverage_bams/high_coverage_bams.txt",
		low_bam_list = "temp/final_low_coverage_bams/low_coverage_bams.txt",
		bed = "temp/cnv_bed/" + panel + "_ROI_b37_CNV.bed",
		bam_indexes = expand("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_final.bam.bai", zip, sample_name=sample_names, sample_number=sample_numbers)
	output:
		expand("temp/exome_depth/" + seqid + "_{sample_name}_{sample_number}_final_cnv_fixed.vcf.gz", zip, sample_name=sample_names, sample_number=sample_numbers),
		expand("temp/exome_depth/" + seqid + "_{sample_name}_{sample_number}_final_cnv_fixed.vcf.gz.tbi", zip, sample_name=sample_names, sample_number=sample_numbers),
		"temp/exome_depth/ExomeDepth.log",
		"temp/exome_depth/" + seqid + "_ExomeDepth_Metrics.txt",
	params:
		ref = config["reference"],
		seqid = seqid,
		prefix = "temp/exome_depth/",
		sequence_dict = config["reference_sequence_dict"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; bash scripts/call_cnvs.sh {input.high_bam_list} {params.ref} {input.bed} {params.seqid} {params.prefix} {params.sequence_dict} {input.low_bam_list}"

#-----------------------------------------------------------------------------------------------------------------#
# Panel Specific Rules
#-----------------------------------------------------------------------------------------------------------------#
if panel == "IlluminaTruSightCancer":


	# Generate single bedfile from gene beds in order to enable custom reporting of gaps and coverage for TSC panel
	rule create_hotspots_bed_file:
		input:
			ancient(config["hotspot_bed_dir"])
		output:
			"temp/hotspot_bed/IlluminaTruSightCancer_CustomROI_b37.bed"
		shell:
			"cat {input}*.bed | sort -k1,1 -k2,2n > {output}"

	# Merge the Manta and CNV calls into a single csv file along with QC data
	rule create_combined_sv_report:
		input:
			bed = "temp/hotspot_bed/IlluminaTruSightCancer_CustomROI_b37.bed",
			manta = expand("temp/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz", zip, sample_name=sample_names, sample_number=sample_numbers),
			manta_index = expand("temp/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz.tbi", zip, sample_name=sample_names, sample_number=sample_numbers),
			exome_depth = expand("temp/exome_depth/" + seqid + "_{sample_name}_{sample_number}_final_cnv_fixed.vcf.gz.tbi", zip, sample_name=sample_names, sample_number=sample_numbers),
			high_coverage_bams = "temp/final_high_coverage_bams/high_coverage_bams.txt",
			exome_depth_metrics = "temp/exome_depth/" + seqid + "_ExomeDepth_Metrics.txt",
		output:
			seqid + "_cnvReport.csv"
		params:
			coverage_dir = "output/depth/depth_of_coverage/",
			exome_dir = "temp/exome_depth/",
			manta_dir = "temp/manta/",
			run = seqid			
		shell:
			"python scripts/generateCNVReport.py "
			"--runid {params.run} "
			"--output {output} "
			"--bed {input.bed} "
			"--exome_metrics {input.exome_depth_metrics} "
			"--manta_dir {params.manta_dir} "
			"--exome_dir {params.exome_dir} "
			"--coverage_dir {params.coverage_dir} "

	# Create custom coverage data for each sample
	rule get_custom_coverage:
		input:
			depth = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage.gz",
			index = "{sample_name}/" + seqid + "_{sample_name}_{sample_number}_DepthOfCoverage.gz.tbi",
			bed = config["hotspot_bed_dir"] + "{bedfile}.bed",
		output:
			"{sample_name}/hotspot_coverage/{sample_name}_{sample_number}_{bedfile}.coverage",
			"{sample_name}/hotspot_coverage/{sample_name}_{sample_number}_{bedfile}.gaps",
			"{sample_name}/hotspot_coverage/{sample_name}_{sample_number}_{bedfile}.missing",
			"{sample_name}/hotspot_coverage/{sample_name}_{sample_number}_{bedfile}.totalCoverage"
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
			"--outdir {wildcards.sample_name}/hotspot_coverage/ "
			"--outname {wildcards.sample_name}_{wildcards.sample_number}_{wildcards.bedfile} "

	rule collect_custom_coverage:
		input:
			get_all_custom_coverage
		output:
			"temp/hotspot_coverage/custom.finished"
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
			"temp/validated_vcf/{seqid}.validated",
			expand("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_Gaps.bed", zip, sample_name=sample_names, sample_number=sample_numbers),
			"variant_reports/{seqid}_finished.txt",
			"output/jointvcf_all_variants_filtered_genotype_roi_meta_nomt/{seqid}_all_variants_filtered_genotype_roi_meta_nomt.vcf",
			expand("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_recalibration.csv", zip, sample_name=sample_names, sample_number=sample_numbers),
			seqid + "_combined_QC.txt",
			"{seqid}.relatedness2",
			"{seqid}_CollectVariantCallingMetrics.variant_calling_detail_metrics"
		output:
			"{seqid}.finished"
		shell:
			"touch {output}"

else:

	# Collect custom coverage and SV combined report
	if panel == "IlluminaTruSightCancer":

		rule final:
			input:
				"temp/validated_vcf/{seqid}.validated",
				expand("temp/manta/{sample_name}_{sample_number}_diploidSV.vcf.gz", zip, sample_name=sample_names, sample_number=sample_numbers),
				expand("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_Gaps.bed",zip, sample_name=sample_names, sample_number=sample_numbers),
				"temp/hotspot_coverage/custom.finished",
				"variant_reports/{seqid}_finished.txt",
				seqid + "_cnvReport.csv",
				"{seqid}_all_variants_filtered_genotype_roi_meta_nomt.vcf",
				seqid + "_combined_QC.txt",
				"{seqid}.relatedness2",
				"{seqid}_CollectVariantCallingMetrics.variant_calling_detail_metrics"
			output:
				"{seqid}.finished"
			shell:
				"touch {output}"

	else:


		rule final:
			input:
				"temp/validated_vcf/{seqid}.validated",
				"variant_reports/{seqid}_finished.txt",
				expand("{sample_name}/" + seqid + "_{sample_name}_{sample_number}_Gaps.bed",zip, sample_name=sample_names, sample_number=sample_numbers),
				"output/jointvcf_all_variants_filtered_genotype_roi_meta_nomt/{seqid}_all_variants_filtered_genotype_roi_meta_nomt.vcf",
				seqid + "_combined_QC.txt",
				"{seqid}.relatedness2",
				"{seqid}_CollectVariantCallingMetrics.variant_calling_detail_metrics"
			output:
				"{seqid}.finished"
			shell:
				"touch {output}"




















