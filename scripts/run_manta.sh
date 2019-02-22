#!/bin/bash
set -euo pipefail

#Bash script for launching manta. Only actually runs manta if the bam is in the high coverage list
#Otherwise we just create empty files

HIGH_COV_BAMS=$1
BAM=$2
FASTA=$3
RUN_DIR=$4
THREADS=$5
SAMPLE=$6
FINAL_DIR=$7
LOW_COV_BAMS=$8

# Loop through high coverage bams
for i in $(cat $HIGH_COV_BAMS); do

	# If we find our bam in the list
	if [ "$i" == "$BAM" ]
	then

		configManta.py \
		--bam $BAM \
		--referenceFasta $FASTA \
		--exome \
		--runDir $RUN_DIR \

		"$RUN_DIR"/runWorkflow.py \
		--quiet \
		-m local \
		-j $THREADS \

		# Copy results to directory
		cp "$RUN_DIR"/results/variants/diploidSV.vcf.gz "$FINAL_DIR"/"$SAMPLE"_diploidSV.vcf.gz
		cp "$RUN_DIR"/results/variants/diploidSV.vcf.gz.tbi "$FINAL_DIR"/"$SAMPLE"_diploidSV.vcf.gz.tbi

		#rm -r "$RUN_DIR"

	fi

done

# Now check whether our bam is low coverage - if so create empty files.
for i in $(cat $LOW_COV_BAMS); do

	if [ "$i" == "$BAM" ]
	then

		echo "NO_MANTA" > "$FINAL_DIR"/"$SAMPLE"_diploidSV.vcf.gz
		echo "NO_MANTA" > "$FINAL_DIR"/"$SAMPLE"_diploidSV.vcf.gz.tbi

	fi

done