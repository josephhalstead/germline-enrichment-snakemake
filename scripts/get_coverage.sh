#!/bin/bash
set -euo pipefail

#Legacy script from the old GermlineEnrichment Pipline
#Only to be used as part of the snakemake script

ROI_BED=$1
PANEL=$2
HOTSPOTS=$3
DEPTH=$4
SEQID=$5
SAMPLEID=$6
MINIMUM_COVERAGE=$7
REF_SEQ=$8
PREFIX=$9

#Make BED file of all genes overlapping ROI
bedtools intersect -wa \
-a $REF_SEQ \
-b $ROI_BED | \
awk -F "\t" '$3 == "gene" { print $1"\t"$4-1"\t"$5 }' | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > "$PREFIX"/"$PANEL"_TargetGenes.bed

#make target bed
bedtools intersect \
-a $REF_SEQ \
-b "$PREFIX"/"$PANEL"_TargetGenes.bed | \
grep "NM_[0-9]*\.[0-9]*" | \
awk -F "\t" '$3 == "exon" { print $1"\t"$4-1"\t"$5 }' | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > "$PREFIX"/"$PANEL"_Targets.bed

awk '$1 ~ /^MT/' $ROI_BED >> "$PREFIX"/"$PANEL"_Targets.bed

#Intersect CDS for all genes, pad by p=n and merge coordinates by gene
bedtools intersect \
-a $REF_SEQ \
-b "$PREFIX"/"$PANEL"_Targets.bed | \
awk -F'[\t|;|=]' -v p=5 '$3 == "CDS" { gene="null"; for (i=9;i<NF;i++) if ($i=="gene"){gene=$(i+1); break}; genes[gene] = genes[gene]$1"\t"($4-1)-p"\t"$5+p"\t"gene";" } END { for (gene in genes) print genes[gene] }' | \
while read line; do
    echo "$line" | \
    tr ';' '\n'| \
    sort -k1,1V -k2,2n -k3,3n | \
    bedtools merge -c 4 -o distinct;
done | \
sort -k1,1V -k2,2n -k3,3n > "$PREFIX"/"$PANEL"_ClinicalCoverageTargets.bed

cat $HOTSPOTS "$PREFIX"/"$PANEL"_ClinicalCoverageTargets.bed | \
sort -k1,1V -k2,2n -k3,3n > "$PREFIX"/"$PANEL"_ClinicalCoverageTargetsHotspots.bed

#Make PASS BED
tabix -R "$PREFIX"/"$PANEL"_ClinicalCoverageTargetsHotspots.bed \
$DEPTH | \
awk -v minimumCoverage="$MINIMUM_COVERAGE" '$3 >= minimumCoverage { print $1"\t"$2-1"\t"$2 }' | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > "$PREFIX"/"$SEQID"_"$SAMPLEID"_PASS.bed

#Calculate overlap between PASS BED and ClinicalCoverageTargets
bedtools coverage \
-a "$PREFIX"/"$PANEL"_ClinicalCoverageTargetsHotspots.bed \
-b "$PREFIX"/"$SEQID"_"$SAMPLEID"_PASS.bed | \
tee "$PREFIX"/"$SEQID"_"$SAMPLEID"_ClinicalCoverageTargetMetrics.txt | \
awk '{pass[$4]+=$6; len[$4]+=$7} END { for(i in pass) printf "%s\t %.2f%\n", i, (pass[i]/len[i]) * 100 }' | \
sort -k1,1 > "$PREFIX"/"$SEQID"_"$SAMPLEID"_ClinicalCoverageGeneCoverage.txt

#Make GAP BED
bedtools subtract \
-a "$PREFIX"/"$PANEL"_ClinicalCoverageTargetsHotspots.bed \
-b "$PREFIX"/"$SEQID"_"$SAMPLEID"_PASS.bed | \
sort -k1,1V -k2,2n -k3,3n \
> "$PREFIX"/"$SEQID"_"$SAMPLEID"_Gaps.bed

# Remove unneeded files
rm "$PREFIX"/"$PANEL"_TargetGenes.bed
rm "$PREFIX"/"$PANEL"_Targets.bed
rm "$PREFIX"/"$PANEL"_ClinicalCoverageTargets.bed
rm "$PREFIX"/"$PANEL"_ClinicalCoverageTargetsHotspots.bed
rm "$PREFIX"/"$SEQID"_"$SAMPLEID"_PASS.bed
