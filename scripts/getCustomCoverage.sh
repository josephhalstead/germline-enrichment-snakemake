#!/bin/bash
set -euo pipefail

# script calculates gaps and coverage using original ROI files. This has been requested for
# truSightCancer as the method currently used assumes all exons in a given gene are included
# in the panel - which is often not the case.

# USE: bash inside run folder - will iterate over all samples. Depends on R script
# calculateTargetCoverage.R (located in same folder).

CUSTOM_BED=$1
DEPTH_FILE=$2
MIN_COV=$3
PASS_BED=$4
GAP_BED=$5
METRICS=$6
CLIN_COVERAGE_GENE=$7
CLIN_COVERAGE_TARGETS=$8


echo "making PASS file"

#Make PASS BED
tabix -R $CUSTOM_BED \
$DEPTH_FILE | \
awk -v minimumCoverage="$MIN_COV" '$3 >= minimumCoverage { print $1"\t"$2-1"\t"$2 }' | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > $PASS_BED

echo "making Gaps file"

#Make GAP BED
bedtools subtract \
-a $CUSTOM_BED \
-b $PASS_BED | \
sort -k1,1V -k2,2n -k3,3n \
> $GAP_BED

echo "making clincoverage file"

bedtools coverage \
-a $CUSTOM_BED \
-b $PASS_BED | \
tee $METRICS | \
awk '{pass[$4]+=$6; len[$4]+=$7} END { for(i in pass) printf "%s\t %.2f%\n", i, (pass[i]/len[i]) * 100 }' | \
sort -k1,1 > $CLIN_COVERAGE_GENE

echo "start R"

Rscript --vanilla scripts/calculateTargetCoverage.R \
$METRICS \
$CLIN_COVERAGE_TARGETS \
$GAP_BED

