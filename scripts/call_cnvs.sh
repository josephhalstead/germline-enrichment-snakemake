
BAM_LIST=$1
REF=$2
BED=$3
SEQID=$4
PREFIX=$5
DICT=$6

if [[ -e $BAM_LIST ]] && [[ $(wc -l $BAM_LIST | awk '{print $1}') -gt 4 ]]; then


	 #call CNVs using read depth
    Rscript scripts/ExomeDepth.R \
    -b $BAM_LIST\
    -f $REF \
    -r $BED \
    -p $PREFIX \
    2>&1 | tee $PREFIX/ExomeDepth.log

    #print ExomeDepth metrics
    echo -e "BamPath\tFragments\tCorrelation" > $PREFIX/"$SEQID"_ExomeDepth_Metrics.txt
    paste $BAM_LIST \
    <(grep "Number of counted fragments" $PREFIX/ExomeDepth.log | cut -d' ' -f6) \
    <(grep "Correlation between reference and tests count" $PREFIX/ExomeDepth.log | cut -d' ' -f8) >> $PREFIX/"$SEQID"_ExomeDepth_Metrics.txt

 	#add CNV vcf headers and move to sample folder
    for vcf in $(ls $PREFIX/*_cnv.vcf); do

        sampleId=$(basename ${vcf%.*})

        #add VCF headers
        picard UpdateVcfSequenceDictionary \
        I="$vcf" \
        O="$PREFIX"/"$sampleId"_fixed.vcf \
        SD=$DICT

        #gzip and tabix
        bgzip "$PREFIX"/"$sampleId"_fixed.vcf
        tabix -p vcf "$PREFIX"/"$sampleId"_fixed.vcf.gz

    done

else

	touch $PREFIX/no_cnvs_final_cnv.vcf;
	touch $PREFIX/no_cnvs_final_cnv.txt;


fi
