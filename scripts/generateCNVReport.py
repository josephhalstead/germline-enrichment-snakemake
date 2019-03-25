import pandas
import allel
import argparse
from pathlib import Path

# Get the script arguments
parser = argparse.ArgumentParser()
parser.add_argument('--runid')
parser.add_argument('--output')
parser.add_argument('--bed')
parser.add_argument('--exome_metrics')
parser.add_argument('--manta_dir')
parser.add_argument('--exome_dir')
parser.add_argument('--coverage_dir')

args = parser.parse_args()
runid = args.runid
output_file = args.output
bed = args.bed
manta_dir = args.manta_dir
exome_dir = args.exome_dir
coverage_dir = args.coverage_dir
exome_metrics = args.exome_metrics

# Utility Functions
def fix_end_position(df):
	"""
	Fixes the END columns:

	If the END column is equal to -1 then change the end columns \
	to the start position + 1

	"""

	if df['END'] != 'NA':

		if int(df['END']) == -1:

			return int(df['POS']) + 1

	return df['END']

def assign_qc(df):
	"""
	Assign a QC status e.g. is r squared is less than a certain value

	"""

	qc_comment = []

	if df['Correlation'] < 0.98:

		qc_comment.append('R2<0.98')

	if df['depth'] < 160:

		qc_comment.append('Depth<160')

	if len(qc_comment) == 0:

		return 'PASS'

	else:

		return '|'.join(qc_comment)

def get_overlapping_genes(df, grouped_bed):
	"""
	Get the custom roi file and check to see if any of our SVs overlap with the regions.

	"""

	sv_start = df['POS']

	if sv_start == 'NA':

		return 'NA'

	sv_chrom = str(df['CHROM'])
	sv_start = int(df['POS'])
	sv_end = int(df['END'])

	gene_list = []

	for gene in grouped_bed.itertuples():

		gene_chrom = str(gene.Chrom)
		gene_start = int(gene.Start)
		gene_end = int(gene.End)

		# variant start is within the gene
		if gene_chrom == sv_chrom and gene_start <= sv_start <= gene_end:

			gene_list.append(gene.gene)

		# variant end is within the gene
		elif gene_chrom == sv_chrom and gene_start <= sv_end <= gene_end:

			gene_list.append(gene.gene)

		# variant start is before gene and variant end is after gene e.g. variant engulfs gene.
		elif gene_chrom == sv_chrom and (sv_start <= gene_start and sv_end >= gene_end):

			gene_list.append(gene.gene)

	if len(gene_list) > 0:

		return '|'.join(list(set(gene_list)))

	else:

		return 'NA'


# Create main dataframe
dataFrame = pandas.DataFrame()

# Read exome depth metrics file to get high coverage samples
exome_depth_metrics = pandas.read_table(exome_metrics) 
num_rows_exome_depth_metrics = exome_depth_metrics.shape[0]

# Add sampleid column to exome depth
exome_depth_metrics['sampleid'] = exome_depth_metrics['BamPath'].apply(lambda x: Path(x).stem.replace('_final', ''))

samples_list = list(exome_depth_metrics.groupby('sampleid').count().index)


# Read BED file
bed_file = pandas.read_csv(bed, names=['Chrom', 'Start', 'End', 'Comment'], sep='\t', header=None)

# Extract gene to new column
bed_file['gene'] = bed_file['Comment'].apply(lambda x: x.split('.')[0])

# Groupby gene and get the smallest start and largest end for each gene.
grouped_bed = bed_file.groupby(['Chrom', 'gene']).agg({'Start': 'min', 'End': 'max'})

grouped_bed = grouped_bed.reset_index()


#Loop through all samples and read in the relevant files
for sample in samples_list:

	sample_name = sample.split('_')[4]
	sample_number = sample.split('_')[5]
	
	# Get the mean depth from the summary file
	depth_of_coverage_summary = pandas.read_table(sample_name + '/' +  sample + '_DepthOfCoverage.sample_summary')
	depth_of_coverage_summary_mean = depth_of_coverage_summary['mean'][1]
	
	#Create dataframe from manta vcf
	manta_df = allel.vcf_to_dataframe(manta_dir + '/' + sample_name + '_' + sample_number + '_diploidSV.vcf.gz' , fields=['*'])

	if manta_df is not None:

		manta_df_2 = manta_df[['CHROM', 'POS', 'END', 'REF', 'SVTYPE']]
		manta_df_2.columns = ['CHROM', 'POS', 'END', 'REF', 'SVTYPE']
		manta_df_2['Regions'] = '-'
		manta_df_2['QC'] = 'PASS'
		manta_df_2['method'] = 'Manta'
		manta_df_2['ALT_1'] = manta_df_2['SVTYPE']

	else:


		# Make a dataframe with a single row with everything set to NA
		data = [{'CHROM':'NA', 'POS':'NA', 'END':'NA', 'REF':'NA', 'SVTYPE':'NA', 'Regions':'-'}]
		manta_df_2 = pandas.DataFrame(data)
		manta_df_2 = manta_df_2[['CHROM', 'POS', 'END', 'REF', 'SVTYPE', 'Regions']]
		manta_df_2['QC'] = 'PASS'
		manta_df_2['method'] = 'Manta'
		manta_df_2['ALT_1'] = manta_df_2['SVTYPE']

	
	#Create dataframe from Exome depth vcf
	exome_df = allel.vcf_to_dataframe(exome_dir + '/' + sample + '_final_cnv_fixed.vcf.gz' , fields=['*'])

	if exome_df is not None:

		exome_df_2 = exome_df[['CHROM', 'POS', 'END', 'REF', 'ALT_1', 'Regions']]
		exome_df_2['QC'] = 'PASS'
		exome_df_2['method'] = 'exomeDepth'

	else:
		exome_data = [{'CHROM':'NA', 'POS':'NA', 'END':'NA', 'REF':'NA', 'ALT_1':'NA', 'Regions':'NA'}]
		exome_df_2 = pandas.DataFrame(exome_data)
		exome_df_2 = exome_df_2[['CHROM', 'POS', 'END', 'REF', 'ALT_1', 'Regions']]
		exome_df_2['QC'] = 'PASS'
		exome_df_2['method'] = 'exomeDepth'


	#Combine both Manta and ExomeDepth dataframes and add depth and sample ID
	if (exome_df_2 is not None) and (manta_df_2 is not None):

		manta_exome= exome_df_2.append(manta_df_2)
		manta_exome['depth'] = depth_of_coverage_summary_mean
		manta_exome['sampleid'] = sample

	elif (exome_df_2 is None) and (manta_df_2 is not None):

		manta_exome= manta_df_2
		manta_exome['depth']= depth_of_coverage_summary_mean
		manta_exome['sampleid']= sample

	elif (exome_df_2 is not None) and (manta_df_2 is None):

		manta_exome= exome_df_2
		manta_exome['depth'] = depth_of_coverage_summary_mean
		manta_exome['sampleid'] = sample

	else:

		manta_exome = None

	if manta_exome is not None:

		manta_exome['END'] = manta_exome.apply(fix_end_position,axis=1)

		# Annotate with data from the exome depth logs
		manta_exome = pandas.merge(manta_exome, exome_depth_metrics, left_on = 'sampleid', right_on = 'sampleid')
		manta_exome['QC'] = manta_exome.apply(assign_qc, axis=1)

		# Add gene annotation
		manta_exome['Gene'] = manta_exome.apply(get_overlapping_genes, axis=1, args=(grouped_bed,))

		#Add sample dataframe to main dataframe containing all samples
		dataFrame = dataFrame.append(manta_exome)

#output the dataframe containing the CNVs for all sample as a TSV
if manta_exome is not None:

	final_df = dataFrame[['sampleid', 'method', 'CHROM', 'POS', 'END', 'REF', 'ALT_1', 'Regions', 'QC', 'depth', 'Gene']]

	# Remove < and > from <DUP> and <DEL>
	final_df['ALT_1'] = final_df['ALT_1'].apply(lambda x: x.strip('<').strip('>'))

	# Create new columns with different names
	final_df['Sample'] = final_df['sampleid']
	final_df['Method'] = final_df['method']
	final_df['Chr'] = final_df['CHROM']
	final_df['Start'] = final_df['POS']
	final_df['End'] = final_df['END']
	final_df['Ref'] = final_df['REF']
	final_df['Type'] = final_df['ALT_1']
	final_df['Depth'] = final_df['depth']

	# Write to TSV
	final_df[['Sample', 'Method', 'Chr', 'Start', 'End', 'Ref', 'Type', 'Regions', 'QC', 'Depth', 'Gene']].to_csv(output_file, index=False, sep=',')













