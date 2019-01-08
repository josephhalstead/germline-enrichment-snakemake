"""
Takes the <sample>_QC.txt files nested within each sample folder and merges them into a file in the run folder.
Also includes the sex from the <sample>.variables file (taken from samplesheet).
Usage: takes run folder as argument:
python merge_qc_files.py <path/to/run/folder>
Author: Erik Waskiewicz 
"""


import os
import csv
import sys
import glob

# Get input filepath as first argument
filepath = sys.argv[1]

allfiles = os.listdir(filepath)


# open empty file to store data and add headers
output_file = open(os.path.join(filepath, 'combined_QC.txt'), 'w')
output_file.write('Sample\tRawSequenceQuality\tTotalReads\tEstimatedContamination\tGender\tSamplesheetSex\n')


for f in allfiles:
    sex = ""
    # Only loop through directories
    if os.path.isdir(os.path.join(filepath, f)):
        # open .variables file within each sample directory
        var_path = os.path.join(filepath, f, (f + '.variables'))
        #print (var_path)
        if os.path.isfile(var_path):
            var_file = open(var_path, 'r')
            output_file.write(f + '\t')
            for row in var_file:
                # look for sequence identifier and use it to open QC.txt file
                if row.startswith('seqId'):
                    qc_path = os.path.join(filepath, f, f + '*_QC.txt')
                    if os.path.isfile(qc_path) is False:
                        qc_path = os.path.join(filepath, f, f + '*_QC.txt')
                    qc_path = glob.glob(qc_path)

                    if len(qc_path) == 1:

                        qc_path = qc_path[0]

                    else:

                        raise

                    if os.path.isfile(qc_path):
                        with open(qc_path, 'r') as qc_file:
                            reader = csv.DictReader(qc_file, delimiter='\t')
                            # extract data from QC.txt file
                            for r in reader:
                                output_file.write(r['RawSequenceQuality'] + '\t') if r.get('RawSequenceQuality') \
                                    else output_file.write('\t')
                                output_file.write(r['TotalReads'] + '\t') if r.get('TotalReads') \
                                    else output_file.write('\t')
                                output_file.write(r['EstimatedContamination'] + '\t') \
                                    if r.get('EstimatedContamination') else output_file.write('\t')
                                output_file.write(r['Sex'] + '\t') if r.get('Sex') else output_file.write('\t')

                # look for sex in .variables file and append to data, else add empty field
                if row.startswith('sex'):
                    sex = row.rstrip()
            output_file.write(sex + '\n')


output_file.close()