""""

Purpose:

	The data about the samples in a run are stored in text file i.e. .variable files.

	These need to be parsed into a yaml config to be used with snakemake.	

	Takes a template config file and adds information from the variables files to it.


Usage:

	python scripts/make_config.py \
		--config configs/IlluminaTruSightCancer/IlluminaTruSightCancer_development_vm_template.yaml \
	 	--input input/ \
		--output test.yaml 

Requires:

	- Python 3
	- pyyaml

"""
import yaml
import argparse
import glob

# Get parser arguments
parser = argparse.ArgumentParser(description='Create a config file for a pipline run - requires a template config and variable files.')
parser.add_argument('--config', type=str, nargs=1, required=True,
					help='The location of the YAML config template')
parser.add_argument('--input', type=str, nargs=1, required=True,
					help='The folder with the sample data and variable files in.')
parser.add_argument('--output', type=str, nargs=1, required=True,
					help='The file to put the final output yaml in.')
args = parser.parse_args()
config = args.config[0]

# Parse YAML
config_dict = yaml.load(open(config))


# Get variable files

variable_files = glob.glob(f'{args.input[0]}/*/*.variables')

sample_info = {}
run_info = {}

for variable_file in variable_files:


	# Get sample id from file name
	sample_id = variable_file.split('/')[-1].split('.')[0]

	# Create a empty dict to store the sample specific information in
	sample_info[sample_id] = {}

	fh = open(variable_file, 'r')

	file_lines =  fh.readlines()

	for line in file_lines:

		processed_line = line.strip()

		# Get each characteristic from the variable file and add to dictionary

		# Sample level information

		if processed_line.split('=')[0] == 'worklistId':

			sample_info[sample_id]['worklistId'] = processed_line.split('=')[1].strip('"')

		elif processed_line.split('=')[0] == 'sex':

			sample_info[sample_id]['sex'] = processed_line.split('=')[1].strip('"')

		elif processed_line.split('=')[0] == 'familyId':

			sample_info[sample_id]['familyId'] = processed_line.split('=')[1].strip('"')

		elif processed_line.split('=')[0] == 'paternalId':

			sample_info[sample_id]['paternalId'] = processed_line.split('=')[1].strip('"')

		elif processed_line.split('=')[0] == 'maternalId':

			sample_info[sample_id]['maternalId'] = processed_line.split('=')[1].strip('"')		
			
		elif processed_line.split('=')[0] == 'phenotype':

			sample_info[sample_id]['phenotype'] = processed_line.split('=')[1].strip('"')

		# Worksheet level information

		elif processed_line.split('=')[0] == 'seqId':

			run_info['seqId'] = 	processed_line.split('=')[1].strip('"')

		elif processed_line.split('=')[0] == 'platform':

			run_info['platform'] = 	processed_line.split('=')[1].strip('"')

		elif processed_line.split('=')[0] == 'pipelineName':

			run_info['pipelineName'] = 	processed_line.split('=')[1].strip('"')

		elif processed_line.split('=')[0] == 'pipelineVersion':

			run_info['pipelineVersion'] = processed_line.split('=')[1].strip('"')

		elif processed_line.split('=')[0] == 'panel':

			run_info['panel'] = processed_line.split('=')[1].strip('"')	


	# If the variable file is null for any characteristic then add None as the value
	for characteristic in ['sex', 'familyId', 'paternalId', 'maternalId','phenotype']:

		if characteristic not in sample_info[sample_id]:

			sample_info[sample_id][characteristic] = None


# Add sample info to config dict
config_dict['sample_info'] = sample_info

# Add run info to config_dict
config_dict = {**config_dict, **run_info}

# Write config dict
output_yaml = args.output[0]

with open(output_yaml, 'w') as outfile:
    yaml.dump(config_dict, outfile, default_flow_style=False)






