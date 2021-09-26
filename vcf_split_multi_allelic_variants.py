import logging
import os
import re

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S %Z")
logger = logging.getLogger(__name__)

VERSION = 'VERSION: 1.0.0'
logger.info(VERSION)

# Function to run shell command
def run_job(cmd):
    logger.info(f'cmd: {cmd}')
    try:
        subprocess.check_call(cmd, shell=True)
        logger.info(f'Command: {cmd} ran successfully.\n\n')
    except subprocess.CalledProcessError as e:
        logger.error(f'Failed to run: {cmd}')
        logger.error(f'Error: {e}')
        sys.exit(1)


# If empty vcf file is provided, generate an empty output vcf file
def isvcf_empty(invcf, outvcf):
	if os.path.getsize(invcf) == 0:
		logger.info('Empty input vcf file is passed. Empty output vcf file will be generated')
		run_job(f'touch {outvcf}')
		exit(0)

# These two dict will store header metadata for info and format tags, mainly tagtype whether A/R/./\d+
info_header_dict = {}
format_header_dict = {}

# Function to parse header metadata for info and format field
def fetch_metadata_for_info_and_format(line):
	res = re.search('^##INFO=<ID=(.*),Number=(.*),Type=(.*),', line)
	if res:
		tempval = res.groups()
		info_header_dict[tempval[0]] = [tempval[1], tempval[2]]
	res = re.search('##FORMAT=<ID=(.*),Number=(.*),Type=(.*),', line)
	if res:
		tempval = res.groups()
		format_header_dict[tempval[0]] = [tempval[1], tempval[2]]

# Function to split genotype for multi-allelic variants
def process_genotype(alt, genotype):
	processed_genotype_list = ['./.' for i in range(len(alt))]

	# Keep genotype as is
	if genotype == "." or genotype == './.' or genotype == '.|.':
		for i in range(0, len(alt)):
			processed_genotype_list[i] = f'{genotype}'
	elif re.search('[0-9.]+\/[0-9.]+', genotype) or re.search('/[0-9.]+\|[0-9.]+/', genotype):
		(gt1, gt2) = re.split('\/|\|', genotype)
		res = re.search('[0-9.]+(.)[0-9.]+', genotype)
		phase = res.groups()[0]

		# Homozygous
		if gt1 == gt2:
			if gt1 == '0':
				# If 0/0, then when split, genotype would be 0/0 for all alt
				processed_genotype_list = [f'0{phase}0' for i in range(len(alt))]
			else:
				# If 1/1, then when split, for first alt allele (0th index of processed_genotype_list) put 1/1. rest ./.
				# If 2/2, then when split, for second alt allele (1th index of processed_genotype_list) put 1/1. rest ./.
				processed_genotype_list[int(gt1)-1] = f'1{phase}1'
		# Heterozygous
		else:
			# If 0/1 then one is ref, other is 1st alt allele. so for 1st alt allele (0th index) vcf row it is 0/1 else ./.
			# If 0/2 then one is ref, other is 2nd alt allele. so for 2nd alt allele (1st index) vcf row it is 0/1 else ./.
			# If 0/3 then one is ref, other is 3rd alt allele. so for 3rd alt allele (1st index) vcf row it is 0/1 else ./.
			if gt1 == '0':
				processed_genotype_list[int(gt2)-1] = '0{phase}1'
			elif gt2 == '0':
				# If 1/0 then first alt allele / ref allele. so for 1st alt allele (0th index) vcf row it is 0/1 else ./.
				# If 2/0 then second alt allele / ref allele. so for 2nd alt allele (1st index) vcf row it is 0/1 else ./.
				processed_genotype_list[int(gt1)-1] = '1{phase}0'
			elif gt1!='0' and gt2 != '0':
				# 1/3: first allele / third allele. Since we do w.r.t reference and referene allele is not genotyped 
				# first allele vcf row will be: 1/.
				# third allele vcf row will be: ./1
				processed_genotype_list[int(gt1)-1] = '1{phase}.'
				processed_genotype_list[int(gt2)-1] = '.{phase}1'
	return processed_genotype_list

# Function to split format values
def process_format(alt, format_field, sample_field):
	processed_format_list = [[] for i in range(len(alt))]
	
	format_tags = format_field.split(':')
	format_values = sample_field.split(':')
	for i in range(0, len(format_tags)):
		tag = format_tags[i]
		value = format_values[i]
		values = value.split(',')
		if tag != 'GT':
			if tag not in format_header_dict:
				logger.error(f'Header information is not present for the tag {tag}')
				exit(1)
			tagtype = format_header_dict[tag][0]
			try:
				if tagtype == 'A':
					for i in range(0, len(alt)):
						processed_format_list[i].append(f'{values[i]}')
				elif tagtype == 'R':
					for i in range(0, len(alt)):
						processed_format_list[i].append(f'{values[0]},{values[i+1]}')
				elif tagtype == '.' or re.search('\d+', tagtype):
					# 'A'
					if len(values) == len(alt):
						for i in range(0, len(alt)):
							processed_format_list[i].append(f'{values[i]}')
					# 'R'
					elif len(values) == len(alt) + 1:
						for i in range(0, len(alt)):
							processed_format_list[i].append(f'{values[0]},{values[i+1]}')
					# If any other cases, then don't split the tag and retain information as is.
					else:
						for i in range(0, len(alt)):
							processed_format_list[i].append(f'{value}')
				else:
					for i in range(0, len(alt)):
						processed_format_list[i].append(f'{value}')
			except:
				logger.warn(f'For the FORMAT tag {tag}, the tagtype is {tagtype}. Total alt record is {len(alt)} but some issues observed with its value: {value}\nMaybe #Alt alleles != #values. So not spliting this tag')
				for i in range(0, len(alt)):
					processed_format_list[i].append(f'{value}')
		else:
			# Tag is GT and process GT differently
			processed_genotype_list = process_genotype(alt, value)
			for i in range(0, len(alt)):
				processed_format_list[i].append(f'{processed_genotype_list[i]}')
	return processed_format_list

# Function to split INFO field
def process_info(alt, info_field):
	# e.g depending on the tagtype (A/R/.), values will be processed
	
	# Suppose tagtype is A
	# chr1	24152721	.	G	GA,GAA	3184.65	.	AF=0.5,0.6
	# 2 alt allele and 2 values per alt
	# When split, GA will get 0.5 and GAA will get 0.6

	# Suppose tagtype is R
	# chr1	24152721	.	G	GA,GAA	3184.65	.	AF=0.5,0.6,0.7
	# 2 alt allele and 3 values so 1st one is ref and next two is for alt
	# When split, GA will get 0.5,0.6 and GAA will get 0.5,0.7

	# For tagtype as dot or tagtype as number, duduce A/R based on counting
	# If len(alt) == len(values), Type A
	# If len(alt)+1 ==  len(values), Type R
	# Any other case means keep value as is. dont split.`

	processed_info_list = [[] for i in range(len(alt))]

	# Loop over each info tag
	for tag in info_field.split(';'):
		if ',' in tag:
			# AF=0.5,0.6
			(key, value) = tag.split('=')
			values = value.split(',')
			if key not in info_header_dict:
				logger.error(f'Header information is not present for the tag {key}')
				exit(1)
			# tagtype 'A' or 'R' or '.'
			tagtype = info_header_dict[key][0]
			try:
				if tagtype == 'A':
					for i in range(0, len(alt)):
						processed_info_list[i].append(f'{key}={values[i]}')
				elif tagtype == 'R':
					for i in range(0, len(alt)):
						processed_info_list[i].append(f'{key}={values[0]},{values[i+1]}')
				elif tagtype == '.' or re.search('\d+', tagtype):
					# 'A'
					if len(values) == len(alt):
						for i in range(0, len(alt)):
							processed_info_list[i].append(f'{key}={values[i]}')
					# 'R'
					elif len(values) == len(alt) + 1:
						for i in range(0, len(alt)):
							processed_info_list[i].append(f'{key}={values[0]},{values[i+1]}')
					# If any other cases, then don't split the tag and retain information as is.
					else:
						for i in range(0, len(alt)):
							processed_info_list[i].append(f'{key}={value}')
				else:
					for i in range(0, len(alt)):
						processed_info_list[i].append(f'{key}={value}')
			except:
				logger.warn(f'For the INFO tag {key}, the tagtype is {tagtype}, but some issues observed with its value: {value}\nMaybe #Alt alleles != #values. So not spliting this tag')
				for i in range(0, len(alt)):
					processed_info_list[i].append(f'{key}={value}')
		else:
			# AN=2 or Flag type of tag such as HS
			for i in range(len(alt)):
				processed_info_list[i].append(tag)
	return(processed_info_list)


def main(args):
	logger.info(f'Input vcf file: {args.invcf}')
	logger.info(f'Output vcf file: {args.outvcf}')
	isvcf_empty(args.invcf, args.outvcf)

	format_column_present = None
	multi_allelic_variant_count = 0
	with open(args.invcf) as fp, open(args.outvcf, 'w') as fw:
		for line in fp.readlines():
			line = line.strip()
			if line.startswith('##'):
				fw.write(line+'\n')
				fetch_metadata_for_info_and_format(line)
			elif line.startswith('#CHROM'):
				fw.write(line+'\n')
				total_columns = len(line.split('\t'))
				if total_columns < 8:
					logger.info(f'Total column in VCF file observed {total_columns} is < 8. Minimum 8 columns required')
					exit(1)
				elif total_columns == 8:
					logger.info('VCF file has no FORMAT and SAMPLE column')
					format_column_present = False
				elif total_columns == 10:
					logger.info('VCF file is standard VCF file with single sample')
					format_column_present = True
				elif total_columns > 10:
					logger.info('Multi sample VCF file is observed. Only first sample will be processed. Please split the vcf file sample wise')
					format_column_present = True
				else:
					logger.info(f'Incorrect number of columns {total_columns} found in the VCF file')
					exit(1)
			else:
				fields = line.split('\t')
				info_field = fields[7]
				if re.search(',', fields[4]):
					multi_allelic_variant_count += 1
					alt = fields[4].split(',')
					logger.info(f'Multiallelic variant: {fields[0]}/{fields[1]}/{fields[3]}/{fields[4]}\tTotal alleles: {len(alt)}')
					processed_info_list = process_info(alt, info_field)
					for tempind in range(0, len(alt)):
						fw.write('\t'.join(fields[0:4])+f'\t{alt[tempind]}\t{fields[5]}\t{fields[6]}\t'+';'.join(processed_info_list[tempind])+'\t')
						if format_column_present:
							processed_format_list = process_format(alt, fields[8], fields[9])
							fw.write(f'{fields[8]}\t{":".join(processed_format_list[tempind])}\n')
						else:
							fw.write(f'\n')
				else:
					if format_column_present:
						fw.write('\t'.join(fields[0:10])+'\n')
					else:
						fw.write('\t'.join(fields[0:8])+'\n')
	logger.info('-\n')
	logger.info(f'Total multi_allelic_variants observed in input vcf file is {multi_allelic_variant_count}');
	logger.info('Done')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Script to split multi allelic variants to individual rows in output vcf file")
    parser.add_argument("--invcf", type=str, required=True, help="Input VCF file")
    parser.add_argument("--outvcf", type=str, required=False, help="Output processed vcf file", default = "single_allelic_out.vcf")
    args = parser.parse_args()
    main(args)


