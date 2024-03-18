#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
	import os
	import sys
	import argparse
	import pandas as pd
except ImportError:
	sys.exit("One of the required modules can't be found...")


#########################
## Input - parse arguments
#########################

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chromo_cov_tsv_list", help="List of samtools coverage output tsv paths", action = "store", nargs='+', default=[], required=True)
parser.add_argument("-i", "--individual", help="Individual to which coverage stats belong", required=True)
parser.add_argument("-o", "--out_dir", help="Directory to store the output file, if unspecified assumes current work dir", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
chromo_cov_fn_list = args.chromo_cov_tsv_list
indiv = args.individual
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()


#########################
## Automated Analysis
#########################

frames = []

## Loop over consensus seqs and estimate missing data
for chromo_tsv in chromo_cov_fn_list:
	chromo_df = pd.read_csv(chromo_tsv, sep='\t', header=0)
	frames.append(chromo_df)

indiv_df = pd.concat(frames)

## Write out file
tsv_fn = os.path.join(output_path, (indiv + '_coverage.tsv'))
indiv_df.to_csv(tsv_fn, sep='\t', index=False)
