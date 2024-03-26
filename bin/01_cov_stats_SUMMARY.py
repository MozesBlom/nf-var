#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
	import os
	import sys
	import argparse
	import pandas as pd
	import numpy as np
	import seaborn as sns
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches
except ImportError:
	sys.exit("One of the required modules can't be found...")


#########################
## Input - parse arguments
#########################

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--indiv_cov_tsv_list", help="Samtools coverage summary per individual, fn list", action = "store", nargs='+', default=[], required=True)
parser.add_argument("-s", "--sex_chr_list", help="List of sex chromosomes if any", action = "store", nargs='+', default=[], required=False)
parser.add_argument("-p", "--plots_per_panel", help="How many indivs to plot next to each other?", action = "store", default=10, required=False)
parser.add_argument("-o", "--out_dir", help="Directory to store the output file, if unspecified assumes current work dir", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
indiv_cov_fn_list = args.indiv_cov_tsv_list
sex_chr_list = args.sex_chr_list
number_per_panel = int(args.plots_per_panel)
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()


#########################
## Automated Analysis
#########################

df = pd.DataFrame(columns=['#rname', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq', 'individual'])

## Loop over consensus seqs and estimate missing data
for indiv_tsv in indiv_cov_fn_list:
	indiv_id = os.path.basename(indiv_tsv).split('_coverage')[0]
	indiv_df = pd.read_csv(indiv_tsv, sep='\t', header=0)
	indiv_df['individual'] = str(indiv_id)
	df = pd.concat([df, indiv_df], ignore_index=True)

# Finally, add another column where we specify if the chromosome is a sex chromo or autosome
# wacky hack to remove weird string characters from input
sex_chr_list_updated = []

for chromo in sex_chr_list:
	x = chromo.replace('[', '')
	y = x.replace(']', '')
	z = y.replace(',', '')
	sex_chr_list_updated.append(z)

def find_chromo_type(chromo):
	if chromo in sex_chr_list_updated:
		return 'sex chromosome'
	else:
		return 'autosome'

df['chromo_type'] = df['#rname'].apply(find_chromo_type)

# List of individuals in df
indivs = df['individual'].unique()

# Summary with averages per individual and chromo type
frames = []
# Looping another time over the df to calculate averages per individual. This second round is needed because we included chromo type now
for indiv in indivs:
	df_indiv = df.loc[df['individual'] == indiv]
	if len(df['chromo_type'].unique()) > 1:
		df_indiv_grouped = df_indiv.groupby('chromo_type', as_index=False).mean(numeric_only=True)
		df_indiv_grouped['individual'] = str(indiv)
		frames.append(df_indiv_grouped[['individual', 'chromo_type', 'coverage', 'meandepth']])
	else:
		chromo_type_df = df_indiv.mean(numeric_only=True).to_frame().T
		chromo_type_df['individual'] = str(indiv)
		chromo_type_df['chromo_type'] = df['chromo_type'].unique()
		frames.append(chromo_type_df[['individual', 'chromo_type', 'coverage', 'meandepth']])

df_summary = pd.concat(frames)

## Write out file
tsv_fn = os.path.join(output_path, ('all_indivs_summary.tsv'))
df_summary.to_csv(tsv_fn, sep='\t', index=False)


#########################
##
## Plot Coverage distributions
##
#########################

# Yield successive n-sized 
# chunks from l. 
def divide_chunks(l, n):
    # looping till length l 
	for i in range(0, len(l), n):
		yield l[i:i + n] 


## - Specify plot, if need be multiple tiers
sns.set_style("whitegrid", {'axes.grid' : False})
plot_nrows = float(len(indivs)/int(number_per_panel))
if plot_nrows > 1:
	if int(plot_nrows) < plot_nrows:
		fig, axes = plt.subplots(nrows=int(plot_nrows)+1, sharex=False, figsize=(19.6,19.6))
	else:
		fig, axes = plt.subplots(nrows=int(plot_nrows), sharex=False, figsize=(19.6,19.6))
	indivs_panel_list = list(divide_chunks(indivs, number_per_panel))
else:
	fig, axes = plt.subplots(nrows=1, sharex=False, figsize=(19.6,19.6))
	indivs_panel_list = [indivs]
fig.subplots_adjust(hspace=0.2)


## Identify the max value for depth
upper_limit_depth = int(df['meandepth'].max())

## Now loop over indiv lists, subset the dataframe and create (multi-tiered) violin plots
ax_counter = 0
for indivs_list in indivs_panel_list:
	df_panel = df[df['individual'].isin(indivs_list)]
	if len(indivs_panel_list) > 1:
		plot = sns.violinplot(data=df_panel, x="individual", y="meandepth", hue="chromo_type",
						split=True, inner="quart", fill=False, palette = {"autosome": "#9b59b6", "sex chromosome": "#3498db"},
						ax=axes[ax_counter], legend=False)
		axes[ax_counter].set(ylim=(0, upper_limit_depth))
		xtick_labels = plot.get_xticklabels()
		ytick_labels = plot.get_yticklabels()
		axes[ax_counter].set_xticklabels(xtick_labels, fontweight="bold", fontsize=12, rotation=45)
		axes[ax_counter].set_yticklabels(ytick_labels, fontweight="bold", fontsize=12)
		axes[ax_counter].set_xlabel('', fontsize=1, labelpad=20)
		axes[ax_counter].set_ylabel('Mean depth', fontweight="bold", fontsize=18, labelpad=20)
	else:
		plot = sns.violinplot(data=df_panel, x="individual", y="meandepth", hue="chromo_type",
						split=True, inner="quart", fill=False, palette = {"autosome": "#9b59b6", "sex chromosome": "#3498db"},
						legend=False)
		axes.set(ylim=(0, upper_limit_depth))
		xtick_labels = plot.get_xticklabels()
		ytick_labels = plot.get_yticklabels()
		axes.set_xticklabels(xtick_labels, fontweight="bold", fontsize=12, rotation=45)
		axes.set_yticklabels(ytick_labels, fontweight="bold", fontsize=12)
		axes.set_xlabel('Individuals', fontweight="bold", fontsize=18, labelpad=20)
		axes.set_ylabel('Mean depth', fontweight="bold", fontsize=18, labelpad=20)
	ax_counter += 1

## Now format a bit
sns.despine()
if ax_counter > 1:
	axes[ax_counter-1].set_xlabel('Individuals', fontweight="bold", fontsize=18, labelpad=20)
else:
	pass

# Include custom legend
autosome = mpatches.Patch(color='#9b59b6', label="autosome")
sex_chromosome = mpatches.Patch(color='#3498db', label="sex chromosome")
if plot_nrows > 1:
	axes[0].legend(handles=[autosome, sex_chromosome], bbox_to_anchor=(1.05, 0.8), fontsize=14, frameon=False)
else:
	axes.legend(handles=[autosome, sex_chromosome], bbox_to_anchor=(1.05, 0.8), fontsize=14, frameon=False)

plt.tight_layout()

# Save output
plt.savefig(os.path.join(output_path, ('all_indivs_summary.pdf')), format="pdf", bbox_inches = "tight", dpi=300)
