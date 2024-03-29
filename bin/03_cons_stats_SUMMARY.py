#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
	import os
	import sys
	import argparse
	import pandas as pd
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
except ImportError:
	sys.exit("One of the required modules can't be found...")


#########################
## Input - parse arguments
#########################

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cons_stats_fn_list", help="List of consensus stats files --> paths", action = "store", nargs='+', default=[], required=True)
parser.add_argument("-o", "--out_dir", help="Directory to store the overview file, if unspecified assumes current work dir", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
consensus_stats_list = args.cons_stats_fn_list
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()


#########################
## Automated Analysis
#########################

df = pd.DataFrame(columns=['Indiv', 'Missing_data', 'Missing_data_PERCENTAGE'])

## Loop over consensus stats files, parse data and store
for indiv in consensus_stats_list:
	indiv_id = os.path.basename(indiv).split('_')[0]
	indiv_df = pd.read_csv(indiv, sep='\t', index_col=0)
	indiv_sum_df = pd.DataFrame([[indiv_id, indiv_df.loc['Total']['Missing_data'], indiv_df.loc['Total']['Missing_data_PERCENTAGE']]], columns=['Indiv', 'Missing_data', 'Missing_data_PERCENTAGE'])
	df = pd.concat([df, indiv_sum_df], ignore_index=True)

## Write out file
tsv_fn = os.path.join(output_path, ('all_indivs_missing_data.tsv'))
df.to_csv(tsv_fn, sep='\t', index=False)

# ## Now let's plot as well a basic graph to highlight missing data per indiv
# fig, axs = plt.subplots(figsize=(4, 12)) 
# df.plot(x="Indiv", y="Missing_data_PERCENTAGE", kind="barh", fontsize=4, legend=False, ax=axs)
# axs.set_xlim(0, 100)
# axs.set_title("Ratio missing data per individual")
# fig.savefig(os.path.join(output_path, "all_indivs_missing_data.pdf"), dpi=300)

import seaborn as sns
import matplotlib.pyplot as plt

## - Specify plot, if need be multiple tiers
sns.set_theme(style="whitegrid")

# Make the PairGrid
g = sns.PairGrid(df.sort_values("Missing_data_PERCENTAGE", ascending=False),
                 x_vars=df.columns[-1], y_vars=["Indiv"],
                 height=10)

# Draw a dot plot using the stripplot function
g.map(sns.stripplot, size=10, orient="h", jitter=False,
      palette="flare_r", linewidth=1, edgecolor="w")

# Use the same x axis limits on all columns and add better labels
g.set(xlim=(0, 100))

titles = ["Percentage missing data"]

for ax, title in zip(g.axes.flat, titles):

	ax.xaxis.grid(False)
	ax.yaxis.grid(True)

	ax.set_ylabel('Individuals', fontweight="bold", fontsize=18, labelpad=20)
	ax.set_xlabel('Missing data percentage', fontweight="bold", fontsize=18, labelpad=20)

sns.despine(left=True, bottom=True)

# Save output
plt.savefig(os.path.join(output_path, "all_indivs_missing_data.pdf"), bbox_inches='tight', dpi=300)
