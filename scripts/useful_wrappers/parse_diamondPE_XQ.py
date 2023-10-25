#!/usr/bin/env python
# this script is used to merge the diamond blast output for paired-end reads
# the rational is, it will count the R2 taxon into account only if there was nohit in R1 taxon.

__author__ = "xiaoqiong Gu"
__version__ = "0.1"

import argparse
import pandas as pd

# parse the arguments
parser = argparse.ArgumentParser(description='Parse Diamond R1 & R2 blast outputs to a count table in CSV format')
parser.add_argument('-1', '--reads1', help='Diamond output for R1 reads')
parser.add_argument('-2', '--reads2', help='Diamond output for R2 reads')
parser.add_argument('-k','--kraken', help='kraken2 output for R1&R2 reads')
parser.add_argument('-o', '--outfile', default='arg_taxon_table.csv', help='Count table')
args = parser.parse_args()


PE1 = pd.read_csv(args.reads1, sep='\t', index_col=0, header=None)
PE2 = pd.read_csv(args.reads2, sep='\t', index_col=0, header=None)
kraken = pd.read_csv(args.kraken,sep='\t')

df_inner = PE1.merge(PE2, left_index=True, right_index=True, how="inner")
df_outer = PE1.merge(PE2, left_index=True, right_index=True, how="outer")
df_merge=pd.DataFrame()
df_merge['ARG']=df_inner['1_x']
df_added = df_outer[df_outer['1_x'].isnull()]['1_y'].to_frame()
df_added = df_added.rename(columns={'1_y':'ARG'})
frame = [df_merge, df_added]
arg = pd.concat(frame)
arg.index.name = 'ID'
arg.reset_index()

kraken = kraken.iloc[:, 0:3]
kraken.columns = ['level', 'ID', 'taxon']
result = arg.merge(kraken, left_on='ID',right_on='ID',how='inner')
result = result.groupby(['ARG','taxon']).size().to_frame()
result = result.reset_index()
result.columns = ['ARG','taxon','count']
result = result.sort_values(ascending=False, by='count')
result.to_csv(args.outfile, index=False)

