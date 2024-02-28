#!/usr/bin/env python
import bioinfokit
import pandas
import argparse
from bioinfokit.analys import norm


# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument

parser.add_argument('-i', required=True)
parser.add_argument('-o', required=True)

# Read arguments from command line
args = parser.parse_args()




# load dataset
df = pandas.read_csv(args.i, sep='\t')

# make contig column as index column (contig can be gene,contig,etc)
df = df.set_index('contig')

# now, normalize raw counts using RPKM method
# gene length must be in bp
nm = norm()
nm.rpkm(df=df, gl='length')

# get RPKM normalized dataframe
rpkm_df = nm.rpkm_norm

rpkm_df.to_csv(args.o, encoding='utf-8', sep='\t')
