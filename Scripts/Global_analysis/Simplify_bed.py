# script to simplify the promoter bed file for use in bedtools

import pandas as pd

file = pd.read_csv('data/refseqv1.0/refseqv1.0_promoters.bed', sep='\t', header=None, usecols=[0,1,2,5,7,9], names=['Chrom', 'Start', 'Stop', 'Strand', 'Type','ID'])

file = file.loc[file['Type'] == 'mRNA']  # Take annotation of only the genes
file.loc[:, 'Name'] = file['ID'].explode().str.rsplit('.').str[0].str.lstrip('ID=')
file.loc[:, 'Number'] = file['ID'].explode().str.rsplit('.').str[1].str[0]
file = file.loc[file['Number'] == '1']  # Only use the first mRNA
file['Score'] = 0
file = file.drop(columns = ['Number', 'Type'])
file.loc[:, 'Start'] = file.loc[:, 'Start'].astype('int')
file.loc[:, 'Stop'] = file.loc[:, 'Stop'].astype('int')

column_names = ['Chrom', 'Start', 'Stop', 'Name', 'Score', 'Strand']
file = file.reindex(columns = column_names)

file.to_csv('data/refseqv1.0/refseqv1.0_promoters_simple.bed', sep='\t', header=False, index=False)