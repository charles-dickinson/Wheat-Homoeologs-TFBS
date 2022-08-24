# change bed file produced by first part of get_promoters bash script into 
# a simple bed file (required for bedtools getfasta)

import pandas as pd

file = pd.read_csv('data/refseqv1.1/refseqv1.1_promoters.bed', sep = '\t', header=None, usecols = [0,1,2,5,7,9], names = ['Chrom', 'Start', 'Stop', 'Strand', 'Type',
                                                                              'ID'])                                                                             
#two_one_module_genes = pd.read_csv('data/circadian/produced/two_one_module_genes/two_one_modules_phase_per_gene.tsv',
                                    #sep = '\t', header=0)

#two_module_gene_list = two_one_module_genes['gene']

all_phase_genes = pd.read_csv('data/circadian/produced/all_module_genes/all_module_genes_list.csv',
                              sep = '\t', header=0)

all_phase_gene_list = all_phase_genes['gene']

file = file.loc[file['Type'] == 'mRNA']  # Take annotation of only the genes
file.loc[:, 'Name'] = file['ID'].explode().str.rsplit('.').str[0].str.lstrip('ID=')
file.loc[:, 'Number'] = file['ID'].explode().str.rsplit('.').str[1].str[0]
file = file.loc[file['Number'] == '1']  # Only use the first mRNA
file['Score'] = 0 # Bed files want scores - 0 just as placemarker
file = file.drop(columns = ['Number', 'Type'])
file.loc[:, 'Start'] = file.loc[:, 'Start'].astype('int')
file.loc[:, 'Stop'] = file.loc[:, 'Stop'].astype('int')


#file = file[file.Name.isin(two_module_gene_list)]
file = file[file.Name.isin(all_phase_gene_list)]
column_names = ['Chrom', 'Start', 'Stop', 'Name', 'Score', 'Strand']
file = file.reindex(columns = column_names)
file = file.drop_duplicates(subset = 'Name', keep='first')
file.to_csv('data/refseqv1.1/refseqv1.1_promoters_all_module_genes.bed',
            sep = '\t', header = False, index = False)
#file.to_csv('data/refseqv1.1/refseqv1.1_promoters_two_one_module_genes.bed,
#            sep = '\t', header = False, index = False)
