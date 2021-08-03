#!/usr/bin/env python3

import datetime
import re
import pandas as pd
import skbio


df = pd.read_csv("$summary", delimiter="\t",  names=['seq_id', 'tax_id', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
tax_table = df['class'].value_counts()
tax_table_class = pd.DataFrame(list(zip(tax_table.index.tolist(), tax_table.tolist())), columns =['tax_id', 'read_count'])
tax_table_class.to_csv("blast_report_class.csv", header=1, columns=["tax_id", "read_count"], index=False)

tax_table = df['order'].value_counts()
tax_table_order = pd.DataFrame(list(zip(tax_table.index.tolist(), tax_table.tolist())), columns =['tax_id', 'read_count'])
tax_table_order.to_csv("blast_report_order.csv", header=1, columns=["tax_id", "read_count"], index=False)

tax_table = df['family'].value_counts()
tax_table_family = pd.DataFrame(list(zip(tax_table.index.tolist(), tax_table.tolist())), columns =['tax_id', 'read_count'])
tax_table_family.to_csv("blast_report_family.csv", header=1, columns=["tax_id", "read_count"], index=False)

tax_table = df['genus'].value_counts()
tax_table_genus = pd.DataFrame(list(zip(tax_table.index.tolist(), tax_table.tolist())), columns =['tax_id', 'read_count'])
tax_table_genus.to_csv("blast_report_genus.csv", header=1, columns=["tax_id", "read_count"], index=False)

tax_table = df['species'].value_counts()
tax_table_species = pd.DataFrame(list(zip(tax_table.index.tolist(), tax_table.tolist())), columns =['tax_id', 'read_count'])
tax_table_species.to_csv("blast_report_species.csv", header=1, columns=["tax_id", "read_count"], index=False)

shannon = skbio.diversity.alpha.shannon(tax_table_class['read_count'])
simpson = skbio.diversity.alpha.simpson(tax_table_class['read_count'])
file1 = open("blast_diversity_class.csv","w") 
file1.write(str(tax_table_class['read_count'].sum()) + "," + str(round(shannon,3)) + "," + str(round(simpson,3)) + "\\n") 
file1.close()

shannon = skbio.diversity.alpha.shannon(tax_table_order['read_count'])
simpson = skbio.diversity.alpha.simpson(tax_table_order['read_count'])
file1 = open("blast_diversity_order.csv","w") 
file1.write(str(tax_table_order['read_count'].sum()) + "," + str(round(shannon,3)) + "," + str(round(simpson,3)) + "\\n") 
file1.close()

shannon = skbio.diversity.alpha.shannon(tax_table_family['read_count'])
simpson = skbio.diversity.alpha.simpson(tax_table_family['read_count'])
file1 = open("blast_diversity_family.csv","w") 
file1.write(str(tax_table_family['read_count'].sum()) + "," + str(round(shannon,3)) + "," + str(round(simpson,3)) + "\\n") 
file1.close()

shannon = skbio.diversity.alpha.shannon(tax_table_genus['read_count'])
simpson = skbio.diversity.alpha.simpson(tax_table_genus['read_count'])
file1 = open("blast_diversity_genus.csv","w") 
file1.write(str(tax_table_genus['read_count'].sum()) + "," + str(round(shannon,3)) + "," + str(round(simpson,3)) + "\\n") 
file1.close()

shannon = skbio.diversity.alpha.shannon(tax_table_species['read_count'])
simpson = skbio.diversity.alpha.simpson(tax_table_species['read_count'])
file1 = open("blast_diversity_species.csv","w") 
file1.write(str(tax_table_species['read_count'].sum()) + "," + str(round(shannon,3)) + "," + str(round(simpson,3)) + "\\n") 
file1.close()
