
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 17:15:37 2016

@author: lpinello
"""

from flask import *
import pandas as pd
import requests 
import json
import os
import re


import flask_excel as excel
excel.init_excel(app)

#import sys
#reload(sys)
#sys.setdefaultencoding('utf8')

dir_path = os.path.abspath(os.curdir)
dir_path_ex = os.path.dirname(os.path.realpath(__file__))


#global settings

N_MAX_GENES=100
#GENES_COMMON_PATH = 2
pd.set_option('max_colwidth',2000)
omim_api_key='jrnmZyB4QS6R5ANxZGswmA'

_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
		return os.path.join(_ROOT, 'data', path)


#load local data
df_omim_gene=pd.read_table(get_data('mim2gene.txt'),skiprows=4)
df_omim_gene=df_omim_gene.dropna(subset=['Approved Gene Symbol (HGNC)'])
df_exac=pd.read_table(get_data('forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz'),compression='gzip').set_index('gene')

df_broad=pd.read_table(get_data('Broad_CTD2_compound.txt'),sep='\t')
df_broad_names = df_broad['Broad druggable gene'].str.strip()
df_broad_drugs = df_broad['Broad CTD2 440 compound plates']
df_broad_targets = df_broad['TARGET_NODES']

df_rxgene=pd.read_table(get_data('Screened_Compounds_cancerrxgene.txt'),sep='\t')
df_rxgene_targets = df_rxgene['Drug Name']
df_rxgene_drugs = df_rxgene['Target']

df_hgnc=pd.read_table(get_data('hgnc_id_tab.txt'),sep='\t')
df_hgnc_id=df_hgnc.hgnc_id
df_hgnc_symbol=df_hgnc.symbol
df_hgnc_name=df_hgnc.name
del df_hgnc
df_kegg=pd.read_table(get_data('CPDB_pathways_genes.tab'),sep='\t')
df_kegg_pathway = df_kegg.pathway
df_kegg_gene_list = df_kegg.hgnc_symbol_ids
del df_kegg
#broad_druggable_genes=set([g.strip() for g in pd.read_table(get_data('Druggable_gene_list_broad.txt')).set_index('Broad druggable gene').index])

#supporting functions
def omim_has_variant(omim_id,apiKey=omim_api_key):
	rest_url_omim_variants='http://api.omim.org/api/entry?mimNumber=%d&include=allelicVariantList&format=json&apiKey=%s' % (omim_id,apiKey)
	response = requests.get(rest_url_omim_variants)
	if response.ok:
		json_data = json.loads(response.content)

	has_variant_table=False
	for entry in json_data['omim']['entryList']:
		if entry['entry'].has_key('allelicVariantList'):
			has_variant_table=True
			break

	return has_variant_table

def get_links_from_gene(gene_name):
	df_hits=df_omim_gene.ix[df_omim_gene['Approved Gene Symbol (HGNC)']==gene_name]
	exac_links=[]
	omim_main_links=[]
	omim_variants_links=[]
	clinvar_links=[]
	gnomad_links=[]
	pharos_links = []

	OMIM_ID=''
	ENS_ID=''

	for idx,row in df_hits.iterrows():
		OMIM_ID=row['# MIM Number']
		ENS_ID=row['Ensembl Gene ID (Ensembl)'].split(',')[0]
		omim_main_links.append('http://www.omim.org/entry/{0} '.format(OMIM_ID))

		if omim_has_variant(OMIM_ID):
			omim_variants_links.append('http://omim.org/allelicVariant/{0} '.format(OMIM_ID))
		else:
			omim_variants_links.append('')

		exac_links.append('http://exac.broadinstitute.org/gene/{0}'.format(ENS_ID))
		clinvar_links.append('http://www.ncbi.nlm.nih.gov/clinvar?term=%d%%5BMIM%%5D' % OMIM_ID)
		gnomad_links.append('http://gnomad.broadinstitute.org/gene/{0}'.format(ENS_ID))
		pharos_links.append('https://pharos.nih.gov/idg/search?q={0}'.format(ENS_ID))

	return OMIM_ID,ENS_ID,','.join(omim_main_links),','.join(omim_variants_links),','.join(exac_links),','.join(clinvar_links),','.join(gnomad_links),','.join(pharos_links)

def get_exac_lof_hom_count(ensembl_gene_id):
	rest_url_exac='http://exac.hms.harvard.edu/rest/gene/variants_in_gene/%s' % ensembl_gene_id
	response = requests.get(rest_url_exac)


	if response.ok:

		try:
			json_data = json.loads(response.content)
		except:
			return 'N/A','N/A'

	hom_count=0
	lof_count=0
	for hit in json_data:

		if hit['filter']=='PASS':
			if hit['hom_count']>0:
				hom_count+=1

			if hit['category']=='lof_variant':
				lof_count+=1


	return hom_count,lof_count

def get_gene_drug_interactions_dgifb(gene_names):

	gene_names_string=','.join(gene_names)
	rest_url_dgifb='http://dgidb.genome.wustl.edu/api/v1/interactions.json?genes=%s' % gene_names_string
	response = requests.get(rest_url_dgifb)
	if response.ok:
		 json_data = json.loads(response.content)

	genes_drug_interactions={}

	for matchedTerm in json_data['matchedTerms']:
		genes_drug_interactions[matchedTerm['geneName']]=sorted(set([interaction['drugName'] for interaction in matchedTerm['interactions']]))

	return genes_drug_interactions

def get_filled_dataframe(list_of_genes):

	if len(list_of_genes)>0:
		df_gene_list=pd.DataFrame(list_of_genes,columns=['Gene Symbol'])
		all_omim_main_links=[]
		all_omim_variants_links=[]
		all_exac_links=[]
		all_clinvar_links=[]
		all_broad_druggable=[]
		all_broad_druggable_names=[]
		all_hgnc_id=[]
		all_kegg_path=[]
		#all_hgnc_id_temp = pd.Series()
		hgnc_id_iterate=pd.Series()
		path_iterate = pd.Series()
		path_iterate_full = pd.Series()
		path_temp_full = pd.Series()
		all_omim_ids=[]
		all_ens_ids=[]
		all_gnomad_links=[]
		all_pharos_links = []
		gene_pathways_counts=pd.Series()
		dgi_counts = pd.Series()
		all_exac_n_lof=[]
		all_exac_pLI=[]
		all_exac_mis_z=[]
		all_dgifb_count=[]
		all_dgifb_interactions=[]
		all_dgifb_html_links=[]
		all_int_maps=[]
		int_maps_api_str_full = pd.Series()
		output_counts = pd.Series()
		all_genes_path_full = pd.Series()
		all_genes_path_count = pd.Series()
		genes_drug_interactions=get_gene_drug_interactions_dgifb(df_gene_list['Gene Symbol'])

		counter = 0
		path_analysis_dup = pd.Series()
		all_broad_idx_drugs = pd.Series()
		broad_idx_drugs = pd.Series()


		for gene_name in df_gene_list['Gene Symbol']:
			counter = counter+1
			path_iterate = pd.Series()
			path_temp_full = pd.Series()

			OMIM_ID,ENS_ID,omim_main_links, omim_variants_links,exac_links,clinvar_links,gnomad_links,pharos_links=get_links_from_gene(gene_name)
			all_omim_main_links.append(omim_main_links if omim_main_links else 'N/A' )
			all_omim_variants_links.append(omim_variants_links if omim_variants_links else 'N/A')

			int_maps_api_str1 = 'http://string-db.org/api/image/network?identifier='
			int_maps_api_str2 ='&species=9606'
			int_maps_api_str = pd.Series(''.join([int_maps_api_str1+gene_name+int_maps_api_str2]))
			int_maps_api_str_full = int_maps_api_str_full.append(int_maps_api_str.reset_index(drop=True), ignore_index=True)

			linker1 = ","
			linker2 = ", "
			temp_var1 = ''.join([linker2+gene_name+linker1])
			temp_var2 = ''.join([gene_name + linker1])
			temp_var3 = ''.join([linker2 + gene_name])
			gene_name_path_idx1 = pd.Series(df_kegg_gene_list.str.contains(temp_var1), name='bools')
			del temp_var1
			gene_name_path_idx2 = pd.Series(df_kegg_gene_list.str.contains(temp_var2), name='bools')
			del temp_var2
			gene_name_path_idx3 = pd.Series(df_kegg_gene_list.str.contains(temp_var3), name='bools')
			del temp_var3
			path_hits1 = gene_name_path_idx1[gene_name_path_idx1.values]
			path_hits2 = gene_name_path_idx2[gene_name_path_idx2.values]
			path_hits3 = gene_name_path_idx3[gene_name_path_idx3.values]
			path_temp1 = df_kegg_pathway.iloc[path_hits1.index]
			path_temp2 = df_kegg_pathway.iloc[path_hits2.index]
			path_temp3 = df_kegg_pathway.iloc[path_hits3.index]
			path_temp_full = path_temp_full.append(path_temp1)
			path_temp_full = path_temp_full.append(path_temp2)
			path_temp_full = path_temp_full.append(path_temp3)


			path_temp = pd.Series(path_temp_full.unique())

			path_analysis_dup = path_analysis_dup.append(path_temp.reset_index(drop=True), ignore_index=True)
			path_analysis_temp = pd.Series()
			output_counts = output_counts.append(pd.Series(len(path_temp)).reset_index(drop=True), ignore_index=True)
			path_iterate = pd.Series([path_temp.str.cat(sep=', ')])
			path_iterate_full = path_iterate_full.append(path_iterate.reset_index(drop=True), ignore_index=True)
			all_kegg_path = path_iterate_full
			del path_temp

			gene_name_idx = df_hgnc_symbol.str.find(gene_name, start=0, end=None)
			gene_name_len = len(gene_name)
			df_hits = gene_name_idx.index[gene_name_idx == 0]

			if len(df_hits) == 1:
				hgnc_id_temp = df_hgnc_id.iloc[df_hits]
				hgnc_id_iterate = hgnc_id_iterate.append(hgnc_id_temp.reset_index(drop=True),ignore_index=True)
			else:
				len_of_gene_name = []
				temps = df_hgnc_symbol.iloc[df_hits]
				for idx, temp in enumerate(temps):
					len_of_gene_name.append(len(temp))

				len_of_gene_name.index(gene_name_len)
				df_hits_sub1 = df_hits[len_of_gene_name.index(gene_name_len)]
				df_hits_sub = pd.Series(df_hits_sub1)
				hgnc_id_temp = df_hgnc_id.iloc[df_hits_sub]
				hgnc_id_iterate = hgnc_id_iterate.append(hgnc_id_temp.reset_index(drop=True),ignore_index=True)

			all_hgnc_id = hgnc_id_iterate

			if gene_name in df_exac.index:
				all_exac_n_lof.append(df_exac.ix[gene_name,'n_lof'])
				all_exac_pLI.append('%.2f' % df_exac.ix[gene_name,'pLI'])
				all_exac_mis_z.append('%.2f' % df_exac.ix[gene_name,'mis_z'])
				if exac_links:
					all_exac_links.append(exac_links)
				else:
					all_exac_links.append('http://exac.broadinstitute.org/gene/%s' % df_exac.ix[gene_name,'transcript'].split('.')[0])
			else:
				all_exac_n_lof.append('N/A')
				all_exac_pLI.append('N/A')
				all_exac_mis_z.append('N/A')
				all_exac_links.append('N/A')

			all_clinvar_links.append(clinvar_links if clinvar_links else 'N/A')
			all_gnomad_links.append(gnomad_links if gnomad_links else 'N/A')
			all_pharos_links.append(pharos_links if pharos_links else 'N/A')

			all_omim_ids.append(OMIM_ID)
			all_ens_ids.append(ENS_ID)

			if 'broad_idx_drugs' in locals():
				del broad_idx_drugs
			if 'full_temps_broad_targets' in locals():
				del full_temps_broad_targets
			if 'full_temps_broad_drugs' in locals():
				del full_temps_broad_drugs
			full_temps_broad_targets = pd.Series()
			full_temps_broad_drugs = pd.Series()
			broad_idx_drugs = pd.Series()

			if len(df_broad_names.index[df_broad_names.str.find(gene_name)==0]) > 0:
				simple_method = 0
				skip_step = 0
				broad_idx1 = df_broad_targets.str.contains(gene_name)
				broad_idx = df_broad_targets.values[broad_idx1 == True]
				broad_gene_name_len = len(gene_name)

				if len(broad_idx) == 1:
					if str(broad_idx).find(";") > 0:
						broad_genes_temp = pd.Series(broad_idx).str.split(';',expand=True)
						broad_genes_temp = broad_genes_temp.transpose()
						broad_genes_temp = broad_genes_temp[0]

						len_of_broad_gene_name = []
						broad_genes_temp = broad_genes_temp.str.strip()
						broad_genes_temp_len = len(broad_genes_temp)
						full_temps_broad_targets = full_temps_broad_targets.append(broad_genes_temp.reset_index(drop=True),ignore_index=True)

						for i in range(len(full_temps_broad_targets)):
							full_temps_broad_drugs = full_temps_broad_drugs.append(pd.Series(df_broad_drugs.values[broad_idx1 == True]).reset_index(drop=True),ignore_index=True)
						skip_step = 1
					else:
						broad_genes_temp_len = pd.Series(broad_idx).str.len()

						if broad_gene_name_len == broad_genes_temp_len[0]:
							broad_idx_drugs = broad_idx_drugs.append(pd.Series(df_broad_drugs.values[broad_idx1 == True]).reset_index(drop=True),ignore_index=True)
							simple_method = 1
						else:
							broad_idx_drugs = broad_idx_drugs.append(pd.Series('N/A').reset_index(drop=True), ignore_index=True)

				else:
					len_of_broad_gene_name = []
					broad_idx_temp = pd.Series()
					rel_temps_broad_targets = pd.Series(df_broad_targets.values[broad_idx1 == True])
					rel_temps_broad_drugs = pd.Series(df_broad_drugs.values[broad_idx1 == True])

					if len(df_broad_targets.str.contains(";")==True) > 0:
						broad_target_temp_idx1 = rel_temps_broad_targets.str.contains(";")
						broad_target_temp_idx_withs = pd.Series(broad_target_temp_idx1.index[broad_target_temp_idx1 == True])
						broad_target_temp_idx_withouts = pd.Series(broad_target_temp_idx1.index[broad_target_temp_idx1 == False])

						if len(broad_target_temp_idx_withs) >= 1:
							for broad_target_temp_idx_with in broad_target_temp_idx_withs:
								full_temps_broad_targets1 = pd.Series(rel_temps_broad_targets.iloc[broad_target_temp_idx_with]).str.split(';',expand=True)
								temps_broad_targets = full_temps_broad_targets1.transpose()
								full_temps_broad_targets = full_temps_broad_targets.append(temps_broad_targets.reset_index(drop=True))

								for i in range(len(temps_broad_targets)):
									full_temps_broad_drugs = full_temps_broad_drugs.append(pd.Series(rel_temps_broad_drugs.iloc[broad_target_temp_idx_with]).reset_index(drop=True),ignore_index=True)

						if len(broad_target_temp_idx_withouts) >= 1:
							for broad_target_temp_idx_without in broad_target_temp_idx_withouts:
								full_temps_broad_targets = full_temps_broad_targets.append(pd.Series(rel_temps_broad_targets.iloc[broad_target_temp_idx_without]).reset_index(drop=True),ignore_index=True)
								full_temps_broad_drugs = full_temps_broad_drugs.append(pd.Series(rel_temps_broad_drugs.iloc[broad_target_temp_idx_without]).reset_index(drop=True),ignore_index=True)

					else:
						full_temps_broad_targets = rel_temps_broad_targets
						full_temps_broad_drugs = rel_temps_broad_drugs

				if simple_method == 0:
					if skip_step == 0:
						full_temps_broad_targets = pd.Series(full_temps_broad_targets[0])
					else:
						full_temps_broad_targets = pd.Series(full_temps_broad_targets)

					full_temps_broad_drugs = pd.Series(full_temps_broad_drugs)

					broad_idx_new = full_temps_broad_targets.str.contains(gene_name)
					broad_gene_name_len_new = len(gene_name)
					full_temps_broad_targets_temps = full_temps_broad_targets[broad_idx_new == True]
					full_temps_broad_drugs_temps = full_temps_broad_drugs[broad_idx_new == True]

					for idx, full_temps_broad_targets_temp in enumerate(full_temps_broad_targets_temps):
						len_of_broad_gene_name.append(len(full_temps_broad_targets_temp))

					len_of_broad_gene_name = pd.Series(len_of_broad_gene_name)
					broad_idx_new_final = full_temps_broad_targets_temps.index[len_of_broad_gene_name == broad_gene_name_len_new]
					broad_hits_rel1 = full_temps_broad_drugs_temps[broad_idx_new_final]
					df_broad_rel = pd.Series(broad_hits_rel1)

					if df_broad_rel.empty:
						broad_idx_drugs = broad_idx_drugs.append(pd.Series('N/A').reset_index(drop=True),ignore_index=True)
					else:
						broad_idx_drugs = broad_idx_drugs.append(df_broad_rel.reset_index(drop=True), ignore_index=True)
			else:
				broad_idx_drugs = broad_idx_drugs.append(pd.Series('N/A').reset_index(drop=True), ignore_index=True)
				broad_idx_drugs = pd.Series(broad_idx_drugs).unique()

			if len(broad_idx_drugs) > 1:
				broad_idx_drugs_temp = pd.Series([broad_idx_drugs.str.cat(sep=', ')])
				all_broad_idx_drugs = all_broad_idx_drugs.append(broad_idx_drugs_temp.reset_index(drop=True), ignore_index=True)
			else:
				all_broad_idx_drugs = all_broad_idx_drugs.append(pd.Series(broad_idx_drugs).reset_index(drop=True), ignore_index=True)

			if genes_drug_interactions.has_key(gene_name):
				all_dgifb_count.append(len(genes_drug_interactions[gene_name]))
				all_dgifb_interactions.append(', '.join(genes_drug_interactions[gene_name]))
				all_dgifb_html_links.append('http://dgidb.genome.wustl.edu/interaction_search_results?search_mode=genes&identifiers=%s' % gene_name)
			else:
				all_dgifb_count.append(0)
				all_dgifb_interactions.append('None')
				all_dgifb_html_links.append('N/A')

		df_gene_list['HGNC ID'] = all_hgnc_id
		df_gene_list['KEGG Pathway'] = all_kegg_path
		df_gene_list['DGIdb #Interactions']=all_dgifb_count
		df_gene_list['DGIdb Interactions']=all_dgifb_interactions
		df_gene_list['DGIdb']=all_dgifb_html_links
		#df_gene_list['BROAD #interactions']=all_broad_druggable
		#df_gene_list['BROAD interactions']=all_broad_druggable_names
		df_gene_list['CTD2 Interactions'] = all_broad_idx_drugs
		df_gene_list['OMIM']=all_omim_main_links
		df_gene_list['OMIM variants']=all_omim_variants_links
		df_gene_list['ClinVar']=all_clinvar_links
		df_gene_list['gnomAD'] = all_gnomad_links
		df_gene_list['Pharos'] = all_pharos_links
		df_gene_list['ExAC']=all_exac_links
		df_gene_list['ExAC #LoF']=all_exac_n_lof
		df_gene_list['ExAC Missense z']=all_exac_mis_z
		df_gene_list['ExAC pLI']=all_exac_pLI
		df_gene_list['Interaction Map'] = int_maps_api_str_full


		common_paths_full = path_analysis_dup[path_analysis_dup.duplicated(keep=False)]
		common_paths1 = pd.Series(common_paths_full.unique())
		common_paths = pd.Series()
		common_paths['Common Pathways'] = common_paths1
		common_paths['Counts'] = pd.Series(common_paths_full.value_counts())

		list_common_pathways2 = ()
		gene_pathways = pd.Series()
		GENES_COMMON_PATH = int(request.form['num_common_pathway'])
		no_pathways_found = 0

		if len(common_paths['Counts'].loc[(common_paths['Counts'] >= GENES_COMMON_PATH)]) > 0:

			gene_pathway_temp = pd.Series(common_paths['Counts'].index[(common_paths['Counts'] >= GENES_COMMON_PATH)])
			gene_pathways = gene_pathways.append(gene_pathway_temp)

			df_gene_list_paths = pd.DataFrame(gene_pathways, columns=['KEGG Pathway'])
			df_gene_list_temp_path = pd.DataFrame()

			gene_path_idx = pd.Series()
			hgnc_subset_iterate_full = pd.Series()
			symbol_id_subset_iterate_full = pd.Series()
			dgidb_num_subset_iterate_full = pd.Series()
			dgidb_int_subset_iterate_full = pd.Series()
			dgidb_subset_iterate_full = pd.Series()
			broad_num_subset_iterate_full = pd.Series()
			broad_int_subset_iterate_full = pd.Series()
			omim_subset_iterate_full = pd.Series()
			omim_var_subset_iterate_full = pd.Series()
			clinvar_subset_iterate_full = pd.Series()
			gnomad_subset_iterate_full = pd.Series()
			pharos_subset_iterate_full = pd.Series()
			exac_subset_iterate_full = pd.Series()
			exac_lof_subset_iterate_full = pd.Series()
			exac_mis_subset_iterate_full = pd.Series()
			exac_pli_subset_iterate_full = pd.Series()

			for gene_pathway in gene_pathways:
				all_genes_path1 = df_kegg_pathway.str.find(gene_pathway, start=0, end=None)

				path_temp_full1 = pd.Series()
				path_analysis1_dup = pd.Series()

				linker1 = ","
				linker2 = ", "
				temp_var11 = ''.join([linker2 + gene_pathway + linker1])
				temp_var22 = ''.join([gene_pathway + linker1])
				temp_var33 = ''.join([linker2 + gene_pathway])
				gene_name_path_idx11 = pd.Series(df_gene_list['KEGG Pathway'].str.find(temp_var11))
				del temp_var11
				gene_name_path_idx22 = pd.Series(df_gene_list['KEGG Pathway'].str.find(temp_var22))
				del temp_var22
				gene_name_path_idx33 = pd.Series(df_gene_list['KEGG Pathway'].str.find(temp_var33))
				del temp_var33

				df_gene_list['KEGG Pathway'].to_csv('test.csv')

				path_hits11 = gene_name_path_idx11.index[gene_name_path_idx11 >= 0]
				path_hits22 = gene_name_path_idx22.index[gene_name_path_idx22 >= 0]
				path_hits33 = gene_name_path_idx33.index[gene_name_path_idx33 >= 0]

				path_temp11 = df_gene_list['Gene Symbol'].iloc[path_hits11]
				path_temp22 = df_gene_list['Gene Symbol'].iloc[path_hits22]
				path_temp33 = df_gene_list['Gene Symbol'].iloc[path_hits33]
				path_temp_full1 = path_temp_full1.append(path_temp11)
				path_temp_full1 = path_temp_full1.append(path_temp22)
				path_temp_full1 = path_temp_full1.append(path_temp33)

				path_temp1 = pd.Series(path_temp_full1.unique())
				del path_temp_full1
				path_analysis1_dup = path_analysis1_dup.append(path_temp1.reset_index(drop=True), ignore_index=True)
				unique_path = path_temp1
				del path_temp1

				all_genes_path_temp = all_genes_path1.index[all_genes_path1 == 0]
				num_path_found = len(all_genes_path_temp)
				pathway_len = len(gene_pathway)

				if num_path_found > 1:
					len_of_input_name_paths = []
					temps_paths = df_kegg_pathway.iloc[all_genes_path_temp]
					for idx, temps_path in enumerate(temps_paths):
						len_of_input_name_paths.append(len(temps_path))

					if pathway_len in len_of_input_name_paths:
						len_of_input_name_paths.index(pathway_len)
						path_sub1 = all_genes_path_temp[len_of_input_name_paths.index(pathway_len)]
						path_sub = pd.Series(path_sub1)
						del all_genes_path_temp
						all_genes_path_temp = path_sub

				all_genes_path = df_kegg_gene_list.iloc[all_genes_path_temp]
				all_genes_path_count_temp = pd.Series(all_genes_path).str.count(',')
				all_genes_path_count_temp = all_genes_path_count_temp + 1
				all_genes_path_count = all_genes_path_count.append(all_genes_path_count_temp.reset_index(drop=True),ignore_index=True)
				all_genes_path_full = all_genes_path_full.append(pd.Series(all_genes_path).reset_index(drop=True),ignore_index=True)
				del all_genes_path_temp
				del all_genes_path

			#for gene_pathway in gene_pathways:
				if 'relevant_genes' in locals():
					del relevant_genes
				relevant_genes = []

				symbol_id_subset = pd.Series()
				hgnc_id_subset = pd.Series()
				dgidb_num_subset = pd.Series()
				dgidb_int_subset = pd.Series()
				dgidb_subset = pd.Series()
				broad_num_subset = pd.Series()
				broad_int_subset = pd.Series()
				omim_subset = pd.Series()
				omim_var_subset = pd.Series()
				clinvar_subset = pd.Series()
				gnomad_subset = pd.Series()
				pharos_subset = pd.Series()
				exac_subset = pd.Series()
				exac_lof_subset = pd.Series()
				exac_mis_subset = pd.Series()
				exac_pli_subset = pd.Series()
				#relevant_genes = df_gene_list_temp_path

				relevant_genes = unique_path
				del unique_path

				rel_genes_temp1 = pd.isnull(relevant_genes)
				rel_genes_temp = relevant_genes.values[rel_genes_temp1 == False]
				rel_genes_temp = pd.Series(rel_genes_temp)
				del relevant_genes
				relevant_genes = rel_genes_temp
				loop_max = len(relevant_genes)
				counter = 0

				for relevant_gene in relevant_genes:
					counter = counter + 1
					relevant_gene = str(relevant_gene)
					gene_path_idx1 = df_gene_list['Gene Symbol'].str.find(relevant_gene, start=0, end=None)
					rel_gene_name_len = len(relevant_gene)
					gene_path_idx = gene_path_idx1.index[gene_path_idx1 == 0]

					if len(gene_path_idx) == 1:
						if df_gene_list['DGIdb #Interactions'].iloc[gene_path_idx].values == 0:
							dgidb_int_subset = dgidb_int_subset.append(pd.Series('N/A').reset_index(drop=True),ignore_index=True)
						else:
							dgidb_int_subset = dgidb_int_subset.append(df_gene_list['DGIdb Interactions'].iloc[gene_path_idx].reset_index(drop=True),ignore_index=True)
					else:
						len_of_rel_gene_name = []
						rel_temps = df_gene_list['Gene Symbol'].iloc[gene_path_idx]
						for idx, rel_temp in enumerate(rel_temps):
							len_of_rel_gene_name.append(len(rel_temp))

						len_of_rel_gene_name.index(rel_gene_name_len)
						df_hits_rel1 = gene_path_idx[len_of_rel_gene_name.index(rel_gene_name_len)]
						df_hits_rel = pd.Series(df_hits_rel1)

						if df_gene_list['DGIdb #Interactions'].iloc[df_hits_rel].values == 0:
							dgidb_int_subset = dgidb_int_subset.append(pd.Series('N/A').reset_index(drop=True),ignore_index=True)
						else:
							dgidb_int_subset = dgidb_int_subset.append(df_gene_list['DGIdb Interactions'].iloc[df_hits_rel].reset_index(drop=True),ignore_index=True)

						#hgnc_id_temp = df_hgnc_id.iloc[df_hits_rel]
						#hgnc_id_iterate = hgnc_id_iterate.append(hgnc_id_temp.reset_index(drop=True), ignore_index=True)

					symbol_id_subset = symbol_id_subset.append(df_gene_list['Gene Symbol'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					hgnc_id_subset = hgnc_id_subset.append(df_gene_list['HGNC ID'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					dgidb_num_subset = dgidb_num_subset.append(df_gene_list['DGIdb #Interactions'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					dgidb_subset = dgidb_subset.append(df_gene_list['DGIdb'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					#broad_num_subset = broad_num_subset.append(df_gene_list['BROAD #interactions'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					broad_int_subset = broad_int_subset.append(df_gene_list['CTD2 Interactions'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					omim_subset = omim_subset.append(df_gene_list['OMIM'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					omim_var_subset = omim_var_subset.append(df_gene_list['OMIM variants'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					clinvar_subset = clinvar_subset.append(df_gene_list['ClinVar'].iloc[gene_path_idx].reset_index(drop=True),ignore_index=True)
					gnomad_subset = gnomad_subset.append(df_gene_list['gnomAD'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					pharos_subset = pharos_subset.append(df_gene_list['Pharos'].iloc[gene_path_idx].reset_index(drop=True), ignore_index=True)
					exac_subset = exac_subset.append(df_gene_list['ExAC'].iloc[gene_path_idx].reset_index(drop=True),ignore_index=True)
					exac_lof_subset = exac_lof_subset.append(df_gene_list['ExAC #LoF'].iloc[gene_path_idx].reset_index(drop=True),ignore_index=True)
					exac_mis_subset = exac_mis_subset.append(df_gene_list['ExAC Missense z'].iloc[gene_path_idx].reset_index(drop=True),ignore_index=True)
					exac_pli_subset = exac_pli_subset.append(df_gene_list['ExAC pLI'].iloc[gene_path_idx].reset_index(drop=True),ignore_index=True)

					if counter == loop_max:
						#gene_pathways_counts = gene_pathways_counts.append(pd.Series(loop_max).reset_index(drop=True),ignore_index=True)
						gene_pathways_counts = gene_pathways_counts.append(pd.Series(len(hgnc_id_subset)).reset_index(drop=True),ignore_index=True)
						dgi_counts = dgi_counts.append(pd.Series(dgidb_num_subset.sum()).reset_index(drop=True), ignore_index=True)

						dgidb_num_subset = dgidb_num_subset.apply(str)
						broad_num_subset = broad_num_subset.apply(str)
						exac_lof_subset = exac_lof_subset.apply(str)
						exac_mis_subset = exac_mis_subset.apply(str)
						exac_pli_subset = exac_pli_subset.apply(str)
						symbol_id_subset_iterate = pd.Series([symbol_id_subset.str.cat(sep=', ')])
						hgnc_subset_iterate = pd.Series([hgnc_id_subset.str.cat(sep=', ')])
						dgidb_num_subset_iterate = pd.Series([dgidb_num_subset.str.cat(sep=', ')])
						dgidb_int_subset_iterate = pd.Series([dgidb_int_subset.str.cat(sep=', ')])
						dgidb_subset_iterate = pd.Series([dgidb_subset.str.cat(sep=', ')])
						#broad_num_subset_iterate = pd.Series([broad_num_subset.str.cat(sep=', ')])
						broad_int_subset_iterate = pd.Series([broad_int_subset.str.cat(sep=', ')])
						omim_subset_iterate = pd.Series([omim_subset.str.cat(sep=', ')])
						omim_var_subset_iterate = pd.Series([omim_var_subset.str.cat(sep=', ')])
						clinvar_subset_iterate = pd.Series([clinvar_subset.str.cat(sep=', ')])
						gnomad_subset_iterate = pd.Series([gnomad_subset.str.cat(sep=', ')])
						pharos_subset_iterate = pd.Series([pharos_subset.str.cat(sep=', ')])
						exac_subset_iterate = pd.Series([exac_subset.str.cat(sep=', ')])
						exac_lof_subset_iterate = pd.Series([exac_lof_subset.str.cat(sep=', ')])
						exac_mis_subset_iterate = pd.Series([exac_mis_subset.str.cat(sep=', ')])
						exac_pli_subset_iterate = pd.Series([exac_pli_subset.str.cat(sep=', ')])

						symbol_id_subset_iterate_full = symbol_id_subset_iterate_full.append(symbol_id_subset_iterate.reset_index(drop=True), ignore_index=True)
						hgnc_subset_iterate_full = hgnc_subset_iterate_full.append(hgnc_subset_iterate.reset_index(drop=True), ignore_index=True)
						dgidb_num_subset_iterate_full = dgidb_num_subset_iterate_full.append(dgidb_num_subset_iterate.reset_index(drop=True), ignore_index=True)
						dgidb_int_subset_iterate_full = dgidb_int_subset_iterate_full.append(dgidb_int_subset_iterate.reset_index(drop=True), ignore_index=True)
						dgidb_subset_iterate_full = dgidb_subset_iterate_full.append(dgidb_subset_iterate.reset_index(drop=True), ignore_index=True)
						#broad_num_subset_iterate_full = broad_num_subset_iterate_full.append(broad_num_subset_iterate.reset_index(drop=True), ignore_index=True)
						broad_int_subset_iterate_full = broad_int_subset_iterate_full.append(broad_int_subset_iterate.reset_index(drop=True), ignore_index=True)
						omim_subset_iterate_full = omim_subset_iterate_full.append(omim_subset_iterate.reset_index(drop=True), ignore_index=True)
						omim_var_subset_iterate_full = omim_var_subset_iterate_full.append(omim_var_subset_iterate.reset_index(drop=True), ignore_index=True)
						clinvar_subset_iterate_full = clinvar_subset_iterate_full.append(clinvar_subset_iterate.reset_index(drop=True), ignore_index=True)
						gnomad_subset_iterate_full = gnomad_subset_iterate_full.append(gnomad_subset_iterate.reset_index(drop=True), ignore_index=True)
						pharos_subset_iterate_full = pharos_subset_iterate_full.append(pharos_subset_iterate.reset_index(drop=True), ignore_index=True)
						exac_subset_iterate_full = exac_subset_iterate_full.append(exac_subset_iterate.reset_index(drop=True), ignore_index=True)
						exac_lof_subset_iterate_full = exac_lof_subset_iterate_full.append(exac_lof_subset_iterate.reset_index(drop=True), ignore_index=True)
						exac_mis_subset_iterate_full = exac_mis_subset_iterate_full.append(exac_mis_subset_iterate.reset_index(drop=True), ignore_index=True)
						exac_pli_subset_iterate_full = exac_pli_subset_iterate_full.append(exac_pli_subset_iterate.reset_index(drop=True), ignore_index=True)

			df_gene_list_paths['Gene Symbol'] = symbol_id_subset_iterate_full
			df_gene_list_paths['HGNC ID'] = hgnc_subset_iterate_full
			df_gene_list_paths['# of Input Genes in Pathway'] = gene_pathways_counts
			df_gene_list_paths['All Genes in Pathway'] = all_genes_path_full
			df_gene_list_paths['DGIdb #Interactions'] = dgidb_num_subset_iterate_full
			df_gene_list_paths['DGIdb Interactions'] = dgidb_int_subset_iterate_full
			df_gene_list_paths['DGIdb'] = dgidb_subset_iterate_full
			#df_gene_list_paths['BROAD #interactions'] = broad_num_subset_iterate_full
			df_gene_list_paths['CTD2 Interactions'] = broad_int_subset_iterate_full
			df_gene_list_paths['OMIM'] = omim_subset_iterate_full
			df_gene_list_paths['OMIM variants'] = omim_var_subset_iterate_full
			df_gene_list_paths['ClinVar'] = clinvar_subset_iterate_full
			df_gene_list_paths['gnomAD'] = gnomad_subset_iterate_full
			df_gene_list_paths['Pharos'] = pharos_subset_iterate_full
			df_gene_list_paths['ExAC'] = exac_subset_iterate_full
			df_gene_list_paths['ExAC #LoF'] = exac_lof_subset_iterate_full
			df_gene_list_paths['ExAC Missense z'] = exac_mis_subset_iterate_full
			df_gene_list_paths['ExAC pLI'] = exac_pli_subset_iterate_full
			#df_gene_list_paths['Interaction Map'] = pd.Series("N/A")

		else:
			df_gene_list_paths = pd.DataFrame()
			no_pathways_found = 1
			df_gene_list_paths['Gene Symbol'] = pd.Series("No common pathways identified")
			df_gene_list_paths['HGNC ID'] = pd.Series("N/A")
			df_gene_list_paths['KEGG Pathway'] = pd.Series("N/A")
			df_gene_list_paths['# of Input Genes in Pathway'] = pd.Series("N/A")
			df_gene_list_paths['All Genes in Pathway'] = pd.Series("N/A")
			df_gene_list_paths['DGIdb #Interactions'] = pd.Series("N/A")
			df_gene_list_paths['DGIdb Interactions'] = pd.Series("N/A")
			df_gene_list_paths['DGIdb'] = pd.Series("N/A")
			#df_gene_list_paths['BROAD #interactions'] = pd.Series("N/A")
			df_gene_list_paths['CTD2 Interactions'] = pd.Series("N/A")
			df_gene_list_paths['OMIM'] = pd.Series("N/A")
			df_gene_list_paths['OMIM variants'] = pd.Series("N/A")
			df_gene_list_paths['ClinVar'] = pd.Series("N/A")
			df_gene_list_paths['gnomAD'] = pd.Series("N/A")
			df_gene_list_paths['Pharos'] = pd.Series("N/A")
			df_gene_list_paths['ExAC'] = pd.Series("N/A")
			df_gene_list_paths['ExAC #LoF'] = pd.Series("N/A")
			df_gene_list_paths['ExAC Missense z'] = pd.Series("N/A")
			df_gene_list_paths['ExAC pLI'] = pd.Series("N/A")
			#df_gene_list_paths['Interaction Map'] = pd.Series("N/A")

	else:
		df_gene_list = pd.DataFrame()
		no_pathways_found = 1
		df_gene_list['Gene Symbol'] = pd.Series("No genes identified")
		df_gene_list['HGNC ID'] = pd.Series("N/A")
		df_gene_list['KEGG Pathway'] = pd.Series("N/A")
		df_gene_list['DGIdb #Interactions'] = pd.Series("N/A")
		df_gene_list['DGIdb Interactions'] = pd.Series("N/A")
		df_gene_list['DGIdb'] = pd.Series("N/A")
		#df_gene_list['BROAD #interactions'] = pd.Series("N/A")
		df_gene_list['CTD2 Interactions'] = pd.Series("N/A")
		df_gene_list['OMIM'] = pd.Series("N/A")
		df_gene_list['OMIM variants'] = pd.Series("N/A")
		df_gene_list['ClinVar'] = pd.Series("N/A")
		df_gene_list['gnomAD'] = pd.Series("N/A")
		df_gene_list['Pharos'] = pd.Series("N/A")
		df_gene_list['ExAC'] = pd.Series("N/A")
		df_gene_list['ExAC #LoF'] = pd.Series("N/A")
		df_gene_list['ExAC Missense z'] = pd.Series("N/A")
		df_gene_list['ExAC pLI'] = pd.Series("N/A")
		df_gene_list['Interaction Map'] = pd.Series("N/A")
		#df_gene_list['All Genes in Pathway'] = pd.Series("N/A")


		df_gene_list_paths = pd.DataFrame()
		no_pathways_found = 1
		df_gene_list_paths['Gene Symbol'] = pd.Series("No common pathways identified")
		df_gene_list_paths['HGNC ID'] = pd.Series("N/A")
		df_gene_list_paths['KEGG Pathway'] = pd.Series("N/A")
		df_gene_list_paths['# of Input Genes in Pathway'] = pd.Series("N/A")
		df_gene_list_paths['All Genes in Pathway'] = pd.Series("N/A")
		df_gene_list_paths['DGIdb #Interactions'] = pd.Series("N/A")
		df_gene_list_paths['DGIdb Interactions'] = pd.Series("N/A")
		df_gene_list_paths['DGIdb'] = pd.Series("N/A")
		#df_gene_list_paths['BROAD #interactions'] = pd.Series("N/A")
		df_gene_list_paths['CTD2 Interactions'] = pd.Series("N/A")
		df_gene_list_paths['OMIM'] = pd.Series("N/A")
		df_gene_list_paths['OMIM variants'] = pd.Series("N/A")
		df_gene_list_paths['ClinVar'] = pd.Series("N/A")
		df_gene_list_paths['gnomAD'] = pd.Series("N/A")
		df_gene_list_paths['Pharos'] = pd.Series("N/A")
		df_gene_list_paths['ExAC'] = pd.Series("N/A")
		df_gene_list_paths['ExAC #LoF'] = pd.Series("N/A")
		df_gene_list_paths['ExAC Missense z'] = pd.Series("N/A")
		df_gene_list_paths['ExAC pLI'] = pd.Series("N/A")

	#df_gene_list_paths_final = pd.DataFrame()
	df_gene_list_paths_final = pd.DataFrame(df_gene_list_paths['KEGG Pathway'], columns=['KEGG Pathway'])
	df_gene_list_paths_final['# of Input Genes in Pathway'] = pd.Series(df_gene_list_paths['# of Input Genes in Pathway'].values[pd.isnull(df_gene_list_paths['# of Input Genes in Pathway']) == False])
	df_gene_list_paths_final['Gene Symbol'] = pd.Series(df_gene_list_paths['Gene Symbol'].values[pd.isnull(df_gene_list_paths['Gene Symbol']) == False])
	df_gene_list_paths_final['HGNC ID'] = pd.Series(df_gene_list_paths['HGNC ID'].values[pd.isnull(df_gene_list_paths['HGNC ID']) == False])
	df_gene_list_paths_final['KEGG Pathway'] = pd.Series(df_gene_list_paths['KEGG Pathway'].values[pd.isnull(df_gene_list_paths['KEGG Pathway']) == False])
	df_gene_list_paths_final['All Genes in Pathway'] = pd.Series(df_gene_list_paths['All Genes in Pathway'].values[pd.isnull(df_gene_list_paths['All Genes in Pathway']) == False])
	df_gene_list_paths_final['DGIdb #Interactions'] = pd.Series(df_gene_list_paths['DGIdb #Interactions'].values[pd.isnull(df_gene_list_paths['DGIdb #Interactions']) == False])
	df_gene_list_paths_final['DGIdb Interactions'] = pd.Series(df_gene_list_paths['DGIdb Interactions'].values[pd.isnull(df_gene_list_paths['DGIdb Interactions']) == False])
	#df_gene_list_paths_final['DGIdb'] = pd.Series(df_gene_list_paths['DGIdb'].values[pd.isnull(df_gene_list_paths['DGIdb']) == False])
	#df_gene_list_paths_final['BROAD #interactions'] = pd.Series(df_gene_list_paths['BROAD #interactions'].values[pd.isnull(df_gene_list_paths['BROAD #interactions']) == False])
	df_gene_list_paths_final['CTD2 Interactions'] = pd.Series(df_gene_list_paths['CTD2 Interactions'].values[pd.isnull(df_gene_list_paths['CTD2 Interactions']) == False])
	#df_gene_list_paths_final['OMIM'] = pd.Series(df_gene_list_paths['OMIM'].values[pd.isnull(df_gene_list_paths['OMIM']) == False])
	#df_gene_list_paths_final['OMIM variants'] = pd.Series(df_gene_list_paths['OMIM variants'].values[pd.isnull(df_gene_list_paths['OMIM variants']) == False])
	#df_gene_list_paths_final['ClinVar'] = pd.Series(df_gene_list_paths['ClinVar'].values[pd.isnull(df_gene_list_paths['ClinVar']) == False])
	#df_gene_list_paths_final['gnomAD'] = pd.Series(df_gene_list_paths['gnomAD'].values[pd.isnull(df_gene_list_paths['gnomAD']) == False])
	#df_gene_list_paths_final['Pharos'] = pd.Series(df_gene_list_paths['Pharos'].values[pd.isnull(df_gene_list_paths['Pharos']) == False])
	#df_gene_list_paths_final['ExAC'] = pd.Series(df_gene_list_paths['ExAC'].values[pd.isnull(df_gene_list_paths['ExAC']) == False])
	#df_gene_list_paths_final['ExAC #LoF'] = pd.Series(df_gene_list_paths['ExAC #LoF'].values[pd.isnull(df_gene_list_paths['ExAC #LoF']) == False])
	#df_gene_list_paths_final['ExAC Missense z'] = pd.Series(df_gene_list_paths['ExAC Missense z'].values[pd.isnull(df_gene_list_paths['ExAC Missense z']) == False])
	#df_gene_list_paths_final['ExAC pLI'] = pd.Series(df_gene_list_paths['ExAC pLI'].values[pd.isnull(df_gene_list_paths['ExAC pLI']) == False])
#df_gene_list_paths_final['Interaction Map'] = pd.Series(df_gene_list_paths['Interaction Map'].values[pd.isnull(df_gene_list_paths['Interaction Map']) == False])
	#del df_gene_list_paths_final['index']

	if no_pathways_found == 0:
		df_gene_list_paths_final['# of Input Genes in Pathway'] = df_gene_list_paths['# of Input Genes in Pathway']

	if no_pathways_found == 0:
		df_gene_list_paths_final_sorted = df_gene_list_paths_final.sort_values(by='# of Input Genes in Pathway', ascending=False)
		del df_gene_list_paths
		df_gene_list_paths_final_sorted = df_gene_list_paths_final_sorted.reset_index()
		df_gene_list_paths = df_gene_list_paths_final_sorted.copy()
		del df_gene_list_paths_final_sorted

	del df_gene_list_paths
	df_gene_list_paths = df_gene_list_paths_final.copy()
	df_gene_list_paths_output = df_gene_list_paths.copy()
	del df_gene_list_paths_final
	df_gene_list_hidden_table=df_gene_list.copy()
	df_gene_list_paths_hidden_table = df_gene_list_paths.copy()

	df_gene_list_output = df_gene_list.copy()

	urlify=lambda x: ', '.join(['<a href="{0}" target="_blank">link</a>'.format(l) for l in x.split(', ')]) if x!='N/A' else x
	urlify_w = lambda x: ', '.join(['<a href="#" id="gotohome" class="gotohome" name="{0}">{0}</a>'.format(l) for l in x.split(', ')]) if x!= 'N/A' else x
	urlify_id = lambda x: ', '.join(['<a id="{0}">{0}</a>'.format(l) for l in x.split(', ')]) if x != 'N/A' else x
	urlify_image = lambda x: ', '.join(['<img src="{0}" width="700">'.format(l) for l in x.split(', ')]) if x!='N/A' else x

	df_gene_list['Gene Symbol'] = df_gene_list['Gene Symbol'].apply(urlify_id)
	#df_gene_list['Counts'] = df_gene_list['Counts'].apply(urlify_ex)
	#df_gene_list['KEGG Pathway'] = df_gene_list['KEGG Pathway'].apply(urlify_ex)
	button_var = pd.Series()
	button_var1 = pd.Series()
	button_var2 = pd.Series()
	button_var3 = pd.Series()
	button_var4 = pd.Series()
	button_var5 = pd.Series()
	button_var6 = pd.Series()

	df_gene_list_paths['All Genes in Pathway'] = df_gene_list_paths['All Genes in Pathway'].str.replace(',',', ')

	ctd2_count_temp1 = pd.Series(df_gene_list_paths['CTD2 Interactions']).str.count(',')
	ctd2_count_temp2 = pd.Series(df_gene_list_paths['CTD2 Interactions']).str.count('N/A')
	ctd2_counts = ctd2_count_temp1 + 1 - ctd2_count_temp2


	for i in range(len(df_gene_list['KEGG Pathway'])):
		button_var_temp = pd.Series(['<button type="button" class="btn btn-info" data-toggle="collapse" data-target=#',str(i),'>Pathways (',str(output_counts[i]),')</button><div id="',str(i),'" class="collapse">',df_gene_list['KEGG Pathway'][i],'</div>']).str.cat(sep='')
		button_var = button_var.append(pd.Series(button_var_temp).reset_index(drop=True),ignore_index=True)
	df_gene_list['KEGG Pathway'] = button_var

	for j in range(len(df_gene_list['DGIdb Interactions'])):
		button_var_temp2 = pd.Series(['<button type="button" class="btn btn-info" data-toggle="collapse" data-target=#int',str(j),'>DGIdb Interactions (',str(df_gene_list['DGIdb #Interactions'][j]),')</button><div id="int',str(j),'" class="collapse">',df_gene_list['DGIdb Interactions'][j],'</div>']).str.cat(sep='')
		button_var2 = button_var2.append(pd.Series(button_var_temp2).reset_index(drop=True),ignore_index=True)
	df_gene_list['DGIdb Interactions'] = button_var2
	del df_gene_list['DGIdb #Interactions']

	if no_pathways_found == 0:
		for m in range(len(df_gene_list_paths['All Genes in Pathway'])):
			button_var_temp1 = pd.Series(['<button type="button" class="btn btn-info" data-toggle="collapse" data-target=#pathways', str(m),'>All Genes in Pathway (', str(all_genes_path_count[m]), ')</button><div id="pathways', str(m),'" class="collapse">', df_gene_list_paths['All Genes in Pathway'][m], '</div>']).str.cat(sep='')
			button_var1 = button_var1.append(pd.Series(button_var_temp1).reset_index(drop=True), ignore_index=True)
		df_gene_list_paths['All Genes in Pathway'] = button_var1

		for h in range(len(df_gene_list_paths['DGIdb Interactions'])):
			button_var_temp4 = pd.Series(['<button type="button" class="btn btn-info" data-toggle="collapse" data-target=#intp', str(h),'>DGIdb Interactions (', str(dgi_counts[h]), ')</button><div id="intp', str(h), '" class="collapse">',df_gene_list_paths['DGIdb Interactions'][h], '</div>']).str.cat(sep='')
			button_var4 = button_var4.append(pd.Series(button_var_temp4).reset_index(drop=True), ignore_index=True)
		df_gene_list_paths['DGIdb Interactions'] = button_var4

		for w in range(len(df_gene_list_paths['HGNC ID'])):
			button_var_temp5 = pd.Series(['<button type="button" class="btn btn-info" data-toggle="collapse" data-target=#hgnc',str(w),'>HGNC IDs (',str(df_gene_list_paths['# of Input Genes in Pathway'][w]),')</button><div id="hgnc',str(w),'" class="collapse">',df_gene_list_paths['HGNC ID'][w],'</div>']).str.cat(sep='')
			button_var5 = button_var5.append(pd.Series(button_var_temp5).reset_index(drop=True),ignore_index=True)
		df_gene_list_paths['HGNC ID'] = button_var5

		for q in range(len(df_gene_list_paths['CTD2 Interactions'])):
			button_var_temp6 = pd.Series(['<button type="button" class="btn btn-info" data-toggle="collapse" data-target=#ctd',str(q),'>CTD2 Interactions (',str(ctd2_counts[q]),')</button><div id="ctd',str(q),'" class="collapse">',df_gene_list_paths['CTD2 Interactions'][q],'</div>']).str.cat(sep='')
			button_var6 = button_var6.append(pd.Series(button_var_temp6).reset_index(drop=True),ignore_index=True)
		df_gene_list_paths['CTD2 Interactions'] = button_var6

		del df_gene_list_paths['DGIdb #Interactions']

	df_gene_list['DGIdb']=df_gene_list['DGIdb'].apply(urlify)
	df_gene_list['OMIM'] = df_gene_list['OMIM'].apply(urlify)
	df_gene_list['OMIM variants'] = df_gene_list['OMIM variants'].apply(urlify)
	df_gene_list['ExAC'] = df_gene_list['ExAC'].apply(urlify)
	df_gene_list['ClinVar']=df_gene_list['ClinVar'].apply(urlify)
	df_gene_list['gnomAD'] = df_gene_list['gnomAD'].apply(urlify)
	df_gene_list['Pharos'] = df_gene_list['Pharos'].apply(urlify)
	df_gene_list['Interaction Map'] = df_gene_list['Interaction Map'].apply(urlify_image)

	for k in range(len(df_gene_list['Interaction Map'])):
		button_var_temp3 = pd.Series(['<button type="button" class="btn btn-info" data-toggle="collapse" data-target=#map',str(k),'>Map</button><div id="map',str(k),'" class="collapse">',df_gene_list['Interaction Map'][k],'</div>']).str.cat(sep='')
		button_var3 = button_var3.append(pd.Series(button_var_temp3).reset_index(drop=True),ignore_index=True)
	df_gene_list['Interaction Map'] = button_var3

	if no_pathways_found == 0:
		#df_gene_list_paths['DGIdb'] = df_gene_list_paths['DGIdb'].apply(urlify)
		#df_gene_list_paths['OMIM'] = df_gene_list_paths['OMIM'].apply(urlify)
		#df_gene_list_paths['OMIM variants'] = df_gene_list_paths['OMIM variants'].apply(urlify)
		#df_gene_list_paths['ExAC'] = df_gene_list_paths['ExAC'].apply(urlify)
		#df_gene_list_paths['ClinVar'] = df_gene_list_paths['ClinVar'].apply(urlify)
		#df_gene_list_paths['gnomAD'] = df_gene_list_paths['gnomAD'].apply(urlify)
		#df_gene_list_paths['Pharos'] = df_gene_list_paths['Pharos'].apply(urlify)
		df_gene_list_paths['Gene Symbol'] = df_gene_list_paths['Gene Symbol'].apply(urlify_w)



	df_gene_list.index += 1
	df_gene_list_output.index += 1
	df_gene_list_hidden_table.index += 1
	df_gene_list_paths.index += 1
	df_gene_list_paths_output.index += 1
	df_gene_list_paths_hidden_table.index += 1

	return df_gene_list,df_gene_list_output,df_gene_list_hidden_table,df_gene_list_paths,df_gene_list_paths_output,df_gene_list_paths_hidden_table

import uuid

output_filename = str(uuid.uuid4()).encode('ascii','ignore')



from flask import Flask, make_response, send_file
app = Flask(__name__)

@app.route('/')
def home():
	return render_template('submission.html',N_MAX_GENES=N_MAX_GENES)

@app.route('/help')
def help():
	return render_template('help.html')

@app.route('/help/example_one')
def example_one():
	return send_file(dir_path_ex + '/riger_analysis_rep_1_shalem_et_al_2014.csv''',
					 mimetype='text/csv',
					 attachment_filename='riger_analysis_rep_1_shalem_et_al_2014.csv''',
					 as_attachment=True)

@app.route('/help/example_two')
def example_two():
	return send_file(dir_path_ex + '/riger_analysis_rep_2_shalem_et_al_2014.csv''',
					 mimetype='text/csv',
					 attachment_filename='riger_analysis_rep_2_shalem_et_al_2014.csv''',
					 as_attachment=True)

@app.route('/submit', methods=['GET', 'POST'])
def submit():
	if request.method == 'POST':

		list_of_genes=[]
		GENES_COMMON_PATH = int(request.form['num_common_pathway'])

		for line in request.form['gene_symbols'].splitlines():
			line = line.strip()
			if not line:continue

			genes_to_append=re.split(';|,| |\n',line)
			list_of_genes+=genes_to_append

		list_of_genes=[g.strip().upper() for g in list_of_genes if g.strip()]
		list_of_genes=list(set(list_of_genes)) #remove dups

		initial_list = len(list_of_genes)
		identified_gene = []
		matches = []
		df_matches = pd.Series()
		found_gene_iterate = pd.Series()
		for list_of_gene in list_of_genes:
			input_name_len = len(list_of_gene)
			matches_temp = df_hgnc_symbol.str.find(list_of_gene, start=0, end=None)
			df_matches = matches_temp.index[matches_temp == 0]

			if len(df_matches) == 1:
				found_gene_temp = df_hgnc_symbol.iloc[df_matches]
				found_gene_iterate = found_gene_iterate.append(found_gene_temp.reset_index(drop=True),ignore_index=True)
			elif len(df_matches > 1):
				len_of_input_name = []
				temps = df_hgnc_symbol.iloc[df_matches]
				for idx, temp in enumerate(temps):
					len_of_input_name.append(len(temp))

				if input_name_len in len_of_input_name:
					len_of_input_name.index(input_name_len)
					df_matches_sub1 = df_matches[len_of_input_name.index(input_name_len)]
					df_matches_sub = pd.Series(df_matches_sub1)
					found_gene_temp = df_hgnc_symbol.iloc[df_matches_sub]
					found_gene_iterate = found_gene_iterate.append(found_gene_temp.reset_index(drop=True),ignore_index=True)

		found_genes = pd.Series()
		full_list_of_genes_for_analysis = pd.Series()
		all_input_names_found = found_gene_iterate
		full_list_of_genes_for_analysis = full_list_of_genes_for_analysis.append(all_input_names_found)
		full_list_of_genes_for_analysis = full_list_of_genes_for_analysis.append(pd.Series(list_of_genes).reset_index(drop=True),ignore_index=True)
		pd.DataFrame(full_list_of_genes_for_analysis)
		found_genes = full_list_of_genes_for_analysis[full_list_of_genes_for_analysis.duplicated()]
		found_genes = found_genes.reset_index(drop=True)
		idx1 = pd.Index(all_input_names_found)
		idx2 = pd.Index(pd.Series(list_of_genes))
		missing_temp = idx2.difference(idx1)
		missing_genes = pd.Series(missing_temp.values)
		if len(missing_genes) < 1:
			missing_genes = pd.Series('All genes were identified and included in the analysis')
		list_of_genes = found_genes


		missing_genes = pd.DataFrame(missing_genes, columns=['Missing Genes'])

		df_gene_list, df_gene_list_output, df_gene_list_hidden_table, df_gene_list_paths, df_gene_list_paths_output, df_gene_list_paths_hidden_table = get_filled_dataframe(list_of_genes)

		import uuid

		missing_genes_output = pd.DataFrame()
		missing_genes_output['List of Missing Genes'] = missing_genes
		missing_genes_output.index += 1
		missing_genes.index += 1

		df_gene_list_output.to_csv('dtg_gene_' + output_filename + '.csv',encoding='utf-8')
		df_gene_list_paths_output.to_csv('dtg_pathway_' + output_filename + '.csv',encoding='utf-8')
		missing_genes.to_csv('dtg_missing_' + output_filename + '.csv',encoding='utf-8')

		return render_template('view.html', tables=[df_gene_list.to_html(columns=df_gene_list.columns[:], classes='report_gene', escape=False),df_gene_list_paths.to_html(columns=df_gene_list_paths.columns[:], classes='report_gene', escape=False),missing_genes.to_html(columns=missing_genes.columns[:],classes='report_gene', escape=False)]),
						   #titles=['na', 'Druggable Genes', 'Druggable Pathways','Missing Genes']),
		#return render_template('view.html', table=df_gene_list.to_html(columns=df_gene_list.columns[:-3], classes='report_gene', escape=False),render_template('druggable_gene.html',table=df_gene_list_paths.to_html(columns=df_gene_list.columns[:-3], classes='report_path', escape=False)

	else:
		return redirect('/')

#@app.route('/individual_gene', methods=['GET', 'POST']) # this is a job for GET, not POST
#def individual_gene():
	#individual_output = pd.read_csv(dir_path + '/dtg_gene_' + output_filename + '.csv''')
	#list_of_genes = 'TAL1'

	#return render_template('view2.html', tables=df_gene_list_trun.to_html(columns=df_gene_list_trun.columns[:],classes='report_gene', escape=False)

@app.route('/download_gene') # this is a job for GET, not POST
def download_gene():
	return send_file(dir_path + '/dtg_gene_' + output_filename + '.csv''',
	#return send_file('C:\Users\Matthew Canver\Desktop\DTG_Full\dtg_gene_' + output_filename + '.csv''',
					 mimetype='text/csv',
					 attachment_filename='dtg_gene_'+ output_filename + '.csv''',
					 as_attachment=True)

@app.route('/download_pathway') # this is a job for GET, not POST
def download_pathway():
	return send_file(dir_path + '/dtg_pathway_' + output_filename + '.csv''',
	#return send_file('C:\Users\Matthew Canver\Desktop\DTG_Full\dtg_gene_' + output_filename + '.csv''',
					 mimetype='text/csv',
					 attachment_filename='dtg_pathway_'+ output_filename + '.csv''',
					 as_attachment=True)

@app.route('/download_missing') # this is a job for GET, not POST
def download_missing():
	return send_file(dir_path + '/dtg_missing_' + output_filename + '.csv''',
	#return send_file('C:\Users\Matthew Canver\Desktop\DTG_Full\dtg_gene_' + output_filename + '.csv''',
					 mimetype='text/csv',
					 attachment_filename='dtg_missing_'+ output_filename + '.csv''',
					 as_attachment=True)

@app.errorhandler(404)
def page_not_found(e):
	return redirect('/')


if __name__ == "__main__":
	app.run(debug=True)
	excel.init_excel(app)
	app.run()


