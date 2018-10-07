<b>DrugThatGene</b> is a a web-based application to automate the analysis of potential therapeutic targets and/or therapeutic pathways identified from functional genetic screens. DTG integrates data from human genetic databases and small molecule databases (See 'Databases Utilized by DrugThatGene' heading). A list of "hits" from a pooled CRISPR screen is required in order to use DTG to aid in the identification of drugs/small molecules.

<b>Databases Utilized by DrugThatGene</b>
<br>Here is the list of databases utilied by DTG:
<br>-The Drug Gene Interaction Database (DGIdb): Database of drugs with known targets
<br>-Cancer Target Discovery and Development (CTD2): Database of drugs with known targets with a focus on cancer targets
<br>-Genome Aggregation Database (gnomAD): Database aggregating exome and genome sequencing data. The database consists of >120,000 exomes and >15,000 whole genome sequences.
<br>-ClinVar: Database of relationships between human variation and phenotypes
<br>-Online Mendelian Inheritance in Man (OMIM): A database cataloging human genes and genetic disorders
<br>-Exome Aggregation Consortium (ExAC): Database of exomes from 60,706 unrelated individuals sequenced as part of various disease-specific and population genetic studies
<br>-Pharos: A comprehensive, integrated database of targets within the druggable genome
<br>-STRING: Database of protein-protein interaction networks
<br>-CORUM: The comprehensive resource of mammalian protein complexes
<br>-Kyoto Encyclopedia of Genes and Genomes (KEGG) Pathways: Pathways derived from the KEGG database based on molecular interaction, reaction and relation networks

<b>DrugThatGene Table Outputs</b>

<i>"Druggable Genes" Tab</i>
<br>Here is a description of each of the table outputs:
<br>-Gene Symbol: Official gene symbol for each gene from the input gene list.
<br>-HGNC ID: HGNC ID for the associated official gene symbol to allow for verification that the correct gene has been identified.
<br>-KEGG Pathways: List of KEGG Pathways for the input gene.
<br>-Complexes: List of protein complexes for the input gene.
<br>-DGIdb #interactions: Number of drug-gene interactions identified within the DGIdb database for the input gene.
<br>-DGIdb interactions: List of drugs/small molecules known to interact with the input gene from the DGIdb database.
<br>-DGIdb: Link to DGIdb for the input gene. See above for a description of the DGIdb database.
<br>-CTD2 interactions: List of drugs/small molecules known to interact with the input gene from the CTD2 database.
<br>-OMIM: Link to OMIM for the input gene. See above for a description of the OMIM database.
<br>-OMIM variants: Disease associated-variants in the OMIM database
<br>-ClinVar: Link to ClinVar for the input gene. See above for a description of the ClinVar database.
<br>-gnomAD: Link to gnomAD for the input gene. See above for a description of the gnomAD database.
<br>-Pharos: Link to Pharos for the input gene. See above for a description of the Pharos database.
<br>-ExAC: Link to ExAC for the input gene. See above for a description of the ExAC database.
<br>-ExAC #LoF: Number of loss of function (LoF) variants present within the ExAC database.
<br>-ExAC Missense Z: Z score for the deviation of observed counts from the expected number. Positive Z scores indicate increased constraint (intolerance to variation) and therefore that the gene had fewer variants than expected. Negative Z scores are given to genes that had a more variants than expected.
<br>-ExAC pLI: pLI is a score that indicates the probability that a gene is intolerent to a loss of function mutation. The closer pLI is to one, the more LoF intolerant the gene appears to be. pLI >= 0.9 is considered to be an extremely LoF intolerant gene.
<br>-Interaction Map: Protein-protein interaction network from STRING

<i>"Druggable Pathways" Tab</i>
<br>This tab displays data grouped by common pathways. The output table is sorted by the fraction of input genes in each pathway (# of input genes in pathway/total number of genes in pathway)
<br>-KEGG Pathways: Common KEGG Pathways from the input gene list.
<br>-Gene Symbol: Official gene symbol for each gene from the input gene list.
<br>-HGNC ID: HGNC ID for the associated official gene symbol to allow for verification that the correct gene has been identified.
<br>-All Genes in Pathway: List of all genes within the relevant KEGG pathway.
<br>-DGIdb interactions: List of drugs/small molecules known to interact with the input gene from the DGIdb database.
<br>-CTD2 interactions: List of drugs/small molecules known to interact with the input gene from the CTD2 database.
<br>-# of Input Genes in Pathway: The number of genes from the input gene list in the relevant KEGG pathway. The number is parentheses is the fraction of input genes in each pathway (# of input genes in pathway/total number of genes in pathway)

<i>"Druggable Complexes" Tab</i>
<br>This tab displays the data grouped by common protein complexes. The output table is sorted by the fraction of input genes in each protein complex (# of input genes in protein complex/total number of genes in protein complex)
<br>-Complexes: Common CORUM protein complexes from the input gene list.
<br>-Gene Symbol: Official gene symbol for each gene from the input gene list.
<br>-HGNC ID: HGNC ID for the associated official gene symbol to allow for verification that the correct gene has been identified.
<br>-All Genes in Complex: List of all genes within the relevant protein complex.
<br>-Complex Function: Function of the relevant protein complex (determined by FunCat).
<br>-GO Description: Gene ontology (GO) description for the relevant protein complex.
<br>-DGIdb interactions: List of drugs/small molecules known to interact with the input gene from the DGIdb database.
<br>-CTD2 interactions: List of drugs/small molecules known to interact with the input gene from the CTD2 database.
<br>-# of Input Genes in Complex: The number of genes from the input gene list in the relevant protein complex. The number is parentheses is the fraction of input genes in each protein complex (# of input genes in protein complex/total number of genes in protein complex)

<i>"Missing Genes" Tab</i>
<br>Entered genes that are not Official Gene Symbols (refer to the HGNC database: https://www.genenames.org/) are excluded from analysis. The list of these "missing genes" is presented here.

<b>Performing DTG Analysis</b>
<br>Steps to perform DTG analysis:
<br>(1) Copy and paste a list of official gene symbols (https://www.genenames.org/) for human genes separated by spaces or on individual lines into the gene list field. The input is case insensitive and duplicate entries will be excluded.
<br>(2) Enter a number (for example "2", but not "two") to set the number of genes in a pathway required for the software to identify a "common pathway".
<br>(3) Enter a number (for example "2", but not "two") to set the number of genes in a protein complex required for the software to identify a "common protein complex".
<br>(4) Click the "Submit" button
<br>(5) The analysis will be output with four tables: "Druggable Genes", "Druggable Pathways", "Druggable Complexes", and "Missing Genes"
<br>(6) To download the analyses in .csv format, click "Download Druggable Gene Analysis", "Download Druggable Pathway Analysis", and/or "Download Druggable Complex Analysis"
<br>(7) Entered genes that are not Official Gene Symbols (refer to the HGNC database: https://www.genenames.org/) are excluded from analysis. To download a list of genes that were excluded from the analysis for this reason, click "Download Genes Missing from Analysis"

<b>Run DTG on Your Local Machine</b>
<br>(1) Download and install Docker from this link: https://docs.docker.com/engine/installation/.
<br>(2) After Docker installation, type the following command:
<br>docker run -it pinellolab/drugthatgene -p 9999:9999
<br>(3) The website can be accessed at localhost:9999
<br>

<b>Run DTG on Your Local Machine with Increased Gene List Limit</b>
<br>(1) Download all files within the DrugThatGene folder and 'run.py' from https://github.com/pinellolab/DrugThatGene. Put the DrugThatGene folder and 'run.py' in the same directory.
<br>(2) To increase the limit for the number of genes run by DTG, change the value for the N_MAX_GENES variable within the __init__.py file (default value is 200). Increasing the value above 200 will allow lists >200 genes to be run by DTG. 
<br>(3) Run the file 'run.py' in a terminal.
<br>(4) The website can be accesssed at localhost:9999

<b>Example Data</b>
<br>(1) 130 essential genes for AML cell survival <i>in vitro</i> and <i>in vivo</i> [1] - available for download from this Github site (https://github.com/pinellolab/DrugThatGene/blob/master/DrugThatGene/AML_essential_genes_yamauchi_et_al_2018.csv)
<br>(2) Top 100 genes from Riger analysis from genome-wide CRISPR screen in A375 cells (replicate #1) [2] - available for download from this Github site (https://github.com/pinellolab/DrugThatGene/blob/master/DrugThatGene/riger_analysis_rep_1_shalem_et_al_2014.csv)</li>
<br>(3) Top 100 genes from Riger analysis from genome-wide CRISPR screen in A375 cells (replicate #2) [2] - available for download (replicate #2 gene list) from this Github site (https://github.com/pinellolab/DrugThatGene/blob/master/DrugThatGene/riger_analysis_rep_2_shalem_et_al_2014.csv)
<br><br>[1] Yamauchi, T. et al. (2018). Genome-wide CRISPR-Cas9 Screen Identifies Leukemia-Specific Dependence on a Pre-mRNA Metabolic Pathway Regulated by DCPS. <i>Cancer Cell</i>, 33, 386–400.e5.
<br>[2] Shalem, O. et al. (2014). Genome-scale CRISPR-Cas9 knockout screening in human cells. <i>Science</i>, 343, 84–87.
<br>

<b>How to Cite DrugThatGene</b>
<br>If you use DrugThatGene in your work, please cite:

Canver, MC., Bauer, DE., Maeda, T., Pinello, L. (2018). DrugThatGene: integrative analysis to streamline the identification of druggable genes, pathways, and protein complexes from CRISPR screens.
