DrugThatGene is a a web-based application to automate the analysis of potential therapeutic targets and/or therapeutic pathways identified from functional genetic screens. DTG integrates data from human genetic databases and small molecule databases (See 'Databases Utilized by DrugThatGene' heading). A list of "hits" from a pooled CRISPR screen is required in order to use DTG to aid in the identification of drugs/small molecules.
Databases Utilized by DrugThatGene
Here is the list of databases utilied by DTG:
The Drug Gene Interaction Database (DGIdb): Database of drugs with known targets
Cancer Target Discovery and Development (CTD2): Database of drugs with known targets with a focus on cancer targets
Genome Aggregation Database (gnomAD): Database aggregating exome and genome sequencing data. The database consists of >120,000 exomes and >15,000 whole genome sequences.
ClinVar: Database of relationships between human variation and phenotypes
Online Mendelian Inheritance in Man (OMIM): A database cataloging human genes and genetic disorders
Exome Aggregation Consortium (ExAC): Database of exomes from 60,706 unrelated individuals sequenced as part of various disease-specific and population genetic studies
Pharos: A comprehensive, integrated database of targets within the druggable genome
STRING: Database of protein-protein interaction networks
CORUM: The comprehensive resource of mammalian protein complexes
Kyoto Encyclopedia of Genes and Genomes (KEGG) Pathways: Pathways derived from the KEGG database based on molecular interaction, reaction and relation networks
DrugThatGene Table Outputs
"Druggable Genes" Tab
Here is a description of each of the table outputs:
Gene Symbol: Official gene symbol for each gene from the input gene list.
HGNC ID: HGNC ID for the associated official gene symbol to allow for verification that the correct gene has been identified.
KEGG Pathways: List of KEGG Pathways for the input gene.
Complexes: List of protein complexes for the input gene.
DGIdb #interactions: Number of drug-gene interactions identified within the DGIdb database for the input gene.
DGIdb interactions: List of drugs/small molecules known to interact with the input gene from the DGIdb database.
DGIdb: Link to DGIdb for the input gene. See above for a description of the DGIdb database.
CTD2 interactions: List of drugs/small molecules known to interact with the input gene from the CTD2 database.
OMIM: Link to OMIM for the input gene. See above for a description of the OMIM database.
OMIM variants: Disease associated-variants in the OMIM database
ClinVar: Link to ClinVar for the input gene. See above for a description of the ClinVar database.
gnomAD: Link to gnomAD for the input gene. See above for a description of the gnomAD database.
Pharos: Link to Pharos for the input gene. See above for a description of the Pharos database.
ExAC: Link to ExAC for the input gene. See above for a description of the ExAC database.
ExAC #LoF: Number of loss of function (LoF) variants present within the ExAC database.
ExAC Missense Z: Z score for the deviation of observed counts from the expected number. Positive Z scores indicate increased constraint (intolerance to variation) and therefore that the gene had fewer variants than expected. Negative Z scores are given to genes that had a more variants than expected.
ExAC pLI: pLI is a score that indicates the probability that a gene is intolerent to a loss of function mutation. The closer pLI is to one, the more LoF intolerant the gene appears to be. pLI >= 0.9 is considered to be an extremely LoF intolerant gene.
Interaction Map: Protein-protein interaction network from STRING
"Druggable Pathways" Tab
This tab displays data grouped by common pathways. The output table is sorted by the fraction of input genes in each pathway (# of input genes in pathway/total number of genes in pathway)
KEGG Pathways: Common KEGG Pathways from the input gene list.
Gene Symbol: Official gene symbol for each gene from the input gene list.
HGNC ID: HGNC ID for the associated official gene symbol to allow for verification that the correct gene has been identified.
All Genes in Pathway: List of all genes within the relevant KEGG pathway.
DGIdb interactions: List of drugs/small molecules known to interact with the input gene from the DGIdb database.
CTD2 interactions: List of drugs/small molecules known to interact with the input gene from the CTD2 database.
# of Input Genes in Pathway: The number of genes from the input gene list in the relevant KEGG pathway. The number is parentheses is the fraction of input genes in each pathway (# of input genes in pathway/total number of genes in pathway)
"Druggable Complexes" Tab
This tab displays the data grouped by common protein complexes. The output table is sorted by the fraction of input genes in each protein complex (# of input genes in protein complex/total number of genes in protein complex)
Complexes: Common CORUM protein complexes from the input gene list.
Gene Symbol: Official gene symbol for each gene from the input gene list.
HGNC ID: HGNC ID for the associated official gene symbol to allow for verification that the correct gene has been identified.
All Genes in Complex: List of all genes within the relevant protein complex.
Complex Function: Function of the relevant protein complex (determined by FunCat).
GO Description: Gene ontology (GO) description for the relevant protein complex.
DGIdb interactions: List of drugs/small molecules known to interact with the input gene from the DGIdb database.
CTD2 interactions: List of drugs/small molecules known to interact with the input gene from the CTD2 database.
# of Input Genes in Complex: The number of genes from the input gene list in the relevant protein complex. The number is parentheses is the fraction of input genes in each protein complex (# of input genes in protein complex/total number of genes in protein complex)
"Missing Genes" Tab
Entered genes that are not Official Gene Symbols (refer to the HGNC database: https://www.genenames.org/) are excluded from analysis. The list of these "missing genes" is presented here.
Performing DTG Analysis
Steps to perform DTG analysis:
Copy and paste a list of official gene symbols (https://www.genenames.org/) for human genes separated by spaces or on individual lines into the gene list field. The input is case insensitive and duplicate entries will be excluded.
Enter a number (for example "2", but not "two") to set the number of genes in a pathway required for the software to identify a "common pathway".
Enter a number (for example "2", but not "two") to set the number of genes in a protein complex required for the software to identify a "common protein complex".
Click the "Submit" button
The analysis will be output with four tables: "Druggable Genes", "Druggable Pathways", "Druggable Complexes", and "Missing Genes"
To download the analyses in .csv format, click "Download Druggable Gene Analysis", "Download Druggable Pathway Analysis", and/or "Download Druggable Complex Analysis"
Entered genes that are not Official Gene Symbols (refer to the HGNC database: https://www.genenames.org/) are excluded from analysis. To download a list of genes that were excluded from the analysis for this reason, click "Download Genes Missing from Analysis"
Example Data
Top 100 genes from Riger analysis from genome-wide CRISPR screen in A375 cells (replicate #1):  Download replicate #1 gene list
Top 100 genes from Riger analysis from genome-wide CRISPR screen in A375 cells (replicate #2):  Download replicate #2 gene list

Riger analyses are from: Shalem, O., Sanjana, N.E., Hartenian, E., Shi, X., Scott, D.A., Mikkelsen, T.S., Heckl, D., Ebert, B.L., Root, D.E., Doench, J.G., et al. (2014). Genome-scale CRISPR-Cas9 knockout screening in human cells. Science 343, 84–87.
How to Cite DrugThatGene
If you use DrugThatGene in your work, please cite:

Canver, MC., Bauer, DE., Maeda, T., Pinello, L. (2018). DrugThatGene: integrative analysis to streamline the identification of druggable genes, pathways, and protein complexes from CRISPR screens. Manuscript in preparation.
