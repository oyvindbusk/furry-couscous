# furry-couscous - genepanelbuilder
Simple script to generate genelists with geneID and bed files with exons. Also checks if geneID is valid.

## How to obtain flatfile with all transcripts:
* In the UCSC browser. Select track NCBI RefSeq - table RefSeq Select (get autodl of this). all fields from selected table.
* * All the alternative chromosomes are removed and "chr" is removed before chromosome name.
* In the UCSC browser. Select track NCBI RefSeq - table RefSeq All (get autodl of this). all fields from selected table.
* In the UCSC browser. Select track NCBI RefSeq - table RefSeq Other (get autodl of this). all fields from selected table.
* From genenames.org (ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/) dl protein-coding_gene.txt

## How it is run:
```
python panel_builder.py -h
usage: panel_builder.py [-h] -g GENELIST -n NAME

Takes a file with gene names as inputm and queries each name against
genenames.org, to check if the genename is present

optional arguments:
  -h, --help            show this help message and exit
  -g GENELIST, --genelist GENELIST
                        Input file
  -n NAME, --name NAME  Name of panel

# Example: 
python panel_builder.py -g genepanels/alle_genlister/Vask_bindevev_v5_080321.txt -n Vask_bindevev_v5_080321

``` 

## What is does:
In addition to generating files, the script does the following:
* Searches databases on the internet for the genes in the txt
* Checks for duplicates in the list
* What is not found in refseq-select is searched for in refseq-all, then finally gencode
* What is not found is reported
* Noncoding genes are also reported

## Ouput:
A folder with the name of the genepanel is created with the following files:
*   genpanel_duplicated.regions.txt - duplicated regions (if any)
*   genpanel_full_gene.BED - bedfile with one row pr gene
*   genpanel.BED - bedfile with one row pr exon
*   genpanel.genetikkportalimport.txt - file used to import gene panel to genetikkportalen.no
*   genpanel.txt - The original gene list used as input







