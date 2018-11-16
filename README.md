# Gff3Curator


## Descriptions

These are the scripts I used to curate the gene models predicted by Augustus.
Need to view those gene models by [Integrative Genomics Viewer (IGV)](http://software.broadinstitute.org/software/igv/)

## Procedures

### wgff.20160111.sh

This is the main script

> wgff.20160111.sh ScaffoldNAme

### wheat_annotation_extractor20150914.sh

Configure **wheat_annotation_extractor20150914.sh** so you can extract each evidence (EST, protein and _ab initio_ predictions)

### gff_sum_exons_for_annotation20150914.pl

collect all the exons for your choice

> gff_sum_exons_for_annotation20150914.pl fasta region mincount *.gff[.gz]

### gff_border_20151001.pl

Extract Exon borders so you can see exon coordinates in IGV

> gff_border_20151001.pl in.gff.gz out.gff.gz

## Author:

Fu-Hao Lu

Post-Doctoral Scientist in Micheal Bevan laboratory

Cell and Developmental Department, John Innes Centre

Norwich NR4 7UH, United Kingdom

E-mail: <Fu-Hao.Lu@jic.ac.uk>
