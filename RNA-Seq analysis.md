# RNA-Seq analysis
Using paper from this study: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142432 

" Tissue neutrophil heterogeneity in germ free mice"

This command searches for the project metadata and saves in csv format

``` bash
bio search PRJNA597017 -H --csv -a > PRJNA597017.csv
```
The useful information is in the samples_title column

Using this command, the run info and sample info can be combined

``` bash
cat PRJNA597017.csv | csvcut -c run_accession,sample_title > run_sample_info.txt
```

Want to make a design file to arrange sample names