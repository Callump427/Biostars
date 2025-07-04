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

We're just interested in the lung samples, so I extracted the required runs before we get the data from the SRA: 

```bash
cat run_sample_info.txt | grep lung | cut -f 1 -d , > wanted_runs.txt
```

Then can run fastq-dump with parallel command to do it for all accessions in the wanted.txt file

```bash
cat wanted_runs.txt | parallel fastq-dump --split-files {}
```

We then ran fastq on all samples in the working directory, followed by multiqc to aggregate the results into 1 analysis file 
```bash
fastqc *.fastq
mutiqc .
```
The dot just means do this in the current directory.