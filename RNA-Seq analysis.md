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

I will use Salmon to do RNA Seq quantification against the Mouse transcriptome. Comparing to transcriptome allows us to see changes in transcript level against known transcripts.

We'd need to compare to genome if we want to look at different splicing, discover new mutations etc.

We have downloaded mouse transcriptome and annotations from Gencode:

```bash
# Download GENCODE mouse transcriptome (e.g., release M33 for GRCm39)
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.transcripts.fa.gz

# GTF from GENCODE
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz
```

We first have to build an Index of the transcriptome:

```bash
# Build the index (using 31-mers, default for reads >75bp)
salmon index \
  -t gencode.vM33.transcripts.fa.gz \  # input transcriptome FASTA
  -i salmon_mouse_index \              # name of output index folder
  -k 31                                # k-mer size, 31 is safe
    
```
I removed -- type quasi as this has been depreciated (caused an error initially) 

Then used salmon_quantification scrip to make count files for all the samples we're interested in (salmon_quantification.sh)

Important to make the script executable first:
```bash
chmod +x salmon_quantification.sh
```



** Insert other stuff here**


We can extract trancript and gene ID's from our transcripts file using grep with perl regular expressions and saving the output

```bash
grep -P -o 'ENSMUST\d{11}' gencode.vM33.transcripts.fa > ENMUST.txt
grep -P -o 'ENSMUSG\d{11}' gencode.vM33.transcripts.fa > ENMUSG.txt
```