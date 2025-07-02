# Practice using markdown 
## For Biostars Bioinformatics workthrough

This will demonstrate various bioinformatics practices with the hope to use skills to analyse RNA-seq data


This command uses the accession no for mouse genome to download it and assocaited GFF3 file

```bash
datasets download genome accession GCF_000001635.27 --include gff3,genome
```
This method was very slow... kb/sec download speed. I instead downloaded on the NIH NCBI genome site: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/

Makefiles are useful files that let you target different elements. They require a tab/

``` bash
.RECIPEPREFIX = >

foo:
> echo Hello John!

bar:
> echo Hello Jane!
> echo Hello Everyone!
```
The .RECIPEPREFIX lets you use > instead of tabs.