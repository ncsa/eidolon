# Reading Bed Data
`eidolon` can now read bed data and filter the reads based on regions specified. This comes with some caveats. First is that rust can't handle any header lines or non-header rows in the bed file, because of potential range of possibilities. Second, the names must match what is in the fasta. Third, `eidolon` can only filter, and does not use the rest of the bed information. It treats each record as a region of interest only.

Bed data will have some challenges. For example, if the contig names in the bed file don't match the assumed contig name from the fasta file, how will rust be able to know which is what? To make things work, we are going to assume, for now, that the bed file contig names match the names as derived by `eidolon`. To determine this:

First, in a linux terminal, use your input fasta to find the names of the contigs like this:

```bash
$ grep "^>" input.fasta
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
>NC_001134.8 Saccharomyces cerevisiae S288C chromosome II, complete sequence
>NC_001135.5 Saccharomyces cerevisiae S288C chromosome III, complete sequence
>etc
```
This will give you the full fasta name. `eidolon` determines the short name of the input fasta by skipping the initial '>' character, then taking the string up to the first whitespace delimiter. In the above example, the short names are "NC_001133.9", "NC_001134.8", etc. 

Next, we check the bed. This example awk command will look for unique values in the first column of your bed file and print out what it finds. 

```bash
$ awk '!seen[$1]++ {print $1}' file.bed
1
2
3
etc
```
In this case, the bed file uses a simple numbering scheme to number the chromosomes. We can see a mapping with the roman numeral chromosomes in the name, but it's tricky to handle consistently. So, the following command will allow you, the user, to map these chromosomes and thus instruct `eidolon`. Yours will vary based on the inputs.

```bash
awk -F'\t' '{$1=($1=="1"?"NC_001133.9":$1); $1=($1=="2"?"NC_001134.8":$1); print}' OFS='\t' file.bed > file_renamed.bed
```
You will have to tailor this command to your dataset, but this is one way to map the names from the bed to the fasta in a way `eidolon` will understand.
