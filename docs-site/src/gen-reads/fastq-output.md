# Fastq Output
The fastq output will have a key name that identifies the block where the read was drawn from, for quick comparisons in alignments. The output BAM file will contain the original sequence and cigar string. The name will have the format `EIDOLON_generated_<contig short name>_<fragment_start>_<fragment_end>_<uniq>/1` (or `/2` for the second read in a pair), where start and end are zero-padded to 10 digits and `<uniq>` is a 16-digit hex per-fragment tag that keeps same-position fragments from colliding (see #210).

```bash
@EIDOLON_generated_Chromosome_0000000000_0000000353_000000000000002a/1
CTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTTAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGC
+
>AC7<GDEGGGGEFGA<GFCGG;GGGGGF>GEGGEGGGFGFEFGCEGGGGGGCG:AEFGFFGG>FG;GDGA9$GGAGF=GFG=EFFCGGGFGGGGGC$BFGEFFAGGG9F7E@>?GFGGGG>EBFGFDG)DGC6DEDFA2EG:EGG%FFB
```
The above example comes from the chromosome in the reference (which was simply named "Chromosome") between 0 and 353, which is the size
of a fragment in the simulated DNA. If there is a pair with this read, it will have the same coordinates, though it started at index 352 instead.
