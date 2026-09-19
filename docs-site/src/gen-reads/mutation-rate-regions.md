# Custom Mutation Rate Regions with a BED File
In addition to targeting certain regions, you can use a BED file with a custom field entered in any column after the third. 

```yaml
mutation_regions: /path/to/mutation_regions.bed
```

The text pattern is "mut_rate=0.0001" followed by a delimiter (end of line, space, semicolon, comma, bar) where the number after the equal sign is any float, and will be the mutation rate for the region defined in columns 1-3. 
```text
chr1    39930   39957   Bed_info    mut_rate=0.002
chr2    0   111199  Other_bed_info  mut_rate=0.02
```
Any region outside of the regions defined in the mutation regions BED will be assigned the default mutation rate, set by the mutation model, which can be overridden a custom default defined in the config file. For example, if you only want to have mutations in exomes, but you want reads from the full genome, you could set the `mutation_rate` to 0.0 and use a bed file defining all exomes, with a column appende `mut_rate=0.0011` or whatever mutation rate you desire. This BED works in conjunction with the `target_bed`, and any region outside the target bed will be excluded from reads anyway, and thus will have no variants. 
