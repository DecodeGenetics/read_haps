## read_haps contamination experiments

Experiments for read_haps paper were run using commands in batch and batch_raw.

Symlinks called `genome.fa` and `genome.fa.fai` should be created in the working directory and point to a GRCh38 reference. Symlinks to read_haps and verifyBamID should also be created in working directory. Further, a VCF needs to generated using

```sh
zcat 1000G_EUR.vcf.gz > 1000G_EUR.vcf
```

A directory called GIAB is assumed containing GIAB bam and vcf files. These need to be downloaded, see more in GIAB/README.md .
