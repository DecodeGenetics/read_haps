## Testing read_haps
You can run `make test` in the top-level directory to run a simple test on a small dataset. The small test does not require any external data. We also provide results of genome-wide tests on public data.

### HG002 60x BAM
We ran a test using public GIAB BAM file for [HG002](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam) and compared against its [truth set](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz).

```sh
read_haps -fa ${REFERENCE} ${HG002.BAM} high_quality_markers_deCODE_2015.txt ${HG002.VCF} > test_out
```
where `${REFERENCE}` is a GRCh38 reference, ${HG002.BAM} and ${HG002.VCF} are the GIAB BAM and truth VCF for HG002, respectively.

The running time (h:mm:ss) was:
real 2:05:33
user 2:04:05
sys  0:01:23

We ran it on a machine with Intel(R) Xeon(R) CPU E5-v4 @ 2.60GHz on a single core. The output is shown in `test_expected` file.

```sh
$ cat test_expected
SNP_PAIRS ERROR_PAIRS DOUBLE_ERROR_PAIR_COUNT DOUBLE_ERROR_FRACTION REL_ERROR_FRACTION NONSENSE_FRACTION PASS_FAIL REASON
1169470 1842 58 4.95951e-05 0.000168402 0.00163578 PASS -
```

### Identifying SNP pairs that show evidence of three haplotypes
We also have a test where we identify evidence for three haplotypes

```sh
read_haps --pairs --fa small_genome.fa HG002_chr1_10_20MB.bam small_hq_markers HG002_chr1_10_20MB.vcf.gz > small_pairs
diff small_pairs small_pairs_expected # compare
```

One marker shows limited evidence of three haplotypes:
$ awk '$1 ~ /chr/ && $3 > 0 && $4 > 0' small_pairs 
chr1:17393039 chr1:17393050 1 1 0

No markers show strong evidence of threee haplotypes. 