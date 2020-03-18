# Mapping RNA-seq data

The following steps require the raw RNA-seq fastq files as input to generate _RAWcounts_ and _TPM_
files as output.
The following steps were performed on each of the fastq files.

## Trim adaptors from raw reads

```
bbduk.sh in=<FASTQ> out=<FASTQTRIM> minlen=20 ref=<ADAPTERS> ktrim=14 minkmerhits=7 mink=11 hdist=1 stats=<STATSFILE>
```

## Map to reference genome

```bash
hisat2 -q --dta --rna-strandness R -k 12 --summary-file <OUTSUMFILE> --met-file <METFILE> --no-unal -x <INDEX> -U <FASTQTRIM> -S <OUTsam>

samtools sort -o <OUTsortedbam> -T <TMP> <OUTsam>

samtools view -bq 5 <OUTsortedbam> <OUTuniqbam>
```

## Calculate raw counts

```bash
htseq-count -f <OUTuniqbam> -r pos -s reverse -t exon -i gene_id -m intersection-nonempty $bam <GRCz11GTF> > <RAWcounts>
```

## Calculate TPM per gene

```bash
kallisto quant -i <KALindx> -o <OUTDIR> -b 100 --rf-stranded --single -l 290 -s 20 <FASTQTRIM>
```
