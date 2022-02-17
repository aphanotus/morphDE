# Supplemental Methods

for Angelini, Steele and O'Brien. **Gene expression during development of a dispersal polyphenism identifies regulators and effectors of growth in muscle and germline**. *Integrative & Comparative Biology*. 

> :warning: **Pre-submission DRAFT**.  Some of this text may be moved into the Methods section of the manuscript's main text.

In order to characterize the development of polyphenic traits in the wings, flight muscles and gonads of the soapberry bug, *Jadera haematoloma*, we sequenced transcriptomes from a multidimensional matrix of life stage, tissue, sex, food regime and host population.

## Sampling design

Transcriptome studies proceeded through three phases. Initially, we sampled whole bodies from 12 individuals, including three biological replicates of each sex and morph combination from the Plantation Key, FL population (phase 1). Next we sampled dorsal thorax and gonads from three individuals of each sex and morph combination from the Plantation Key, FL and Aurora, CO populations (phase 2). Samples from these phases were prepared as full-length 150-bp paired-end Illumina libraries and used to assembly a reference transcriptome. 

The final phase utilized 3'-tag sequencing with the exclusive goal of quantifying gene expression. The sampling design at this phase included two stages (fifth instars and nascent adults), two tissues (dorsal thorax and gonads), two sexes, long- and short-wing morphs, high and low food regimes, four populations (collected from Key Largo, FL, Plantation Key, FL, Aurora, CO, and Frederick, MD), and three biological replicates of each unique combination of the preceding factors. Thus the full matrix contained 384 samples. Seventeen samples were lost or did not pass library quality controls.

## Nucleic acid extractions

Individuals were collected from laboratory cultures, anesthetized using CO~2~ exposure, photographed, and flash frozen in liquid nitrogen. Samples were stored at -80˚C until further processing. The thorax was removed from head and abdomen with a clean scalpel blade, and the dorsal thorax was separated above the legs while tissue was still frozen. Gonads were then dissected from the abdomen of adults. For juveniles, the entire abdomen was included. Nucleic acid extractions used the Invitrogen PureLink RNA extraction kit.

## High-throughput RNA sequencing

Samples used for transcriptome assembly (phases 1 and 2) were delivered on dry ice to Beckman Coulter Genomics. :warning: **More details!**

RNA samples were shipped overnight on dry ice to the DNA Technologies Core at the Genome and Biomedical Sciences Facility of the University of California at Davis. Messenger RNA was poly-A selected from these samples, reverse transcribed, and used to construct 90-bp single-end Lexogen libraries. These four samples were run on one Illumina HiSeq lane producing 1,487,663,804 reads. These data were inspected for quality using FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and trimmed using Trimmomatic versions 0.33, as in the example below.

```bash
java -jar trimmomatic-0.33.jar SE \
 raw.data/A1-d25-2_S161_L006_R1_001.fastq \
 filtered.reads/A1-d25-2_S161_L006_R1_001.trimmed.fastq \
 ILLUMINACLIP:/export/local/src/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 \
 HEADCROP:12 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:50
```

The remaining 76 RNA samples were used for 3′-end tag sequencing (3-Seq; (*30*, *31*)). Following the protocol of the QuantSeq 3' mRNA-Seq Library Prep Kit (Lexogen, Greenland, New Hampshire, USA) produced 90-bp single-end libraries, which were sequenced in one Illumina HiSeq lane. This run included an additional 20 samples that were not a part of this study, for a total of 96 multiplexed samples.

## Annotation of the transcriptome

First collapse transcript isoforms using CD-HIT-EST.

```bash
cd /research/drangeli/phase2_BCG_RNAseq/Jhae_GU_trinity_assembly/
cd-hit-est -i GU.Trinity.fa -o lh GU.Trinity.c90.fa -c 0.90 -n 10 -d 0 -M 64000 -T 12
```

This took about 3.5 days on node 26.

Then run EnTAP from node 28 against the transcriptome assembly. This call will utilize EnTAP databases set-up for another project, including the invertebrate RefSeq dataset, SwissProt and Trembl.

```bash
EnTAP --runP --threads 24 \
-i /research/drangeli/phase2_BCG_RNAseq/Jhae_GU_trinity_assembly/GU.Trinity.c90.fa \
-d /export/groups/drangeli/entap_config/bin/invert.dmnd \
-d /export/groups/drangeli/entap_config/bin/uniprot_sprot.dmnd \
-d /export/groups/drangeli/entap_config/bin/uniprot_trembl.dmnd \
--ini /export/groups/drangeli/Jhae.genome/entap_config.ini 
```

The sequences were matched to species that make sense, primarily *Halyomorpha halys*, *Cimex lectularius* and *Nilaparavata lugens* for RefSeq and *Oncopeltus*, *Lygus* and *Riptortus* from TrEMBL. The contaminant sequence species are also as expected, bacteria and yeasts. Just 12 sequences with homology to humans! 

Importantly, EnTAP assigns GO terms to sequences. It will probably be worth running EnTAP again on the masked genome and its gene predictions. But the EnTAP transcriptome annotation may also be useful. Later, I'll want to map transcripts from this transcriptome assembly to the genome assembly. 

## Inspecting the raw 3seq reads

We had (how many?) samples sequenced using 3′-tag sequencing (3seq) at the UC Davis [DNA Technologies & Expression Analysis Core](https://dnatech.genomecenter.ucdavis.edu/).

Working from `/export/groups/drangeli/Jhae.morph.txome/davis.3seq/`.

I spot-checked several samples from each plate. For example: 

```bash
fastqc 200611/DAJHX2_UMI_S187_L008_R1_001.fastq.gz --outdir=fastqc
fastqc 200701/DAJHX2_UMI_S2_L001_R1_001.fastq.gz --outdir=fastqc
```

Reads appear to be high quality. Even samples flagged by the center as "degraded" produced 3' libraries with high quality reads. Read length and quantity were also consistent.

## Trimming 3seq reads

Quality trimming was preformed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) version 0.33. A simple [bash script](https://github.com/aphanotus/Jhae.genome/blob/main/trimming.3seq.reads.sh) was used to apply trimming to all read files and see the output in a separate folder. The parameters were `ILLUMINACLIP:/export/local/src/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:16 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:79`. Fastqc inspection afterward revealed notable improvements.

## Mapping 3seq reads to the transcriptome

[Kallisto](https://pachterlab.github.io/kallisto/about) will quantify the abundance of transcripts from RNA-Seq data using a pseudoalignment method. First, create a Kallisto index from the assembled transcripts, collapsed by cd-hit-est. I intentionally avoided the TransDecoder CDS files, because if there are fragmentary or non-coding RNAs that were sequenced, then some reads will correspond to those RNAs too.

```bash
cd /research/drangeli/phase2_BCG_RNAseq/Jhae_GU_trinity_assembly/

kallisto index --index=GU.Trinity.c90.kallisto.index GU.Trinity.c90.fa
```

In order to run Kallisto on single-end sequence, it wants to know the mean and standard deviation of read lengths. First, use bash to find the total number of filtered reads.

```bash
cd /export/groups/drangeli/Jhae.morph.txome/davis.3seq/filtered.reads

gunzip -c 180824/*.fastq.gz 190917/*.fastq.gz 200611/*.fastq.gz 200701/*.fastq.gz \
| sed -n -e '2~4p' \
| awk '{ print length($0); }' > read.lengths.txt

wc -l read.lengths.txt
```

That’s 1,149,235,725 reads. Next, calculate the mean and SD.

```bash
cat read.lengths.txt | awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++) { printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)} }'
```

- mean = 84.892739  bp
- SD = 0.507056 bp

Next, I wrote a [bash script](https://github.com/aphanotus/Jhae.genome/blob/main/Jhae.kallisto.mapping.sh) to run `kallisto quant` for each sample and to rename the output files with the sample names. The result is a folder `kallisto.vs.JhaeGUtx` that contains `.abundance.tsv` files for each sample, listing the counts for each transcript. 

## Summary of mapping results

3seq reads mapped to `/export/groups/drangeli/Jhae.morph.txome/davis.3seq/`

- genome, via hisat2: `./hisat2.vs.genome`
- transcriptome, via kallisto: `./kallisto.vs.JhaeGU.tx`

:information_source: **Next: Collate kallisto and hisat2 results into a tabular form for analysis by DEseq2.**

