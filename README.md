# Supplemental Methods

for Angelini, Steele and O'Brien. **Gene expression during development of a dispersal polyphenism identifies regulators and effectors of growth in muscle and germline**. *Integrative & Comparative Biology*. 

> :warning: **Pre-submission DRAFT**.  Some of this text may be moved into the Methods section of the manuscript's main text. Details that may be important now, but should likely be removed before publication are given as block quotes.

> ℹ️ While any process is running on nscc, you can monitor it's progress indirectly by checking the CPU usage on [node 26](http://schupflab.colby.edu/mrtg/nscc/nscc-n26.load.html) or [n28](http://schupflab.colby.edu/mrtg/nscc/nscc-n28.load.html).

In order to characterize the development of polyphenic traits in the wings, flight muscles and gonads of the soapberry bug, *Jadera haematoloma*, we sequenced transcriptomes from a multidimensional matrix of life stage, tissue, sex, food regime and host population.

## Sampling design

Transcriptome studies proceeded through three phases. Initially, we sampled whole bodies from 12 individuals, including three biological replicates of each sex and morph combination from the Plantation Key, FL population (phase 1). Next we sampled dorsal thorax and gonads from three individuals of each sex and morph combination from the Plantation Key, FL and Aurora, CO populations (phase 2). Samples from these phases were prepared as full-length 150-bp paired-end Illumina libraries and used to assembly a reference transcriptome. 

The final phase utilized 3'-tag sequencing with the exclusive goal of quantifying gene expression. The sampling design at this phase included two stages (fifth instars and nascent adults), two tissues (dorsal thorax and gonads), two sexes, long- and short-wing morphs, high and low food regimes, four populations (collected from Key Largo, FL, Plantation Key, FL, Aurora, CO, and Frederick, MD), and three biological replicates of each unique combination of the preceding factors. Thus the full matrix contained 384 samples. Seventeen samples were lost or did not pass library quality controls.

## Nucleic acid extractions

Individuals were collected from laboratory cultures, anesthetized using CO~2~ exposure, photographed, and flash frozen in liquid nitrogen. Samples were stored at -80˚C until further processing. The thorax was removed from head and abdomen with a clean scalpel blade, and the dorsal thorax was separated above the legs while tissue was still frozen. Gonads were then dissected from the abdomen of adults. For juveniles, the entire abdomen was included. Nucleic acid extractions used the Invitrogen PureLink RNA extraction kit.

## Sequencing and assembly of a reference transcriptome

Samples used for transcriptome assembly (phases 1 and 2) were delivered on dry ice to Beckman Coulter Genomics (Danvers, Massachusetts) for poly-A selection, preparation of 125 bp paired-end TruSeq libraries, and sequencing using Illumina HiSeq. Libraries from phases 1 and 2 were sequenced separately. Reads from all samples were trimmed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) version 0.33 ([Bolger et al. 2014](https://doi.org/10.1093/bioinformatics/btu170)) with the parameters `ILLUMINACLIP:/export/local/src/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:12 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:50`.  The resulting reads were assembled into transcripts by [Trinity](http://trinityrnaseq.github.io/) version 2.0.6 ([Grabherr et al. 2011](https://www.nature.com/articles/nbt.1883)). Before annotation, potentially redundant transcripts with greater than 90% identity were collapsed using CD-HIT-EST version 4.6 ([Li & Godzik 2006](https://doi.org/10.1093/bioinformatics/btl158)).  

```bash
cd /research/drangeli/phase2_BCG_RNAseq/Jhae_GU_trinity_assembly/
cd-hit-est -i GU.Trinity.fa -o GU.Trinity.c90.fa -c 0.90 -n 10 -d 0 -M 64000 -T 12
```
> This took about 3.5 days on node 26.

[TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki) identified coding sequences. At each stage, [BUSCO](https://busco.ezlab.org/) version 5.2.2 ([Manni et al. 2021](https://doi.org/10.1093/molbev/msab199)) was used to assess transcriptome completeness.

```bash
busco --in GU.Trinity.fa \
--lineage_dataset /research/drangeli/DB/busco.lineages/hemiptera_odb10 \
--out BUSCO.report.for.GU \
--augustus --augustus_species fly --mode transcriptome --cpu 24
```

Consolidation by `cd-hist-est` dramatically reduced the number of sequences flagged as duplicates (Table S1).

**Table S1: Details of transcriptome datasets at each stage of processing**

| dataset               | contigs | BUSCO short summary                           |
|:--------------------- | -------:|:--------------------------------------------- |
| initial assembly      | 549,485 | C:97.0%[S:54.2%,D:42.8%],F:1.3%,M:1.7%,n:2510 |
| initial assembly CDS  | 85,937  |                                               |
| 90% consolidation     | 395,125 | C:96.6%[S:77.7%,D:18.9%],F:1.4%,M:2.0%,n:2510 |
| 90% consolidation CDS | 50,771  | C:94.4%[S:77.8%,D:16.6%],F:2.2%,M:3.4%,n:2510 |


## Annotation of the transcriptome

Transcript annotation was performed by EnTAP ([Hart et al. 2019](https://doi.org/10.1111/1755-0998.13106)) using the NCBI invertebrate RefSeq dataset, SwissProt and Trembl databases. EnTAP assigns sequence names and GO terms based on similarity.

```bash
EnTAP --runP --threads 24 \
-i /research/drangeli/phase2_BCG_RNAseq/Jhae_GU_trinity_assembly/GU.Trinity.c90.fa \
-d /export/groups/drangeli/entap_config/bin/invert.dmnd \
-d /export/groups/drangeli/entap_config/bin/uniprot_sprot.dmnd \
-d /export/groups/drangeli/entap_config/bin/uniprot_trembl.dmnd \
--ini /export/groups/drangeli/Jhae.genome/entap_config.ini 
```

The *J. haematoloma* sequences were matched to other Hemiptera, primarily *Halyomorpha halys*, *Cimex lectularius* and *Nilaparavata lugens* for the RefSeq comparisons and *Oncopeltus*, *Lygus* and *Riptortus* from TrEMBL. The contaminant sequence species are also as expected, coming from bacteria and yeasts. Only 12 potential contaminant sequences had homology to humans. 

## Sequencing for gene expression quantification

For gene expression analysis, RNA samples were shipped overnight on dry ice to the [DNA Technologies Core](https://dnatech.genomecenter.ucdavis.edu/) at the University of California at Davis. Libraries were prepared using the QuantSeq 3' mRNA-Seq Library Prep Kit (Lexogen, Greenland, New Hampshire, USA). Samples were sequenced in batches of up to 96 samples, producing 90-bp single-end reads. These batches were sequenced in three lanes yielding :warning: **HOW MANY???** :warning:reads. These data were inspected for quality using FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), before and after trimmed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) version 0.33, as in the example command below.

```bash
java -jar trimmomatic-0.33.jar SE \
raw.data/A1-d25-2_S161_L006_R1_001.fastq \
filtered.reads/A1-d25-2_S161_L006_R1_001.trimmed.fastq \
ILLUMINACLIP:/export/local/src/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 \
HEADCROP:16 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:79
```

> A simple [bash script](https://github.com/aphanotus/Jhae.genome/blob/main/trimming.3seq.reads.sh) was used to apply trimming to all read files and save the output in a separate folder.
>
> Reads appeared to be high quality. Even samples where RNA was initially flagged by the center as "degraded" produced 3' libraries with high quality reads. Read length and quantity were also consistent.
>
> FastQC inspection after trimming revealed notable improvements.

## Mapping 3seq reads to the transcriptome

[Kallisto](https://pachterlab.github.io/kallisto/about) was used to quantify the abundance of transcripts from 3seq data using a pseudo-alignment method ([Bray et al. 2016](https://www.nature.com/articles/nbt.3519)). First, a Kallisto index was created from the assembled transcripts, collapsed by cd-hit-est. 

> We chose to use this target, rather than the TransDecoder CDS files, because if there are fragmentary or non-coding RNAs that were sequenced, then some reads will correspond to those RNAs too.

```bash
cd /research/drangeli/phase2_BCG_RNAseq/Jhae_GU_trinity_assembly/
kallisto index --index=GU.Trinity.c90.kallisto.index GU.Trinity.c90.fa
```

In order to run Kallisto on single-end sequence, the program requires the mean and standard deviation of read lengths. The total number of reads was found using bash commands.

```bash
cd /export/groups/drangeli/morphDE/davis.3seq/filtered.reads

gunzip -c 180824/*.fastq.gz 190917/*.fastq.gz 200611/*.fastq.gz 200701/*.fastq.gz \
| sed -n -e '2~4p' \
| awk '{ print length($0); }' > read.lengths.txt

wc -l read.lengths.txt
```

Total read number was 1,149,235,725. 

Mean and  standard deviation were calculated using `awk`.

```bash
cat read.lengths.txt | awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++) { printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)} }'
```

- mean = 84.892739  bp
- SD = 0.507056 bp

A [bash script](https://github.com/aphanotus/Jhae.genome/blob/main/Jhae.kallisto.mapping.sh) was used to run `kallisto quant` for each sample and to rename the output files with the sample names. 

> The result is a folder `kallisto.vs.JhaeGUtx` that contains `.abundance.tsv` files for each sample, listing the counts for each transcript. 

The resulting output files were then merged into a single comma-separated values (CSV) file with columns for each sample.

> This file `kallisto.counts.csv` is 309 Mb, which is too large to post to GitHub.

## Differential gene expression analysis

Gene expression differences were examined using the Bioconductor package DESeq2 ([Love et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)) in R version 4.1.1 ([R Core Team 2021](https://www.R-project.org/)). The function `DESeq` fits a generalized linear model where counts are modeled using a negative binomial distribution with a gene-specific dispersion parameter. Means were estimated using approximate posterior estimation for generalized linear model ([Zhu et al. 2018](https://academic.oup.com/bioinformatics/article/35/12/2084/5159452)). We extracted log~2~ fold change and p-values corrected by independent hypothesis weighing ([Ignatiadis et al. 2016](https://www.nature.com/articles/nmeth.3885)) and applied a critical threshold of 0.05. Differential gene expression was visualized by customized plots generated using ggplot2. Variance stabilizing transformation, using the function `vst`,  was applied before ordination using principle component analysis.

> For more details see the R script [`dge.analysis.R`](https://github.com/aphanotus/morphDE/blob/main/dge.analysis.R)

