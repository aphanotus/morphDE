# Supplemental Methods

for Angelini, Steele and O'Brien. **Gene expression during development of a dispersal polyphenism identifies regulators and effectors of growth in muscle and germline**. *Integrative & Comparative Biology*. 

> :warning: **Pre-submission DRAFT**.  Some of this text may be moved into the Methods section of the manuscript's main text. Details that may be important now, but should likely be removed before publication are given as block quotes.

> ℹ️ While any process is running on nscc, you can monitor it's progress indirectly by checking the CPU usage on [node 26](http://schupflab.colby.edu/mrtg/nscc/nscc-n26.load.html) or [n28](http://schupflab.colby.edu/mrtg/nscc/nscc-n28.load.html).

In order to characterize the development of polyphenic traits in the wings, flight muscles and gonads of the soapberry bug, *Jadera haematoloma*, we sequenced transcriptomes from a multidimensional matrix of life stage, tissue, sex, food regime and host population.

## Sampling design

Transcriptome studies proceeded through three phases. Initially, we sampled whole bodies from 12 individuals, including three biological replicates of each sex and morph combination from the Plantation Key, FL population (phase 1). Next we sampled dorsal thorax and gonads from three individuals of each sex and morph combination from the Plantation Key, FL and Aurora, CO populations (phase 2). Samples from these phases were prepared as full-length 150-bp paired-end Illumina libraries and used to assembly a reference transcriptome. 

The final phase utilized 3'-tag sequencing with the exclusive goal of quantifying gene expression. The sampling design at this phase included two stages (fifth instars and nascent adults), two tissues (dorsal thorax and gonads), two sexes, long- and short-wing adult morphs, high and low food regimes, four populations (collected from Key Largo, FL, Plantation Key, FL, Aurora, CO, and Frederick, MD), and three biological replicates of each unique combination of the preceding factors. Thus the full matrix contained 288 samples. Seven samples were lost or did not pass library quality controls, resulting in 281 samples in the analysis.

## Nucleic acid extractions

Individuals were collected from laboratory cultures, anesthetized using CO<sub>2</sub> exposure, photographed, and flash frozen in liquid nitrogen. Samples were stored at -80˚C until further processing. The thorax was removed from head and abdomen with a clean scalpel blade, and the dorsal thorax was separated above the legs while tissue was still frozen. Gonads were then dissected from the abdomen of adults. For juveniles, the entire abdomen was included. Nucleic acid extractions used the Invitrogen PureLink RNA extraction kit.

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

**Table S1**: Details of transcriptome datasets at each stage of processing.

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

For gene expression analysis, RNA samples were shipped overnight on dry ice to the [DNA Technologies Core](https://dnatech.genomecenter.ucdavis.edu/) at the University of California at Davis. Libraries were prepared using the QuantSeq 3' mRNA-Seq Library Prep Kit (Lexogen, Greenland, New Hampshire, USA) in the three batches for sequencing in separate Illumina HiSeq lanes. Batch 3 was sequenced twice and reads counts for each run were added together for each sample after filtering.

These batches were sequenced in three lanes yielding 115.9×10^9^ bp. These data

Short read sequences were inspected for quality using FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), before and after trimmed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) version 0.33, as in the example command below.

```bash
java -jar trimmomatic-0.33.jar SE \
raw.data/A1-d25-2_S161_L006_R1_001.fastq \
filtered.reads/A1-d25-2_S161_L006_R1_001.trimmed.fastq \
ILLUMINACLIP:/export/local/src/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 \
HEADCROP:16 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:79
```

> A simple [bash script](https://github.com/aphanotus/Jhae.genome/blob/main/trimming.3seq.reads.sh) was used to apply trimming to all read files and save the output in a separate folder.FastQC inspection after trimming revealed notable improvements.
>

**Table S2**: Read numbers before and after quality filtering.

| batch   |     raw reads | filtered reads | passing filter |
| :------ | ------------: | -------------: | -------------- |
| 1       |   341,370,316 |    323,188,378 | 94.7%          |
| 2       |   321,564,969 |    313,590,249 | 97.5%          |
| 3       |   484,588,962 |    419,931,564 | 86.7%          |
| overall | 1,147,524,247 |  1,056,710,191 | 92.1%          |

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

Gene expression differences were examined using the Bioconductor package DESeq2 version 1.32.0 ([Love et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)) in R version 4.1.1 ([R Core Team 2021](https://www.R-project.org/)). The function `DESeq` fits a generalized linear model where counts are modeled using a negative binomial distribution with a gene-specific dispersion parameter. Means were estimated using approximate posterior estimation for generalized linear models ([Zhu et al. 2018](https://academic.oup.com/bioinformatics/article/35/12/2084/5159452)). We extracted log<sub>2</sub> fold-change and p-values corrected by independent hypothesis weighing ([Ignatiadis et al. 2016](https://www.nature.com/articles/nmeth.3885)) and applied a critical threshold of 0.05. Contrasts were made for sex, wing morph and food regime, after accounting for batch effects (multiple sequencing runs). Data were analyzed separately for each stage (fifth instar juveniles and nascent adults) and for each tissue (thorax and gonads). Differential gene expression was visualized by customized plots generated using ggplot2. Variance stabilizing transformation, using the function `vst`,  was applied before ordination using principle component analysis.

> For more details see the R script [`analysis.dge.R`](https://github.com/aphanotus/morphDE/blob/main/analysis.dge.R)

Batch effects were introduced by the need to sequence samples in multiple lanes (Figure S1).

![](https://i.imgur.com/MAjpAP7.jpg)

**Figure S1**. Principle components of uncorrected counts, color-coded by biological factors (A) and sequencing batch (B). 

DESeq2 allows variance from confounding factors, such as sequencing batch, to be accounted for through generalized linear modeling. The normalized counts that result from the model then reflect biologically meaningful variance (Figure S2).

<img src="https://i.imgur.com/XFd0Alf.jpg =300x" style="zoom:30%;" />

**Figure S2**. Principle components of transcript counts for each sample, following normalization for sequencing run (batch). Samples are color-coded by biological origin. 

More complex models can be applied to examine specific contrasts of biological interest. These models partition the variance in gene expression sequentially to factors in the model. For more complex models, less power exists to detect differences, but it allows for the influence of multiple factors to be disentangle from one another. 


![](https://i.imgur.com/098suDV.jpg)

**Figure S3**. Volcano plots, showing log<sub>2</sub> fold-change (LFC) and -log<sub>10</sub> FDR-adjusted p-values for each dataset and contrast. Genes are highlighted in red if they have an absolute LFC values > 1 (two-fold difference) and an FDR-adjusted p-value < 10<sup>-3</sup> (panels A-C,F-G) or 0.05 (panels D-E,H-J). These thresholds are indicated by vertical and horizontal lines. Highlighted transcript are annotated where possible. 

> **Results**: Notice that fewer transcripts show differential expression in the juvenile datasets compared adults (e.g. A vs. B and F vs. G). Sex is a much more dominant factor in gene expression in the gonad (G), compared to the thorax (B). Many genes, especially those related to muscle and energy metabolism, are up-regulated in the thorax of long-winged adults (C). However the gonads show relatively few expression differences between morphs, after accounting for the significant differences imposed by sex (H). After modeling other factors, food regime has a comparatively small influence on gene expression (D,E,I,J).

## Enrichment analysis

Enrichment analyses were performed based on annotations with terms from the GO and KEGG ontologies, as well as EggNOG protein domains, applied to transcripts by EnTAP. We also tested for enrichment of xenic or contaminant transcripts among those that were differentially expressed, and neither of these categories was significantly over-represented. For each term, p-values from the hypergeometric test, implemented in the base R function `phyper`, were adjusted for false discovery rates (FDR).  

> For more details on the enrichment analysis see the R script [`enrichment.R`](https://github.com/aphanotus/morphDE/blob/main/enrichment.R)

**Table S3**: Terms from GO and KEGG ontologies and EggNOG protein domains over-represented in differentially expressed genes in each dataset and contrast.

| tissue | stage | contrast | ontology | enriched terms | top terms                                                    |
| ------ | ----- | -------- | -------- | -------------- | ------------------------------------------------------------ |
| thorax | L5    | by sex   | EggNOG   | 15             | dsx/mab-3 DNA-binding domain, 9.13x10<sup>-7</sup>; Complementary sex determiner-like SR domain, 0.0180; proline-rich protein BAT2 N-terminal domain, 0.0180; CD34/Podocalyxin family, 0.0180; DFDF domain, 0.0180; Sex determination protein N terminal domain, 0.0180; Barentsz/CASC3 domain, 0.0215; Meckelin, 0.0215; CDC48, 0.0298 |
| thorax | L5    | by sex   | GO       | 0              |                                                              |
| thorax | L5    | by sex   | KEGG     | 0              |                                                              |
| thorax | adult | by sex   | EggNOG   | 26             | Reverse transcriptase, 2.11x10<sup>-9</sup>; Endonuclease/Exonuclease/phosphatase family, 1.38x10<sup>-8</sup>; Rnase H, 7.70x10<sup>-7</sup>; ZnF C2HC, 5.78x10<sup>-6</sup>; Lipoprotein N-terminal domain, 6.48x10<sup>-5</sup>; Vitellogenin-like N-terminal domain, 6.48x10<sup>-5</sup>; DUF1943, 2.54x10<sup>-4</sup>; Co/Zn superoxide dismutase, 8.63x10<sup>-4</sup>; CUB domain, 1.77x10<sup>-3</sup>; Trypsin, 3.17x10<sup>^</sup>-3 |
| thorax | adult | by sex   | GO       | 267            | GO:0022626(L=5) cytosolic ribosome, 2.98×10^-41; GO:0006412(L=5) translation, 1.07×10<sup>-29</sup>; GO:0003723(L=4) RNA binding, 1.08×10<sup>-27</sup>; GO:0034645(L=4) cellular macromolecule biosynthetic process, 4.73×10<sup>-27</sup>; GO:0009059(L=4) macromolecule biosynthetic process, 1.16×10<sup>-26</sup>; GO:0000022(L=5) mitotic spindle elongation, 4.66×10<sup>-26</sup>; GO:0051231(L=4) spindle elongation, 3.81×10<sup>-25</sup>; GO:0044391(L=5) ribosomal subunit, 5.19×10<sup>-25</sup>; GO:0022625(L=6) cytosolic large ribosomal subunit, 9.80×10<sup>-24</sup>; GO:0044249(L=3) cellular biosynthetic process, 3.9×10<sup>-23</sup> |
| thorax | adult | by sex   | KEGG     | 3              | map04974 Protein digestion and absorption, 6.71x10<sup>-6</sup>; map04972 Pancreatic secretion, 1.04x10<sup>-4</sup>; map00513 Various types of N-glycan biosynthesis,7.60x10<sup>-3</sup> |
| gonad  | L5    | by sex   | EggNOG   | 46             | COIL, 1.62x10<sup>-69</sup>; S/T kinases catalytic domain, 5.50x10<sup>-9</sup>; TRANS, 3.37x10<sup>-8</sup>; WD40 domain, 3.31x10<sup>-7</sup>; RNA recognition motif, 5.60x10<sup>-7</sup>; protein kinase, 1.14x10<sup>-6</sup>; Protein tyrosine kinase, 1.14x10<sup>-6</sup>; RNA recognition motif, 1.28x10<sup>-6</sup>; M17 cytosol aminopeptidase, 2.15x10<sup>-5</sup>; Leucin-rich repeats, 4.30x10<sup>-5</sup> |
| gonad  | L5    | by sex   | GO       | 1299           | GO:0005622(L=3) intracellular anatomical structure, 1.40×10<sup>-94</sup>; GO:0043229(L=3) intracellular organelle, 7.32×10<sup>-84</sup>; GO:0043231(L=4) intracellular membrane-bounded organelle, 1.47×10<sup>-64</sup>; GO:0043227(L=2) membrane-bounded organelle, 5.40×10<sup>-64</sup>; GO:0005737(L=4) cytoplasm, 2.40×10<sup>-55</sup>; GO:0044237(L=2) cellular metabolic process, 7.05×10<sup>-55</sup>; GO:0044238(L=2) primary metabolic process, 1.08×10<sup>-52</sup>; GO:0016043(L=2) cellular component organization, 3.73×1<sup>-51</sup>; GO:0071704(L=2) organic substance metabolic process, 6.15×10<sup>-50</sup>; GO:0043228(L=2) non-membrane-bounded organelle, 1.66×10<sup>-44</sup> |
| gonad  | L5    | by sex   | KEGG     | 31             | map01100 Metabolic pathways, 4.58×10<sup>-9</sup>; map05012 Parkinson disease, 1.22×10<sup>-7</sup>; map00710 Carbon fixation in photosynthetic organisms, 1.83×10<sup>-6</sup>; map05010 Alzheimer disease, 3.77×10<sup>-6</sup>; map00051 Fructose and mannose metabolism, 9.74×10<sup>-5</sup>; map01110 Biosynthesis of secondary metabolites, 1.04×10<sup>-4</sup>; map00030 Pentose phosphate pathway, 9.86×10<sup>-4</sup>; map04340 Hedgehog signaling pathway, 9.86×10<sup>-4</sup>; map05152 Tuberculosis, 9.86×10<sup>-4</sup>; map05110 Vibrio cholerae infection, 1.13×10<sup>-3</sup> |
| gonad  | adult | by sex   | EggNOG   | 0              |                                                              |
| gonad  | adult | by sex   | GO       | 0              |                                                              |
| gonad  | adult | by sex   | KEGG     | 0              |                                                              |
| thorax | adult | by morph | EggNOG   | 48             | COIL, 4.72×10<sup>-41</sup>; TRANS, 5.43×10<sup>-40</sup>; SIGNAL, 4.54×10<sup>-19</sup>; Immunoglobulin I-set domain, 2.84×10<sup>-5</sup>; Immunoglobulin C-2 Type, 3.53×10<sup>-5</sup>; Immunoglobulin V-set domain, 3.77×10<sup>-5</sup>; CUB domain, 4.06×10<sup>-5</sup>; Immunoglobulin domain, 1.15×10<sup>-4</sup>; EF-hand calcium binding motif, 9.41×10<sup>-4</sup>; Allatostatin, 2.30×10<sup>-3</sup> |
| thorax | adult | by morph | GO       | 1251           | GO:0005737(L=4) cytoplasm, 1.67×10<sup>-124</sup>; GO:0005622(L=3) intracellular anatomical structure, 1.68×10<sup>-118</sup>; GO:0005739(L=5) mitochondrion, 1.4×10<sup>-100</sup>; GO:0043229(L=3) intracellular organelle, 4.93×10<sup>-98</sup>; GO:0043227(L=2) membrane-bounded organelle, 4.05×10<sup>-87</sup>; GO:0043231(L=4) intracellular membrane-bounded organelle, 7.23×10<sup>-87</sup>; GO:0044237(L=2) cellular metabolic process, 1.30×10<sup>-68</sup>; GO:0071704(L=2) organic substance metabolic process, 2.37×10<sup>-56</sup>; GO:0044238(L=2) primary metabolic process, 2.97×10<sup>-56</sup>; GO:0005759(L=6) mitochondrial matrix, 6.3×10<sup>-50</sup> |
| thorax | adult | by morph | KEGG     | 16             | map01100 Metabolic pathways, 5.266×10<sup>-32</sup>; map05012 Parkinson disease, 3.35×10<sup>-26</sup>; map05010 Alzheimer disease, 5.08×10<sup>-19</sup>; map01110 Biosynthesis of secondary metabolites, 7.45×10<sup>-13</sup>; map04974 Protein digestion and absorption, 7.02×10<sup>-6</sup>; map00640 Propanoate metabolism, 1.16×10<sup>-3</sup>; map00513 Various types of N-glycan biosynthesis, 1.21×10<sup>-3</sup>; map01120 Microbial metabolism in diverse environments, 2.22×10<sup>-3</sup>; map04260 Cardiac muscle contraction, 2.92×10<sup>-3</sup>; map00020 Citrate cycle (TCA cycle) 4.59×10<sup>-3</sup> |
| gonad  | adult | by morph | EggNOG   | 11             | COIL, 0.00155; PI3-kinase C2 domain, 0.00155; PI3K p85-binding domain, 0.00155; PI3K ras-binding domain, 0.00155; PI3K accessory domain, 0.00155; Cyclin C, 0.00183; Guanylate cyclase catalytic domain, 0.00183; IPT domain, 0.00183; Phosphatidylinositol 3/4-kinase, 0.00183; PI3K catalytic domain, 0.00183 |
| gonad  | adult | by morph | GO       | 370            | GO:0001705(L=5) ectoderm formation, 0.009967; GO:0001712(L=6) ectodermal cell fate commitment, 0.009967; GO:0001715(L=7) ectodermal cell fate specification, 0.009967; GO:0001727(L=5) lipid kinase activity, 0.009967; GO:0004016(L=3) adenylate cyclase activity, 0.009967; GO:0005942(L=3) phosphatidylinositol 3-kinase complex, 0.009967; GO:0005943(L=5) phosphatidylinositol 3-kinase complex, class IA, 0.009967; GO:0006171(L=8) cAMP biosynthetic process, 0.009967; GO:0007154(L=2) cell communication, 0.009967; GO:0007165(L=3) signal transduction, 0.009967 |
| gonad  | adult | by morph | KEGG     | 71             | map04062 Chemokine signaling, 0.000229; map04725 Cholinergic synapse, 0.000229; map04914 Progesterone-mediated oocyte maturation, 0.000229; map04915 Estrogen signaling, 0.000319; map05166 Human T-cell leukemia virus 1 infection, 0.000348; map04020 Calcium signaling, 0.00627; map04070 Phosphatidylinositol signaling, 0.00627; map04150 mTOR signaling, 0.00627; map04210 Apoptosis, 0.00627; map04370 VEGF signaling, 0.00627; map04620 Toll-like receptor signaling, 0.00627; map04630 JAK-STAT signaling, 0.00627; map04911 Insulin secretion, 0.00627; map04913 Ovarian steroidogenesis, 0.00627 |
| thorax | L5    | by food  | EggNOG   | 7              | tetrapeptide XPTX repeats, 0.00144; catecholamine-binding domain, 0.00162; dopamine beta-monooxygenase N-terminal-like domain, 0.00162; SEA domain, 0.00162; EGF2, 0.00360; EGF-like, 0.00960; EGF, 0.00977 |
| thorax | L5    | by food  | GO       | 2              | GO:0016052(L=4) carbohydrate catabolic process, 0.03984; GO:0030246(L=2) carbohydrate binding, 0.03984; GO:0005975(L=3) carbohydrate metabolic process, 0.09147; GO:0009056(L=2) catabolic process, 0.1903; GO:1901575(L=3) organic substance catabolic process, 0.1903; GO:0044238(L=2) primary metabolic process, 0.4692; GO:0071704(L=2) organic substance metabolic process, 0.4692 |
| thorax | L5    | by food  | KEGG     | 0              |                                                              |
| thorax | adult | by food  | EggNOG   | 5              | DUF1943 domain, 0.01779; Lipoprotein N-terminal domain, 0.01779; Co/Zn superoxide dismutase, 0.01779; Vitellogenin N-terminal domain, 0.01779; von Willebrand factor type D domain, 0.01930 |
| thorax | adult | by food  | GO       | 0              |                                                              |
| thorax | adult | by food  | KEGG     | 0              |                                                              |
| gonad  | L5    | by food  | EggNOG   | 4              | DAN domain, 0.00911; Peroxisome biogenesis factor 1, 0.00911; Trypsin-like serine protease, 0.0119; Trypsin, 0.0119 |
| gonad  | L5    | by food  | GO       | 0              |                                                              |
| gonad  | L5    | by food  | KEGG     | 0              |                                                              |
| gonad  | adult | by food  | EggNOG   | 8              | CUB domain, 0.0256; DMPK coiled coil, 0.0256; Trypsin-like serine protease, 0.0256; Trypsin, 0.0256; microtubule-associated protein 7 family, 0.0308; P21-Rho-binding domain, 0.0308; CNH domain, 0.0330; Sodium solute symporter, 0.0375 |
| gonad  | adult | by food  | GO       | 0              |                                                              |
| gonad  | adult | by food  | KEGG     | 2              | map04974 Protein digestion and absorption, 0.00304; map04972 Pancreatic secretion, 0.00420 |

