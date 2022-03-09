#!/bin/bash

# Designed to be run from a folder containing filtered, unpaired (e.g. 3seq) reads

NUMTHREADS=12
EXT=".fastq.gz"
FILENAMES=(*$EXT)
BASENAMES=( "${FILENAMES[@]/$EXT}" )

INDEX="/research/drangeli/phase2_BCG_RNAseq/Jhae_GU_trinity_assembly/GU.Trinity.c90.kallisto.index"
OUTPUTFOLDER="/export/groups/drangeli/Jhae.morph.txome/davis.3seq/kallisto.vs.JhaeGUtx"

echo
echo "Mapping ${#BASENAMES[@]} files."
echo

for i in "${BASENAMES[@]}"; do
  SAMPLENAME=${i%_L00?_R1_001}
  echo "Mapping $SAMPLENAME"
  kallisto quant --index=$INDEX --output-dir=./ --single --single-overhang -l 84.892739 -s 0.507056 --threads=$NUMTHREADS $i$EXT
  mv abundance.tsv $OUTPUTFOLDER/$SAMPLENAME.abundance.tsv
  mv abundance.h5 $OUTPUTFOLDER/$SAMPLENAME.abundance.h5
  mv run_info.json $OUTPUTFOLDER/$SAMPLENAME.run_info.json
done

echo
echo "Done mapping ${#BASENAMES[@]} files."
echo
