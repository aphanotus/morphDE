#!/bin/bash

# Designed to be run from a folder containing read data

EXT=".fastq.gz"
FILENAMES=(*$EXT)
FILENUMBER=${#FILENAMES[@]}
OUTPUTFILENAME="read.counts.tsv"

echo
echo "Read counts from $FILENUMBER files will be recorded in $OUTPUTFILENAME"
READLENGTH=$(zcat ${FILENAMES[1]} | head -n 2 | tail -n 1 | awk '{print length}')
echo "Reads are $READLENGTH bp"
echo

touch $OUTPUTFILENAME
echo "" > $OUTPUTFILENAME

TOTALREADS=0
ITER=1

for i in "${FILENAMES[@]}"; do
  printf "$ITER of $FILENUMBER \t$i:\t"
  NREADS=$(zcat $i | grep -c "^@")
  printf "$NREADS reads\n"
  printf "$i\t$NREADS\n" >> $OUTPUTFILENAME
  TOTALREADS=$((TOTALREADS + NREADS))
  ITER=$((ITER + 1))
done

printf "total\t$TOTALREADS\n" >> $OUTPUTFILENAME

echo
echo "Read counts from $FILENUMBER files have been recorded in $OUTPUTFILENAME"
echo "Total read number: $TOTALREADS"
echo "Reads are $READLENGTH bp"
echo "Total data: $((READLENGTH * TOTALREADS / 1000000)) million bp"
echo
