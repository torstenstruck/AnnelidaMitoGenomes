#!/bin/sh

#read in species names and the size of the mitochondrial genome from a file
while IFS=$'\t' read -r SPECIES SIZE
do
 echo "$SPECIES"
#assign bed-file to a variable
 BEDFILE=$(ls ./02_Annotation/Used/${SPECIES}/${SPECIES}.bed)
# echo "$BEDFILE"
 while read -r LINE
 do
  TESTLINE=$LINE
 done < ${BEDFILE}
#assign alignment-file to a variable
 ALIFILE=$(ls ./01_Data/MitochondrialGenomes/Used/${SPECIES}*)
# echo "$ALIFILE"
 while read -r LINE
 do
  TESTLINE=$LINE
 done < ${ALIFILE}
done < ListSpeciesSize.txt