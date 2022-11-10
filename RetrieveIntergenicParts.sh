#!/bin/sh

#To run this script the following files are needed in the following folders:
#the sequence of each whole mitochondrial genomes as a fasta file with the extension ".fasta" in the folder "./01_Data/MitochondrialGenomes/Used/"
#from the MITOS2 analyses the ".bed"-file of each whole mitochondrial genome in a folder for with the species names (e.g., "Terebratulina_retusa". All these folders need to be in the folder "./02_Annotation/Used/"
#additionally, a text file "ListSpeciesSize.txt" is needed that lists per species the total number of position of each mitochondrial genome separated by a tab

#generate header for output file
echo "Species	Number	Average" > Positive_Intergenic_Number_Average.txt

#read in species names and the size of the mitochondrial genome from a file
while IFS=$'\t' read -r SPECIES SIZE
do
 echo "Working on ${SPECIES}"
#assigne variable for last position of a gene, number of positive intergenic regions, count of of their bases, and count of Ns in fasta-file
 STOP=0
 COUNT=0
 NUMBER=0
 NS=0
#assigne bed-file for retrieving intergenic regions and read in the lines, assign first column to first, second to start, third to third and rest to rest
 BEDFILE=$(ls ./02_Annotation/Used/${SPECIES}/${SPECIES}.bed)
 while IFS=$'\t' read -r FIRST START THIRD REST
 do
#calculate the size of the intergenic region by substrating the stop position of the previous gene from start position of the present one minus 1
  DIFF=$((START-STOP))
  DIFF=$((DIFF-1))
#count only positive intergenic regions, which have at least one position, calculate the total number of bases and the number of such intergenic regions
  if [ $DIFF -gt 0 ]; then
   COUNT=$((COUNT+DIFF))
   NUMBER=$((NUMBER+1))
 fi
#assign the stop position of the present gene, to be available for the next line
 STOP=$THIRD
 done < ${BEDFILE}
#calculate the size of the last intergenic region if present and add only if present
 DIFF=$((SIZE-STOP))
 if [ $DIFF -gt 0 ]; then
  COUNT=$((COUNT+DIFF))
  NUMBER=$((NUMBER+1))
 fi
#assigne fasta-file to count the number of Ns included to put together incomplete genomes, then count the numbers of Ns in the sequence lines 
 ALIFILE=$(ls ./01_Data/MitochondrialGenomes/Used/${SPECIES}*)
 while read -r LINE
 do
  if [[ $LINE != ">"* ]]; then
   NUMBERN=$(echo $LINE| grep -o N | wc -l)
   NS=$((NS+NUMBERN))
  fi
 done < ${ALIFILE}
#substrat the numbers of Ns from the count of intergenic positions and calculate the average size of the intergenic positions
 COUNT=$((COUNT-NS))
 AVERAGE=$((COUNT/NUMBER))
#export species name, number and average size of postive intergenic regions to output file
 echo "${SPECIES}	${NUMBER}	${AVERAGE}" >> Positive_Intergenic_Number_Average.txt
done < ListSpeciesSize.txt