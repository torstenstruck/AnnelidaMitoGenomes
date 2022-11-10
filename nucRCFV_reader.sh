#!/bin/sh

## Generate the RCFV values of the aligned nuc sequences:
for FILE in *.fasta
do
 perl NuclRCFVReader.pl ${FILE} ${FILE}_Results
done

#### submit job as "sh this_script_name" ####
