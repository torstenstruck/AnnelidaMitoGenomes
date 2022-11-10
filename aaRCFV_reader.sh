#!/bin/sh

## Generate the RCFV values of the aligned AA sequences:
for FILE in *.fasta
do
 perl ProtRCFVReader.pl ${FILE} ${FILE}_Results
done

#### submit job as "sh this_script_name" ####
