#!/bin/sh
#(c) Torsten Struck 2019

#To run this script the following files are needed in the following folders:
#the sequence of each whole mitochondrial genomes as a fasta file with the extension ".fasta" in the folder "./01_Data/MitochondrialGenomes/Used/"
#from the MITOS2 analyses the ".fas" and "*.txt"-files containing the nucleotide sequences of each gene and gene order of each whole mitochondrial genomes in a folder for with the species names (e.g., "Terebratulina_retusa". All these folders need to be in the folder "./02_Annotation/Used/"

base=$(pwd)
today=$(date '+%Y_%m_%d_%H_%M')
logfile="$base/log_${today}.txt"
mkdir $base/03_MitochondrialProperties/
mkdir $base/03_MitochondrialProperties/WholeGenome/
mkdir $base/03_MitochondrialProperties/CodingGenes/
mkdir $base/03_MitochondrialProperties/StructuralInformation/

#compile whole mitochondrial genomes
echo "Compiling complete mitochondrial genome sequences into one file"
echo "Compiling complete mitochondrial genome sequences into one file" > $logfile
cd $base/01_Data/MitochondrialGenomes/Used
cat *.fasta | sed "s/_[a-zA-Z][a-zA-Z][0-9][0-9][0-9][0-9][0-9][0-9]$//" | sed "s/_[a-zA-Z][a-zA-Z][0-9][0-9][0-9][0-9]$//" | sed "s/_NC_[0-9][0-9][0-9][0-9][0-9][0-9]$//" | sed "s/_[a-zA-Z][a-zA-Z][0-9][0-9][0-9][0-9]$//" > $base/03_MitochondrialProperties/WholeGenome/All_MitochondrialGenomes.fas
echo "Done\n\n"
echo "Done\n\n" >> $logfile

#compile all mitochondrial gene names
echo "Compiling mitochondrial gene names"
echo "Compiling mitochondrial gene names" >> $logfile
cd $base/02_Annotation/Used/
count1=0
#for loop to visit each directory
for dir in *; do
  if [ -d "$dir" ]; then
#count the number of visited directories as a control number
    count1=$((count1 + 1))
    echo "$count1	$dir"
    echo "$count1	$dir" >> $logfile
    cd $dir
#retrieve name of fasta file with sequence information and assign new names while exchanging windows line breaks for unix ones
    filenamenuc=(*.fas)
    newfilenuc="$filenamenuc.fasta"
    tr -d '\015' < $filenamenuc > $newfilenuc
    echo ">" >> $newfilenuc
#retrieve the sequence names in each file
    grep ">\w" $newfilenuc | cut -f4 -d ";" | tr -d " " >> ../GeneNamesNuc.txt
    cd ..
  fi
done
sort GeneNamesNuc.txt | uniq | sed '/^$/d' > $base/03_MitochondrialProperties/CodingGenes/GeneNamesAllNuc.txt
rm GeneNamesNuc.txt
echo "Done\n\n"
echo "Done\n\n" >> $logfile

#compile all nuclear mitochondrial genes in separate files
echo "Compiling nucleotide sequences of mitochondrial genes into separate files"
echo "Compiling nucleotide sequences of mitochondrial genes into separate files" >> $logfile
cd $base/02_Annotation/Used/
while read gene
do
  echo "Working on $gene"
  echo "\n$gene\nNumber_species\tSpecies_file\tNumber_sequences" >> $logfile
  outfile="$base/03_MitochondrialProperties/CodingGenes/$gene.fas"
  count2=0
#for loop to visit each directory
  for dir in *; do
    if [ -d "$dir" ]; then
#count the number of visited directories as a control number
      count2=$((count2 + 1))
      cd $dir
      infile=(*.fas.fasta)
#retrieve the sequences of the gene
      sed -n "/>*$gene$/,/>/p" < $infile | sed '$d' | sed "s/^>.*$/>$dir/" >> $outfile
      countseq=$(grep '>' $outfile | wc -l) 
      echo "$count2\t$infile\t$countseq" >> $logfile
      cd ..
    fi
  done
done < $base/03_MitochondrialProperties/CodingGenes/GeneNamesAllNuc.txt
echo "Done\n\n"
echo "Done\n\n" >> $logfile

#compile structural information about mitochondrial genomes into one file
echo "Compiling structural information about mitochondrial genomes"
echo "Compiling structural information about mitochondrial genomes" >> $logfile
echo "\nNumber_species\tSpecies_file\tNumber_sequences" >> $logfile
cd $base/02_Annotation/Used/
count4=0
fasoutfile="$base/03_MitochondrialProperties/StructuralInformation/GeneOrder.fas"
tableoutfile="$base/03_MitochondrialProperties/StructuralInformation/GeneOrder.txt"
echo "Species	GeneOrder" > $tableoutfile
#for loop to visit each directory
for dir in *; do
  if [ -d "$dir" ]; then
#count the number of visited directories as a control number
    count4=$((count4 + 1))
    cd $dir
#retrieve name of fasta file with sequence information and assign new names while exchanging windows line breaks for unix ones
    filename=($dir.txt)
    newfile=($dir.text)
    tr -d '\015' < $filename > $newfile
#retrieve the sequences of the gene
    cat $newfile | sed '/^$/d' | sed "s/^>.*$/>$dir/" >> $fasoutfile
    countseq2=$(grep '>' $fasoutfile | wc -l)
    echo "$count4\t$newfile\t$countseq2" >> $logfile
    cd ..
  fi
done
#converting the fasta file format into a tab-delimited table format
sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' < $fasoutfile | sed "s/ /	/g" | tr '>' '\n' | sed '/^$/d' >> $tableoutfile
echo "Done\n\n"
echo "Done\n\n" >> $logfile
