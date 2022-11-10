#!/bin/sh
#the following files are needed:
#Grouped_factors.txt
#different files ending on *.group

COUNT1=0
COUNT2=0

for FILE in *.group; 
do
 COUNT1=$((COUNT1+1))
 while read LINE
 do
  echo "Working on $LINE"
  COUNT2=$((COUNT2+1))
  grep "$LINE	" Grouped_factors.txt  >> $FILE.corr
 done < $FILE
 tr -d '\n' < $FILE.corr > $FILE.corr.one
 echo >> $FILE.corr.one
done
echo "Number of files: $COUNT1"
echo "Number of lines: $COUNT2"

cat *.one > Grouped_Cluster_Corr.txt

rm *.corr *.one
