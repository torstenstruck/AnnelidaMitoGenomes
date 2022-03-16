#!/bin/sh

declare -i FOUND
echo "Pairs not found:" > NotFoundPairs.txt

while read -r LINE
do
 echo "$LINE" > ${LINE}.txt
done < Correlation_HighlyFactors.txt

while IFS=$'\t' read -r not1 V1 V2 not2
do
 echo "Working on: $V1	$V2"
 FOUND=0
 while read -r LINE
 do
  if [ $LINE = $V1 ] && [ $FOUND -eq 0 ]; then
   echo "$V1" >> ${LINE}.txt
   echo "$V2" >> ${LINE}.txt
   FOUND=1; 
   break 
  elif [ $LINE = $V2 ] && [ $FOUND -eq 0 ]; then
   echo "$V1" >> ${LINE}.txt
   echo "$V2" >> ${LINE}.txt
   FOUND=1; 
   break 
  elif [ $FOUND -eq 0 ]; then
   sort ${LINE}.txt | dos2unix | uniq > tmp.txt 
   while read -r ADDLINE
   do
	if [ "$LINE" != "$ADDLINE" ] && [ $FOUND -eq 0 ]; then
     if [ $ADDLINE = $V1 ] && [ $FOUND -eq 0 ]; then
      echo "$V1" >> ${LINE}.txt
      echo "$V2" >> ${LINE}.txt
      FOUND=1; 
	  break 2
     elif [ $ADDLINE = $V2 ] && [ $FOUND -eq 0 ]; then
      echo "$V1" >> ${LINE}.txt
      echo "$V2" >> ${LINE}.txt
      FOUND=1; 
	  break 2
	 fi
	fi
   done < tmp.txt
  fi
 done < Correlation_HighlyFactors.txt
 if [ $FOUND -eq 0 ]; then
  echo "$V1	$V2" >> NotFoundPairs.txt
 fi
done < Correlation_HighPairs.txt

rm -f tmp.txt 
sed -i '1d' NotFoundPairs.txt
echo "Pairs not found second time:" > NotFoundPairs2.txt

while IFS=$'\t' read -r V1 V2
do
 echo "Working again on: $V1	$V2"
 FOUND=0
 while read -r LINE
 do
   sort ${LINE}.txt | dos2unix | uniq > tmp.txt 
   while read -r ADDLINE
   do
	if [ "$LINE" != "$ADDLINE" ] && [ $FOUND -eq 0 ]; then
     if [ $ADDLINE = $V1 ] && [ $FOUND -eq 0 ]; then
      echo "$V1" >> ${LINE}.txt
      echo "$V2" >> ${LINE}.txt
      FOUND=1;
	  break 2
     elif [ $ADDLINE = $V2 ] && [ $FOUND -eq 0 ]; then
      echo "$V1" >> ${LINE}.txt
      echo "$V2" >> ${LINE}.txt
      FOUND=1;
	  break 2
	 fi
	fi
   done < tmp.txt
 done < Correlation_HighlyFactors.txt
 if [ $FOUND -eq 0 ]; then
  echo "$V1	$V2" >> NotFoundPairs2.txt
 fi
done < NotFoundPairs.txt

rm -f tmp.txt #NotFoundPairs.txt
sed -i '1d' NotFoundPairs2.txt

echo "Pairs not found third time:" > NotFoundPairs3.txt

while IFS=$'\t' read -r V1 V2
do
 echo "Working third time on: $V1	$V2"
 FOUND=0
 while read -r LINE
 do
   sort ${LINE}.txt | dos2unix | uniq > tmp.txt 
   while read -r ADDLINE
   do
	if [ "$LINE" != "$ADDLINE" ] && [ $FOUND -eq 0 ]; then
     if [ $ADDLINE = $V1 ] && [ $FOUND -eq 0 ]; then
      echo "$V1" >> ${LINE}.txt
      echo "$V2" >> ${LINE}.txt
      FOUND=1;
	  break 2
     elif [ $ADDLINE = $V2 ] && [ $FOUND -eq 0 ]; then
      echo "$V1" >> ${LINE}.txt
      echo "$V2" >> ${LINE}.txt
      FOUND=1;
	  break 2
	 fi
	fi
   done < tmp.txt
 done < Correlation_HighlyFactors.txt
 if [ $FOUND -eq 0 ]; then
  echo "$V1	$V2" >> NotFoundPairs3.txt
 fi
done < NotFoundPairs2.txt

rm -f tmp.txt #NotFoundPairs2.txt
sed -i '1d' NotFoundPairs3.txt

echo "Pairs not found fourth time:" > NotFoundPairs4.txt

while IFS=$'\t' read -r V1 V2
do
 echo "Working fourth time on: $V1	$V2"
 FOUND=0
 while read -r LINE
 do
   sort ${LINE}.txt | dos2unix | uniq > tmp.txt 
   while read -r ADDLINE
   do
	if [ "$LINE" != "$ADDLINE" ] && [ $FOUND -eq 0 ]; then
     if [ $ADDLINE = $V1 ] && [ $FOUND -eq 0 ]; then
      echo "$V1" >> ${LINE}.txt
      echo "$V2" >> ${LINE}.txt
      FOUND=1;
	  break 2
     elif [ $ADDLINE = $V2 ] && [ $FOUND -eq 0 ]; then
      echo "$V1" >> ${LINE}.txt
      echo "$V2" >> ${LINE}.txt
      FOUND=1;
	  break 2
	 fi
	fi
   done < tmp.txt
 done < Correlation_HighlyFactors.txt
 if [ $FOUND -eq 0 ]; then
  echo "$V1	$V2" >> NotFoundPairs4.txt
 fi
done < NotFoundPairs3.txt

rm -f tmp.txt #NotFoundPairs3.txt
sed -i '1d' NotFoundPairs4.txt

echo "Groups of correlated parameters:" > Grouped_factors.txt

while read -r LINE
do
 sort ${LINE}.txt | dos2unix | uniq > ${LINE}_group.txt
 rm -f ${LINE}.txt
 sed -z 's/\n/\t/g' ${LINE}_group.txt >> Grouped_factors.txt
 echo '' >> Grouped_factors.txt
 rm -f ${LINE}_group.txt
done < Correlation_HighlyFactors.txt

