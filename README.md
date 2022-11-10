# AnnelidaMitoGenomes
This repository contains the scripts used to conduct the analyses of annelid mitochondrial genomes. Below is a step-by-step description on how the scripts can be used as they have been used in the paper itself. All data to be used and to repeat these analyses can be found in the zipped folder "[TestData_Compiling](https://github.com/torstenstruck/AnnelidaMitoGenomes/blob/main/TestData_Compiling.zip)" here at GitHub as well as at DataDryad. The folder contains a folder structure that works for these analytical procedures. It also contains all the scripts used herein. 

## Analytical procedures

1. In addition to our own mitochondrial genomes, retrieve mitochondrial genomes from NBCI using the following settings as of 17.11.2021 (plus recently published Hydroides and Ophryotrocha genomes, which were not part of the blast database):

        ```
        tblastx
        Stygocapitella subterranea COI-barcode:
        >A355_03E	TACTTTGTATTTTATTTTGGGGATCTGGTCTGGGTTGTTGGGTGCTTCTATAAGGGTTATTATTCGGGTTGAGCTTAGTCAACCTGGGTCTTGGTTAGGTAGTGATCAGATTTATAATACGGTTGTTACTGCTCACGCTTTGCTTATGATTTTCTTTTTGGTCATGCCTGTGATGATTGGGGGGTTTGGAAATTGACTTATTCCTTTAATGATTTCTAGTCCTGATATGGCTTTTCCTCGTATAAATAATATGAGTTTTTGGTTATTACCTCCTGCTTTAATTTTGTTATTAATATCATCTATTTTAGAGAAGGGTGTTGGTACTGGTTGGACTATTTATCCTCCTCTGTCTAGGTCTTTAGGTCATAGAGGTTCTTCGGTCGATTTGGCTATTTTTTCTTTACATTTGGCTGGGGTTTCTTCTATTTTGGGGGCTGCAAATTTTATTACTACTATTTTTAATATGCGGGCATTAGGGTTGCGGTTGGAACGTGTTCCTTTGTTTATTTGGTCAGTTGTTATTACTGCTATTTTGTTGTTACTTTCTTTACCTGTGTTGGCTGGTGCTATTACTATGTTATTAACTGATCGTAATATTAATACTTCGTTTTTTGATCCTGCTGGTGGTGGTGATCCTATTCTTTTCCAACATCTCTTT
        limited to records that include: Annelida (taxid:6340), and entrez query: 12000:200000[slen]
        only one mitochondrial genome per species was kept (the larger one)
        ```

        ```
        DATA: The mitochondrial genomes used herein can be found in 01_Data/MitochondrialGenomes/Used.
        ```

2. Retrieve 18S rRNA sequences (at least 1500 bp) from NCBI using MegaBlast and AF412810 to match mitochondrial genomes; the following order was applied:  
     - if possible, from the same individual (mostly own sequences, but also for NGS mitochondrial genomes such as Glyceridae and Aphroditiformia)  
     - from the same species  
     - from the same genus (only if a and b are not already covered by one species from this genus)  
     - Ramisyllys and Trypanosyllis were excluded as they had very divergent sequences making the alignment problematic and introducing unnecessary extreme long branches in the tree

        ```
        DATA: The 18S rRNA sequences used herein can be found in 01_Data/18Sdata.
        ```

3. Annotate mitochondrial genomes using MITOS2 with the following settings:  
     - Reference: RefSeq 63 Metazoa  
     - Genetic Code: 5 Invertebrate  
     - advanced settings: exclude OH and OL search  
     - all mitochondrial genomes were manually investigated to detect problematic issues and possible genes not found by MITOS2  
     - information on overlaps and missing genes saved in SPECIES_NAME_additional_info.txt  

        ```
        DATA: The annotated mitochondrial genomes used herein can be found in 02_Annotation/Used.
        ```

4. Compile mitochondrial datasets  
     - Compiling structural and sequence information for the mitochondrial genomes into different files running the custom-made shell script "CompileDatasets.sh" as sh CompileDatasets.sh from the top-level folder (the folder where the script is in the test data). This will generate a new folder "03_MitochondrialProperties" with all the relevant information.  
     - Compile the information on intergenic regions using the scripts "CheckFileNames.sh" (to check if the file names are correct and occur only once) and "RetrieveIntergenicParts.sh".  

        ```
        NOTES: 
        The sequence of each whole mitochondrial genomes as a fasta file with the extension ".fasta" needs to be in the folder "./01_Data/MitochondrialGenomes/Used/".
        From the MITOS2 analyses the ".bed"-file of each whole mitochondrial genome needs to be in a folder with the species names (e.g., "Terebratulina_retusa". All these folders need to be in the folder "./02_Annotation/Used/".
        Additionally, a text file "ListSpeciesSize.txt" is needed that lists per species the total number of position of each mitochondrial genome separated by a tab. See provided example.
        ```

5. Generate reference tree using 18S sequences  
     - Compile information for constraint tree based on:  
          - backbone Annelida: Struck, 2018, Handbook (Fig. 7.2)  
          - Clitellata: Erseus et al. 2020, Zooligica Scripta & Anderson et al. (2017) BMC Evol. Biol.  
          - Hirudinea: Phillips et al. 2019, GBE  
          - Terebellida: Stiller et al. 2020  
          - Orbiniida: Struck et al., 2015 (Fig. 1)  
          - Eunicida: Struck et al., 2015 (Fig. 1)  
          - Phyllodocida: Tilic et al., 2022 (Fig. 2 ML)  
          - Aphroditiformia: Zhang Y, Sun J, Rouse GW, Wiklund H, Pleijel F, Watanabe HK, Chen C, Qian P-Y, Qiu J-W. 2018. Phylogeny, evolution and mitochondrial gene order rearrangement in scale worms (Aphroditiformia, Annelida). Molecular Phylogenetics and Evolution 125:220-231. (Fig. 2)  
          - Protodriliformia: Struck et al., 2015 (Fig. 1)  
          - Sabellida: Tilic et al. (2020) MPE  
          - Sipuncula: Lemer et al. (2015)  
     - Compared species names using  

          ```
          for both datasets of the mitochondrial genomes and 18S rRNA:
          grep '>' < All_MitochondrialGenomes.fas | sort > SpeciesNames.txt 
          cmp SpeciesNames.txt ../../01_Data/TreeReconstruction/SpeciesNames.txt
          ```
	
     - Differences in species names were fixed  
     - Alignment with MAFFT/7.470, --auto and --reorder  
     - Trimmed sequence ends using AliView, so that less than 10 sequences had no information on the first and last column  
     - Mask potentially non-homologous positions using AliScore and AliCut  

          ```
          perl Aliscore.02.2.pl -i 18S_data_aligned_trimmed.fasta -N -r 1253698754426984236
          perl ALICUT_V2.3.pl -s
          ```

     - Estimate branch length for constraint tree as well as relationships within groups (where appropriate) using IQtree (IQ-TREE/1.6.12-foss-2018b):  

          ```
          iqtree -s ALICUT_18S_data_aligned_trimmed.fasta -m MFP+LM -g ConstraintTree.tre -o Lineus_viridis -pre Masked_18S -nt AUTO
          iqtree -s 18S_data_aligned_trimmed.fasta -m MFP+LM -g ConstraintTree.tre -o Lineus_viridis -pre Unmasked_18S -nt AUTO
          ```

     - Generating ultrametric trees from both trees using chronos of the APE package in R (script Generate_UltrametricTree.R)  

          ```
          Notes:
	  The following files are needed Masked_18S.treefile and Unmasked_18S.treefile in the same folder as the script.
          ```
	  
6. Determine the different properties for different parts of the genome  
     - Complete mitochondrial data in the folder "03_MitochondrialProperties/WholeGenome"  
          - Open AliView and save as fasta format to ensure that all sequences have the same length (adding - at the end)  
          - Compiled information such as frequencies and skew values using BaCoCa.v1.109:  

          ```
          perl BaCoCa.v1.109.pl -i All_MitochondrialGenomes_alignment.fas
          ```

          - Compiled information such as nRCFV values using RCFV_Reader:  

          ```
          perl NuclRCFVReader.pl All_MitochondrialGenomes_alignment.fas WholeGenome
          ```

     - Code structural information in the folder "03_MitochondrialProperties/StructuralInformation"  
          - Let all gene orders start with COX1 at position 1 and align the gene order by including NA for genes in the missing region for incomplete genomes and - for the most likely position in complete genomes as well as by removing duplications and coding them as absent/present in a spearate sheet and count each  
          - Export aligned gene orders with and without tRNAs as tab-delimited text file; then convert it to a fasta file and remove NA and - again  
          - Upload .fas file with and without tRNAs to [CRex](http://pacosy.informatik.uni-leipzig.de/crex/form)  
          - Run the analyses with default settings (except for "remove duplicates) for datasets determining "common interval", "breakpoints" and "reversal distance"  
          - Copy the distance matrices for each dataset and setting into a separate excel sheet (DistanceMatricesGeneOrder.xlsx in "03_MitochondrialProperties")  
          - For the TRex analyses the species were sorted based on their phylogenetic affiliation and the dataset without tRNAs has been modified in the following way to reduce missingness:  
              - Missing (-) and lacking genes (NA) have been included in accordance to the close relatives and then after the analyses removed again as losses or without further consideration  
          - Export aligned and filled gene order without tRNAs as tab-delimited text file and convert to fasta file  
          - Run the TRex-Analyses using the 18S tree (Masked_18S_rerooted.treefile) and GeneOrder_aligned_reducedMissingness.fas (both are in "TreeREx_Analyses"):  

          ```
          trex -f GeneOrder_aligned_reducedMissingness.fas -t Masked_18S_rerooted.treefile -d GeneOrder_aligned_reducedMissingness.dot -v -s -w -W > GeneOrder_aligned_reducedMissingness.out
          dot -Tpdf GeneOrder_aligned_reducedMissingness.dot > GeneOrder_aligned_reducedMissingness.pdf
          ```

          - The number of the different gene orders were counted using "Count_SequenceOrders.R" for both the gene order with and without tRNA using:  

          ```
	  Notes:
	  The following files are needed:
	  GeneOrder_aligned_with_tRNA.txt
	  GeneOrder_aligned_without_tRNA.txt
          ```
	  
	  - Generate fasta files from the files Counts_wotRNAs.txt and Counts_wotRNAs_MatchedSpecies.txt as well as Counts_wtRNAs.txt and Counts_wtRNAs_MatchedSpecies.txt, respectively  
          - A qMGR analysis were conducted at [qGMR](http://qmgr.hnnu.edu.cn/) using the most common gene order as the ground pattern and without the outgroup providing the two fasta files (Counts_Sequence_Species_wotTRNAs.fas and Counts_Sequence_Species_wtTRNAs.fas in "qMGR_Analyses_outgroup_excluded")  
  		
     - protein-coding and rRNA genes  
          - For each protein- and rRNA-coding gene determine the species lacking in the dataset, adjust constraint tree and compare species names from tree:
  
          ```
          while read gene; do grep $gene < GeneOrder_aligned_without_tRNA.txt | awk '{print $1}' | sort > ${gene}_SpeciesNames.txt; done <GeneNamesAllNuc_without_tRNAs.txt
          while read gene; do grep -v $gene < GeneOrder_aligned_without_tRNA.txt | awk '{print $1}' > ${gene}_SpeciesNames_Not.txt; done <GeneNamesAllNuc_without_tRNAs.txt
          
	  adjust tree and also generate a sorted list of species names within tree
          for each gene: 
	  cmp atp6_SpeciesNames.txt atp6_aligned.fasta_ConstraintTree.tre_SpeciesNames.txt
	  ```
	  
	  ```
	  DATA: These files can be found in the folder "01_Data/TreeReconstruction/ConstraintTrees_Genes"
	  ```

          - For each gene (including the rRNA genes, but not tRNAs): Compared species names between dataset and constraint tree using:  

          ```
	  for file in *.fas; do grep '>' < $file | sort | sed 's/>//' > ${file}_SpeciesNames.txt; done 
	  for each dataset: cmp ../../01_Data/TreeReconstruction/ConstraintTrees_Genes/atp6_aligned.fasta_ConstraintTree.tre_SpeciesNames.txt atp6.fas_SpeciesNames.txt
	  ```

          - fix differences in species names  
          - move to tRNA genes to new folder (e.g., Genes_NotUsed) #Not part of test data provided at DataDryad  
           
          ```
	  DATA: Gene files are in the folder "03_MitochondrialProperties/CodingGenes/unaligned_data/"
	  ```

          - Alignment with Mega 11.0.10 using MUSCLE with codons, invertebrate mitochondrial code 5 and default settings to obtain both an amino acid alignment and a nucleotide alignment based on the amino acid one  
          - Alignment with Mega 11.0.10 using MUSCLE without codons and default settings for rRNA genes  

          ```
	  DATA: Alignment files are in the folders "03_MitochondrialProperties/CodingGenes" and "03_MitochondrialProperties/ProteinCodingGenes"			
	  ```

          - For both nucleotides and amino acids:  
              - Generate different supermatrices (only protein_coding, only rRNA, both together) using FASconCAT:  
              
          ```
	  perl FASconCAT-G_v1.05.pl -s -l
	  ```

          - For each gene & supermatrix:  
              - Reconstruct a tree to obtain tree-based measurements using  

          ```
          for file in *.fasta; 
          do 
              iqtree -s ${file} -m MFP -g ${file}_ConstraintTree.tre -pre ${file} -nt AUTO
              or
              iqtree -s ${file} -m MFP+MERGE -g ${file}_ConstraintTree.tre -spp ${file}_partitions.txt -pre ${file} -nt AUTO
          done
          ```

          ```
          DATA: Files are in the folder "03_MitochondrialProperties/CodingGenes/" and "03_MitochondrialProperties/ProteinCodingGenes/"
          ```

          - Compile information such as frequencies and skew values using BaCoCa.v1.109:  

          ```
          for FILE in *.fasta
          do
              tr '[:lower:]' '[:upper:]' < ${FILE} > ${FILE}.fas
	      perl BaCoCa.v1.109.pl -i ${FILE}.fas
	      mv BaCoCa_Results ${FILE}_BaCoCa_Results
          done
          ```

          ```
          DATA: Files are in the folder "03_MitochondrialProperties/CodingGenes/" and "03_MitochondrialProperties/ProteinCodingGenes/"
          ```

          - Compiled information such as ntsRCFV values using RCFV_Reader:  

          ```
          perl NuclRCFVReader.pl|ProtRCFVReader.pl FILE_NAME FILE_NAME_Results
          ```

          ```
          DATA: Files are in the folder "03_MitochondrialProperties/CodingGenes/" and "03_MitochondrialProperties/ProteinCodingGenes/"
          ```

          - Calculate evolutionary distances, LB scores and tip-to-root distance using  

          ```
          ls *.treefile > TreeNames.txt
          perl TreSpEx.v1.2.pl -fun e -ipt TreeNames.txt -tf SpeciesNames.txt
          ```

          ```
          DATA: Files are in the folder "03_MitochondrialProperties/CodingGenes/" and "03_MitochondrialProperties/ProteinCodingGenes/"
          ```

          - Compile the different sequence-based properties generated by BaCoCA, RCFVReader and TreSpEx into a single xlsx-file (see CompiledProperties_Sequence.xlsx in "03_MitochondrialProperties")  
  
7. Correla tion analyses of the mitochondrial data  
     - Gene order distance data with and without tRNA included: Determine the correlation of the three matrices to each other using CorrelationsGeneOrder.r  
```
		Notes:
		#the following files are needed:
		#Breakpoints_wtRNA.csv
		#CommonInterval_wtRNA.csv
		#ReverseDistance_wtRNA.csv
```
```
		DATA: Files are in the folder "04_Macroevolution/01_CorrelationAnalyses/GeneOrder/with_tRNAs" and "04_Macroevolution/01_CorrelationAnalyses/GeneOrder/without_tRNAs"
```

     - Sequence property data: determine correlations between parameters and reduce highly positively and negatively ones to one using CorrelationsSeqProperties_Final.R:  
```
		Notes:
		#files needed:
		#CompiledProperties.csv
```
```
		DATA: Files are in the folder "04_Macroevolution/01_CorrelationAnalyses/GeneOrder/SequenceProperties"
```

     - Retrieve all groups of correlated parameters "RetrieveCorrelatedGroups.sh":  
          - The file "Correlation_HighPairs.txt" includes the highly correlated pairs and the quotations marks need to be removed from it as well as the first line  
          - Make a list of all parameters left over after cleaning by generating a file "Correlation_HighlyFactors.txt" containing all the column names of the datamatrix reduced to the correlated characters "CompiledSeqProperties_Corr_Clean.txt"  
```
	        Notes:
		#the following files are needed:
		#Correlation_HighlyFactors.txt
		#Correlation_HighPairs.txt
```
```
		DATA: Files are in the folder "04_Macroevolution/01_CorrelationAnalyses/GeneOrder/SequenceProperties/RetrieveCorrelatedGroups"
```

     - Generate a heatmap using "GenerateHeatmap.R"  
          - Add up the different categories in excel (see "Grouped_factors.xlsx") and export as "PropertiesCorrelatedFactors.txt"  
```
		Notes:
		#the following files are needed:
		#PropertiesCorrelatedFactors.txt
```
```
		DATA: Files are in the folder "04_Macroevolution/01_CorrelationAnalyses/GeneOrder/SequenceProperties/RetrieveCorrelatedGroups"
```

     - Retrieve all clustered groups together using "CompileBothGroupsTogether.sh" (in the folder "RetrieveClusteredGroups"):  
          - Besides the "Grouped_factors.txt" generate a file for each grouped cluster ending on ".group" and containing all the categories of the collapsed clusters  
```
		Notes:
		#the following files are needed:
		#Grouped_factors.txt
		#different files ending on *.group
```
```
		DATA: Files are in the folder "04_Macroevolution/01_CorrelationAnalyses/GeneOrder/SequenceProperties/RetrieveClusteredGroups"
```

     - Repeat step "Generate a heatmap" above using clustered groups instead of correlated groups  
```
		Notes: 
		#the following files are needed:
		#PropertiesClusteredFactors.txt
```
```
		DATA: Files are in the folder "04_Macroevolution/01_CorrelationAnalyses/GeneOrder/SequenceProperties/RetrieveClusteredGroups"
```

     - Run "ExplorationDataAnalyses.R" to explore the different molecular properties in more detail (in the folder "ExplorationData")  
```
		Notes:
		#the following file is needed:
		#CompiledProperties.csv
```
```
		DATA: Files are in the folder "04_Macroevolution/01_CorrelationAnalyses/GeneOrder/SequenceProperties/ExplorationData"
```

8. Phylogenetic Least Square Regression analyses on the molecular properties with the gene orders with and without tRNAs as responses  
     - Calculate the mean RD value for each species for both the gene orders with and without tRNAs from the files ReverseDistance_wtRNA.csv and ReverseDistance_wotRNA.csv  
     - Add these mean values as columns to the "CompiledProperties.csv"  
     - Run "pgls_all.R" using additionally "GeneOrder_aligned_reducedMissingness.tsv", "GeneOrder_aligned_with_tRNA.tsv" and "Masked_18S_ultrametric_rerooted.treefile"  
          - Before the analyses of combined molecular properties against the RD values with and without tRNA determine the best-fitting properties shared and independently based on the three major categories frequency-based, rate-based and structure-based (see Combined_sign_predictors.xlsx)  
```
		Notes:
		#the following files are needed:
		#Masked_18S_ultrametric_rerooted.treefile
		#CompiledProperties.csv
		#GeneOrder_aligned_with_tRNA.tsv
		#GeneOrder_aligned_reducedMissingness.tsv
```
```
		DATA: Files are in the folder "04_Macroevolution/02_PhylogeneticLeastSquareRegression/Molecular_data"
```

9. Phylogenetic Least Square Regression analyses on the life history data with the gene orders with and without tRNAs as responses  
     - From "LifeHistoryEcology.xlsx", exclude the species from the variables with multiple coding or missing data, add two columns for mean reverse distances for gene orders with and without tRNAs and save each variable as separate csv-files; the remaining ones with all species save as one csv-file including two columns for mean reverse distances for gene orders with and without tRNAs  
     - Run "pgls_with_tRNA.R" and "pgls_without_tRNA.R":   
```
	Notes:
	#files needed:
	#Masked_18S_ultrametric_rerooted.treefile
	#the different *.csv files
```
```
	DATA: Files are in the folder "04_Macroevolution/02_PhylogeneticLeastSquareRegression/Life_history"
```

10. Analyse the correlation of the Reverse Distance to each other as well as to the phylogenetic placement of the corresponding taxa  
     - Generate a comma-separated list of families assigned to species (see Family_Species_list.csv)  
     - Transform the treefiles of the analyses of the Nuc_supermatrix_partition, AA_supermatrix_partition, PB_Nuc_chain2 and PB_AA_chain1 to equally spaced cladograms where each branch has a length of 1 and save them as .treefile files  
     - Save the same treefiles as a cladogram without any branch length as .txt files  
```
	Notes:
	#the following files are needed
	#without_tRNAs/ReverseDistance_wotRNA.csv
	#with_tRNAs/ReverseDistance_wtRNA.csv
	#Family_Species_list.csv
	#Nuc_supermatrix_partition_rooted_cladogram.txt.treefile
	#PB_Nuc_chain2_rooted_cladogram.treefile
	#Masked_18S_rooted_cladogram.treefile
	#AA_supermatrix_partition_rooted_cladogram.txt.treefile
	#PB_AA_chain1_rooted_cladogram.treefile"
	#Nuc_supermatrix_partition_rooted_cladogram.txt
	#PB_Nuc_chain2_rooted_cladogram.txt
	#Masked_18S_rooted_cladogram.txt
	#AA_supermatrix_partition_rooted_cladogram.txt
	#PB_AA_chain1_rooted_cladogram.txt
```
```
	DATA: Files are in the folder "04_Macroevolution/01_CorrelationAnalyses/GeneOrder"
```

## Phylogenetic reconstruction:
1. Masking and concatenation running the following commands in the folders "05_Phylogeny/00_Data/Nucleotides" and "05_Phylogeny/00_Data/AminoAcids":  
```
	for FILE in *.fasta; 
	do
		perl Aliscore.02.2.pl -N -r 200000000000 -i $FILE; 
		mv ${FILE}* $backdir
		cd $backdir
	done
	perl ALICUT_V2.31.pl -s

	mkdir Masked

	for FILE in ALICUT*.fasta;
	do
		mv $FILE Masked/
	done
	cd Masked/
	perl FASconCAT-G_v1.05.pl -s -p -p -l > Concatation_log.txt
	cd ..
```

2. Phylogenetic reconstruction using IQtree:  
     - Generate trees  
```
	iqtree -s AA_supermatrix.phy -spp AA_supermatrix_partition.txt -m MFP+MERGE -bb 1000
	iqtree -s Nuc_supermatrix.phy -spp Nuc_supermatrix_partition.txt -m MFP+MERGE -bb 1000
```

     - Generate ultrametric trees from these treefiles:  
          - Reroot the trees at one of the outgroups, double the branch length leading to that outgroup and shorten the one leading to all other taxa to 0.00000000001 to mimic an unrooted tree despite being rooted; save rooted treefile  
          - Generate ultrametric trees using "AA_Generate_UltrametricTree.R" and "Nuc_Generate_UltrametricTree.R":  
```
		Notes:
		#the following files are needed:
		#Nuc_supermatrix_partition_rooted.txt.treefile
		#AA_supermatrix_partition_rooted.txt.treefile
```
```
		DATA: Files are in the folder "05_Phylogeny/01_Unconstraint"
```

3. Constraint phylogenetic reconstruction using IQtree:  
     - Generate trees  
```
	iqtree -s AA_supermatrix.phy -spp AA_supermatrix_partition.txt -m MFP+MERGE -g ConstraintTree_unrooted.tre -o Owenia_fusiformis -pre Constraint_AA
	iqtree -s Nuc_supermatrix.phy -spp Nuc_supermatrix_partition.txt -m MFP+MERGE -g ConstraintTree_unrooted.tre -o Owenia_fusiformis -pre Constraint_Nuc
```

     - Generate ultrametric trees from these treefiles:  
          - Reroot the trees at one of the outgroups, double the branch length leading to that outgroup and shorten the one leading to all other taxa to 0.00000000001 to mimic an unrooted tree despite being rooted; save rooted treefile  
          - Generate ultrametric trees using "AA_Generate_UltrametricTree.R" and "Nuc_Generate_UltrametricTree.R":  
```
		Notes:
		#the following files are needed:
		#Constraint_Nuc_rooted.treefile
		#Constraint_AA_rooted.treefile
```
```
		DATA: Files are in the folder "05_Phylogeny/02_Constraint"
```

4. Phylogenetic reconstruction using PhyloBayes with two chains per dataset (AA and Nuc):  
     - Generate trees  
        - Example of command lines to be used:  
```
	pb_mpi -d AA_supermatrix.phy -cat -gtr -x 1 40000 PB_AA_chain2
	tracecomp -x 10000 PB_AA_chain1 PB_AA_chain2
	bpcomp -x 10000 10 PB_AA_chain1 PB_AA_chain2
```

     - Generate ultrametric trees from these treefiles:  
          - Reroot the trees at one of the outgroups, double the branch length leading to that outgroup and shorten the one leading to all other taxa to 0.00000000001 to mimic an unrooted tree despite being rooted; save rooted treefile  
          - Generate ultrametric trees using "AA_Generate_UltrametricTree.R" and "Nuc_Generate_UltrametricTree.R":  
```
		Notes:
		#the following files are needed:
		#PB_Nuc_chain2.con.rooted.tre
		#PB_AA_chain1.con.rooted.tre
```
```
		DATA: Files are in the folder "05_Phylogeny/03_PhyloBayes"
```
