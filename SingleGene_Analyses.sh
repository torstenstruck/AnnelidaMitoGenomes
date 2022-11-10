#!/bin/bash

for FILE in *.fasta
do
echo "#!/bin/sh" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh
echo "#SBATCH --job-name=iqtree_${FILE}" >> iqtree_${FILE}.sh
echo "#SBATCH --output=iqtree_${FILE}_slurm-%j.txt" >> iqtree_${FILE}.sh
echo "#SBATCH --account=nn9408k" >> iqtree_${FILE}.sh
echo "#SBATCH --time=2-00:00:00" >> iqtree_${FILE}.sh
echo "#SBATCH --nodes=1" >> iqtree_${FILE}.sh
echo "#SBATCH --cpus-per-task=16 " >> iqtree_${FILE}.sh
echo "#SBATCH --mem-per-cpu=4G" >> iqtree_${FILE}.sh 

echo "set -o errexit # exit on errors" >> iqtree_${FILE}.sh
echo "set -o nounset  # Treat any unset variables as an error" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "module --quiet purge  # Reset the modules to the system default" >> iqtree_${FILE}.sh
echo "module load IQ-TREE/1.6.12-foss-2018b" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "## Create work dir" >> iqtree_${FILE}.sh
echo 'backdir=$(pwd)' >> iqtree_${FILE}.sh
echo 'workdir=$USERWORK/$SLURM_JOB_ID' >> iqtree_${FILE}.sh
echo 'mkdir -p $workdir' >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "## Copy input files to the work directory and move to work dir:" >> iqtree_${FILE}.sh
echo "cp ${FILE} ${FILE}_ConstraintTree.tre \$workdir" >> iqtree_${FILE}.sh
echo 'cd $workdir' >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "## Generate the phylogenetic tree of the aligned sequences:" >> iqtree_${FILE}.sh
echo "iqtree -s ${FILE} -m MFP -g ${FILE}_ConstraintTree.tre -pre ${FILE} -nt AUTO" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "## Make sure the results are copied back to the submit directory:" >> iqtree_${FILE}.sh
echo 'cp -r * ${backdir}' >> iqtree_${FILE}.sh
echo 'cd ${backdir}' >> iqtree_${FILE}.sh
echo 'rm -rf ${workdir}' >> iqtree_${FILE}.sh

sbatch iqtree_${FILE}.sh
done