#!/bin/sh

#SBATCH --job-name=TreeREx_without_tRNA
#SBATCH --account=nn9408k
#SBATCH --output=TreeREx_without_tRNA_slurm-%j.txt
#SBATCH --mem-per-cpu=3G
#SBATCH --ntasks=16
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=t.h.struck@nhm.uio.no

set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

## Create and move to work dir
backdir=$(pwd)
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory:
cp GeneOrder_aligned_reducedMissingness.fas  Masked_18S_rerooted.treefile  $workdir

cd $workdir

# place TreeREx commands after this line

/cluster/projects/nn9408k/Programs/TreeREx/trex-1.85-64bit/trex -f GeneOrder_aligned_reducedMissingness.fas -t Masked_18S_rerooted.treefile -d GeneOrder_aligned_reducedMissingness.dot -v -s -w -W > GeneOrder_aligned_reducedMissingness.out
dot -Tpdf GeneOrder_aligned_reducedMissingness.dot > GeneOrder_aligned_reducedMissingness.pdf

## Make sure the results are copied back to the submit directory:
cp -r * $backdir
cd $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####