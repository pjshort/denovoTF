# simulate de novos and split into lots of chunks
cd ~/software/DenovoSim/scripts/
sim_dir="$HOME/experiments/simulated_data/diagnosed_sim_chunks"
logs="$HOME/experiments/simulated_data/logs"
data="$HOME/experiments/simulated_data/data"
iterations=10000
chunks=1000
mkdir $sim_dir # where simulated files will be stored
mkdir $logs # where farm logs will be written
mkdir $data # where combined simulation data will be written

bsub -R'select[mem>850] rusage[mem=850]' -M850 -o $logs/diagnosed_simulation.out /software/R-3.1.2/bin/Rscript \
simulate.R --verbose --n_snps=125 --n_probands=119 --iterations=10000 --n_chunks=1000 --base_name=$sim_dir/sim \
--regions=~/software/denovoTF/data/DDD_well_cov_regions.txt

cd $HOME/software/denovoTF/scripts

# run denovoLOBGOB on each of the chunks
bsub -J "denovoLOBGOBsimulation[1-\$chunks]" -R'select[mem>450] rusage[mem=450]' -M450 -o \
$logs/denovoTF.%I.out /software/R-3.1.2/bin/Rscript denovoLOBGOB.R \
--verbose --de_novos=$sim_dir/sim.\$LSB_JOBINDEX.txt \
--out=$sim_dir/simulated_LOBGOB.\$LSB_JOBINDEX.txt --min_score=95%

# combine the results into a single large dataframe
bsub -R'select[mem>500] rusage[mem=500]' -M500 -o $logs/combine_data.out \
/software/R-3.1.2/bin/Rscript combine.R --n_chunks=$chunks --base_name=$sim_dir/simulated_LOBGOB \
--out=$data/simulated_de_novos_JASPAR_tfbs_annotation.txt

