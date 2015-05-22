# simulate de novos and split into lots of chunks
tmp_dir="/tmp/denovoLOBGOB_sim_chunks"
logs="~/experiments/denovoLOBGOB/simulations"
n_chunks=1000
iterations=1000
mkdir $tmp_dir # where the simulation data files will be temporarily stored

bsub -R'select[mem>500] rusage[mem=500]' -M500 -o $logs/generate_data.out \
/software/R-3.1.2/bin/Rscript simulateDN.R --verbose --base_name=$tmp_dir/sim_data \
--regions="../data/DDD_well_cov_regions.txt"  --n_snps=453 --n_probands=425 --iterations=$iterations \
--n_chunks=$n_chunks


# run denovoLOBGOB on each of the chunks
bsub -J "denovoLOBGOBsimulation[1-\$n_chunks]" -R'select[mem>450] rusage[mem=450]' -M450 -o \
$logs/test.%I.out /software/R-3.1.2/bin/Rscript denovoLOBGOB.R \
--verbose --de_novos=$tmp_dir/sim_data.\$LSB_JOBINDEX.txt \
--out=$tmp_dir/simulated_LOBGOB.\$LSB_JOBINDEX.txt --min_score=95%

# combine the results into a single large dataframe
bsub -R'select[mem>500] rusage[mem=500]' -M500 -o $logs/combine_data.out \
/software/R-3.1.2/bin/Rscript combine.R --n_chunks=$n_chunks --base_name=$tmp_dir/simulated_LOBGOB \
--out="../results/simulated_de_novos_JASPAR_tfbs_annotation.txt"

