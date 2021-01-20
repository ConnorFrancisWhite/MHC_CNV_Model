#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=2012MB
#SBATCH --time=0-04:00:00
#int population = atoi(argv[2]);
#int generations = atoi(argv[3]);
#int exon_fitness_change = generations/10;
#int exons = atoi(argv[4]);
#int loci = atoi(argv[5]);
#int initial_loci = atoi(argv[6]);


#double recomb_prob = atof(argv[7]);
#double exon_mut_prob = atof(argv[8]);
#double allele_exon_ratio = atof(argv[9]);
#int max_flag = atoi(argv[10]);
#int min_flag = atoi(argv[11]);
#int mean_flag = atoi(argv[12]);


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/coolbeans/arma/armadillo-8.300.2 
./Exon ${NUMBERARG} 5000 5000 1 41 8 0.0001 0.000001 0 1 0 0 0 1


