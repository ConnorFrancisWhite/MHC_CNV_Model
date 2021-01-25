#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
//#include <random>
#include <algorithm> 
#include <armadillo>
#include <ctime>
//Fitness calculations
double fitness_calculation_simple(int A1, int A2, double F1, double F2);
double fitness_calculation(int A1, int A2, double F1, double F2,  double p1, double p2);
double fitness_calculation_mutliple_loci(arma::imat& genes, arma::mat& allele_fitness, arma::mat& gene_frequency);
double fitness_calculation_mutliple_loci_mean(arma::irowvec genes, arma::mat allele_fitness, arma::mat gene_frequency);
double fitness_calculation_max(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency);
double fitness_calculation_min(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency);
double fitness_calculation_mean(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency);
double fitness_calculation_mean_chrom(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency);
double fitness_calculation_mean_no_freq(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency);

double fitness_calculation_max_Lcost(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency);
double fitness_calculation_mean_Lcost(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency);

double fitness_calculation_same(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency);

arma::irowvec Recombination(arma::icube& genes, int loci, int parent, int gen, int inherited_chromosome);
arma::irowvec Recombination_2(arma::imat& genes, int loci, int parent, int inherited_chromosome);
arma::irowvec Allele_Recombination_2(arma::imat& genes, int loci, int parent, int inherited_chromosome);
arma::irowvec Exon_Mutation(arma::irowvec genes, int & exons, arma::mat& allele_fitness, arma::mat& gene_frequency);
arma::irowvec Exon_Mutation_Recycle(arma::irowvec genes, int & exons, arma::mat& allele_fitness, arma::mat& gene_frequency);
arma::irowvec Gene_Duplication(arma::imat& genes, int loci, int parent, int inherited_chromosome);

void Change_Allele_Fitness(arma::mat& allele_fitness);
void Change_Proportion_Allele_Fitness(arma::mat& allele_fitness,double prop);
void Change_Proportion_Allele_Fitness_Normal_Dist(arma::mat& allele_fitness,double prop);

void gene_freq_calc(arma::imat& genes, arma::mat& gene_frequency);
void gene_freq_calc_all_alleles(arma::imat& genes, arma::mat& gene_frequency);

long seedgen(void);


