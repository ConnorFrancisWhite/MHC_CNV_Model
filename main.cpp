#include <fstream>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <ctime>
#include <string>
#include <armadillo>

#include "Datum.h"
#include "BinaryTree.h"
#include "functions.h"


using namespace std;
//using namespace arma;

double R = (double)RAND_MAX;


int main(int argc, char **argv){
	string number = argv[1];
	
	//seed by time or by run number or by seedgen function.
	//srand(atoi(argv[1]));
	//srand (time(0));
	long seed = seedgen();
	srand(seed);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int population = atoi(argv[2]);
	int generations = atoi(argv[3]);
	int exon_fitness_change = generations/10;
	int exons = atoi(argv[4]);
	int exon_space = 100;
	int loci = atoi(argv[5]);
 	int initial_loci = atoi(argv[6]);
 	
 	
 	double recomb_prob = atof(argv[7]);
 	double exon_mut_prob = atof(argv[8]);
	double allele_exon_ratio = atof(argv[9]);
	
 	int max_flag = atoi(argv[10]);
 	int min_flag = atoi(argv[11]);
 	int mean_flag = atoi(argv[12]);
 	double Env_flag_Uni = atof(argv[13]);
 	double Env_flag_Norm = atof(argv[14]);
 	double Sigma = atof(argv[15]);
 	
 	//cout << Sigma << endl;
 	int Mult_Gen_flag = atoi(argv[16]);
 	int Detailed_End = atoi(argv[17]);
 	
 	int skip_gen = Mult_Gen_flag;
 	
 	cout << "Population " << population << endl;
 	cout << "Generations " << generations << endl;
 	cout << "Exons " << exons << endl;
 	cout << "Loci " << loci << endl;
 	cout << "Initial Loci No " << initial_loci << endl;
 	cout << "Recombination Prob " << recomb_prob << endl;
 	cout << "Mutation Prob " << exon_mut_prob << endl;
 	cout << "Exon or Allele Recombination " << allele_exon_ratio << endl;
	
	cout << "Max Flag " << max_flag << endl;
	cout << "Min Flag " << min_flag << endl;
	cout << "Mean Flag " << mean_flag << endl;
	cout << "Env Flag Uni " << Env_flag_Uni << endl;
	cout << "Env Flag Norm " << Env_flag_Norm << endl;
	
	if(Env_flag_Norm > 0.0){
		cout << "Sigma " << Sigma << endl;
	}
	
	cout << "Miltiple Generations " << Mult_Gen_flag << endl;
	
	if(Mult_Gen_flag > 0){
		cout << "Detailed End Flag " << Detailed_End << endl;
	}
	//Defining the fitness of each allele
	arma::mat allele_fitness(exon_space, exon_space);
	
	for(int exon1 = 0; exon1<exon_space; exon1++){
		for(int exon2 = 0; exon2<exon_space; exon2++){
		
			allele_fitness(exon1,exon2) = rand()/R;
		}
	}
	
	
	
	//Defining the intital genes and gene frequencies
	//arma::icube genes(population*2,loci*2,generations);
	
	arma::imat genes_1(population*2,loci*2);
	arma::imat genes_2(population*2,loci*2);
	genes_1.fill(0);
	genes_2.fill(0);
	arma::mat gene_frequency(exon_space, exon_space);
	gene_frequency.fill(0);
	
	
	//Initialising where every exon is in intial population
	/*
	for(int person = 0; person < population*2; person++){
		for(int locus = 0; locus < initial_loci; locus++){
				
			genes_1(person,0 + locus*2) = rand()%exons+1;
			genes_1(person,1 + locus*2) = rand()%exons+1;		
		}
	}
	*/
	
	//This makes sure that initial exons only appear on one side.
	for(int person = 0; person < population*2; person++){
		for(int locus = 0; locus < initial_loci; locus++){
				
			//genes_1(person,0 + locus*2) = rand()%(exons/2)+1;
			//genes_1(person,1 + locus*2) = rand()%(exons/2)+1 + exons/2;		
			genes_1(person,0 + locus*2) = (rand()%(exons/2)+1)*2 - 1;
			genes_1(person,1 + locus*2) = (rand()%(exons/2)+1)*2;	
		}
	}
	
	
	exons = exon_space;
	
	//GENE FREQUENCY CALULATER/////////////////////////////////////////////////////////////
	gene_freq_calc(genes_1, gene_frequency);
	

	//defining intiital populaton fitness
	arma::vec population_fitness(population);
	population_fitness.fill(0);
	std::vector<Datum> popf;
 
 	popf.resize(population);
	BinaryTree *T;
	
	double (*fitness_function)(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency);
	
	if(max_flag == true){
		//fitness_function =  &fitness_calculation_max;
		fitness_function =  &fitness_calculation_max_soft_boundary;
	}
	if(min_flag == true){
	
		fitness_function =  &fitness_calculation_min;
	}
	if(mean_flag == true){
	
		fitness_function =  &fitness_calculation_mean_chrom;
	}
	
	//fitness_function =  &fitness_calculation_same;
	
	for(int person = 0; person < population; person++){
	
	
		arma::imat genome = genes_1.submat(person*2,0,person*2 + 1,loci*2-1);
	//HERE IS WHERE YOU CHANGE THE FITNESSS///////////////////////////////////////////////////////////////////////////////////////////////////

		//popf[person] =  Datum(person + 1, fitness_calculation_max(genome, allele_fitness, gene_frequency ));
		//population_fitness(person) =  fitness_calculation_max(genome, allele_fitness, gene_frequency);
		popf[person] = Datum(person + 1, (*fitness_function)(genome, allele_fitness, gene_frequency ));
		population_fitness(person) =  (*fitness_function)(genome, allele_fitness, gene_frequency);
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//defining integers for who the parents will be.
	int mother;
	int father;
	
	//arma::imat Adjacency_Matrix(1, 1);
	//Adjacency_Matrix.fill(0);
	
	for(int gen = 1; gen <= generations/2; gen++){
	
		//building fenwick tree for each generation
		T = new BinaryTree();
   		T->build(T,popf);
   			
   		//resseting the gene frequency for the next gen
   		gene_frequency.fill(0);
   		
   		//looping over ech person in the population 
   		for(int person = 0; person < population; person++){
   		
   			double range = sum(population_fitness);	
			double motherf = range*rand()/R;
			double fatherf = range*rand()/R;	
		
			//choosing the parents
			if(range == 0){ 
			
				mother = rand()%population;
				father = rand()%population;
								
			}
			else{
				mother = T->search(T,motherf,1).n-1;
				father = T->search(T,fatherf,1).n-1;
			}
			
			
			
			int chromosomem = rand()%2;
			int chromosomef = rand()%2;
			
			//finding which chromsome is the smallest and governing that to be the length that determins how likeley recombination is to happen.
			int chromosomem_other = abs(chromosomem - 1);
			int chromosomef_other = abs(chromosomef - 1);
					
			
			arma::uvec mother_blank_space = find(genes_1.submat(chromosomem+2*mother, 0,    chromosomem+2*mother, loci*2 - 1) == 0);
			arma::uvec father_blank_space = find(genes_1.submat(chromosomef+2*father, 0,    chromosomef+2*father, loci*2 - 1) == 0);
			
			arma::uvec mother_other_blank_space = find(genes_1.submat(chromosomem_other+2*mother, 0,    chromosomem_other+2*mother, loci*2 - 1) == 0);
			arma::uvec father_other_blank_space = find(genes_1.submat(chromosomef_other+2*father, 0,    chromosomef_other+2*father, loci*2 - 1) == 0);
			
			
			int mother_length;
			int father_length;
			
			int mother_recomb_length;
			int father_recomb_length;

			
			mother_length = loci - mother_blank_space.size()/2;
			father_length = loci - father_blank_space.size()/2;
			
			if(mother_blank_space.size() > mother_other_blank_space.size()){
			//if(mother_blank_space.size() < mother_other_blank_space.size()){
				mother_recomb_length = mother_length;
			}
			else{
				mother_recomb_length = loci - mother_other_blank_space.size()/2;
			}
			
			if(father_blank_space.size() > father_other_blank_space.size()){
			//if(father_blank_space.size() < father_other_blank_space.size()){
				father_recomb_length = father_length;
			}
			else{
				father_recomb_length = loci - father_other_blank_space.size()/2;
			}
			
			//recombination or not
			//mother
			if(rand()/R < recomb_prob*mother_recomb_length){
				if(rand()/R < allele_exon_ratio){
					genes_2.row(0 + 2*person) = Recombination_2(genes_1,loci,mother,chromosomem);
				}
				else{
					//genes_2.row(0 + 2*person) = Gene_Duplication(genes_1,loci,mother,chromosomem);
				
					genes_2.row(0 + 2*person) = Allele_Recombination_2(genes_1,loci,mother,chromosomem);
				}
			}
			else{
			
				genes_2.submat(0+2*person, 0,    0+2*person, loci*2 - 1) = genes_1.submat(chromosomem+2*mother, 0,    chromosomem+2*mother, loci*2 - 1);
			}
			//father
			if(rand()/R < recomb_prob*father_recomb_length){
				if(rand()/R < allele_exon_ratio){
					genes_2.row(1 + 2*person)  = Recombination_2(genes_1,loci,father,chromosomef);
				}
				else{
					//genes_2.row(1 + 2*person)  = Gene_Duplication(genes_1,loci,father,chromosomef);
				
					genes_2.row(1 + 2*person)  = Allele_Recombination_2(genes_1,loci,father,chromosomef);
				}
			}
			else{
				
				genes_2.submat(1+2*person, 0,    1+2*person, loci*2 - 1) = genes_1.submat(chromosomef+2*father, 0,    chromosomef+2*father, loci*2 - 1);
			}
			
			//Adding in new exon due to mutation
			if(rand()/R < exon_mut_prob*mother_length){
				
				genes_2.row(0 + 2*person) = Exon_Mutation_Recycle(genes_1.row(chromosomem+2*mother), exons, allele_fitness,  gene_frequency);
			}
			if(rand()/R < exon_mut_prob*father_length){
				
				genes_2.row(1 + 2*person) = Exon_Mutation_Recycle(genes_1.row(chromosomef+2*father), exons, allele_fitness,  gene_frequency);
			}
   
   		}
   		
   		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   		//ENVIRONMENTAL CHANGE
   		//This part changes allele_fitness remember to turn off when running
		if(Env_flag_Uni > 0.0){
	  
			Change_Proportion_Allele_Fitness(allele_fitness, Env_flag_Uni);
		}
		if(Env_flag_Norm > 0.0){
	  
			Change_Proportion_Allele_Fitness_Normal_Dist_Boundary_Force(allele_fitness, Env_flag_Norm,Sigma);
		}
   		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   		
		//HERE YOU HAVE MADE NEW FUNCTION FOR CALCULATING GENE FREQUENCY YOU HAVE TWO TYPES/////////////////////
		gene_freq_calc(genes_2, gene_frequency);
		
		
   		for(int person = 0; person < population; person++){
		
		
			arma::imat genome = genes_2.submat(person*2,0,person*2 + 1,loci*2-1);
	
			//HERE IS WHERE YOU CHANGE THE FITNESSS///////////////////////////////////////////////////////////////////////////////////////////////////
			//popf[person] = Datum(person + 1, fitness_calculation_max(genome, allele_fitness, gene_frequency ));
			//population_fitness(person) =  fitness_calculation_max(genome, allele_fitness, gene_frequency);
			popf[person] = Datum(person + 1, (*fitness_function)(genome, allele_fitness, gene_frequency ) );
			population_fitness(person) =  (*fitness_function)(genome, allele_fitness, gene_frequency );
		
		}
		delete T;
		
		/*
		if((gen+1)*2 - 1 >= generations - 10){
			string gene_file = "data/gene_gen_" + to_string((gen+1)*2-1) + "_job_"  + number + ".mat";
			genes_2.save(gene_file,arma::hdf5_binary);
			string freq_file = "data/allele_frequency_gen_" + to_string((gen+1)*2-1) + "_job_"  + number + ".mat";
			gene_frequency.save(freq_file,arma::hdf5_binary);
			string Pop_Fitness_file = "data/Population_Fitness_gen_" + to_string((gen+1)*2-1) + "_job_"  + number + ".mat";
			population_fitness.save(Pop_Fitness_file,arma::hdf5_binary);
	
		}
		*/
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//building fenwick tree for each generation
		T = new BinaryTree();
   		T->build(T,popf);
   			
   		//resseting the gene frequency for the next gen
   		gene_frequency.fill(0);
   		
   		//looping over ech person in the population 
   		for(int person = 0; person < population; person++){
   		
   			double range = sum(population_fitness);	
			double motherf = range*rand()/R;
			double fatherf = range*rand()/R;			
			//choosing the parents
			
			//This checks if there is only 1 allele and how that is selected.
			if(range == 0){ 
				mother = rand()%population;
				father = rand()%population;
			}
			else{
				mother = T->search(T,motherf,1).n-1;
				father = T->search(T,fatherf,1).n-1;
			}
			
			int chromosomem = rand()%2;
			int chromosomef = rand()%2;
			
			//finding which chromsome is the smallest and governing that to be the length that determins how likeley recombination is to happen.
			int chromosomem_other = abs(chromosomem - 1);
			int chromosomef_other = abs(chromosomef - 1);
					
			
			arma::uvec mother_blank_space = find(genes_2.submat(chromosomem+2*mother, 0,    chromosomem+2*mother, loci*2 - 1) == 0);
			arma::uvec father_blank_space = find(genes_2.submat(chromosomef+2*father, 0,    chromosomef+2*father, loci*2 - 1) == 0);
			
			arma::uvec mother_other_blank_space = find(genes_2.submat(chromosomem_other+2*mother, 0,    chromosomem_other+2*mother, loci*2 - 1) == 0);
			arma::uvec father_other_blank_space = find(genes_2.submat(chromosomef_other+2*father, 0,    chromosomef_other+2*father, loci*2 - 1) == 0);
			
			int mother_length;
			int father_length;
			
			int mother_recomb_length;
			int father_recomb_length;
			
			mother_length = loci - mother_blank_space.size()/2;
			father_length = loci - father_blank_space.size()/2;
			
			
			if(mother_blank_space.size() > mother_other_blank_space.size()){
			//if(mother_blank_space.size() < mother_other_blank_space.size()){
				mother_recomb_length = mother_length;
			}
			else{
				mother_recomb_length = loci - mother_other_blank_space.size()/2;
			}
			
			if(father_blank_space.size() > father_other_blank_space.size()){
			//if(father_blank_space.size() < father_other_blank_space.size()){
				father_recomb_length = father_length;
			}
			else{
				father_recomb_length = loci - father_other_blank_space.size()/2;
			}
			/*
			cout << "Length of chromosome 1 = " << loci - mother_blank_space.size()/2 << endl;
			cout << "Length of chromosome 2 = " << loci - mother_other_blank_space.size()/2 << endl;
			cout << "mother length = " << mother_length << endl;
			cout << "mother recomb length = "<< mother_recomb_length << endl;
			*/
			//recombination or not
			//mother
			if(rand()/R < recomb_prob*mother_recomb_length){
				if(rand()/R < allele_exon_ratio){
					genes_1.row(0 + 2*person) = Recombination_2(genes_2,loci,mother,chromosomem);
				}
				else{
					//genes_1.row(0 + 2*person) = Gene_Duplication(genes_2,loci,mother,chromosomem);
				
					genes_1.row(0 + 2*person) = Allele_Recombination_2(genes_2,loci,mother,chromosomem);
				}
			}
			else{
			
				genes_1.submat(0+2*person, 0,    0+2*person, loci*2 - 1) = genes_2.submat(chromosomem+2*mother, 0,    chromosomem+2*mother, loci*2 - 1);
			}
			//father
			if(rand()/R < recomb_prob*father_recomb_length){
				if(rand()/R < allele_exon_ratio){
					genes_1.row(1 + 2*person)  = Recombination_2(genes_2,loci,father,chromosomef);
				}
				else{
					//genes_1.row(1 + 2*person)  = Gene_Duplication(genes_2,loci,father,chromosomef);
				
					genes_1.row(1 + 2*person)  = Allele_Recombination_2(genes_2,loci,father,chromosomef);
				}
			}
			else{
				
				genes_1.submat(1+2*person, 0,    1+2*person, loci*2 - 1) = genes_2.submat(chromosomef+2*father, 0,    chromosomef+2*father, loci*2 - 1);
			}
			
			//Adding in new exon due to mutation
			if(rand()/R < exon_mut_prob*mother_length){
														
				genes_1.row(0 + 2*person) = Exon_Mutation_Recycle(genes_2.row(chromosomem+2*mother), exons, allele_fitness,  gene_frequency);
			}
			if(rand()/R < exon_mut_prob*father_length){
				
				genes_1.row(1 + 2*person) = Exon_Mutation_Recycle(genes_2.row(chromosomef+2*father), exons, allele_fitness,  gene_frequency);
			}
			
   
   		}
   		
   		
   		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   		//ENVIRONMENTAL CHANGE
   		//This part changes allele_fitness remember to turn off when running
		if(Env_flag_Uni > 0.0){
	  		
			Change_Proportion_Allele_Fitness(allele_fitness, Env_flag_Uni);
		}
		if(Env_flag_Norm > 0.0){
	  		
			Change_Proportion_Allele_Fitness_Normal_Dist_Boundary_Force(allele_fitness, Env_flag_Norm,Sigma);
		}
   		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   		
		//HERE YOU HAVE MADE NEW FUNCTION FOR CALCULATING GENE FREQUENCY YOU HAVE TWO TYPES/////////////////////
		gene_freq_calc(genes_1, gene_frequency);
		
   		for(int person = 0; person < population; person++){
		
		
			arma::imat genome = genes_1.submat(person*2,0,person*2 + 1,loci*2-1);
	
			//HERE IS WHERE YOU CHANGE THE FITNESSS///////////////////////////////////////////////////////////////////////////////////////////////////
			popf[person] = Datum(person + 1, (*fitness_function)(genome, allele_fitness, gene_frequency ) );
			population_fitness(person) =  (*fitness_function)(genome, allele_fitness, gene_frequency );
		
		}
		
		delete T;
		
		
		/*
		if((gen+1)*2 >= generations - 10){
			string gene_file = "data/gene_gen_" + to_string((gen+1)*2) + "_job_"  + number + ".mat";
			genes_1.save(gene_file,arma::hdf5_binary);
			string freq_file = "data/allele_frequency_gen_" + to_string((gen+1)*2) + "_job_"  + number + ".mat";
			gene_frequency.save(freq_file,arma::hdf5_binary);
			string Pop_Fitness_file = "data/Population_Fitness_gen_" + to_string((gen+1)*2) + "_job_"  + number + ".mat";
			population_fitness.save(Pop_Fitness_file,arma::hdf5_binary);
		
		}
		*/
		//Printing results part
		//This is one gen before environment change
		
		
		
		if(Mult_Gen_flag){
			
			//This triggers for more frequenctly recorded data in the last 1/4 timesteps. 
			if(Detailed_End){
				if((gen+1)*2 == generations - generations/4){
			
					skip_gen = skip_gen/10;
				}
			}
			if((gen+1)*2%(skip_gen) == 0){
				string gene_file = "data/gene_gen_" + to_string((gen+1)*2) + "_job_"  + number + ".mat";
				genes_1.save(gene_file,arma::hdf5_binary);
				
				
				string freq_file = "data/allele_frequency_gen_" + to_string((gen+1)*2) + "_job_"  + number + ".mat";
				gene_frequency.save(freq_file,arma::hdf5_binary);
				
				string fitness_file = "data/allele_fitness_gen_" + to_string((gen+1)*2) + "_job_"  + number + ".mat";
				allele_fitness.save(fitness_file,arma::hdf5_binary);
				/*
				string Pop_Fitness_file = "data/Population_Fitness_gen_" + to_string((gen+1)*2) + "_job_"  + number + ".mat";
				population_fitness.save(Pop_Fitness_file,arma::hdf5_binary);
				*/
			}
			
			
		}
		
	}
	
	//cerr<< "Time elapsed = " << elapsedTime << " secs"  <<endl;
  	arma::imat genome = genes_2;
  	//genes.save("data/genes.mat",arma::hdf5_binary);
  	//genome.save(argv[8],arma::hdf5_binary);
  	
	string genome_file = "data/genome_"  + number + ".mat";
	string allele_fitness_file = "data/allele_fitness_"  + number + ".mat";
	string allele_frequency_file = "data/allele_frequency_"  + number + ".mat";
	
	genome.save(genome_file,arma::hdf5_binary);
  	allele_fitness.save(allele_fitness_file,arma::hdf5_binary);
  	//gene_frequency.save(allele_frequency_file,arma::hdf5_binary);
}
