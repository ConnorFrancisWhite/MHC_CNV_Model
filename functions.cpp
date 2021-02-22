#include "functions.h"

using namespace std;


double fitness_calculation_simple(int A1, int A2, double F1, double F2){

	double output;
	//Homozygote
	if(A1 == A2){
	
		output = F1;
		return output;
	
	}
	//heterozygote
	else{
	
		output = F1 + F2;
		return output;
		
	}
}

double fitness_calculation(int A1, int A2, double F1, double F2,  double p1, double p2){
	
	double output;
	//Homozygote
	if(A1 == A2){
	
		output = F1*(1-p1);
		return output;
	
	}
	//Heterozygote
	else{
	
		output = F1*(1-p1) + F2*(1-p2);
		return output;
		
	}
}

double fitness_calculation_mutliple_loci(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency){

	double output = 0;
	
	//Finding the unique alleles for the given person, as not to count the fitness from a allele multiple times. 
	
	
	for(int loci = 0; loci < genome.n_cols/2;loci++){
		
		if(genome(0,loci*2) != 0){
			output = output + allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1));
		}
		if(genome(1,loci*2) != 0){
			output = output + allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1));
		}
	}
	return output;
}
double fitness_calculation_mutliple_loci_mean(arma::irowvec genes, arma::mat allele_fitness, arma::mat gene_frequency){

	double output = 0;
	int zero = 0;
	//Finding the unique alleles for the given person, as not to count the fitness from a allele multiple times. 
	arma::irowvec U = unique(genes);
	
	for(int gene = 0; gene < U.size(); gene++){
		//checks for empty alleles so it does not contribute to the fitness
		if(U(gene)==0){
			zero = 1;
			
		}
		else{
		
			output = output + allele_fitness(U(gene)-1)*(1-gene_frequency(U(gene)-1));
		}
	}
	
	output = output/((double)U.size() - zero);
	return output;
}

double fitness_calculation_max(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency){
	
	double output = 0;

	if(genome(0,genome.n_cols/2-1) != 0 || genome(1,genome.n_cols/2-1) != 0){
	
		output = 0;
	}
	else{
		for(int loci = 0; loci < genome.n_cols/2;loci++){
		
			//maternal
			if(genome(0,loci*2) != 0){
				if(output < allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1))){
					output = allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1));
				}

			}
			//paternal
			if(genome(1,loci*2) != 0){
				if(output < allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1))){
					output = allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1));
				}

			}
	
		}
	}
	
	
	return output;
}

double fitness_calculation_min(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency){
	
	double output = 1;

	if(genome(0,genome.n_cols/2-1) != 0 || genome(1,genome.n_cols/2-1) != 0 || genome(0,0) == 0 || genome(1,0) == 0){
	
		output = 0;
	}
	else{
		for(int loci = 0; loci < genome.n_cols/2;loci++){
		
			//maternal
			if(genome(0,loci*2) != 0){
				if(output > allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1))){
					output = allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1));
					
				}

			}
			//paternal
			if(genome(1,loci*2) != 0){
				if(output > allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1))){
					output = allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1));
					
				}

			}
	
		}
	}
	
	
	return output;
	
}


double fitness_calculation_mean(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency){


	double output = 0;
	int tot_loci = 0;
	
	
	if(genome(0,genome.n_cols/2) != 0 || genome(1,genome.n_cols/2) != 0){
	
		output = 0;
	}
	else{
		for(int loci = 0; loci < genome.n_cols/2;loci++){
		
			//maternal
			if(genome(0,loci*2) != 0){
			
				output = output + allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1));
				tot_loci++;

			}
			//paternal
			if(genome(1,loci*2) != 0){
			
				output = output + allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1));
				tot_loci++;

			}
	
		}
	
		if(tot_loci != 0){	
			output = output/((double)tot_loci);
			//cout << "Individual Fitness: " <<output << endl;
		}
	}
	return output;

}

double fitness_calculation_mean_chrom(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency){


	double output = 0;
	double output_1 = 0;
	double output_2 = 0;
	int tot_loci_1 = 0;
	int tot_loci_2 = 0;
	
	
	if(genome(0,genome.n_cols/2) != 0 || genome(1,genome.n_cols/2) != 0){
	
		output = 0;
	}
	else{
		for(int loci = 0; loci < genome.n_cols/2;loci++){
		
			//maternal
			if(genome(0,loci*2) != 0){
			
				output_1 = output_1 + allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1));
				tot_loci_1++;

			}
			//paternal
			if(genome(1,loci*2) != 0){
			
				output_2 = output_2 + allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1));
				tot_loci_2++;

			}
	
		}
	
		if(tot_loci_1 != 0 && tot_loci_2 != 0){	
			output_1 = output_1/((double)tot_loci_1);
			output_2 = output_2/((double)tot_loci_2);
			output = (output_1 + output_2)/2.0;
			
			
			/*
			cout << "Maternal Fitness: " <<output_1 << endl;
			cout << "Paternal Fitness: " <<output_2 << endl;
			cout << "Individual Fitness: " <<output << endl;
			*/
		}
	}
	
	
	return output;

}

double fitness_calculation_mean_no_freq(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency){


	double output = 0;
	int tot_loci = 0;
	
	
	if(genome(0,genome.n_cols/2) != 0 || genome(1,genome.n_cols/2) != 0){
	
		output = 0;
	}
	else{
		for(int loci = 0; loci < genome.n_cols/2;loci++){
		
			//maternal
			if(genome(0,loci*2) != 0){
			
				output = output + allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1);
				tot_loci++;

			}
			//paternal
			if(genome(1,loci*2) != 0){
			
				output = output + allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1);
				tot_loci++;

			}
	
		}
	
		if(tot_loci != 0){	
			output = output/((double)tot_loci);
		}
	}
	return output;

}

double fitness_calculation_same(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency){
	
	/*
	double output = 0;
	
	
	if(genome(0,genome.n_cols/2-1) != 0 || genome(1,genome.n_cols/2-1) != 0)
	{
	
		output = 0;
	}
	else
	{
	
		
		if(genome(0,0) == 0 || genome(1,0) == 0)
		{
	
			output = 0;
		}
		else
		{
	
			output = 1;
	
		}
		
	}
	*/
	
	
	double output = 1;
	
	if(genome(0,genome.n_cols/2-1) != 0 || genome(1,genome.n_cols/2-1) != 0){
	
		output = 0;
	}

	if(genome(0,0) == 0 || genome(1,0) == 0){

		output = 0;
	}
	
	
	
	return output;
}



double fitness_calculation_max_Lcost(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency){
	
	double output = 0;
	int NML = 0;
	int NFL = 0;
	
	double NL = 0;
	double NLmax = 20;
	
	for(int loci = 0; loci < genome.n_cols/2;loci++){
		
		//maternal
		if(genome(0,loci*2) != 0){
			NML++;
			if(output < allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1))){
				output = allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1));
			}

		}
		//paternal
		if(genome(1,loci*2) != 0){
			NFL++;
			if(output < allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1))){
				output = allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1));
			}

		}
	
	}
	NL = (NML+NFL)/2.0;
	output = output*(1 - NL/NLmax);
	return output;
}
double fitness_calculation_mean_Lcost(arma::imat& genome, arma::mat& allele_fitness, arma::mat& gene_frequency){

	double output = 0;
	int tot_loci = 0;
	
	int NML = 0;
	int NFL = 0;
	
	double NL = 0;
	double NLmax = 20;
	for(int loci = 0; loci < genome.n_cols/2;loci++){
		
		//maternal
		if(genome(0,loci*2) != 0){
			NML++;
			output = output + allele_fitness(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1)*(1 - gene_frequency(genome(0,loci*2)-1,genome(0,loci*2 + 1)-1));
			tot_loci++;

		}
		//paternal
		if(genome(1,loci*2) != 0){
			NFL++;
			output = output + allele_fitness(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1)*(1 - gene_frequency(genome(1,loci*2)-1,genome(1,loci*2 + 1)-1));
			tot_loci++;

		}
	
	}
	
	output = output/((double)tot_loci);
	NL = (NML+NFL)/2.0;
	output = output*(1 - NL/NLmax);
	return output;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//function for recombination of exons
arma::irowvec Recombination(arma::icube& genes, int loci, int parent, int gen, int inherited_chromosome){


	arma::uvec maternal_blank_space = find(genes.subcube(0+2*parent, 0, gen-1,    0+2*parent, loci*2 - 1, gen-1) == 0);
	arma::uvec paternal_blank_space = find(genes.subcube(1+2*parent, 0, gen-1,    1+2*parent, loci*2 - 1, gen-1) == 0);
	
	int maternal_break = rand()%((loci*2 - maternal_blank_space.size())/2);
	int paternal_break = rand()%((loci*2 - paternal_blank_space.size())/2);
	
	arma::irowvec new_maternal(2*loci);
	new_maternal.fill(0);
	arma::irowvec new_paternal(2*loci);
	new_paternal.fill(0);
	
	
	new_maternal.subvec(0, maternal_break*2) = ((genes.slice(gen-1)).row(0+2*parent)).subvec(0,maternal_break*2);
	
	new_maternal.subvec(maternal_break*2+1,maternal_break*2+1 + (loci*2 - paternal_blank_space.size() - paternal_break*2) - 2) = ((genes.slice(gen-1)).row(1+2*parent)).subvec(paternal_break*2+1, loci*2 - paternal_blank_space.size() - 1);
	
	new_paternal.subvec(0, paternal_break*2) = ((genes.slice(gen-1)).row(1+2*parent)).subvec(0,paternal_break*2);
	new_paternal.subvec(paternal_break*2+1,paternal_break*2+1 + (loci*2 - maternal_blank_space.size() - maternal_break*2) - 2) = ((genes.slice(gen-1)).row(0+2*parent)).subvec(maternal_break*2+1, loci*2 - maternal_blank_space.size() - 1);
	
	if(inherited_chromosome == 0){
		return new_maternal;
	}
	else{
		return new_paternal;
	}
	
	/*
	genes.subcube(0+2*parent, 0, gen-1,    1+2*parent, loci*2 - 1, gen-1).print();
	
	cout << maternal_break << endl;
	cout << paternal_break << endl;
	new_maternal.print();
	new_paternal.print();
	
	*/

}

arma::irowvec Recombination_2(arma::imat& genes_1,int loci, int parent, int inherited_chromosome){


	arma::uvec maternal_blank_space = find(genes_1.submat(0+2*parent, 0,    0+2*parent, loci*2 - 1) == 0);
	arma::uvec paternal_blank_space = find(genes_1.submat(1+2*parent, 0,    1+2*parent, loci*2 - 1) == 0);
	

	arma::irowvec new_maternal(2*loci);
	new_maternal.fill(0);
	arma::irowvec new_paternal(2*loci);
	new_paternal.fill(0);
	
	//Only goes through this if the chromosome has at least 1 locus
	if(loci*2 - maternal_blank_space.size() != 0 && loci*2 - paternal_blank_space.size()){
	
		int maternal_break = rand()%((loci*2 - maternal_blank_space.size())/2);
		int paternal_break = rand()%((loci*2 - paternal_blank_space.size())/2);
	
		new_maternal.subvec(0, maternal_break*2) = (genes_1.row(0+2*parent)).subvec(0,maternal_break*2);
	
		new_maternal.subvec(maternal_break*2+1,maternal_break*2+1 + (loci*2 - paternal_blank_space.size() - paternal_break*2) - 2) = (genes_1.row(1+2*parent)).subvec(paternal_break*2+1, loci*2 -paternal_blank_space.size() - 1);
	
		new_paternal.subvec(0, paternal_break*2) = (genes_1.row(1+2*parent)).subvec(0,paternal_break*2);
		new_paternal.subvec(paternal_break*2+1,paternal_break*2+1 + (loci*2 - maternal_blank_space.size() - maternal_break*2) - 2) = (genes_1.row(0+2*parent)).subvec(maternal_break*2+1, loci*2 - maternal_blank_space.size() - 1);
	
	}
	
	
	if(inherited_chromosome == 0){
		return new_maternal;
	}
	else{
		return new_paternal;
	}
	
	/*
	genes.subcube(0+2*parent, 0, gen-1,    1+2*parent, loci*2 - 1, gen-1).print();
	
	cout << maternal_break << endl;
	cout << paternal_break << endl;
	new_maternal.print();
	new_paternal.print();
	
	*/

}


arma::irowvec Allele_Recombination_2(arma::imat& genes_1,int loci, int parent, int inherited_chromosome){

	//cout << "Here 1" << endl;
	arma::uvec maternal_blank_space = find(genes_1.submat(0+2*parent, 0,    0+2*parent, loci*2 - 1) == 0);
	arma::uvec paternal_blank_space = find(genes_1.submat(1+2*parent, 0,    1+2*parent, loci*2 - 1) == 0);
	
	int maternal_break = rand()%((loci*2 - maternal_blank_space.size())/2 + 1);
	int paternal_break = rand()%((loci*2 - paternal_blank_space.size())/2 + 1);
	
	
	arma::irowvec new_maternal(2*loci);
	new_maternal.fill(0);
	arma::irowvec new_paternal(2*loci);
	new_paternal.fill(0);
	

	
	if(maternal_break != 0){
	
		new_maternal.subvec(0, maternal_break*2 - 1) = (genes_1.row(0+2*parent)).subvec(0,maternal_break*2 - 1);
	}
	
	new_maternal.subvec(maternal_break*2,maternal_break*2 + (loci*2 - paternal_blank_space.size() - paternal_break*2)+1) = (genes_1.row(1+2*parent)).subvec(paternal_break*2, loci*2 -paternal_blank_space.size() + 1);
	
	if(paternal_break != 0){
		new_paternal.subvec(0, paternal_break*2 - 1) = (genes_1.row(1+2*parent)).subvec(0,paternal_break*2 - 1);
	}
	
	new_paternal.subvec(paternal_break*2,paternal_break*2 + (loci*2 - maternal_blank_space.size() - maternal_break*2)+1) = (genes_1.row(0+2*parent)).subvec(maternal_break*2, loci*2 - maternal_blank_space.size() + 1);
	
	/*
	cout << maternal_break << endl;
	cout << paternal_break << endl;
	genes_1.submat(0+2*parent, 0,    0+2*parent, loci*2 - 1).print();
	genes_1.submat(1+2*parent, 0,    1+2*parent, loci*2 - 1).print();
	new_maternal.print();
	new_paternal.print();
	
	
	//exit(0);
	*/
	
	//cout << "Here 2" << endl;
	
	if(inherited_chromosome == 0){
		return new_maternal;
	}
	else{
		return new_paternal;
	}
	
	
	
	
	
	
	

}

arma::irowvec Exon_Mutation(arma::irowvec genes, int & exons, arma::mat& allele_fitness, arma::mat& gene_frequency){
	
	
	arma::irowvec output(genes.size());
	
	output = genes; 
	if(genes(0) != 0){
		//increasing the total amount of exons
		exons++;
	
		//resizing allele fitness and gene frequency matrix
	
		allele_fitness.resize(exons,exons);
		gene_frequency.resize(exons,exons);
	
		for(int i = 0; i < exons; i++){
		
			allele_fitness(i,exons-1) = rand()/((double)RAND_MAX);	
			allele_fitness(exons-1,i) = rand()/((double)RAND_MAX);	
	
		}
	
		arma::uvec gene_length = find(output != 0);
	
		int mutated_gene = rand()%(gene_length.size());
	
	
		output(mutated_gene) = exons;
	
	}
	return output;
	
}

arma::irowvec Exon_Mutation_Recycle(arma::irowvec genes, int & exons, arma::mat& allele_fitness, arma::mat& gene_frequency){
	
	//cout <<"Here1" << endl;
	arma::irowvec output(genes.size());
	
	output = genes;
	//output.print(); 
	if(genes(0) != 0){
	
		arma::uvec gene_length = find(output != 0);
	
		int mutated_gene = rand()%(gene_length.size());
	
	
		output(mutated_gene) = rand()%(exons)+1;
	
	}
	
	//output.print();
	
	//exit(0);
	
	//cout <<"Here2" << endl;
	return output;
	
}

arma::irowvec Gene_Duplication(arma::imat& genes, int loci, int parent, int inherited_chromosome){

	
	
	arma::uvec blank_space = find(genes.submat(inherited_chromosome +2*parent, 0,    inherited_chromosome+2*parent, loci*2 - 1) == 0);
	
	int duplicated_gene = rand()%((loci*2 - blank_space.size())/2);
	
	arma::irowvec new_chrom(2*loci);
	new_chrom = genes.submat(inherited_chromosome +2*parent, 0,    inherited_chromosome+2*parent, loci*2 - 1);
	
	new_chrom.subvec(loci*2 - blank_space.size(), loci*2 - blank_space.size() + 1) = (genes.row(inherited_chromosome+2*parent)).subvec(duplicated_gene*2 + 0,duplicated_gene*2 + 1);
	
	
	return new_chrom;
}


void Change_Allele_Fitness(arma::mat& allele_fitness){

	for(int exon1 = 0; exon1<allele_fitness.n_rows; exon1++){
		for(int exon2 = 0; exon2<allele_fitness.n_cols; exon2++){
		
			allele_fitness(exon1,exon2) = rand()/((double)RAND_MAX);
		}
	}



}

void Change_Proportion_Allele_Fitness(arma::mat& allele_fitness,double prop){
		
	
	for(int exon1 = 0; exon1<allele_fitness.n_rows; exon1++){
		for(int exon2 = 0; exon2<allele_fitness.n_cols; exon2++){
		
			if(rand()/((double)RAND_MAX) < prop){
				allele_fitness(exon1,exon2) = rand()/((double)RAND_MAX);
			}
		}
	}

}

void Change_Proportion_Allele_Fitness_Normal_Dist(arma::mat& allele_fitness,double prop,double sigma){
		
					
	for(int exon1 = 0; exon1<allele_fitness.n_rows; exon1++){
		for(int exon2 = 0; exon2<allele_fitness.n_cols; exon2++){
		
			if(rand()/((double)RAND_MAX) < prop){
				mt19937 mt(seedgen());
				std::normal_distribution<double> distribution(allele_fitness(exon1,exon2),sigma);
				
				double temp = 2;			

				while(temp < 0 | temp >1){
				
					temp = distribution(mt);
			
				}
				
				allele_fitness(exon1,exon2) = temp;
			}
		}
	}

}

void Change_Proportion_Allele_Fitness_Normal_Dist_Boundary_Force(arma::mat& allele_fitness,double prop,double sigma){
		
					
	for(int exon1 = 0; exon1<allele_fitness.n_rows; exon1++){
		for(int exon2 = 0; exon2<allele_fitness.n_cols; exon2++){
		
			if(rand()/((double)RAND_MAX) < prop){
				mt19937 mt(seedgen());
				std::normal_distribution<double> distribution(allele_fitness(exon1,exon2),sigma);
				
				double temp = 2;			
				temp = distribution(mt);
				
				
				if(temp > 1){
					temp = 1;
				}
				if(temp < 0){
					temp = 0;
				}
				
				allele_fitness(exon1,exon2) = temp;
			}
		}
	}

}

void gene_freq_calc(arma::imat& genes, arma::mat& gene_frequency){

	
	//cout << "Here 1" << endl;
	for(int person = 0; person < genes.n_rows/2; person++){
		map<pair<int,int>, int> allele_check;
		for(int locus = 0; locus < genes.n_cols/2; locus++){
			if(genes(person*2 + 0, locus*2 + 0) != 0){
				//chromosome 1 
				if(allele_check[make_pair(genes(person*2 + 0, locus*2 + 0),genes(person*2 + 0, locus*2 + 1))]){			
				}
				else{
					allele_check[make_pair(genes(person*2 + 0, locus*2 + 0),genes(person*2 + 0, locus*2 + 1))] = 1;
					gene_frequency(genes(person*2 + 0,0 + locus*2)-1,genes(person*2 + 0,1 + locus*2)-1) = gene_frequency(genes(person*2 + 0,0 + locus*2)-1,genes(person*2 + 0,1 + locus*2)-1) + 1.0;
				}
			}
			if(genes(person*2 + 1, locus*2 + 0) != 0){
				//chromosome 2
				if(allele_check[make_pair(genes(person*2 + 1, locus*2 + 0),genes(person*2 + 1, locus*2 + 1))]){			
				}
				else{
					allele_check[make_pair(genes(person*2 + 1, locus*2 + 0),genes(person*2 + 1, locus*2 + 1))] = 1;
					gene_frequency(genes(person*2 + 1,0 + locus*2)-1,genes(person*2 + 1,1 + locus*2)-1) = gene_frequency(genes(person*2 + 1,0 + locus*2)-1,genes(person*2 + 1,1 + locus*2)-1) + 1.0;
				}
			}
		}
	}
	
	//add one in the denominator to allow some some fitness if everyone has the same allele. 
	gene_frequency = gene_frequency/((double)genes.n_rows/2+1);
	
	//cout << "Here 2" << endl;

}



void gene_freq_calc_all_alleles(arma::imat& genes, arma::mat& gene_frequency){

	
	for(int person = 0; person < genes.n_rows/2; person++){
		for(int locus = 0; locus < genes.n_cols/2; locus++){
			if(genes(0+2*person,0+locus*2) != 0){
				gene_frequency(genes(0+2*person,0+locus*2)-1,genes(0+2*person,1+locus*2)-1) = gene_frequency(genes(0+2*person,0+locus*2)-1,genes(0+2*person,1+locus*2)-1) + 1.0;
			}
			if(genes(1+2*person,0+locus*2) != 0){
				gene_frequency(genes(1+2*person,0+locus*2)-1,genes(1+2*person,1+locus*2)-1) = gene_frequency(genes(1+2*person,0+locus*2)-1,genes(1+2*person,1+locus*2)-1) + 1.0;
			}
		}
	}
	arma::uvec allele_total1 = find(genes > 0);
	gene_frequency = 2*gene_frequency/((double)allele_total1.size());
}


long seedgen(void){

	long s, seed, pid;

	pid = getpid(); /* get process ID */
	//s = time ( &seconds ); /* get CPU seconds since 01/01/1970 */
	s = time (0);
	seed = abs(((s*181)*((pid-83)*359))%104729);	
	return seed;

}


















