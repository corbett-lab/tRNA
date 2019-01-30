#ifndef __FUNCTION_TO_FITNESS_H
#define __FUNCTION_TO_FITNESS_H

////// OHTA FITNESS FUNCTION:
// 1 for ki >= kb
// 1 - s(kb - ki) for ki < kb
//// where ki is the number of functional alleles,
//// and kb is the average number of functional alleles in a random gamete
//////// ^ this is what we're using

//// TRUNCATION SELECTION NOTES (no longer using):
// There is no obvious way to decide how large d may become.
// Is it meaningful for d to be greater than 3-5 standard deviations?
// In most cases there is not information enough to decide.
// However, a natural extreme alternative to truncate selection
// is for fitness to be simply proportional to the character value,
// X. 
// c and d should be units of standard deviation
// c can be zero and d a nuisance parameter from 0 to 5
// 1 or 2 is often a good fit
// for us, X = the total function 
// need mean fitness m (is this desired mean or actual mean?)
// relative efficiency increases when the standard dev is smaller rel to the range of X
// truncation selection is instead 0 below a threshold and 1 above (similar to what we originally used)
// the idea is if w(x) < c or w(x) > c+d, then fitness = 0, else fitness = 1

void function_to_fitness ( double total_function_vector[], double fitness[], double mean_function, vector<individual> &population, cmd_line &options ) {


	//int fitness_count = 0 ;

	for ( int i = 0 ; i < population.size() ; i ++ ) {
		if (total_function_vector[i] >= mean_function) {
			fitness[i] = 1 ;
		}
		else{
			fitness[i] = 1 - ( options.selection_coefficient * ( mean_function - total_function_vector[i] ) ) ;
			if (fitness[i] < 0){
				fitness[i] = 0 ;
			}
		}
	}
	//cout << "fitness_count\t" << fitness_count << endl ;
}

#endif