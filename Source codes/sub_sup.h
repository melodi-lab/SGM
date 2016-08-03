//  1. READ INPUT
//  2. GENERATE BIPARTITE GRAPHS FOR EACH (PEPTIDE, OBSERVED SPECTRUM) PAIR
//  3. GET MAX WEIGHTED BIPARTITE MATCHING FOR EACH BIPARTITE GRAPH
//  4. GENERATE SCORES FOR EACH BIPARTITE GRAPH
//  5. OUTPUT SCORES TO FILE

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <numeric> 
#include <functional>
#include <algorithm>
#include <cmath>
#include "../graphs/graph_generator.h"

#include <cstdlib>
//#include "../fast_code/max_weight_bipartite_match_fast.h"
//#include "../fast_code/max_weight_bipartite_match_complete.h"
//#include "../fast_code/essentials.h"
//#include "../fast_code/somefunction.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;
class sub_sup{
	private:
		vector<int> A1;
		vector<int> A2;
		vector<int> A3;
		vector<float> W;
		int n;
		int k1;
		int k2;
		int n_A1;
		int n_A2;
		int n_all;
		float lambda_1;
		float lambda_2;
		float lambda_power;
		float ground_score1;
		float ground_score2;
		float s_2_0;
		vector<int>current_matching;
		vector<float>modular_bound1;
		vector<float>modular_bound2;
		vector<float>modular_bound3;
		vector<float>gain_empty;
		vector<float>gain_ground;
		vector<float> *test_weight1;
		vector<float> *test_weight2;
		int *test_pow1;
		int *test_pow2;
		int test_t;
		int test_n;
		int test;	
	public:
		float optimal_scores();
		float random_scores(int number_test=1);
		float diversity();
		void calulate_modular_bound(int modular);
		void calulate_modular_bound1();
		void calulate_modular_bound2();
		void calulate_modular_bound3();
		float submodualr_function(vector<int>matching,int modular=0,float projection_t=1);
		float submodualr_function1(vector<int>matching,float projection_t=1);
		float submodualr_function2(vector<int>matching,float projection_t=1);
		float scores(int modular);
		float submodualr_projection(float projection_t);
		float optimal_score();
		void output_average_modular_gain();

		float submodualar_function2_with_modular(vector<int>matching,int modular);
		float GreedyMax(float projection_t=1.0);
		float GreedyMax_lazy(int modular,float projection_t=1.0);
		float sub_sup_process(int modular);
		float sub_sup_process_max1();
		void print_current_matching();
		sub_sup(vector<int>A10, vector<int>A20, vector<float>W0, int k10, int k20,int n_A10, int n_A20, float lambda0,float lambda1,float lambda2);
};


