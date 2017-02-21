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
#include "sub_sup.h"
#include <float.h>

using namespace std;
struct three_vaule_sub_sup {
	int a1;
	int a2;
	float w;
	float s;
	int index;
};
bool sortBys_sub_sup(const three_vaule_sub_sup &lhs, const three_vaule_sub_sup &rhs) {return lhs.s > rhs.s;}
int float_to_int(float x) {
	int i = (int) x;
	if ((x - (float)i) * (x - (float)i) > 0.0001) {
		cout << "erro float to int";
	}
	return i;
}
float sub_sup::concave_function(float x, float type, float lambda, bool inverse) {
	int t = float_to_int(type);
	if (t == 0) {
		if ( inverse == false) {
			return pow(x, lambda);
		}
		else {
			return pow(x, 1 / lambda);
		}
	} else if (t == 1) {
		if ( inverse == false) {
			float l = lambda / 10.0;
			return log(x / l + 1.0) * l;
		}
		else {
			float l = lambda / 10.0;
			if (x / l > log(FLT_MAX)) return FLT_MAX / 1000;

			float result1 = (exp(x / l) - 1.0) * l;
			if (result1 > FLT_MAX / 1000) return FLT_MAX / 1000;
			return result1;
		}
	} else if (t == 2) {
		if ( inverse == false) {

			return atan(x / lambda + 1.0) * lambda;
		}
		else {
			return (tan(x / lambda) - 1.0) * lambda;
		}

	}
	else if (t == 3) {

		if ( inverse == false) {
			float r1 = x;
			float d = lambda_vector[6];

			float r2 = pow(x, lambda) * power_d1;
			if (r1 > r2) return r2;
			else return r1;
		} else {
			float r1 = x;
			float d = lambda_vector[7];
			float r2 = 0;//pow(x,1/lambda)*power_d2;


			if (x < pow(FLT_MAX, lambda) &&  pow(x, 1 / lambda) < FLT_MAX / power_d1) {
				r2 = pow(x, 1 / lambda) * power_d2;
			} else {
				return FLT_MAX / 1000;
			}

			if (r2 > FLT_MAX / 1000) return FLT_MAX / 1000;

			if (r1 > r2) return r1;
			else return r2;

		}

	}
	cout << "erro no type";
	return 0;
}

sub_sup::sub_sup(vector<int>A10, vector<int>A20, vector<float>W0, vector<float> obs_mz_A10, vector<float> obs_mz_A20, vector<float> obs_intensity_A10, int isHighConf0,	int k10, int k20, int n_A10, int n_A20, vector<float>lambda_vector0) {
//sub_sup::sub_sup(vector<int>A10, vector<int>A20, vector<float>W0, int k10, int k20,int n_A10, int n_A20, vector<float>lambda_vector0){

	//lambda_vector0[0] lambda0
	//lambda_vector0[1] lambda1
	//lambda_vector0[2] concave_function
	//lambda_vector0[3] concave parameter
	//lambda_vector0[4] convex_function
	//lambda_vector0[5] convex_parameter
	lambda_vector = lambda_vector0;
	power_d1 = pow(lambda_vector[6], 1 - lambda_vector[3]);
	power_d2 = pow(lambda_vector[7], 1 - 1 / lambda_vector[5]);
	k1 = k10;
	k2 = k20;
	n_all = 1000000;
	//n_all=5;
	n = W0.size();
//creat test
	test = 0;
	test_t = 5;
	test_n = 15;
	if (test == 1) {
		k1 = 1000000;
		k2 = 1000000;
		n_all = 5;
		n = test_n;
		test_pow1 = new int[test_t];
		test_pow2 = new int[test_t];
		for (int i = 0; i < test_t; i++) {
			test_pow1[i] = ((float)rand()) / ((float)RAND_MAX);
			test_pow2[i] = ((float)rand()) / ((float)RAND_MAX) + 1.0;
		}
		test_weight1 = new vector<float>[test_t];
		test_weight2 = new vector<float>[test_t];
		for (int i = 0; i < test_t; i++) {
			for (int j = 0; j < test_n; j++) {
				test_weight1[i].push_back(((float)rand()) / ((float)RAND_MAX));
				test_weight2[i].push_back(((float)rand()) / ((float)RAND_MAX));
			}
		}
	}

	A1 = A10;
	A2 = A20;
	W = W0;


	obs_mz_A1 = obs_mz_A10;
	obs_mz_A2 = obs_mz_A20;
	obs_intensity_A1 = obs_intensity_A10;
	isHighConf = isHighConf0;
	n_A1 = n_A10;
	n_A2 = n_A20;
//	int *n_A1_side=new int[n_A1];
//
//	int **taken;
//	taken=new int*[n_A1];
//	for(int i=0;i<n_A1;i++){
//		n_A1_side[i]=0;
//		taken[i]=new int[n_A2];
//		for(int j=0;j<n_A2;j++){
//		taken[i][j]=0;
//		}
//	}
//	vector<int>A1_temp;
//	vector<int>A2_temp;
//	vector<float>W_temp;
//	for(int i=0;i<W.size();i++){
//	int	a1=A1[i];
//	int	a2=A2[i];
//		if (taken[a1][a2]==1) continue;
////		if (n_A1_side[a1]>=2) continue;
//
//		A1_temp.push_back(a1);
//		A2_temp.push_back(a2);
//		W_temp.push_back(W[i]);
//		taken[a1][a2]=1;
//		n_A1_side[a1]++;
//	}
//	delete taken;
//	delete n_A1_side;
//	A1=A1_temp;
//	A2=A2_temp;
//	W=W_temp;
//	n=W.size();
//	A3=vector<int>();
//	for(int i=0;i<A2.size();i++){
//		A3.push_back(rand()%n_A1);
//	}
//	cout<<"asdf\n";
	A3 = A2;
//	std::random_shuffle(A3.begin(),A3.end());
//	n=test_n;
//	cout<<n<<endl;
//if(n>30)cout<<n<<endl;
	lambda_1 = lambda_vector[0];

	lambda_2 = lambda_vector[1];
	lambda_power = 0;
//	lambda=0.5;
//	lambda=0.9;
//lambda=1.0;
//
	vector<int> ground_set;
	for (int i = 0; i < n; i++) {
		ground_set.push_back(i);
	}
	ground_score1 = submodualr_function1(ground_set);
	ground_score2 = submodualr_function2(ground_set);

	current_matching = std::vector<int>();
	s_2_0 = 0;
	for (int i = 0; i < n; i++) {
		s_2_0 += W[i];
	}
	if (0 == 1) {
		modular_bound1.assign(n, 0);
		modular_bound2.assign(n, 0);
		modular_bound3.assign(n, 0);
		gain_empty.assign(n, 0);
		gain_ground.assign(n, 0);
		for (int i = 0; i < n; i++) {
			vector<int> temp_matching;
			temp_matching.push_back(i);
			gain_empty[i] = submodualr_function2(temp_matching);
		}
		vector<int> ground_set;
		for (int i = 0; i < n; i++) {
			ground_set.push_back(i);
		}
		float score_ground = submodualr_function2(ground_set);
		for (int i = 0; i < n; i++) {
			vector<int> temp_matching = ground_set;
			temp_matching.erase(temp_matching.begin() + i);
			gain_ground[i] = score_ground - submodualr_function2(temp_matching);
		}
	}

}
void sub_sup::print_current_matching() {
	for (int i = 0; i < current_matching.size(); i++) {
		cout << current_matching[i] << " ";
	}
	cout << endl;
}
void sub_sup::calulate_modular_bound1() {
	float current_score = submodualr_function2(current_matching);
	vector<int>matching_order;
	matching_order.assign(n, 0);
	for (int i = 0; i < current_matching.size(); i++) {
		matching_order[i] = 1;
	}
	for (int i = 0; i < n; i++) {
		if (matching_order[i] == 1) continue;
		vector<int>temp_matching = current_matching;
		temp_matching.push_back(i);
		modular_bound1[i] = submodualr_function2(temp_matching) - current_score;
	}
	for (int i = 0; i < current_matching.size(); i++) {
		int j = current_matching[i];
		modular_bound1[j] = gain_ground[j];
	}
}
void sub_sup::calulate_modular_bound2() {
	modular_bound2 = gain_empty;
	float current_score = submodualr_function2(current_matching);
	for (int i = 0; i < current_matching.size(); i++) {
		int j = current_matching[i];
		vector<int>temp_matching = current_matching;
		temp_matching.erase(temp_matching.begin() + i);
		modular_bound2[j] = current_score - submodualr_function2(temp_matching);
	}

}

void sub_sup::calulate_modular_bound3() {
	modular_bound3 = gain_empty;
	for (int i = 0; i < current_matching.size(); i++) {
		int j = current_matching[i];
		modular_bound3[j] = gain_ground[j];
	}
}


float sub_sup::submodualar_function2_with_modular(vector<int>matching, int modular) {
	vector<float>modular_bound;
	if (modular == 1) {
		modular_bound = modular_bound1;
	} else if (modular == 2) {
		modular_bound = modular_bound2;
	} else if (modular == 3) {
		modular_bound = modular_bound3;
	} else {
		cout << "--------------ERROR\n";
	}
	float s = 0;
	for (int i = 0; i < matching.size(); i++) {
		s += modular_bound[matching[i]];
	}
	return s;

}
float sub_sup::submodualr_function1(vector<int>matching, float projection_t) {
//define submodular1
	if (test == 1) {
		float s = 0;
		for (int i = 0; i < test_t; i++) {
			float s1 = 0;

			for (int j = 0; j < matching.size(); j++) {
				int k = matching[j];
				s1 += test_weight1[i][k];
			}
			s += pow(s1, test_pow1[i] / projection_t);
		}
		return s;
	}
	float s1 = 0;
	float pow1 = 2.0;
	vector<float> s_a1;
	s_a1.assign(n_A1, 0);
	for (int j = 0; j < matching.size(); j++) {
		int i = matching[j];
//		if(i<=n/2) continue;
//		s_a1[A1[i]]+=pow(W[i],pow1);
//		s_a1[A1[i]]+=concave_function(W[i],lambda_vector[2],lambda_vector[3],true);
		s_a1[A1[i]] += W[i];
		//	s_a1[A1[i]]+=W[i];
//		s1+=W[i]*W[i];
		//	s1+=pow(W[i],2);

	}
//	return s1;
	float pow2 = 1.0 / pow1;
	for (int i = 0; i < n_A1; i++) {
//		s1+=concave_function(s_a1[i],lambda_vector[2],lambda_vector[3],false);
		s1 += s_a1[i];
	}

	return s1;
}
float sub_sup::submodualr_function2(vector<int>matching, float projection_t) {
	//define submodualr2
	if (test == 1) {
		//		cout<<"test\n";
		float s = 0;
		for (int i = 0; i < test_t; i++) {
			float s1 = 0;

			for (int j = 0; j < matching.size(); j++) {
				int k = matching[j];
				s1 += test_weight2[i][k];
			}
			s += pow(s1, test_pow2[i] / projection_t);
		}
		return s;
	}
	float lambda4 = lambda_power;
//	lambda4=2;
	float s3 = 0;
//	vector<float> s_a2;
	vector<float> s_a3;

	//s_a2.assign(n_A2/2,0);
	s_a3.assign(n_A2, 0);
	int found = 0;
	for (int j = 0; j < matching.size(); j++) {

		int i = matching[j];
		if (i == 22) found = 1;
		//	if(i>n/2) continue;
		int a3 = A2[i];
		a3 = A3[i];
		//	s_a3[a3]+=pow(W[i],1.0/lambda4);
		//s_a3[a3]+=concave_function(W[i],lambda_vector[4],lambda_vector[5],false);
//	s_a3[a3]+=sqrt(W[i]);
		s_a3[a3] += (W[i]);
//		s_a3[a3]+=W[i];
	}
	//for(int i = 0; i<matching.size();i++){
	//	cout<<matching[i]<<"--";
	//}
	//cout<<endl;
	//for(int j = 0; j<n_A2;j++){
	//	cout<<s_a3[j]<<" ";

	//}
	//cout<<endl;
	for (int i = 0; i < n_A2 / 2; i++) {
		//	s3+=pow(s_a3[i]+s_a3[n_A2/2+i],lambda4/projection_t);

//		s3+=concave_function(s_a3[i]+s_a3[n_A2/2+i],lambda_vector[4],lambda_vector[5],true);
//		s3+=concave_function(s_a3[i]+s_a3[n_A2/2+i],lambda_vector[4],lambda_vector[5],true);
		s3 += (s_a3[i] + s_a3[n_A2 / 2 + i]) * (s_a3[i] + s_a3[n_A2 / 2 + i]);
		//s3+=exp((s_a3[i]+s_a3[n_A2/2+i]));
		//s3+=exp(s_a3[i]+s_a3[n_A2/2+i])-1;
	}
	//if (found) cout<<s3<<endl;
	return s3;
}

float sub_sup::submodualr_function(vector<int>matching, int modular, float projection_t) {
	//if (modular==0){

	//return lambda*submodualr_function1(matching,projection_t)+(1-lambda)*submodualr_function2(matching,projection_t);
	float f1 = submodualr_function1(matching, projection_t);
	float f2 = submodualr_function2(matching, projection_t);
//	return f2;
//	return f1;
//	if (lambda_2 < 0.0001) return f2;
//	if (lambda_2 > 0.9999) return f1


	//if (lambda_1>0.95 &&lambda_2>0.95) return f1;
	if (ground_score1 < 0.000001) {
		return 0;
	}
	if (ground_score2 < 0.000001) {
		return 0;
	}
	float ratio1 = f1 / ground_score1;
	float ratio2 = f2 / ground_score2;
	float average1 = 1.10;
	float average2 = 1.28;
	float ss1 = lambda_1 * ground_score1 / average1 + (1 - lambda_1) * ground_score2 / average2;

	float ss2 = lambda_2 * ratio1 + (1 - lambda_2) * ratio2;
	return ss2 * ss1;
//	}
//	return lambda_1*submodualr_function1(matching)+(1-lambda_1)*submodualar_function2_with_modular(matching,modular);



}
float sub_sup::GreedyMax_lazy(int modular, float projection_t) {
	vector<three_vaule_sub_sup> all_W;
	current_matching.assign(0, 0);
	vector<int>matching_order;
	matching_order.assign(n, 0);

	for (int i = 0; i < W.size(); i++) {
		three_vaule_sub_sup one_W;
		one_W.a1 = A1[i];
		one_W.a2 = A2[i];
		one_W.w = W[i];
		one_W.index = i;
		vector<int> matching;
		matching.push_back(i);
		one_W.s = submodualr_function(matching, modular, projection_t);
		all_W.push_back(one_W);
	}
	sort(all_W.begin(), all_W.end(), sortBys_sub_sup);
	if (all_W.size() > 2) {
		if (all_W[0].w < all_W[1].w) {
			cout << "errors";
		}
	}
	float current_score = 0;
	vector<int> b1;
	vector<int> b2;
	b1.assign(n_A1, 0);
	b2.assign(n_A2, 0);
	while (true) {
		if (all_W.size() == 0) break;

		int i = all_W[0].index;

		if (matching_order[i] == 1 || (b1[A1[i]] >= k1) || (b2[A2[i]] >= k2) ) {
			all_W.erase(all_W.begin());
			continue;
		}
		if (all_W.size() == 1) {
			current_matching.push_back(i);
			break;
		}
		if (all_W[0].s < all_W[1].s) {
			for (int i = 0; i < all_W.size(); i++) {
				vector<int>temp_matching = current_matching;
				temp_matching.push_back(all_W[i].index);
				float s0 = submodualr_function(temp_matching, modular, projection_t) - current_score;
				all_W[i].s = s0;
			}

			sort(all_W.begin(), all_W.end(), sortBys_sub_sup);
			continue;
		}
		int temp_i = i;
		current_matching.push_back(temp_i);
		if (current_matching.size() >= n_all) break;
		if (current_matching.size() >= k2 * n_A2) break;
		matching_order[temp_i] = 1;
		b1[A1[temp_i]]++;
		b2[A2[temp_i]]++;
		all_W.erase(all_W.begin());
	}
	return submodualr_function(current_matching, modular);

}
float sub_sup::GreedyMax(float projection_t, std::vector<int> intial_set) {
	//vector<three_vaule_sub_sup> all_W;
	//current_matching.assign(0,0);
	//for(int i=0;i<W.size();i++){
	//	three_vaule_sub_sup one_W;
	//	one_W.a1=A1[i];
	//	one_W.a2=A2[i];
	//	one_W.w=W[i];
	//	one_W.index=i;
	//	vector<int> matching;
	//	matching.push_back(i);
	//	one_W.s=submodualr_function(matching);
	//	all_W.push_back(one_W);
	//}
	//sort(all_W.begin(),all_W.end(),sortBys_sub_sup);
	//float current_score=0;
	//vector<int> b1;
	//vector<int> b2;
	//b1.assign(n_A1,0);
	//b2.assign(n_A2,0);
	//while(true){
	//	if(all_W[0].s<all_W[1].s){
	//		for(int i=0;i<all_W.size();i++){
	//			vector<int>temp_matching=current_matching;
	//			temp_matching.push_back(all_W[i].index);
	//			float s0=submodualr_function(temp_matching)-current_score;
	//			all_W[i].s=s0;
	//		}
	//
	//	sort(all_W.begin(),all_W.end(),sortBys_sub_sup);
	//	continue;
	//	}
	//
	//
	//}

	//current_matching.assign(0,0);
	current_matching = intial_set;
	vector<float> greedy_scores;
	vector<float> intial_gain;
	vector<int>matching_order;
	matching_order.assign(n, 0);
	float current_score = 0;
	vector<int> b1;
	vector<int> b2;
	b1.assign(n_A1, 0);
	b2.assign(n_A2 / 2, 0);
	int max_i = n_A1 * k1;
	if (max_i > n_A2 * k2 / 2) max_i = n_A2 * k2 / 2;
	max_i = 200;

	//initialize
	for (int i = 0; i < current_matching.size(); i++) {
		int j = current_matching[i];
		matching_order[j] = 1;
		b1[A1[j]] += 1;
		int a_2 = A2[j];
		if (a_2 % 2 == 0) {a_2 = a_2 / 2;}
		else {
			a_2 = a_2 - 1;
			a_2 = a_2 / 2;
		}
		//if(b2[A2[i]]>=k2) continue;
		b2[a_2] += 1;
		max_i -= 1;


	}

	for (int kk = 0; kk < max_i; kk++) {
		//	cout<<kk<<endl;
//		greedy_scores.push_back(0);
//		intial_gain.push_back(0);
//
		float temp_s = -10000000;
		float temp_i = 0;
		bool seleted = 0;
		for (int i = 0; i < n; i++) {

			if (matching_order[i] == 1) continue;
			if (b1[A1[i]] >= k1) continue;
			int a_2 = A2[i];
			if (a_2 % 2 == 0) {a_2 = a_2 / 2;}
			else {
				a_2 = a_2 - 1;
				a_2 = a_2 / 2;
			}
			//if(b2[A2[i]]>=k2) continue;
			if (b2[a_2] >= k2) continue;

			if (W[i] > 2000) cout << i << " " << W[i] << endl;
			vector<int>temp_matching = current_matching;
			temp_matching.push_back(i);
			float s0 = submodualr_function(temp_matching, 0, projection_t);
//	float		s0 = W[i];
			if (s0 > temp_s) {
				temp_s = s0;
				temp_i = i;
				seleted = 1;
				//			greedy_scores[greedy_scores.size()-1]=s0;
				vector<int>one_matching;
				one_matching.push_back(temp_i);
//				intial_gain[intial_gain.size()-1]=submodualr_function(one_matching,0,projection_t);
			}
		}
		if (seleted == 0) break;

		current_matching.push_back(temp_i);
		matching_order[temp_i] = 1;
		b1[A1[temp_i]]++;
		int a_2 = A2[temp_i];
		if (a_2 % 2 == 0) {a_2 = a_2 / 2;}
		else {
			a_2 = a_2 - 1;
			a_2 = a_2 / 2;
		}

		b2[a_2]++;
		if (kk >= n_all) break;
	}
//	if(greedy_scores.size()>40){
//	for(int i=0;i<greedy_scores.size();i++){
//		cout<<greedy_scores[i]<<" ";
//
//	}
//	cout<<endl;
//	for(int i=0;i<greedy_scores.size();i++){
//		cout<<intial_gain[i]<<" ";
//
//	}
//	cout<<endl;
//	cout<<"----------------"<<endl;
//	}
	return submodualr_function(current_matching, 0, 1.0);

}
double sub_sup::postprocessing3() {
	double current_score = GreedyMax();
	for (int i = 0; i < n; i++) {
		std::vector<int> temp_current_matching_backup(0, 0);
		temp_current_matching_backup.push_back(i);

		double s = GreedyMax(1.0, temp_current_matching_backup);
		//double s = 0;
		if (s > current_score) {
			sub_sup::counter ++;
			cout << sub_sup::counter << " " << s << " " << current_score << endl;
			return s;
		}

	}
	return GreedyMax();

}
double sub_sup::postprocessing2() {
	double current_score = GreedyMax();
	if (current_matching.size() < 6) return current_score;
	std::vector<int> current_matching_backup = current_matching;

	for (int i = 0; i < 100; i++) {
		std::vector<int> temp_current_matching_backup = current_matching;

		temp_current_matching_backup.erase(temp_current_matching_backup.begin() + rand() % temp_current_matching_backup.size());
		temp_current_matching_backup.erase(temp_current_matching_backup.begin() + rand() % temp_current_matching_backup.size());
		temp_current_matching_backup.erase(temp_current_matching_backup.begin() + rand() % temp_current_matching_backup.size());
		temp_current_matching_backup.erase(temp_current_matching_backup.begin() + rand() % temp_current_matching_backup.size());
		temp_current_matching_backup.erase(temp_current_matching_backup.begin() + rand() % temp_current_matching_backup.size());
		double s = GreedyMax(1.0, temp_current_matching_backup);
		//double s = 0;
		if (s > current_score) {
			sub_sup::counter ++;
			cout << sub_sup::counter << " " << s << " " << current_score << endl;
			return s;
		}
		current_matching = current_matching_backup;
	}
	return current_score;

}

double sub_sup::postprocessing() {
	double current_score = GreedyMax();
	int current_size = current_matching.size();
	int number_matched = number_matched_byion_pair();
	std::vector<int> b2;
	b2.assign(n_A2 / 2, 0);
	std::vector<int> current_selected(n, 0);
	for (int i = 0; i < current_matching.size(); i++) {
		int j = current_matching[i];
		current_selected[current_matching[i]] = 1;
		int a_2 = A2[j];

		if (a_2 % 2 == 0) {a_2 = a_2 / 2;}
		else {
			a_2 = a_2 - 1;
			a_2 = a_2 / 2;
		}
		b2[a_2] += 1;
	}
	int nums = 0;
	std::vector<int> temp_current_matching;
	std::vector<int> temp_current_matching_backup = current_matching;
	double temp_current_max = current_score;
	for (int i = 0; i < b2.size() / 2; i++) {
		if (b2[i] == 0 || b2[b2.size() / 2 + i] == 0) {
			std::vector<int> temp_m = current_matching;
			//which means b-i and y-i are not selected.
			//we can try to select b-i and y-i, then reject two edges.
			//four steps. select one b-i, select one y-i, reject one, reject one.
			//1
			int cout = 0;
			double temp_max = -10000;
			int temp_j1 = -1;
			int temp_j2 = -1;
			int temp_j3 = -1;
			for (int j1 = 0; j1 < n; j1++) {
				if (current_selected[j1] == 1) continue;
				//if (2 * i != A2[j1] && 2 * i + 1 != A2[j1] && 2 * i + b2.size() != A2[j1] && 2 * i + 1 + b2.size() != A2[j1]) continue;
				if (A2[j1] / 2 != i && A2[j1] / 2 != i + b2.size() / 2) continue;
				if (b2[A2[j1] / 2] >= k2) continue;
				for (int j2 = j1 + 1; j2 < n; j2++) {
					if (current_selected[j2] == 1) continue;
					if (A2[j2] / 2 != i && A2[j2] / 2 != i + b2.size() / 2) continue;
					//if (b2[A2[j2]/2] + 1 == k2 + 1) continue;
					if (b2[A2[j2] / 2] + (A2[j2] / 2 == A2[j1] / 2) >= k2) continue;
					//if (A2[j2]/2 == A2[j1]/2) continue;
					//if (2 * i != A2[j2] && 2 * i + 1 != A2[j2] && 2 * i + b2.size() != A2[j2] && 2 * i + 1 + b2.size() != A2[j2]) continue;
					for (int j3 = j2 + 1; j3 < n; j3++) {
						if (current_selected[j3] == 1) continue;
						if (A2[j3] / 2 != i && A2[j3] / 2 != i + b2.size() / 2) continue;
						//if (2 * i != A2[j3] && 2 * i + 1 != A2[j3] && 2 * i + b2.size() != A2[j3] && 2 * i + 1 + b2.size() != A2[j3]) continue;
						if (b2[A2[j3] / 2] + (A2[j3] / 2 == A2[j1] / 2) + (A2[j3] / 2 == A2[j2] / 2) >= k2) continue;
						cout ++;
						vector<int> temp_matching = current_matching;
						temp_matching.push_back(j1);
						temp_matching.push_back(j2);
						temp_matching.push_back(j3);
						double temp_s = submodualr_function(temp_matching, 0, 1.0);
						if (temp_s > temp_max) {
							temp_max = temp_s;
							temp_j1 = j1;
							temp_j2 = j2;
							temp_j3 = j3;
						}
					}
				}
			}
			if (cout == 0) { current_matching = temp_m; continue;}
			cout = 0;
			// if (temp_max < current_score) return current_score;
			current_matching.push_back(temp_j1);
			current_matching.push_back(temp_j2);
			current_matching.push_back(temp_j3);

			// 3
			temp_max = -10000;
			temp_j1 = -1;
			temp_j2 = -1;
			temp_j3 = -1;
			for (int j1 = 0; j1 < current_matching.size(); j1++) {
				for (int j2 = j1 + 1; j2 < current_matching.size(); j2++) {
					for (int j3 = j2 + 1; j3 < current_matching.size(); j3++) {
						vector<int> temp_matching = current_matching;
						temp_matching.erase(temp_matching.begin() + j3);
						temp_matching.erase(temp_matching.begin() + j2);
						temp_matching.erase(temp_matching.begin() + j1);
						double temp_s = submodualr_function(temp_matching, 0, 1.0);
						if (temp_s > temp_max) {
							temp_max = temp_s;
							temp_j1 = j1;
							temp_j2 = j2;
							temp_j3 = j3;
						}
					}
				}
			}
			// if (temp_max < current_score) { current_matching = temp_m;continue;}
			current_matching.erase(current_matching.begin() + temp_j3);
			current_matching.erase(current_matching.begin() + temp_j2);
			current_matching.erase(current_matching.begin() + temp_j1);

			double temp_d = submodualr_function(current_matching, 0, 1.0);
			if (temp_current_max < temp_d) {
				// cout<<"improved"<<current_matching.size()<<" "<<current_size<<endl;
				// cout<<number_matched_byion_pair()<<" "<< number_matched<< endl;
				temp_current_max = temp_d;
				temp_current_matching = current_matching;
			}
			current_matching = temp_m;
		}
	}
	current_matching = temp_current_matching;

	b2.assign(n_A2 / 2, 0);
	int max_i = n_A1 * k1;
	if (max_i > n_A2 * k2 / 2) max_i = n_A2 * k2 / 2;
	max_i = 20;
	if (current_matching.size() > current_size) cout << "error--------\n";
	for (int i = 0; i < current_matching.size(); i++) {
		int j = current_matching[i];
		//matching_order[j] = 1;
		//b1[A1[j]] += 1;
		int a_2 = A2[j];
		if (a_2 % 2 == 0) {a_2 = a_2 / 2;}
		else {
			a_2 = a_2 - 1;
			a_2 = a_2 / 2;
		}
		//if(b2[A2[i]]>=k2) continue;
		b2[a_2] += 1;
		max_i -= 1;
		if (b2[a_2] > k2) cout << "error--------\n";


	}


	if (temp_current_max > current_score) {
		cout << temp_current_max << " " << current_score << endl;
		cout << current_matching.size() << " " << current_size << endl;
		for (int i = 0; i < b2.size(); i++) {
			cout << b2[i] << " ";
		}
		cout << endl;
		for (int i = 0; i < temp_current_matching_backup.size(); i++) {
			cout << temp_current_matching_backup[i] << " ";
		}
		cout << endl;
		for (int i = 0; i < current_matching.size(); i++) {
			cout << current_matching[i] << " ";
		}
		cout << endl;
		cout << number_matched_byion_pair() << " " << number_matched << endl;
		return temp_current_max;
	}
	else return current_score;

	// vector<int> b1;


	//return submodualr_function(current_matching,0,1.0);
}



void sub_sup::calulate_modular_bound(int modular) {
	if (modular == 1) {
		calulate_modular_bound1();
	} else if (modular == 2) {
		calulate_modular_bound2();
	} else if (modular == 3) {
		calulate_modular_bound3();
	} else {
		cout << "--------------ERROR\n";
	}

}
float sub_sup::submodualr_projection(float projection_t) {
	return GreedyMax(projection_t);
//	return GreedyMax_lazy(0,projection_t);


}
float sub_sup::sub_sup_process_max1() {
	float a1 = sub_sup_process(1);
	float a2 = sub_sup_process(2);
	float a3 = sub_sup_process(3);
	float r = 0;
	if (a1 > r) r = a1;
	if (a2 > r) r = a2;
	if (a3 > r) r = a3;
	return r;

}

float sub_sup::sub_sup_process(int modular) {
	vector<int>current_matching_back = current_matching;
	vector<int>current_matching_back2 = current_matching;
	vector<int>current_matching_back3 = current_matching;
	float s_max = 0;

	vector<int>empty_set;
//	current_matching=empty_set;
	int	modular2 = modular;
	modular2 = 0;
	float current_score1 = submodualr_function(current_matching, modular2);
	float current_score2 = 0;
//	return current_score1;
	int i = 0;
	while (true) {
		//	cout<<"iteration "<<i<<endl;
		i = i + 1;
		current_matching_back2 = current_matching;
		if (modular <= 3) {
			calulate_modular_bound(modular);
		} else {

			calulate_modular_bound(1);
			calulate_modular_bound(2);
			calulate_modular_bound(3);
		}
		float temp_score = 0;
		current_score2 = 0;
		if (modular <= 3) {
			temp_score = GreedyMax_lazy(modular, 1.0);
			current_score2 = submodualr_function(current_matching, modular2);
			current_matching_back3 = current_matching;
		}
		else {
			temp_score = GreedyMax_lazy(1, 1.0);
			float a1 = submodualr_function(current_matching, modular2);
			if (a1 > current_score2) {
				current_score2 = a1;
				current_matching_back3 = current_matching;
			}

			temp_score = GreedyMax_lazy(2, 1.0);
			float a2 = submodualr_function(current_matching, modular2);
			if (a2 > current_score2) {
				current_score2 = a2;
				current_matching_back3 = current_matching;
			}

			temp_score = GreedyMax_lazy(3, 1.0);
			float a3 = submodualr_function(current_matching, modular2);
			if (a3 > current_score2) {
				current_score2 = a3;
				current_matching_back3 = current_matching;
			}

		}
		//	current_score2=submodualr_function(current_matching,modular2);
//		if (temp_score<current_score-1e-7) cout<<"error--------\n";
//		if (temp_score<current_score+1e-7) break;
		if (current_score2 <= current_score1 + 1e-5) break;
		current_score1 = current_score2;
		current_matching = current_matching_back3;

	}

	float s = submodualr_function(current_matching_back2, 0);
//	print_current_matching();
	current_matching = current_matching_back;
	return s;
//	return submodualr_function(current_matching,0);
}
float sub_sup::scores(int modular) {

	return GreedyMax();
	modular = 0;
	if (modular == 0) return GreedyMax_lazy(0);
	GreedyMax_lazy(0);
	return sub_sup_process(2);
}
float sub_sup::optimal_scores() {
	float max_s = -10000;

	int m = pow(2, n);
	int *taken = new int[m];
	for (int i = 0; i < m; i++) taken[i] = 0;
	cout << "n: " << n << endl;
	for (int i = 0; i < pow(2, n); i++) {

		if (i % (m / 10) == 0) cout << i << endl;
		if (i >= 1) {
			double a = i;

			a = floor(log(a) / log(2.0));
			int a2 = (int)round(pow(2.0, a));
			//	cout<<i<<" "<<a2<<endl;
			if (taken[i - a2] == 1) {
				taken[i] = 1;

				continue;
			}
		}


		vector<float>b1(n_A1, 0);
		vector<float>b2(n_A2 / 2, 0);

		vector<int>matching;
		int t = i;
		for (int j = 0; j < n; j++) {
			if (t % 2 == 1) matching.push_back(j);
			t = t / 2;
			if (t == 0) break;
		}
		int p = 0;
		if (matching.size() > n_all) {taken[i] = 1; continue;}
		for (int k = 0; k < matching.size(); k++) {
			int j = matching[k];
			if (b1[A1[j]] >= k1) {p = 1; continue;}
			//if(b2[A2[j]]>=k2) {p=1;continue;}
			int a_2 = A2[j];
			if (a_2 % 2 == 0) {a_2 = a_2 / 2;}
			else {
				a_2 = a_2 - 1;
				a_2 = a_2 / 2;
			}

			if (b2[a_2] >= k2) {p = 1; continue;}
			b1[A1[j]]++;
			//b2[A2[j]]++;
			b2[a_2]++;
		}
		if (p == 1) {taken[i] = 1; continue;}
		float s = submodualr_function(matching);
		if (s > max_s) max_s = s;
	}
	return max_s;
}
float sub_sup::random_scores(int number_test) {
	int N = number_test;
	vector<float>result;
	for (int i = 0; i < N; i++) {
		vector<int>matching;
		vector<float>b1(n_A1, 0);
		vector<float>b2(n_A2, 0);
		int count = 0;
		vector<int>chosen(n, 0);
		while (true) {
			int j = rand() % n;
			if (chosen[j] == 1) continue;
			count++;
			if (count > n * 100) break;
			if (matching.size() >= n_all)continue;
			if (b1[A1[j]] >= k1) continue;
			if (b2[A2[j]] >= k2) continue;
			b1[A1[j]]++;
			b2[A2[j]]++;
			matching.push_back(j);
			count = 0;
			chosen[j] = 1;

		}
		result.push_back(submodualr_function(matching));
	}
	float mean = 0;
	for (int i = 0; i < N; i++) {
		mean += result[i] / ((float)N);
	}
	return mean;
}
float sub_sup::diversity() {
	vector<float> v = W;
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double mean = sum / v.size();

	double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
	return stdev;
}
void sub_sup::output_average_modular_gain() {
	float s1 = 0;
	float s2 = 0;
	vector<int> matching;
	for (int i = 0; i < n; i++) {
		matching = vector<int>();
		matching.push_back(i);
		s1 += submodualr_function1(matching);
		s2 += submodualr_function2(matching);
	}
	cout << s1 << " " << s2 << endl;


}
void sub_sup::write_all_bipartite(int sid) {
//	if (isHighConf == 0) return;
	vector<int>matching_order2;
	matching_order2.assign(n, 0);
	for (int i = 0; i < current_matching.size(); i++) {
		matching_order2[current_matching[i]] = 1;
	}

	std::ofstream myfile;
	myfile.open("/n/trombone/s1/wrbai/codes/matlab_codes/2016_12_30/matchings_one.txt", std::ios::app);
	//cout<<A1.size()<<" "<<A2.size()<<" "<<W.size()<<" "<<obs_mz_A1.size()<<" "<<obs_intensity_A1.size()<<" "<<obs_mz_A2.size()<<" "<<matching_order2.size()<<endl;
	//
	myfile << -n_A2 << "\t" << sid << "\t-1\t-1\t-1\t-1\t-1\n";
	for (int i = 0; i < A1.size(); i++) {
		myfile << A1[i] << "\t" << A2[i] << "\t" << W[i] << "\t" << obs_mz_A1[i] << "\t" << obs_intensity_A1[i] << "\t" << obs_mz_A2[i] << "\t" << matching_order2[i] << endl;
	}

}

void sub_sup::write_all_parameters(int sid) {

	// vector<int> A10(A1);
//   	vector<int> A20(A2);
	// vector<float> W0(W);
//   	vector<float> obs_mz_A10(obs_mz_A1);
//   	vector<float> obs_mz_A20(obs_mz_A2);
//   	vector<float> obs_intensity_A10(obs_intensity_A1);
//   	int isHighConf0(isHighConf);
	// int k10 = k1;
//   	int k20 = k2;
	// int n_A10 = n_A1;
//   	int n_A20 = n_A2;
	// vector<float>lambda_vector0(lambda_vector);
	std::ofstream myfile;
	myfile.open("/n/trombone/s1/wrbai/codes/matlab_codes/2016_12_30/input_for_subsup.txt", std::ios::app);
	myfile << "sid: " << sid << endl;

	myfile << "static const int arr1[] = {";
	for (int i = 0; i < A1.size() - 1; i++) {
		myfile << A1[i] << ",";
	}
	myfile << A1[A1.size() - 1] << "};\n";


	myfile << "static const int arr2[] = {";
	for (int i = 0; i < A2.size() - 1; i++) {
		myfile << A2[i] << ",";
	}
	myfile << A2[A2.size() - 1] << "};\n";


	myfile << "static const float arr3[] = {";
	for (int i = 0; i < W.size() - 1; i++) {
		myfile << W[i] << ",";
	}
	myfile << W[W.size() - 1] << "};\n";



	myfile << "static const float arr4[] = {";
	for (int i = 0; i < obs_mz_A1.size() - 1; i++) {
		myfile << obs_mz_A1[i] << ",";
	}
	myfile << obs_mz_A1[obs_mz_A1.size() - 1] << "};\n";



	myfile << "static const float arr5[] = {";
	for (int i = 0; i < obs_mz_A2.size() - 1; i++) {
		myfile << obs_mz_A2[i] << ",";
	}
	myfile << obs_mz_A2[obs_mz_A2.size() - 1] << "};\n";



	myfile << "static const float arr6[] = {";
	for (int i = 0; i < obs_intensity_A1.size() - 1; i++) {
		myfile << obs_intensity_A1[i] << ",";
	}
	myfile << obs_intensity_A1[obs_intensity_A1.size() - 1] << "};\n";



	myfile << "static const float arr7[] = {";
	for (int i = 0; i < lambda_vector.size() - 1; i++) {
		myfile << lambda_vector[i] << ",";
	}
	myfile << lambda_vector[lambda_vector.size() - 1] << "};\n";

	myfile << "int isHighConf0 = " << isHighConf << ";\n";
	myfile << "int k10 = " << k1 << ";\n";
	myfile << "int k20 = " << k2 << ";\n";
	myfile << "int n_A10 = " << n_A1 << ";\n";
	myfile << "int n_A20 = " << n_A2 << ";\n";
//	myfile<<isHighConf<<endl<<k1<<endl<<k2<<endl<<n_A1<<endl<<n_A2<<endl;
}


float sub_sup::GreedyMax_initialByOne(int sid) {
	float scores = GreedyMax();
	float scores0 = scores;
	int number_matched = number_matched_byion_pair();
	int number_matched0 = 0;
	int kk = 0;
	//std::std::vector<int> current_matching0 = current_matching;
	std::vector<int> M(n, 0);
	for (int i = 0; i < current_matching.size(); i++) {
		M[current_matching[i]] = 1;
	}

	float temp_score = 0;
	int changed = 0;
	for (int i = 0; i < M.size(); i++) {
		if (M[i] == 0) {
			std::vector<int> intial_set(1, i);
			temp_score = GreedyMax(1.0, intial_set);
			if (temp_score > scores) {
				//if (number_matched < number_matched_byion_pair()){
				//	cout<<sid<<" "<<number_matched <<" "<< number_matched_byion_pair() <<endl;
				//}
				kk = i;
				scores = temp_score;
				number_matched0 = number_matched_byion_pair();
				changed = 1;
			}


		}

	}
	if (changed == 1) {
		if (number_matched < number_matched0) {

			cout << sid << " " << kk << " " << scores0 << " " << scores << endl;
			cout << number_matched << " " << number_matched0 << endl;
		}

		// cout<<number_matched<<endl;
		// std::vector<int> b2;
		// b2.assign(n_A2/2,0);
		// for(int i = 0; i < current_matching.size(); i++){
		// int j = current_matching[i];
		// int a_2=A2[j];
		// if (a_2%2==0) {a_2=a_2/2;}
		// else{
		// 		a_2=a_2-1;
		// 		a_2=a_2/2;
		// }
		// b2[a_2] += 1;
		// }
		// for(int i = 0; i <b2.size()/2; i++){
		// 	cout<<b2[i] <<" " << b2[b2.size()/2 + i]<<endl;
		// }

	}



}
int sub_sup::number_matched_byion_pair() {
	std::vector<int> b2;
	b2.assign(n_A2 / 2, 0);
	for (int i = 0; i < current_matching.size(); i++) {
		int j = current_matching[i];
		int a_2 = A2[j];
		if (a_2 % 2 == 0) {a_2 = a_2 / 2;}
		else {
			a_2 = a_2 - 1;
			a_2 = a_2 / 2;
		}
		b2[a_2] += 1;
	}
	int nums = 0;
	for (int i = 0; i < b2.size() / 2; i++) {
		if (b2[i] >= 1 && b2[b2.size() / 2 + i] >= 1) {
			nums++;
		}

	}
	return nums;



}
int sub_sup::counter = 0;
// sub_sup* read_all_input(){
// 	vector<int> A10;
//    	vector<int> A20;
// 	vector<float> W0;
//    	vector<float> obs_mz_A10;
//    	vector<float> obs_mz_A20;
//    	vector<float> obs_intensity_A10;
//    	int isHighConf0;
// 	int k10 ;
//    	int k20 ;
// 	int n_A10 ;
//    	int n_A20 ;
// 	vector<float>lambda_vector0;
// 	std::ifstream myfile;
// 	myfile.open("/n/trombone/s1/wrbai/codes/matlab_codes/input_for_subsup.txt",std::ios::app);

// 	std::string line;
// 	i = 0;
// while (std::getline(myfile, line))
// {
//     std::istringstream iss(line);
//     if (i == 0){
//     	do
//     	{
//         	int sub;
//         	iss >> sub;
//        	 	A10.push_back(sub);
//     	} while (iss);
// 	}
//     if (i == 1){
//     	do
//     	{
//         	int sub;
//         	iss >> sub;
//        	 	A20.push_back(sub);
//     	} while (iss);
// 	}
//     if (i == 2){
//     	do
//     	{
//         	int sub;
//         	iss >> sub;
//        	 	A10.push_back(sub);
//     	} while (iss);
// 	}
//     if (i == 3){
//     	do
//     	{
//         	int sub;
//         	iss >> sub;
//        	 	A10.push_back(sub);
//     	} while (iss);
// 	}
//     if (i == 4){
//     	do
//     	{
//         	int sub;
//         	iss >> sub;
//        	 	A10.push_back(sub);
//     	} while (iss);
// 	}
//     if (i == 5){
//     	do
//     	{
//         	int sub;
//         	iss >> sub;
//        	 	A10.push_back(sub);
//     	} while (iss);
// 	}
//     if (i == 6){
//     	do
//     	{
//         	int sub;
//         	iss >> sub;
//        	 	A10.push_back(sub);
//     	} while (iss);
// 	}
// 	i++;
// }


// }

//int main()
//{
//	vector<int> A1;
//	A1.assign(10,0);
//	vector<int> A2;
//	A2.assign(10,0);
//	vector<float> W;
//	W.push_back(1);
//	W.push_back(10);
//	W.push_back(3);
//	W.push_back(100000);
//	W.push_back(100);
//	W.push_back(1.1231);
//	W.push_back(10);
//	W.push_back(1);
//	W.push_back(100);
//	W.push_back(10);
//
//
//	int k1=4;
//	int k2=4;
//	int n_A1=1;
//	int n_A2=2;
//	sub_sup A(A1,A2,W,k1,k2,n_A1,n_A2,0.5);
//	cout<<"----\n";
//	A.GreedyMax();
//	A.print_current_matching();
//	A.sub_sup_process(1);
//	A.sub_sup_process(2);
//	A.sub_sup_process(3);
//
//	return 0;
//
//}


