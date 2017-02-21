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
#include <unistd.h>
#include <cstdlib>
//#include "../fast_code/max_weight_bipartite_match_fast.h"
//#include "../fast_code/max_weight_bipartite_match_complete.h"
//#include "../fast_code/essentials.h"
//#include "../fast_code/somefunction.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "sub_sup.h"
#include "input_para.h"
using namespace std;

// Define struct variables
struct three_vaule {
	int a1;
	int a2;
	float w;
	float modular_bound;
	float s;
	int ion;
	int type;//1:center 2-4 neutral 5:iosto;
};
struct float_int_s {
	float s;
	int a;
};
bool sortBys_float_int(const float_int_s &lhs, const float_int_s &rhs) {return lhs.s > rhs.s;}
bool sortBys(const three_vaule &lhs, const three_vaule &rhs) {return lhs.s > rhs.s;}
bool sortByw(const three_vaule &lhs, const three_vaule &rhs) {return lhs.w > rhs.w;}
struct output_scores {
	int sid;
	int num_peptides;
	vector<string> peptide_names;
	vector<float> scores;
	vector<float> parameter;
	vector<float> *edges;
};

struct spectrum_obs {
	int num_obs_spec;
	vector<float> mz_obs;
	vector<float> intensity_obs;
};

struct spectrum_theo {
	string name_peptide;
	int num_theo_spec;
	vector<float> mz_theo;
	vector<float> intensity_theo;
	vector<float> bin_left;
	vector<float> bin_right;
};

float find_max(float* a, int n, int k) {
	vector <float> b;
	for (int i = 0; i < n; i++) b.push_back(a[i]);
	std::sort(b.begin(), b.end(), std::greater<float>());
	int c = k;
	if (c > n) c = n;
	float r = 0;
	for (int i = 0; i < c; i++) r += b.at(i);
	return r;
}

void charged_ion_factor(float *factor, int charge, int length_peptide, int length_ion, float s) {
	float x = ( (float)length_ion) / ((float)length_peptide);
	if (charge == 3) {

//		factor[1]=s*(1-x);
//		factor[2]=s*x;
		if (length_ion < length_peptide / 2) {

			factor[1] = s;
			factor[2] = 1 / s;
		} else {
			factor[2] = s;
			factor[1] = 1 / s;
		}
		return;
	}
	return;
}
// ---------------- Main function -------------------- //

//float calculate_abs_ranking_ch2(vector<float> weight_b, vector<float> weight_y, vector<float>other_p, int k2_input, int kind, int index, int charge, vector<string>filenames, int n0, int nn0, vector<float> lambda_vector)
float calculate_abs_ranking_ch2(input_para paras){

	std::vector<float> lambda_vector = paras.lambda_vector;
	int charge = paras.charge;
	float l_s_factor=1.1;

	vector<std::string> high_confidence(15000, "AAA");

	ofstream outs;
	outs.open("/s1/wrbai/codes/output/outs.txt", ios::trunc);
	if (true) {
		ifstream high_confidence_file;
//	high_confidence_file.open("/n/trombone/s1/wrbai/codes/matlab_codes/2016_11_1/plasm-10-ch3-high.txt");
		high_confidence_file.open("/n/trombone/s1/wrbai/codes/matlab_codes/2016_11_1/plasm-10-ch3-high1.txt");
//	high_confidence_file.open("/n/trombone/s1/wrbai/codes/matlab_codes/2016_11_1/plasm-10-ch3-high.txt");
		while (true) {
			int x;
			string y;
			if (high_confidence_file.eof()) break;
			high_confidence_file >> x;
			high_confidence_file >> y;
			//		cout<<x<<"============"<<y<<endl;
			high_confidence[x] = y;
		}
	}
	

	bool TMT_label = 1;	
	
	int start = 1;
	int ii = start - 1;
	bool add_ch3 = 1;

	int crux = 0;

	//load weight for high_res
	vector<float> weight_high_b_ion;
	weight_high_b_ion.assign(3001 + 3000, 0); //wight[0]:mz=-30.00; wight[3000]:mz=0;wight[1000]:mz=10;
	vector<float> weight_high_y_ion;
	weight_high_y_ion.assign(3001 + 3000, 0);
	ifstream weight_b_file;
	weight_b_file.open("/n/trombone/s1/wrbai/codes/pipeline/bion_weight_malaria.txt");
	for (int kk = 0; kk < weight_high_b_ion.size(); kk++) {
		float temp = 0;

		weight_b_file >> temp;

		if (kk > 3000 && kk < 3120) temp = temp;
		weight_high_b_ion[kk] = (temp * 1);
	}

	weight_b_file.close();
	ifstream weight_y_file;
	weight_y_file.open("/n/trombone/s1/wrbai/codes/pipeline/yion_weight_malaria.txt");
	for (int kk = 0; kk < weight_high_y_ion.size(); kk++) {
		float temp = 0;

		weight_y_file >> temp;
		if (kk > 3000 && kk < 3120) temp = temp;
		weight_high_y_ion[kk] = (temp * 1);
	}
	weight_y_file.close();
	//

	int n = paras.n; // Corresponds to the Peptide with a differing score from xcorr

	int nn = paras.nn;

	
//if (n>1000) n=1000;

	cout << n << " " << nn << endl;


	


	string line;
	ifstream myfile (paras.input_ms2_file.c_str());
	vector<int>sid_bin;
	getline (myfile, line);
	int Sid_bin;
	stringstream iss(line);
	string temp_string_bin;
	iss >> temp_string_bin;
	iss >> Sid_bin;
	sid_bin.push_back(Sid_bin);


	//cout<<temp_string_bin<<"--"<<Sid_bin<<endl;
	vector<float>obs_mz_bin[nn];
	vector<float>obs_int_bin[nn];
	float current_bin_mz = 0;
	float current_bin_int = 0;
	int i_bin = 0;
	float max_mz_global = 0;
	while (getline (myfile, line) )
	{
		//cout<<i_bin<<endl;
		if (line[0] == 'Z') continue;
		if (line[0] == 'S') {

			current_bin_mz = 0;
			current_bin_int = 0;
			int Sid_bin;
			stringstream iss(line);
			string temp_string_bin;
			iss >> temp_string_bin;
			iss >> Sid_bin;
			sid_bin.push_back(Sid_bin);
			//cout<<obs_mz_bin[i_bin].size()<<"\t";
			//	cout<<"i_bin: "<<i_bin<<" size: "<<obs_mz_bin[i_bin].size()<<" Sid "<<Sid_bin<<endl;
			i_bin = i_bin + 1;
//				cout<<endl;
			continue;
		}

		stringstream iss(line);
		float temp_1 = 0;
		float temp_2 = 0;
		iss >> temp_1;
		iss >> temp_2;
		//	cout<<i_bin<<"---\n";
		if (temp_1 > max_mz_global) {
			//	cout<<max_mz_global<<endl;
			max_mz_global = temp_1;
		}
		float bin = floor(temp_1 / 0.001);
		if (abs(bin - current_bin_mz) > 0.01) {
//			cout<<" a "<<current_bin_mz<<" "<<bin;
			current_bin_mz = bin;
			current_bin_int = temp_2;
			obs_mz_bin[i_bin].push_back(temp_1);

			obs_int_bin[i_bin].push_back(log(1 + temp_2));
		} else if (temp_2 >= current_bin_int) {
//			cout<<" 1 ";
			obs_mz_bin[i_bin][obs_mz_bin[i_bin].size() - 1] = (temp_1);

			obs_int_bin[i_bin][obs_int_bin[i_bin].size() - 1] = (log(1 + temp_2));
			current_bin_int = temp_2;
		} else {
//	cout<<" 2 ";
		}
		//obs_int_bin[i_bin].push_back(pow(temp_2,0.5));
		//		obs_int_bin[i_bin].push_back(log(1+temp_2));
	}
	//cout<<max_mz_global<<endl;
	for (int kk = 0; kk < sid_bin.size(); kk++) {
		for (int jj = 0; jj < obs_mz_bin[kk].size(); jj++) {
			//		cout<<obs_mz_bin[kk][jj]<<"\t";
		}
		//	cout<<"\n====================\n";
		for (int jj = 0; jj < obs_mz_bin[kk].size(); jj++) {
			//	cout<<obs_int_bin[kk][jj]<<"\t";
		}
		//	cout<<"\n";
//cout<<"---------------------------\n";
	}
	float max_mz = 0.0001;

	if (1 == 1) {
		for (int kk = 0; kk < sid_bin.size(); kk++) {
			float m1 = obs_mz_bin[kk][0];
			float m2 = obs_mz_bin[kk][obs_mz_bin[kk].size() - 1];
			bool nonzero1 = 1;
			for (int cc = 0; cc < obs_mz_bin[kk].size(); cc++) {
				if (nonzero1 && obs_int_bin[kk][cc] > 0.0001) {
					m1 = obs_mz_bin[kk][cc];
					nonzero1 = 0;
				}
				if (obs_int_bin[kk][cc] > 0.00001) {
					m2 = obs_mz_bin[kk][cc];
				}

			}
			float max[10];
			int n_bin[10];
			for (int jj = 0; jj < 10; jj++) max[jj] = 0;
			for (int jj = 0; jj < 10; jj++) n_bin[jj] = 0;
			for (int jj = 0; jj < obs_mz_bin[kk].size(); jj++) {
				int m = ((obs_mz_bin[kk][jj] - m1) / (m2 - m1) * 10.0);
				if (m <= 0) m = 0;
				if (m >= 10) m = 9;
				n_bin[m]++;
				if (obs_int_bin[kk][jj] > max[m]) max[m] = obs_int_bin[kk][jj];

			}


			for (int jj = 0; jj < obs_mz_bin[kk].size(); jj++) {
//				float m1=obs_mz_bin[kk][0];
//float m2=obs_mz_bin[kk][sid_bin.size()-1];

				int m = ((obs_mz_bin[kk][jj] - m1) / (m2 - m1) * 10.0);
				if (m <= 0) m = 0;
				if (m >= 10) m = 9;
				if (max[m] > 0.01) obs_int_bin[kk][jj] = 50.0 * obs_int_bin[kk][jj] / max[m];
				if (obs_int_bin[kk][jj] > 50.0) obs_int_bin[kk][jj] = 50.0;
				//if(obs_int_bin[kk][jj]>20.0) obs_int_bin[kk][jj]=20.0;


			}
		}

	}

	for (int kk = 0; kk < sid_bin.size(); kk++) {
		float max = 0;
		for (int jj = 0; jj < obs_mz_bin[kk].size(); jj++) {
			if (obs_int_bin[kk][jj] > max) max = obs_int_bin[kk][jj];
		}
		//	cout<<"\n====================\n";
		for (int jj = 0; jj < obs_mz_bin[kk].size(); jj++) {
			if (obs_int_bin[kk][jj] < 0.1 * max) obs_int_bin[kk][jj] = 0;
		}
		//	cout<<"\n";
//cout<<"---------------------------\n";
	}



	float bin_width = 1.0005079;
	float bin_offset = 0.68;

// Output scores

	output_scores output[n - start + 1]; // Create an array of output scores struct

// Create input/output streams

	ofstream scores_file2;
	ofstream all_scores_file;//contains all scores of target and decoy PSMs
	ofstream all_scores_file_other_decoy;//contains all scores of larger decoy set PSMs for calibration propose.
	ifstream infile;	
	
	int currentfile = 0;
	infile.open(paras.input_file.c_str());
	scores_file2.open(paras.output_file.c_str());
	int temp_int; float temp_float;
	while (ii < n)
	{

		if (ii % 100 == 0) {
			char hostname[1024];
			hostname[1023] = '\n';
			gethostname(hostname, 1023);
			//cout << hostname << "***" << filenames[1].c_str() << " " << filenames[4].c_str() << " " << " ii=" << ii << " ";
			for (int kkk = 0; kkk < lambda_vector.size(); kkk++) {
				cout << lambda_vector[kkk] << " ";
			}
			cout << endl;
// ------------------ READ OBSERVED SPECTRUM DATA ---------------------- //
		}

		spectrum_obs first;
		infile >> output[ii].sid;

// Input spectrum id number

// Input number of observed spectrum
		infile >> first.num_obs_spec;
// Reserve space for vector variables
		(first.mz_obs).reserve(first.num_obs_spec);
		(first.intensity_obs).reserve(first.num_obs_spec);

		/*
		float mz_obs_spec[first.num_obs_spec];
		float intensity_obs_spec[first.num_obs_spec];

		// Assign pointers to address of arrays

		first.mz_obs = mz_obs_spec;
		first.intensity_obs = intensity_obs_spec;
		*/

		int i = 0;

// Get observed m/z ratio for spectrum 1
		while (i < first.num_obs_spec)
		{
//infile>>temp_float;
			(first.mz_obs).push_back(temp_float);
			i++;
		}

		i = 0;

// Get observed intensity for spectrum 1
		while (i < first.num_obs_spec)
		{
//infile>>temp_float;
			(first.intensity_obs).push_back(temp_float);
			i++;
		}
// -------------------------------- READ THEORETICAL SPECTRUM DATA -------------------------------- //

// Get peptide names, m/z and intensity for the spectra of peptides


		int num_peptides;
		infile >> num_peptides; //Get number of peptides for spectrum 1
		output[ii].num_peptides = num_peptides;
		i = 0;
// print
//cout<<"# peptides = "<<num_peptides<<endl;

// Define array of structs

		spectrum_theo peptides[num_peptides]; // Define an array of theoretical spectrum structs
		bipartite_graph bipartite[num_peptides]; // Define an array of bipartite graphs
//output[ii].peptide_name_ptr = new string[num_peptides]; // Dynamic allocation of memory for storing strings
//output[ii].scores = new float[num_peptides];

// Reserve space for vector variables
		(output[ii].peptide_names).reserve(num_peptides * 2);
		(output[ii].scores).reserve(num_peptides * 2);
		(output[ii].scores).assign(num_peptides, 0);
		(output[ii].parameter).reserve(num_peptides * 2);
		(output[ii].parameter).assign(num_peptides, 0);
		(output[ii].edges) = new vector<float>[num_peptides];



//while (i < num_peptides)
		while (i < num_peptides) //For debugging purposes
		{


// Get name of peptide
			infile >> peptides[i].name_peptide;

			/*
			if (peptides[i].name_peptide == "LEGSYHWYMEK")
			{
			break;
			}
			*/

			(output[ii].peptide_names).push_back(peptides[i].name_peptide);

// Get number of theoretical spectrum
			infile >> peptides[i].num_theo_spec;

			/*
			// Dynamic memory allocation
			peptides[i].mz_theo = new float[peptides[i].num_theo_spec];
			peptides[i].intensity_theo = new float[peptides[i].num_theo_spec];
			*/


// Reserve memory for vector variables

			(peptides[i].mz_theo).reserve(peptides[i].num_theo_spec);
			(peptides[i].intensity_theo).reserve(peptides[i].num_theo_spec);
			(peptides[i].bin_left).reserve(peptides[i].num_theo_spec);
			(peptides[i].bin_right).reserve(peptides[i].num_theo_spec);



			int j = 0;
//cout<<"j = "<<j<<endl;

// Read off the theoretical m/z locations for current peptide
			while (j < peptides[i].num_theo_spec)
			{
				infile >> temp_float;
				(peptides[i].mz_theo).push_back(temp_float);
				j++;
			}

			j = 0;

// Read off the theoretical intensity locations for current peptide
			while (j < peptides[i].num_theo_spec)
			{
				infile >> temp_float;
				(peptides[i].intensity_theo).push_back(temp_float);
//(output[ii].scores).at(i) = (output[ii].scores).at(i) + (peptides[i].intensity_theo).at(j);
				j++;
			}

			i++;
		}



		i = 0;
		float local_max = 0;

		int ii_sid = 0;
		for (int kk = 0; kk < sid_bin.size(); kk++) {
			if (sid_bin[kk] == output[ii].sid) { ii_sid = kk; break;}

		}
// cout<<"num_peptides: "<< num_peptides <<endl;
		while (i < num_peptides)
		{
			int is_top = 1;
//for(int kk=0;kk<60;kk++){if(top_index_i[ii][kk]==i) is_top=1;}
			if (is_top == 0) {
				output[ii].scores.at(i) == -100000;
				i++;
				continue;
			}

			string temp_name = peptides[i].name_peptide;

			vector<float> b_ion;
			vector<float> y_ion;

			int size_theo = peptides[i].mz_theo.size();
			for (int kk = 0; kk < size_theo; kk++) {
				float mze = (peptides[i].mz_theo.at(kk));
				if (TMT_label == 1) {
					if (mze < 10000) {b_ion.push_back(mze);}
					else {
						if (temp_name.at(temp_name.length() - 1) == 'K') {

							y_ion.push_back(mze - 10000.0);
						} else {

							y_ion.push_back(mze - 10000.0);
						}
					}
				}
				else {
					if (mze < 10000) {b_ion.push_back(mze);}
					else y_ion.push_back(mze - 10000);
				}
			}

			std::sort (b_ion.begin(), b_ion.end());
			std::sort(y_ion.begin(), y_ion.end(), std::greater<int>());
		

			vector<int>A1;
			vector<int>A2;
			vector<float>obs_mz_A1;
			vector<float>obs_intensity_A1;
			vector<float>obs_mz_A2;
			vector<float>W2;
			vector<float>We;
			vector<float>We2;
			vector<int>B1;
			vector<int>B2;
			vector<int>ion;
			vector<int>type;
			B1.assign(2100, 0);
			B2.assign((int)(2 * (b_ion.size() + y_ion.size()) + 1000), 0);
			int nb=b_ion.size();
	
			float flanking_factor1=1.0;

			for (int kk = 0; kk < b_ion.size(); kk++) {

				float temp1 = kk;
				float temp2 = b_ion.size();
				//float ion_length_factor=4.5-4.0*pow(2.0*((temp1+1-temp2/2.0)/temp2),2);
				float ion_length_factor = 0.9 + 0.2 * temp1 / temp2;
				ion_length_factor = 1.0;
				float mze=b_ion.at(kk);
				if(mze>=first.intensity_obs.size()) continue;


				int length_ion = kk + 1;
				int length_peptide = b_ion.size() + 1;
				float charge_factor[10];
				for (int jj = 0; jj < 10; jj++) charge_factor[10] = 1.0;
				charged_ion_factor(charge_factor, charge, length_peptide, length_ion, l_s_factor);


				for (int jj = 0; jj < obs_mz_bin[ii_sid].size(); jj++) {
					if (abs(obs_mz_bin[ii_sid][jj] - mze) < 30.0) {
						int diff = (int)((obs_mz_bin[ii_sid][jj] - mze) * 100) + 3000;
						if (diff < 0) continue;

						We2.push_back(ion_length_factor * charge_factor[1]*obs_int_bin[ii_sid][jj]*flanking_factor1 * weight_high_b_ion[diff]);
						if (weight_high_b_ion[diff] < 0.01)continue;
						We.push_back(ion_length_factor * charge_factor[1]*obs_int_bin[ii_sid][jj]*flanking_factor1 * weight_high_b_ion[diff]);
						int n_or_c = 0;
						if (diff < 2800) n_or_c = 1;
//A1.push_back(jj);

						A1.push_back((int)(obs_mz_bin[ii_sid][jj] - bin_offset));
						A2.push_back(2 * kk);
						obs_mz_A1.push_back(obs_mz_bin[ii_sid][jj]);
						obs_mz_A2.push_back(mze);
						obs_intensity_A1.push_back(obs_int_bin[ii_sid][jj]);


						ion.push_back(0);
						type.push_back(diff - 3000);

					}


				}

//add charge3
				if (add_ch3 == 1 && charge >= 3) {
					for (int charge_ion = 2; charge_ion < charge; charge_ion++) {
						float mz_charge_ion = mze;
						mz_charge_ion = (mz_charge_ion + (float)charge_ion - 1) / ((float)charge_ion);
						for (int jj = 0; jj < obs_mz_bin[ii_sid].size(); jj++) {
							if (abs(obs_mz_bin[ii_sid][jj] - mz_charge_ion) < 30.0) {
								int diff = (int)((obs_mz_bin[ii_sid][jj] - mz_charge_ion) * 100) + 3000;

								diff = (int)((obs_mz_bin[ii_sid][jj] * (float)(charge_ion) + 1.0 - (float)(charge_ion) - mze) * 100.0) + 3000;
								if (diff < 0) continue;

								if (diff >= weight_high_y_ion.size()) continue;

								We2.push_back(ion_length_factor * charge_factor[charge_ion]*obs_int_bin[ii_sid][jj]*flanking_factor1 * weight_high_b_ion[diff]);
								if (weight_high_b_ion[diff] < 0.01)continue;
								We.push_back(ion_length_factor * charge_factor[charge_ion]*obs_int_bin[ii_sid][jj]*flanking_factor1 * weight_high_b_ion[diff]);
								int n_or_c = 0;
								if (diff < 2800) n_or_c = 1;
								//A1.push_back(jj);

								A1.push_back((int)(obs_mz_bin[ii_sid][jj] - bin_offset));
								A2.push_back(2 * kk + 1);
								obs_mz_A1.push_back(obs_mz_bin[ii_sid][jj]);
								obs_mz_A2.push_back(mze);
								obs_intensity_A1.push_back(obs_int_bin[ii_sid][jj]);




								ion.push_back(0);
								type.push_back(diff - 3000);

							}

						}
					}
				}
//add charge3 end
			}
			for (int kk = 0; kk < y_ion.size(); kk++) {
				float temp1 = y_ion.size() - kk;
				float temp2 = y_ion.size();
				//float ion_length_factor=1.5-2.0*abs(temp1+1-temp2/2.0)/temp2;

//	float ion_length_factor=4.5-4.0*abs(2.0*((temp1+1-temp2/2.0)/temp2));
//
				float ion_length_factor = 0.9 + 0.2 * temp1 / temp2;

				ion_length_factor = 1.0;
				//if (kk==y_ion.size()-1) ion_length_factor=0.1;


				//decide if y has some loss
//	int n_NH3_loss=0;
//	int n_H2O_loss=1;
//	for (int jj=kk+1;jj<y_ion.size()+1;jj++){
//		for (int j2=0;j2<4;j2++){
//			if (NH3_loss[j2]==amino_acid[jj]) n_NH3_loss++;
//			if (H2O_loss[j2]==amino_acid[jj]) n_H2O_loss++;
//		}
//
//	}
//	cout<<"y"<<n_NH3_loss<<" "<<n_H2O_loss<<endl;
				//end

				float mze = y_ion.at(kk);
		
//iso[1]=n_NH3_loss*factor_Nh3*iso[1];
//iso[0]=n_H2O_loss*factor_Nh3*iso[0];
//if(n_NH3_loss==0) iso[1]=0;
//if(n_H2O_loss==0) iso[0]=0;

//float iso_sum=0;
//for(int jj=0;jj<6;jj++) iso_sum+=iso[jj];
//cout<<iso_sum<<endl;
//for low res
//for(int jj=0;jj<n_iso;jj++){
//if(iso[jj]<0.001) continue;
//if(mze+iso_mz[jj]<first.intensity_obs.size()-1&&mze+iso_mz[jj]>=0&&first.intensity_obs.at(mze+iso_mz[jj])>0){
//We.push_back(iso[jj]*first.intensity_obs.at(mze+iso_mz[jj]));
////daishan
////float temp_50=iso[jj]*first.intensity_obs.at(mze+iso_mz[jj]);
////if(temp_50>49){
////cout<<"ii: "<<ii<<" i: "<<i<<" mze:"<<mze<<" intensity: "<<temp_50<<" "<<y_ion.size()-kk<<endl;
////}
////daishan end
//A1.push_back(mze+iso_mz[jj]);
//A2.push_back(kk+nb);
//ion.push_back(1);
//type.push_back(jj-3);
//}
//}
//end for low res
				int length_ion = y_ion.size() - kk;
				int length_peptide = y_ion.size() + 1;
				float charge_factor[10];
				for (int jj = 0; jj < 10; jj++) charge_factor[10] = 1.0;
				charged_ion_factor(charge_factor, charge, length_peptide, length_ion, l_s_factor);


				for (int jj = 0; jj < obs_mz_bin[ii_sid].size(); jj++) {
					if (abs(obs_mz_bin[ii_sid][jj] - mze) < 30.0) {
//	int diff=(int)((obs_mz_bin[ii_sid][jj]-mze)*100)+3000;
						int diff = (int)((obs_mz_bin[ii_sid][jj] - mze) * 100) + 3000;
						if (diff < 0) continue;

						int n_or_c = 0;
						if (diff < 2800) n_or_c = 1;

						We2.push_back(ion_length_factor * charge_factor[1]*obs_int_bin[ii_sid][jj]*flanking_factor1 * weight_high_y_ion[diff]);
						if (weight_high_y_ion[diff] < 0.01)continue;
						We.push_back(ion_length_factor * charge_factor[1]*obs_int_bin[ii_sid][jj]*flanking_factor1 * weight_high_y_ion[diff]);
						A1.push_back((int)(obs_mz_bin[ii_sid][jj] - bin_offset));
						A2.push_back(2 * (kk + nb) + 1);

						obs_mz_A1.push_back(obs_mz_bin[ii_sid][jj]);
						obs_mz_A2.push_back(mze);
						obs_intensity_A1.push_back(obs_int_bin[ii_sid][jj]);




						ion.push_back(1);
						type.push_back(diff - 3000);

					}


				}

//add charge3
				if (add_ch3 == 1 && charge >= 3) {
					for (int charge_ion = 2; charge_ion < charge; charge_ion++) {
						float mz_charge_ion = mze;
						mz_charge_ion = (mz_charge_ion + (float)charge_ion - 1) / ((float)charge_ion);
						for (int jj = 0; jj < obs_mz_bin[ii_sid].size(); jj++) {
							if (abs(obs_mz_bin[ii_sid][jj] - mz_charge_ion) < 30.0) {
								int diff = (int)((obs_mz_bin[ii_sid][jj] - mz_charge_ion) * 100) + 3000;
								diff = (int)((obs_mz_bin[ii_sid][jj] * (float)(charge_ion) + 1.0 - (float)(charge_ion) - (mze)) * 100.0) + 3000;
								if (diff < 0) continue;
								if (diff >= weight_high_y_ion.size()) continue;

								We2.push_back(ion_length_factor * charge_factor[charge_ion]*obs_int_bin[ii_sid][jj]*flanking_factor1 * weight_high_y_ion[diff]);
								if (weight_high_y_ion[diff] < 0.01)continue;
								We.push_back(ion_length_factor * charge_factor[charge_ion]*obs_int_bin[ii_sid][jj]*flanking_factor1 * weight_high_y_ion[diff]);
								int n_or_c = 0;
								if (diff < 2800) n_or_c = 1;
								//A1.push_back(jj);

								A1.push_back((int)(obs_mz_bin[ii_sid][jj] - bin_offset));
								A2.push_back(2 * (kk + nb));
								obs_mz_A1.push_back(obs_mz_bin[ii_sid][jj]);
								obs_mz_A2.push_back(mze);
								obs_intensity_A1.push_back(obs_int_bin[ii_sid][jj]);



								ion.push_back(0);
								type.push_back(diff - 3000);

							}

						}
					}
				}
//add charge3 end
			}


			float lambda3 = 0.9;
			lambda3 = 0.9;
//lambda3=1.0;


//sub_sup all_graph(A1,A2,We,10,2,2000,2*b_ion.size()+2*y_ion.size(),other_p[3],other_p[6],other_p[7]);
			sub_sup all_graph(A1, A2, We, obs_mz_A1, obs_mz_A2, obs_intensity_A1, 0, 10, 2, 2000, 2 * b_ion.size() + 2 * y_ion.size(), lambda_vector);

			if (0) {

				cout << ii << " " << i << endl;

				outs << all_graph.GreedyMax() << "\t" << all_graph.submodualr_projection(2.0) << "\t" << all_graph.sub_sup_process(1) << "\t" << all_graph.sub_sup_process(2) << "\t" << all_graph.sub_sup_process(3) << "\t" << all_graph.random_scores() << endl;
			}




			output[ii].edges[i] = We;

			(output[ii].scores).at(i) = all_graph.GreedyMax();


			//if (high_confidence[output[ii].sid] == peptides[i].name_peptide && 1 == 1) {
			//	cout << "find a high PSMs " << output[ii].sid << " " << peptides[i].name_peptide << endl;
			//	all_graph.write_all_bipartite(output[ii].sid);
			//	all_graph.write_all_parameters(output[ii].sid);
//all_graph.GreedyMax_initialByOne(output[ii].sid);
			//}
			


			i++;
		}
//dealing with calbration
//should be uncommented

		vector<float>scoreset2 = output[ii].scores;
		vector<float>scoreset;
		for (i = 0; i < scoreset2.size(); i++) scoreset.push_back(scoreset2.at(i));
		float sum = std::accumulate(scoreset.begin(), scoreset.end(), 0.0);
		float mean = sum / scoreset.size();

		float sq_sum = std::inner_product(scoreset.begin(), scoreset.end(), scoreset.begin(), 0.0);
		float stdev = std::sqrt((sq_sum / scoreset.size() - mean * mean));
		for (int i = 0; i < scoreset2.size(); i++) {
			float t = (scoreset2.at(i));

			(output[ii].scores).at(i) = (t - mean);

		}
		i = 0;
		ii = ii + 1;
	} // End spectra while loop

	double number_cand = 0;

	for (int kk = 0; kk < n; kk++) {
		number_cand += (float)(output[kk].scores.size());

	}
	cout << "average:" << number_cand / ((float)(n));



	vector<float_int_s>all_scores;


	ofstream para_file2;
	ofstream para_file3;
	para_file2.open("output/para_temp.txt");
	para_file3.open("output/para2.txt");
	ii = 0;
//cout<<"ii = "<<ii<<endl;
	while (ii < n)
	{

		
		int j = 0;
//****************2014-10-26*****************************//

		int lens = output[ii].scores.size();
		vector<float> target(output[ii].scores.begin(), output[ii].scores.begin() + lens / 2);
		vector<float> decoy(output[ii].scores.begin() + lens / 2, output[ii].scores.end());
		float maxt = -10000;


//if(target.size()<10) cout<<lens<<" "<<target.size()<<" " <<decoy.size()<<endl;


		int maxt_i = 0;
		for (int kk = 0; kk < target.size(); kk++) {
			if (target.at(kk) > maxt) {
				maxt = target.at(kk);
				maxt_i = kk;
			}
		}
		float maxd = -10000;
		int maxd_i = 0;
		for (int kk = 0; kk < decoy.size(); kk++) {
			if (decoy.at(kk) > maxd) {
				maxd = decoy.at(kk);
				maxd_i = kk;
			}
		}
		if (maxt >= maxd) maxd = -1000000000;
		if (maxt < maxd) maxt = -1000000000;
		scores_file2 << "t\t" << output[ii].sid << "\t" << (output[ii].peptide_names).at(maxt_i) << "\t" << maxt << endl;
		scores_file2 << "d\t" << output[ii].sid << "\t" << (output[ii].peptide_names).at(maxd_i + lens / 2) << "\t" << maxd << endl;

		if (maxt >= maxd) {

//para_file2<<output[ii].parameter[maxt_i]<<"\t"<<-10<<"\t"<<maxt<<endl;
			para_file2 << ii << "\t" << maxt_i << "\t" << maxt << endl;
			if (maxt > 162 && maxt < 232) {

				para_file3 << "t\t" << output[ii].edges[maxt_i].size() << "\t" << maxt << endl;
				for (int kk = 0; kk < output[ii].edges[maxt_i].size(); kk++) {
					para_file3 << output[ii].edges[maxt_i][kk] << "\t";
				}
				para_file3 << endl;
			}
		}
		if (maxt < maxd) {

// para_file2<<-10<<"\t"<<output[ii].parameter[maxd_i+lens/2]<<"\t"<<maxd<<endl;
			if (maxd > 162) {
				para_file3 << "d\t" << output[ii].edges[maxd_i + lens / 2].size() << "\t" << maxd << endl;
				for (int kk = 0; kk < output[ii].edges[maxd_i + lens / 2].size(); kk++) {
					para_file3 << output[ii].edges[maxd_i + lens / 2][kk] << "\t";
				}
				para_file3 << endl;
			}

		}
		float_int_s temp_float_int1;
		temp_float_int1.s = maxt;
		temp_float_int1.a = 1;
		all_scores.push_back(temp_float_int1);
		float_int_s temp_float_int2;
		temp_float_int2.s = maxd;
		temp_float_int2.a = 2;
		all_scores.push_back(temp_float_int2);



		ii++;
	}


	sort(all_scores.begin(), all_scores.end(), sortBys_float_int);
//for(int kk=0;kk<100;kk++) cout<<all_scores.at(kk).s<<" "<<all_scores.at(kk).a<<endl;
	float q[11];
	for (int kk = 0; kk < 11; kk++) q[kk] = float(kk) / 100;
	float area = 0;
	float high_t = 0;
	float high_d = 0;
	int j = 0;
	for (int kk = 0; kk < all_scores.size(); kk++) {
		float s1 = all_scores.at(kk).s;
		float a1 = all_scores.at(kk).a;
		if (a1 == 1) high_t += 1;
		else high_d += 1;
		if (high_d / (float(kk)) > q[j]) {j = j + 1; area += high_t; if (j > 10) break;}



	}
// Close any remaining files
//scores_file.close();
	scores_file2.close();
	infile.close();
//edges_file.close();
	ofstream output_para;
	output_para.open("parameter.txt");
	for (int para_i = 0; para_i < lambda_vector.size(); para_i++) {
		output_para << lambda_vector[para_i] << endl;
	}

	output_para.close();
	return area;
}
