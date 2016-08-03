#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <numeric> 
#include <functional>
#include <algorithm>
#include <cmath>
#include <unistd.h>
#include <cstdlib>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "sub_sup.h"
using namespace std;
// Define struct variables
struct three_vaule{
int a1;
int a2;
float w;
float modular_bound;
float s;
int ion;
int type;//1:center 2-4 neutral 5:iosto;
};
struct float_int_s{
float s;
int a;
};
bool sortBys_float_int(const float_int_s &lhs,const float_int_s &rhs){return lhs.s>rhs.s;}



bool sortBys(const three_vaule &lhs,const three_vaule &rhs){return lhs.s>rhs.s;}
bool sortByw(const three_vaule &lhs,const three_vaule &rhs){return lhs.w>rhs.w;}
struct output_scores{
int sid;
int num_peptides;
vector<string> peptide_names;
vector<float> scores;
vector<float> parameter;
vector<float> *edges;
};

struct spectrum_obs{
int num_obs_spec;
vector<float> mz_obs;
vector<float> intensity_obs;
};

struct spectrum_theo{
string name_peptide;
int num_theo_spec;
vector<float> mz_theo;
vector<float> intensity_theo;
vector<float> bin_left;
vector<float> bin_right;
};
float find_max(float* a,int n, int k){
vector <float> b;
for(int i=0;i<n;i++) b.push_back(a[i]);
std::sort(b.begin(), b.end(), std::greater<float>());
int c=k;
if(c>n) c=n;
float r=0;
for(int i=0;i<c;i++) r+=b.at(i);
return r;
}

float max_2(float *a,int n,int *b){

int m1=0;
int m2=0;
int m3=0;
float s1=-1000;
float s2=-1000;
float s3=-1000;
for(int i=0;i<n;i++){
if (a[i]>s1) {
m1=i;
s1=a[i];
}
}
for (int i=0;i<n;i++){
if (a[i]>s2 && i!=m1){
m2=i;
s2=a[i];
}
}
for (int i=0;i<n;i++){
if(a[i]>s3 && i!=m1 && i!=m2){
m3=i;
s3=a[i];
}
}
b[0]=m1;
b[1]=m2;
b[2]=m3;
if(s3<-100) cout<<"error\n";
//if(s1<0) s1=0;
//if(s2<0) s2=0;
return s1+s2+s3;
}

void charged_ion_factor(float *factor,int charge,int length_peptide, int length_ion, float s){
	float x=( (float)length_ion)/((float)length_peptide);
	if (charge==3) {
	
//		factor[1]=s*(1-x);
//		factor[2]=s*x;
		if (length_ion<length_peptide/2) {

		factor[1]=s;
		factor[2]=1/s;
	}else{
		factor[2]=s;
		factor[1]=1/s;
	}
	return;
	}
	return;
}
// ---------------- Main function -------------------- //

float calculate_abs_ranking_ch2(vector<float> weight_b,vector<float> weight_y,vector<float>other_p,int k2_input,int kind,int index,int charge,vector<string>filenames,int n0, int nn0)
{
	kind=-1;
	float sum_all_scores=0;
	float number_all_scores=0;
	vector<int>high_confidence_i;
	vector<int>high_confidence_ii;
	vector<float>high_confidence_w;
	ofstream outs;
	outs.open("/s1/wrbai/codes/output/outs.txt",ios::trunc);
	if(true){
	ifstream high_confidence_file;
	high_confidence_file.open("output/para.txt");
	int counting_3=0;
	while(true){
		int x;
		float y;
		if(high_confidence_file.eof()) break;
		if(counting_3==0){
			high_confidence_file>>x;
			high_confidence_ii.push_back(x);
			counting_3=1;
		}else if(counting_3==1){
			high_confidence_file>>x;
			high_confidence_i.push_back(x);
			counting_3=2;
		}else{
			high_confidence_file>>y;
			high_confidence_w.push_back(y);
			counting_3=0;
		}
	}
	cout<<high_confidence_i.size()<<high_confidence_ii.size()<<high_confidence_w.size()<<endl;
	cout<<high_confidence_i[0]<<" "<<high_confidence_i[1]<<endl;
	}
	int num_high_confidence=high_confidence_ii.size();
	bool all_score_output_flag=0;
	bool all_score_larger_decoy_flag=0;
	if(all_score_output_flag) all_score_larger_decoy_flag=0;
	int index_larger_decoy_output=2;
	srand(index_larger_decoy_output);

	bool TMT_label=1;
	float l_s_factor=1.1;
int k2_matroid=k2_input;
int start = 1;
int ii = start - 1;
bool add_ch3=1;
//int charge=3;
int crux=0;

	//load weight for high_res
	vector<float> weight_high_b_ion;
	weight_high_b_ion.assign(3001+3000,0);//wight[0]:mz=-30.00; wight[3000]:mz=0;wight[1000]:mz=10;
    vector<float> weight_high_y_ion;
	weight_high_y_ion.assign(3001+3000,0);
	ifstream weight_b_file;
   weight_b_file.open("bion_weight_malaria.txt");//working fo malaria
	for(int kk=0;kk<weight_high_b_ion.size();kk++){
	float temp=0;
	
	 weight_b_file>>temp;
	 
	 weight_high_b_ion[kk]=(temp);
	}
	
	weight_b_file.close();
	ifstream weight_y_file;
	weight_y_file.open("yion_weight_malaria.txt");
	for(int kk=0;kk<weight_high_y_ion.size();kk++){
	float temp=0;
	
	 weight_y_file>>temp;
	 weight_high_y_ion[kk]=(temp);
	}
	weight_y_file.close();
int n =34499; // Corresponds to the Peptide with a differing score from xcorr

int nn = 0;





n= n0;



nn=nn0;



//ifstream myfile (stream.str().c_str());
cout<<filenames[1].c_str()<<endl;
ifstream myfile (filenames[1].c_str());
	vector<int>sid_bin;
	getline (myfile,line);
	int Sid_bin;
	stringstream iss(line);	
	string temp_string_bin;
	iss>>temp_string_bin;
	iss>>Sid_bin;
	sid_bin.push_back(Sid_bin);
	

	//cout<<temp_string_bin<<"--"<<Sid_bin<<endl;
		vector<float>obs_mz_bin[nn];
		vector<float>obs_int_bin[nn];
	float current_bin_mz=0;
	float current_bin_int=0;	
		int i_bin=0;
		float max_mz_global=0;
	    while (getline (myfile,line) )
		{
			//cout<<i_bin<<endl;
			if(line[0]=='Z') continue;
			if(line[0]=='S'){
				
		 current_bin_mz=0;
		 current_bin_int=0;	
				int Sid_bin;
				stringstream iss(line);	
				string temp_string_bin;
				iss>>temp_string_bin;
				iss>>Sid_bin;
				sid_bin.push_back(Sid_bin);
				//cout<<obs_mz_bin[i_bin].size()<<"\t";
			//	cout<<"i_bin: "<<i_bin<<" size: "<<obs_mz_bin[i_bin].size()<<" Sid "<<Sid_bin<<endl;
				i_bin=i_bin+1;
//				cout<<endl;
				continue;	
			}	
			
			stringstream iss(line);	
			float temp_1=0;
			float temp_2=0;
			iss>>temp_1;
			iss>>temp_2;
		//	cout<<i_bin<<"---\n";
			if (temp_1>max_mz_global){
		//	cout<<max_mz_global<<endl;
			 max_mz_global=temp_1;
			}
		float bin=floor(temp_1/0.001);
		if (abs(bin-current_bin_mz)>0.01){
//			cout<<" a "<<current_bin_mz<<" "<<bin;
			current_bin_mz=bin;
			current_bin_int=temp_2;
			obs_mz_bin[i_bin].push_back(temp_1);
		
			obs_int_bin[i_bin].push_back(log(1+temp_2));
		}else if(temp_2>=current_bin_int){
//			cout<<" 1 ";
			obs_mz_bin[i_bin][obs_mz_bin[i_bin].size()-1]=(temp_1);
		
			obs_int_bin[i_bin][obs_int_bin[i_bin].size()-1]=(log(1+temp_2));
			current_bin_int=temp_2;
		}else{
//	cout<<" 2 ";
		}
			//obs_int_bin[i_bin].push_back(pow(temp_2,0.5));
	//		obs_int_bin[i_bin].push_back(log(1+temp_2));
		}
		//cout<<max_mz_global<<endl;
		for(int kk=0;kk<sid_bin.size();kk++){
			for (int jj=0;jj<obs_mz_bin[kk].size();jj++){
		//		cout<<obs_mz_bin[kk][jj]<<"\t";
			}
		//	cout<<"\n====================\n";
			for (int jj=0;jj<obs_mz_bin[kk].size();jj++){
			//	cout<<obs_int_bin[kk][jj]<<"\t";
			}
		//	cout<<"\n";
//cout<<"---------------------------\n";
		}
		float max_mz=0.0001;
		
		if(1==1){
	for(int kk=0;kk<sid_bin.size();kk++){
	float m1=obs_mz_bin[kk][0];
	float m2=obs_mz_bin[kk][obs_mz_bin[kk].size()-1];	
		bool nonzero1=1;
		for(int cc=0; cc<obs_mz_bin[kk].size();cc++){
			if(nonzero1 && obs_int_bin[kk][cc]>0.0001){
				m1=obs_mz_bin[kk][cc];
				nonzero1=0;
			}		
			if(obs_int_bin[kk][cc]>0.00001){
				m2=obs_mz_bin[kk][cc];
			}
				
		}
		float max[10];
		int n_bin[10];
		for(int jj=0;jj<10;jj++) max[jj]=0;
		for(int jj=0;jj<10;jj++) n_bin[jj]=0;
		for(int jj=0;jj<obs_mz_bin[kk].size();jj++){
			int m=((obs_mz_bin[kk][jj]-m1)/(m2-m1)*10.0);
			if(m<=0) m=0;
			if(m>=10) m=9;
			n_bin[m]++;
			if(obs_int_bin[kk][jj]>max[m]) max[m]=obs_int_bin[kk][jj];

		}
		
		
		for(int jj=0;jj<obs_mz_bin[kk].size();jj++){
//				float m1=obs_mz_bin[kk][0];
//float m2=obs_mz_bin[kk][sid_bin.size()-1];	

			int m=((obs_mz_bin[kk][jj]-m1)/(m2-m1)*10.0);
			if(m<=0) m=0;
			if(m>=10) m=9;
			if(max[m]>0.01) obs_int_bin[kk][jj]=50.0*obs_int_bin[kk][jj]/max[m];
			if(obs_int_bin[kk][jj]>50.0) obs_int_bin[kk][jj]=50.0;
			//if(obs_int_bin[kk][jj]>20.0) obs_int_bin[kk][jj]=20.0;
		

		}
	}

		}

		for(int kk=0;kk<sid_bin.size();kk++){
			float max=0;
			for (int jj=0;jj<obs_mz_bin[kk].size();jj++){
				if(obs_int_bin[kk][jj]>max) max=obs_int_bin[kk][jj];
			}
		//	cout<<"\n====================\n";
			for (int jj=0;jj<obs_mz_bin[kk].size();jj++){
				if(obs_int_bin[kk][jj]<0.1*max) obs_int_bin[kk][jj]=0;
			}
		//	cout<<"\n";
//cout<<"---------------------------\n";
		}

		if(0==1){
	for(int kk=0;kk<sid_bin.size();kk++){
			float sum=0;
			float max_int=0;
			for (int jj=0;jj<obs_mz_bin[kk].size();jj++){
				sum+=obs_int_bin[kk][jj];
				if(max_int<obs_int_bin[kk][jj]) max_int=obs_int_bin[kk][jj];
			}
			for (int jj=0;jj<obs_mz_bin[kk].size();jj++){
			//	obs_int_bin[kk][jj]=20.0*obs_int_bin[kk][jj]/sum;
		//		obs_int_bin[kk][jj]=50.0*obs_int_bin[kk][jj]/max_int;
				
			}

		
		}
		}
		if(0==1){

			vector<float>temp_int_v;
			temp_int_v.assign(20000,0);
		for(int kk=0;kk<sid_bin.size();kk++){


			for (int jj=0;jj<obs_mz_bin[kk].size();jj++){
				
				float back=0;
				for (int ii=0;ii<obs_mz_bin[kk].size();ii++){
				if(abs(obs_mz_bin[kk][jj]-obs_mz_bin[kk][ii])<75.0) back+=obs_int_bin[kk][ii]/151.0;
				}
				
			temp_int_v[jj]=((obs_int_bin[kk][jj]-back));
			if (temp_int_v[jj]<=0) temp_int_v[jj]=0;
			}
			
			for (int jj=0;jj<obs_mz_bin[kk].size();jj++){obs_int_bin[kk][jj]=temp_int_v[jj];}
		//	cout<<"\n====================\n";
		//	cout<<"\n";
//cout<<"---------------------------\n";
		}

		}

//end read 
//


//


int calculate_by_part=(int)other_p.at(5);

int number_egde_print=0;
float daishan=0;
float n_daishan=0;
//srand(time(NULL));
// Number of spectra

int part_start=0;
int part_end=n;
if(calculate_by_part>0){
	part_start=0+(calculate_by_part-1)*(n/3);
	part_end=(calculate_by_part)*(n/3);
}
if(calculate_by_part==3){
	part_end=n;
}
	//n=20000;
//n=13390;
//n=10000;
//cout<<"n: "<<n<<endl;
//n=3000;
// Bin parameters for bipartite graph construction
float bin_width = 1.0005079;
float bin_offset = 0.68;

// Output scores

output_scores output[n - start + 1]; // Create an array of output scores struct
//daishan
/*
ofstream out;
ifstream in;
//in.open()
in.open("./input/yeast-01-ch2-precursorFilt-20130321/msms-fastSequestTransform-bw1.0005-bo0.68-tol3Da-yeast-01-2-ch2.txt");
out.open("./Wenruo/spectrum.txt",ios::trunc);
for(int h=0;h<100;h++){
string abc;
getline(in,abc);
out<<abc<<endl;
}
out.close();
in.close();
*/
// Create input/output streams
ofstream scores_file2;
ofstream all_scores_file;//contains all scores of target and decoy PSMs
ofstream all_scores_file_other_decoy;//contains all scores of larger decoy set PSMs for calibration propose.
ifstream infile;
string filename1[10];
int number_NH3[10];
for(int kk=0;kk<10;kk++) number_NH3[kk]=0;
int currentfile=0;

filename1[0]=filenames[0];
cout<<filename1[0].c_str()<<endl;
infile.open(filename1[0].c_str());
scores_file2.open(filenames[3].c_str());

ifstream possible_i_file;
possible_i_file.open("/s1/wrbai/codes/output/top_index.txt");
int **top_index_i;
top_index_i=new int*[n];
for(int kk=0;kk<n;kk++) top_index_i[kk]=new int[60];
for(int kk=0;kk<n;kk++){
for(int jj=0;jj<60;jj++)
	possible_i_file>>top_index_i[kk][jj];
}
int temp_int; float temp_float;




while (ii < n)
{
//coutall
if(ii%100==0){
char hostname[1024];
hostname[1023]='\n';
gethostname(hostname,1023);
cout<<hostname<<"***"<<filenames[1].c_str()<<" ii="<<ii<<" lambda1="<<other_p[3]<<" lambda2="<<other_p[6]<<" power="<<other_p[7]<<endl;
// ------------------ READ OBSERVED SPECTRUM DATA ---------------------- //
}

spectrum_obs first;
infile>>output[ii].sid;
//cout<<output[ii].sid<<endl;
// Input spectrum id number
//infile>>output[ii].sid; linshishandiao
// Input number of observed spectrum
infile>>first.num_obs_spec;
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
infile>>num_peptides; //Get number of peptides for spectrum 1
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
(output[ii].peptide_names).reserve(num_peptides*2);
(output[ii].scores).reserve(num_peptides*2);
(output[ii].scores).assign(num_peptides,0);
(output[ii].parameter).reserve(num_peptides*2);
(output[ii].parameter).assign(num_peptides,0);
(output[ii].edges)=new vector<float>[num_peptides];



//while (i < num_peptides)
while (i < num_peptides) //For debugging purposes
{


// Get name of peptide
infile>>peptides[i].name_peptide;

/*
if (peptides[i].name_peptide == "LEGSYHWYMEK")
{
break;
}
*/

(output[ii].peptide_names).push_back(peptides[i].name_peptide);

// Get number of theoretical spectrum
infile>>peptides[i].num_theo_spec;

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
infile>>temp_float;
(peptides[i].mz_theo).push_back(temp_float);
j++;
}

j = 0;

// Read off the theoretical intensity locations for current peptide
while (j < peptides[i].num_theo_spec)
{
infile>>temp_float;
(peptides[i].intensity_theo).push_back(temp_float);
//(output[ii].scores).at(i) = (output[ii].scores).at(i) + (peptides[i].intensity_theo).at(j);
j++;
}

i++;
}




// ------------------------------------------- END SANITY CHECK ------------------------------------------ //

if(ii<part_start || ii>=part_end){ii++; continue;}
// ------------------- Generate bipartite graphs for each observed spectrum - peptide pair ----------------------------- //

if(crux==1)
{
 vector<float>tempvector(first.intensity_obs);
    float sum;
    int num=first.intensity_obs.size();
    for(int mm=0;mm<num;mm++){
        sum = 0;
        for(int l=max(0,mm-75);l<min(num,mm+76);l++){
            sum+=tempvector.at(l);
        }
        first.intensity_obs.at(mm)=first.intensity_obs.at(mm)-sum/151.0;
    }
}
i = 0;


// ---------------------- Get scores (max weighted matchings) for each of the bipartite graphs ---------------------------------- //

i = 0;
float local_max=0;
// Make a transformation of indices and then feed the indices to Max_weight_bipartite_matching

//float temp_float4=0;
//vector<float>spectrum_foreground=first.intensity_obs;
//vector<float>spectrum_background;
//for(int jj=0;jj<spectrum_foreground.size();jj++){
//int m1=jj-75;
//int m2=jj+75;
//if(m1<0) m1=0;
//if(m2>=spectrum_foreground.size()) m2=spectrum_foreground.size()-1;
//float temp=0;
//for(int kk=m1;kk<=m2;kk++){
//temp+=spectrum_foreground.at(kk)/151.0;
//}
//spectrum_background.push_back(spectrum_foreground.at(jj)-temp);
//}
//vector<int> b_ion0;
//vector<int> y_ion0;
//int size_theo=peptides[0].mz_theo.size();
//for(int kk=0;kk<size_theo;kk++){
//int mze=(int)(peptides[0].mz_theo.at(kk));
//if(mze<10000) {b_ion0.push_back(mze);}
//else y_ion0.push_back(mze-10000);
//}
//
//std::sort (b_ion0.begin(), b_ion0.end());
//std::sort(y_ion0.begin(), y_ion0.end(), std::greater<int>());
//
//
//int mz_whole_by=b_ion0.at(0)+y_ion0.at(0);
//
int ii_sid=0;
for (int kk=0;kk<sid_bin.size();kk++){
	if(sid_bin[kk]==output[ii].sid){ ii_sid=kk;break;}

}

while (i < num_peptides)
{
int is_top=1;
//for(int kk=0;kk<60;kk++){if(top_index_i[ii][kk]==i) is_top=1;}
if(is_top==0){
output[ii].scores.at(i)==-100000;
i++;
continue;

}

int number_of_NH3_loss=0;//A b- or y-ion fragment can exhibit ammonia loss (-17 Da, -NH3) if it contains the residues K, R, Q or N.)
string temp_name=peptides[i].name_peptide;
//if(ii==31585&& i==0){
//ofstream spectrum_file;
//spectrum_file.open("output/spectrum_file.txt",ios::trunc);
//for(int kk=0;kk<first.intensity_obs.size();kk++){

//spectrum_file<<first.mz_obs.at(kk)<<"\t"<<first.intensity_obs.at(kk)<<endl;
//}
//}else{
//i=i+1;
//continue;
//}
vector<float> b_ion;
vector<float> y_ion;
vector<int> amino_acid;
vector<int> flanking_1;
vector<int> flanking_2;
for(int kk=0;kk<temp_name.length()-1;kk++){
amino_acid.push_back(temp_name.at(kk));
flanking_1.push_back(temp_name.at(kk));
flanking_2.push_back(temp_name.at(kk+1));
}
//cout<<temp_name<<endl;
for(int kk=0;kk<amino_acid.size();kk++){
int aa=amino_acid.at(kk);
if(aa=='K'||aa=='R'||aa=='Q'||aa=='N') number_of_NH3_loss++;
}
//cout<<number_of_NH3_loss<<endl;
if(i<num_peptides/2) number_NH3[number_of_NH3_loss]++;
double sum_nh3=0;
//for(int kk=0;kk<10;kk++) sum_nh3+=(double)number_NH3[kk];
//for(int kk=0;kk<10;kk++) cout<<(double)number_NH3[kk]/sum_nh3<<" ";
//cout<<endl;
int size_theo=peptides[i].mz_theo.size();
for(int kk=0;kk<size_theo;kk++){
float mze=(peptides[i].mz_theo.at(kk));
if(TMT_label==1){
if(mze<10000) {b_ion.push_back(mze);}
else{
	if (temp_name.at(temp_name.length()-1)=='K'){
	
   	y_ion.push_back(mze-10000.0);
	}else{
	
   	y_ion.push_back(mze-10000.0);
	}
}
}
else{
if(mze<10000) {b_ion.push_back(mze);}
else y_ion.push_back(mze-10000);
}
}



std::sort (b_ion.begin(), b_ion.end());
std::sort(y_ion.begin(), y_ion.end(), std::greater<int>());


if(all_score_larger_decoy_flag){
vector<float>amino_acid_seq;
	float temp_b_ion=0;
	for(int kk=0;kk<b_ion.size();kk++){
	amino_acid_seq.push_back(b_ion[kk]-temp_b_ion);
	temp_b_ion=b_ion[kk];
	}
amino_acid_seq.push_back(y_ion[y_ion.size()-1]);

std::random_shuffle ( amino_acid_seq.begin()+1, amino_acid_seq.end()-1 );

temp_b_ion=0;
for(int kk=0;kk<amino_acid_seq.size()-1;kk++){
b_ion[kk]=temp_b_ion+amino_acid_seq[kk];
temp_b_ion=b_ion[kk];
}

temp_b_ion=0;
for(int kk=0;kk<amino_acid_seq.size()-1;kk++){
y_ion[y_ion.size()-kk-1]=temp_b_ion+amino_acid_seq[amino_acid_seq.size()-kk-1];
temp_b_ion=y_ion[y_ion.size()-kk-1];
}

}

vector<int>A1;
vector<int>A2;
vector<float>W2;
vector<float>We;
vector<float>We2;
vector<int>B1;
vector<int>B2;
vector<int>ion;
vector<int>type;
B1.assign(2100,0);
B2.assign((int)(2*(b_ion.size()+y_ion.size())+1000),0);
//float factor_Nh3=other_p.at(3);
float factor_t=0.5;
float factor_t2=factor_t*factor_t;
float factor_t3=factor_t*factor_t2;

//float iso1[]={0.1101*factor_t,0.0225*factor_t2,0.0121*factor_t2,0.3128*factor_t,0.2364*factor_t,0.0784*factor_t2,0.0112*factor_t3,0.0107*factor_t2,0.0481*factor_t,0.6122,0.2514*factor_t,0.0511*factor_t2};

//float iso2[]={0.1364*factor_t,0.1179*factor_t,0.0345*factor_t2,1,0.4253*factor_t,0.0741*factor_t2};
float iso1[12];
float iso2[6];
for(int kk=0;kk<12;kk++) iso1[kk]=weight_b.at(kk);
for(int kk=0;kk<6;kk++) iso2[kk]=weight_y.at(kk);

float back_b=find_max(iso1,12,k2_matroid);
float back_y=find_max(iso2,6,k2_matroid);

float background=0;
int NH3_loss[]={'K','R','Q','N'};
int H2O_loss[]={'S','T','E','D'};
int Enlarge_bion[]={'K','R'};
float temp_back_b=1;
float temp_back_y=1.0*temp_back_b;
for(int kk=0;kk<b_ion.size();kk++){
	//decide if b has some loss
//	int n_NH3_loss=0;
//	int n_H2O_loss=0;
//	float n_Enlarge_bion=0;
//	for (int jj=0;jj<=kk;jj++){
//		for (int j2=0;j2<4;j2++){
//			if (NH3_loss[j2]==amino_acid[jj]) n_NH3_loss++;
//			if (H2O_loss[j2]==amino_acid[jj]) n_H2O_loss++;
//		}
//		for (int j2=0;j2<2;j2++){
//			if (Enlarge_bion[j2]==amino_acid[jj]) n_Enlarge_bion++; 
//		}
//	}
	int length_ion=kk+1;
	int length_peptide=b_ion.size()+1;
	float charge_factor[10];
	for(int jj=0;jj<10;jj++) charge_factor[10]=1.0;
	charged_ion_factor(charge_factor,charge,length_peptide, length_ion, l_s_factor);


for(int jj=0;jj<12;jj++) iso1[jj]=weight_b.at(jj);
//iso1[4]=n_NH3_loss*factor_Nh3*iso1[4];
//iso1[3]=n_H2O_loss*factor_Nh3*iso1[3];
//if(n_NH3_loss==0) iso1[4]=0;
//if(n_H2O_loss==0) iso1[3]=0;

float flanking_factor1=1.0;
back_b=flanking_factor1*find_max(iso1,12,k2_matroid);
back_b=temp_back_b*flanking_factor1;
int m1=b_ion.at(kk)-75;
int m2=b_ion.at(kk)+75;
if(m1<=0) m1=0;
if(m2>=first.intensity_obs.size()-1) m2=first.intensity_obs.size()-1;
for(int jj=m1;jj<=m2;jj++){
//background+=back_b*1/151.0*(first.intensity_obs.at(jj));
}
for(int jj=0;jj<obs_mz_bin[ii_sid].size();jj++){
if(abs(obs_mz_bin[ii_sid][jj]-b_ion[kk])<75) background+=charge_factor[1]*back_b/151.0*obs_int_bin[ii_sid][jj];


}
//add charge3;
if(add_ch3==1 && charge>=3){
	for(int charge_ion=2;charge_ion<charge;charge_ion++){
		float mz_charge_ion=(float)b_ion[kk];
		mz_charge_ion=(mz_charge_ion+(float)charge_ion-1)/((float)charge_ion);
		for(int jj=0;jj<obs_mz_bin[ii_sid].size();jj++){
		if(abs(obs_mz_bin[ii_sid][jj]-mz_charge_ion)<75) background+=charge_factor[charge_ion]*back_b/151.0*obs_int_bin[ii_sid][jj];
		}
	}	
}	
//add charge3 end;
}
for(int kk=0;kk<y_ion.size();kk++){
//	//decide if y has some loss
//	int n_NH3_loss=0;
//	int n_H2O_loss=1;
//	for (int jj=kk+1;jj<y_ion.size()+1;jj++){
//		for (int j2=0;j2<4;j2++){
//			if (NH3_loss[j2]==amino_acid[jj]) n_NH3_loss++;
//			if (H2O_loss[j2]==amino_acid[jj]) n_H2O_loss++;
//		}
//		
//	}
////	cout<<"y"<<n_NH3_loss<<" "<<n_H2O_loss<<endl;
//	//end
	int length_ion=y_ion.size()-kk;
	int length_peptide=y_ion.size()+1;
	float charge_factor[10];
	for(int jj=0;jj<10;jj++) charge_factor[10]=1.0;
	charged_ion_factor(charge_factor,charge,length_peptide, length_ion, l_s_factor);


for(int kk2=0;kk2<6;kk2++) iso2[kk2]=weight_y.at(kk2);

//iso2[1]=n_NH3_loss*factor_Nh3*iso2[1];
//iso2[0]=n_H2O_loss*factor_Nh3*iso2[0];
//if(n_NH3_loss==0) iso2[1]=0;
//if(n_H2O_loss==0) iso2[0]=0;


float flanking_factor1=1.0;
back_y=flanking_factor1*find_max(iso2,6,k2_matroid);
back_y=flanking_factor1*temp_back_y;
int m1=y_ion.at(kk)-75;
int m2=y_ion.at(kk)+75;
if(m1<=0) m1=0;
if(m2>=first.intensity_obs.size()-1) m2=first.intensity_obs.size()-1;
for(int jj=m1;jj<=m2;jj++){
//background+=back_y*1/151.0*(first.intensity_obs.at(jj));
}
for(int jj=0;jj<obs_mz_bin[ii_sid].size();jj++){
if(abs(obs_mz_bin[ii_sid][jj]-y_ion[kk])<75) background+=charge_factor[1]*back_y*1.0/151.0*obs_int_bin[ii_sid][jj];
}
//add charge3;
if(add_ch3==1 && charge>=3){
	for(int charge_ion=2;charge_ion<charge;charge_ion++){
		float mz_charge_ion=(float)y_ion[kk];
		mz_charge_ion=(mz_charge_ion+(float)charge_ion-1)/((float)charge_ion);
		for(int jj=0;jj<obs_mz_bin[ii_sid].size();jj++){
		if(abs(obs_mz_bin[ii_sid][jj]-mz_charge_ion)<75) background+=charge_factor[charge_ion]*back_y/151.0*obs_int_bin[ii_sid][jj];
		}
	}	
}
//add charge3 end;
}

for (int kk=0;kk<b_ion.size();kk++){

float temp1=kk;
float temp2=b_ion.size();
	//float ion_length_factor=4.5-4.0*pow(2.0*((temp1+1-temp2/2.0)/temp2),2);
	float ion_length_factor=0.9+0.2*temp1/temp2;
	ion_length_factor=1.0;

//	if (kk==0 || kk==y_ion.size()-1) ion_length_factor=0.1;

		
	//decide if b has some loss
//	int n_NH3_loss=0;
//	int n_H2O_loss=0;
//	float n_Enlarge_bion=0;
//	for (int jj=0;jj<=kk;jj++){
//		for (int j2=0;j2<4;j2++){
//			if (NH3_loss[j2]==amino_acid[jj]) n_NH3_loss++;
//			if (H2O_loss[j2]==amino_acid[jj]) n_H2O_loss++;
//		}
//		for (int j2=0;j2<2;j2++){
//			if (Enlarge_bion[j2]==amino_acid[jj]) n_Enlarge_bion++;
//		}
//		
//	}
//	cout<<"b"<<n_NH3_loss<<" "<<n_H2O_loss<<endl;
	//end
float mze=b_ion.at(kk);
if(mze>=first.intensity_obs.size()) continue;
int n_iso=12;
int iso_mz[]={-28,-27,-19,-18,-17,-16,-15,-12,-1,0,1,2};

float iso[]={0.1101*factor_t,0.0225*factor_t2,0.0121*factor_t2,0.3128*factor_t,0.2364*factor_t,0.0784*factor_t2,0.0112*factor_t3,0.0107*factor_t2,0.0481*factor_t,0.6122,0.2514*factor_t,0.0511*factor_t2};

float flanking_factor1=1.0;
for(int jj=0;jj<12;jj++) iso[jj]=flanking_factor1*weight_b.at(jj);
if(crux==1){
float iso_xcorr[]={0.2,0,0,0.2,0.2,0,0,0,0,1,0,0};
for(int jj=0;jj<n_iso;jj++) iso[jj]=iso_xcorr[jj];
}
//iso[4]=n_NH3_loss*factor_Nh3*iso[4];
//iso[3]=n_H2O_loss*factor_Nh3*iso[3];
//if(n_NH3_loss==0) iso[4]=0;
//if(n_H2O_loss==0) iso[3]=0;
//float iso_sum=0;
//for(int jj=0;jj<12;jj++) iso_sum+=iso[jj];
//cout<<"b"<<iso_sum<<endl;
//for low res
//for(int jj=0;jj<n_iso;jj++){
//if(iso[jj]<0.001) continue;
//if(mze+iso_mz[jj]<first.intensity_obs.size()-1 && mze+iso_mz[jj]>=0&&first.intensity_obs.at(mze+iso_mz[jj])>0){
//We.push_back(iso[jj]*first.intensity_obs.at(mze+iso_mz[jj]));
//A1.push_back(mze+iso_mz[jj]);
//A2.push_back(kk);
//ion.push_back(0);
//type.push_back(jj-9);
//}
//}
//end for low res

	int length_ion=kk+1;
	int length_peptide=b_ion.size()+1;
	float charge_factor[10];
	for(int jj=0;jj<10;jj++) charge_factor[10]=1.0;
	charged_ion_factor(charge_factor,charge,length_peptide, length_ion, l_s_factor);


for(int jj=0;jj<obs_mz_bin[ii_sid].size();jj++){
if(abs(obs_mz_bin[ii_sid][jj]-mze)<30.0){
	int diff=(int)((obs_mz_bin[ii_sid][jj]-mze)*100)+3000;
	if (diff<0) continue;

We2.push_back(ion_length_factor*charge_factor[1]*obs_int_bin[ii_sid][jj]*flanking_factor1*weight_high_b_ion[diff]);
	if(weight_high_b_ion[diff]<0.01)continue;
We.push_back(ion_length_factor*charge_factor[1]*obs_int_bin[ii_sid][jj]*flanking_factor1*weight_high_b_ion[diff]);
int n_or_c=0;
if (diff<2800) n_or_c=1;
//A1.push_back(jj);

A1.push_back((int)(obs_mz_bin[ii_sid][jj]-bin_offset));
A2.push_back(2*kk);
ion.push_back(0);
type.push_back(diff-3000);

}


}

//add charge3
if(add_ch3==1 && charge>=3){
	for(int charge_ion=2;charge_ion<charge;charge_ion++){
		float mz_charge_ion=mze;
		mz_charge_ion=(mz_charge_ion+(float)charge_ion-1)/((float)charge_ion);
		for(int jj=0;jj<obs_mz_bin[ii_sid].size();jj++){
			if(abs(obs_mz_bin[ii_sid][jj]-mz_charge_ion)<30.0){
				int diff=(int)((obs_mz_bin[ii_sid][jj]-mz_charge_ion)*100)+3000;

			diff=(int)((obs_mz_bin[ii_sid][jj]*(float)(charge_ion)+1.0-(float)(charge_ion)-mze)*100.0)+3000;
				if (diff<0) continue;

				if (diff>=weight_high_y_ion.size()) continue;

				We2.push_back(ion_length_factor*charge_factor[charge_ion]*obs_int_bin[ii_sid][jj]*flanking_factor1*weight_high_b_ion[diff]);
				if(weight_high_b_ion[diff]<0.01)continue;
				We.push_back(ion_length_factor*charge_factor[charge_ion]*obs_int_bin[ii_sid][jj]*flanking_factor1*weight_high_b_ion[diff]);
				int n_or_c=0;
				if (diff<2800) n_or_c=1;
				//A1.push_back(jj);

				A1.push_back((int)(obs_mz_bin[ii_sid][jj]-bin_offset));
				A2.push_back(2*kk+1);
				ion.push_back(0);
				type.push_back(diff-3000);

			}

		}
	}
}
//add charge3 end
}
for (int kk=0;kk<y_ion.size();kk++){
float temp1=y_ion.size()-kk;
float temp2=y_ion.size();
	//float ion_length_factor=1.5-2.0*abs(temp1+1-temp2/2.0)/temp2;

//	float ion_length_factor=4.5-4.0*abs(2.0*((temp1+1-temp2/2.0)/temp2));
//
	float ion_length_factor=0.9+0.2*temp1/temp2;

	ion_length_factor=1.0;	
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

float mze=y_ion.at(kk);
if(mze>=first.intensity_obs.size()) continue;
int nb=b_ion.size();
int n_iso=6;
int iso_mz[]={-18,-17,-16,0,1,2};
float iso[]={0.1364*factor_t,0.1179*factor_t,0.0345*factor_t2,1,0.4253*factor_t,0.0741*factor_t2};
float flanking_factor1=1.0;
for(int jj=0;jj<6;jj++) iso[jj]=flanking_factor1*weight_y.at(jj);
if(crux==1){
float iso_xcorr[]={0.2,0.2,0,1,0,0};
for(int jj=0;jj<n_iso;jj++) iso[jj]=iso_xcorr[jj];
}
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
	int length_ion=y_ion.size()-kk;
	int length_peptide=y_ion.size()+1;
	float charge_factor[10];
	for(int jj=0;jj<10;jj++) charge_factor[10]=1.0;
	charged_ion_factor(charge_factor,charge,length_peptide, length_ion, l_s_factor);


for(int jj=0;jj<obs_mz_bin[ii_sid].size();jj++){
if(abs(obs_mz_bin[ii_sid][jj]-mze)<30.0){
//	int diff=(int)((obs_mz_bin[ii_sid][jj]-mze)*100)+3000;
	int diff=(int)((obs_mz_bin[ii_sid][jj]-mze)*100)+3000;
	if (diff<0) continue;

int n_or_c=0;
if (diff<2800) n_or_c=1;

We2.push_back(ion_length_factor*charge_factor[1]*obs_int_bin[ii_sid][jj]*flanking_factor1*weight_high_y_ion[diff]);
	if(weight_high_y_ion[diff]<0.01)continue;
We.push_back(ion_length_factor*charge_factor[1]*obs_int_bin[ii_sid][jj]*flanking_factor1*weight_high_y_ion[diff]);
A1.push_back((int)(obs_mz_bin[ii_sid][jj]-bin_offset));
A2.push_back(2*(kk+nb)+1);
ion.push_back(1);
type.push_back(diff-3000);

}


}

//add charge3
if(add_ch3==1 && charge>=3){
	for(int charge_ion=2;charge_ion<charge;charge_ion++){
		float mz_charge_ion=mze;
		mz_charge_ion=(mz_charge_ion+(float)charge_ion-1)/((float)charge_ion);
		for(int jj=0;jj<obs_mz_bin[ii_sid].size();jj++){
			if(abs(obs_mz_bin[ii_sid][jj]-mz_charge_ion)<30.0){
				int diff=(int)((obs_mz_bin[ii_sid][jj]-mz_charge_ion)*100)+3000;
				diff=(int)((obs_mz_bin[ii_sid][jj]*(float)(charge_ion)+1.0-(float)(charge_ion)-(mze))*100.0)+3000;
				if (diff<0) continue;
				if (diff>=weight_high_y_ion.size()) continue;

				We2.push_back(ion_length_factor*charge_factor[charge_ion]*obs_int_bin[ii_sid][jj]*flanking_factor1*weight_high_y_ion[diff]);
				if(weight_high_y_ion[diff]<0.01)continue;
				We.push_back(ion_length_factor*charge_factor[charge_ion]*obs_int_bin[ii_sid][jj]*flanking_factor1*weight_high_y_ion[diff]);
				int n_or_c=0;
				if (diff<2800) n_or_c=1;
				//A1.push_back(jj);

				A1.push_back((int)(obs_mz_bin[ii_sid][jj]-bin_offset));
				A2.push_back(2*(kk+nb));
				ion.push_back(0);
				type.push_back(diff-3000);

			}

		}
	}
}
//add charge3 end
}
//print all edgeweight
//for(int kk=0;kk<We.size();kk++){
//cout<<We.at(kk)<<" ";

//}
//
//count degree//////////
//for(int kk=1;kk<16;kk++) cout<<degree_a2[kk]<<" ";
//for(int kk=1;kk<16;kk++) sum_degree+=(float)degree_a2[kk];
//cout<<endl;
//for(int kk=1;kk<16;kk++) cout<<((float)(degree_a2[kk]))/sum_degree<<" ";
//cout<<endl;
//count degree end//////
float lambda3=0.9;
lambda3=0.9;
//lambda3=1.0;

sub_sup all_graph(A1,A2,We,10,2,2000,2*b_ion.size()+2*y_ion.size(),other_p[3],other_p[6],other_p[7]);
if(0){

	cout<<ii<<" "<<i<<endl;
//cout<<"------------\n";
// cout<<"scores: "<<all_graph.GreedyMax()<<endl;
//A.print_current_matching();
//	cout<<"scores: "<<all_graph.sub_sup_process(1)<<endl;
//	cout<<"scores: "<<all_graph.sub_sup_process(2)<<endl;
//	cout<<"scores: "<<all_graph.sub_sup_process(3)<<endl;
//	cout<<"scores: "<<all_graph.submodualr_projection(2.0)<<endl;
//outs<<all_graph.optimal_scores();
	outs<<all_graph.GreedyMax()<<"\t"<<all_graph.submodualr_projection(2.0)<<"\t"<<all_graph.sub_sup_process(1)<<"\t"<<all_graph.sub_sup_process(2)<<"\t"<<all_graph.sub_sup_process(3)<<"\t"<<all_graph.random_scores()<<endl;
}


//cout<<temp_float3<<" "<<score_ground<<" "<<background<<" "<<output[ii].scores.at(i)<<endl;
//(output[ii].scores).at(i)=4.5*(temp_float3/0.72-background)+0.5*(temp_float5-background);
//(output[ii].scores).at(i)=0*(temp_float3/0.72-other_p.at(0)*background)+other_p.at(1)*(temp_float5-other_p.at(0)*background);
float lambda=other_p.at(3);
//lambda=1.0;
lambda=1;
//(output[ii].scores).at(i)=2.0*pow(pow(f_groundset2/0.78,-1)+pow(f_groundset3/1.05,-1),-1)*(temp_float3/f_groundset2*lambda+(1-lambda)*temp_float5/f_groundset3)-background;
//(output[ii].scores).at(i)=temp_float3/0.78*lambda+temp_float5/1.05*(1-lambda)-0.5*background;
//(output[ii].scores).at(i)=temp_float3/0.78-0*background+modular_score_empty_2;

//(output[ii].scores).at(i)=all_graph.scores(0);
//
//

//(output[ii].parameter).at(i)=all_graph.diversity();
output[ii].edges[i]=We;

(output[ii].scores).at(i)=all_graph.GreedyMax();
sum_all_scores+=output[ii].scores[i];
number_all_scores+=1;
//if(ii%100==0 && i==10) cout<<sum_all_scores/number_all_scores<<"----"<<number_all_scores<<endl;


//(output[ii].scores).at(i)=10;
//if (A1.size()<25 &&A1.size()>20){
//	cout<<"---------\n";
//	cout<<A1.size()<<endl;
//	cout<<"=========\n";
//cout<<"greed:"<<all_graph.GreedyMax()<<endl;
//cout<<"opt"<<all_graph.optimal_scores()<<endl;
//}
//(output[ii].scores).at(i)=all_graph.submodualr_projection(other_p[6]);


//cout<<temp_float3<<" "<<modular_score_empty_2<<endl;
//(output[ii].scores).at(i)=temp_float3*lambda+temp_float5/0.9*(1-lambda)-background;
//(output[ii].scores).at(i)=temp_float3*lambda+temp_float5/0.9*(1-lambda);
//cout<<"-------\n";
//cout<<temp_float3<<" "<<temp_float5<<endl;
//cout<<(output[ii].scores).at(i)<<" "<<background<<endl;
//(output[ii].scores).at(i)=temp_float3-background;


//(output[ii].scores).at(i)=temp_float3*lambda+temp_float5/1.05*(1-lambda)-1.0*background;
//output[ii].scores.at(i)=temp_float3*sqrt(all_sum)*107.0/161.0/0.78-background;
//output[ii].scores.at(i)=temp_float3*sqrt(all_sum)*107.0/165.0/0.78-background;
//output[ii].scores.at(i)=temp_float3;
//output[ii].scores.at(i)=temp_float3*107.0/53.0/0.78-background;


//(output[ii].scores).at(i)=(lambda2*f_groundset2/0.72+(1-lambda1)*f_groundset3/1.2)*(temp_float3/f_groundset2*lambda+(1-lambda)*temp_float5/f_groundset3)-background;
//cout<<temp_float3<<" "<<f_groundset2<<" "<<endl;
//cout<<temp_float5<<" "<<f_groundset3<<" "<<endl;
//cout<<"---------------\n";
//(output[ii].scores).at(i)=(temp_float3-background);



i++;
}
//dealing with calbration
//should be uncommented

vector<float>scoreset2=output[ii].scores;
vector<float>scoreset;
for(i=0;i<scoreset2.size();i++) scoreset.push_back(scoreset2.at(i));
float sum = std::accumulate(scoreset.begin(), scoreset.end(), 0.0);
float mean = sum / scoreset.size();

float sq_sum = std::inner_product(scoreset.begin(), scoreset.end(), scoreset.begin(), 0.0);
float stdev = std::sqrt((sq_sum / scoreset.size() - mean * mean));
for(int i=0; i<scoreset2.size();i++){
//	(output[ii].scores).at(i)=(scoreset2.at(i)-mean)/stdev;
float t=(scoreset2.at(i));	
//cout<<t<<endl;
//float t2=(float)scoreset2.size();
//float t3=(0.5-0.5*erf(t/sqrt(2.0)));
//(output[ii].scores).at(i)=-t2*t3+(t2*(t2-1))*t3*t3/2.0;
(output[ii].scores).at(i)=(t-mean);
//(output[ii].scores).at(i)=(t);
}
i = 0;
ii=ii+1;
} // End spectra while loop

double number_cand=0;

for(int kk=0;kk<n;kk++){
	number_cand+=(float)(output[kk].scores.size());

}
cout<<"average:"<<number_cand/((float)(n));

//ofstream twodecoy1;
//ofstream twodecoy2;
//twodecoy1.open("output/decoy1",ios::trunc);
//twodecoy2.open("output/decoy2",ios::trunc);


 // End outer if - open infile

// Print output scores to file

//cout<<"check"<<endl;
/*
ofstream top_index;
top_index.open("/s1/wrbai/codes/output/top_index.txt",ios::trunc);
for(int kk=0;kk<n;kk++){
	vector<float_int_s> target_s;
	int lens=output[kk].scores.size();
	for(int jj=0;jj<lens/2;jj++){
		float_int_s b;
		b.s=output[kk].scores.at(jj);
		b.a=jj;
		target_s.push_back(b);
	}
	vector<float_int_s> decoy_s;
	for(int jj=lens/2;jj<lens;jj++){
		float_int_s b;
		b.s=output[kk].scores.at(jj);
		b.a=jj;
		decoy_s.push_back(b);
	}

	
	sort(target_s.begin(),target_s.end(),sortBys_float_int);
	sort(decoy_s.begin(),decoy_s.end(),sortBys_float_int);
	for(int jj=0;jj<30;jj++) top_index<<target_s.at(jj).a<<" ";
	for(int jj=0;jj<30;jj++) top_index<<decoy_s.at(jj).a<<" ";
	top_index<<endl;
}
top_index.close();
*/

if(all_score_output_flag){
int ii=0;

	all_scores_file<<"mean\tstdev\tn"<<endl;
while(ii < n){

	vector<float>scoreset=output[ii].scores;
	float sum = std::accumulate(scoreset.begin(), scoreset.end(), 0.0);
	float mean = sum / scoreset.size();
	
	float sq_sum = std::inner_product(scoreset.begin(), scoreset.end(), scoreset.begin(), 0.0);
	float stdev = std::sqrt((sq_sum / scoreset.size() - mean * mean));
	all_scores_file<<mean<<"\t"<<stdev<<"\t"<<output[ii].scores.size()<<endl;
	ii++;
}
}

if(all_score_larger_decoy_flag){
int ii=0;

	all_scores_file_other_decoy<<"mean\tstdev\tn"<<endl;
while(ii < n){

	vector<float>scoreset=output[ii].scores;
	float sum = std::accumulate(scoreset.begin(), scoreset.end(), 0.0);
	float mean = sum / scoreset.size();
	
	float sq_sum = std::inner_product(scoreset.begin(), scoreset.end(), scoreset.begin(), 0.0);
	float stdev = std::sqrt((sq_sum / scoreset.size() - mean * mean));
	all_scores_file_other_decoy<<mean<<"\t"<<stdev<<"\t"<<output[ii].scores.size()<<endl;
	ii++;
}
}

vector<float_int_s>all_scores;
//if (scores_file.is_open())
//{
//scores_file<<"SID"<<"	"<<"Peptide"<<"	"<<"Score"<<"	"<<endl;

if(calculate_by_part<=1) scores_file2<<"Kind\tSid\tPeptide\tScore\n";

ofstream para_file2;
ofstream para_file3;
para_file2.open("output/para_temp.txt");
para_file3.open("output/para2.txt");
ii = 0;
//cout<<"ii = "<<ii<<endl;
while (ii < n)
{

	if(ii<part_start || ii>=part_end){ii++; continue;}
int j = 0;
//****************2014-10-26*****************************//

int lens=output[ii].scores.size();
vector<float> target(output[ii].scores.begin(),output[ii].scores.begin()+lens/2);
vector<float> decoy(output[ii].scores.begin()+lens/2,output[ii].scores.end());
float maxt=-10000;


//if(target.size()<10) cout<<lens<<" "<<target.size()<<" " <<decoy.size()<<endl;

 
int maxt_i=0;
for(int kk=0;kk<target.size();kk++){
if(target.at(kk)>maxt){ 
maxt=target.at(kk);
maxt_i=kk;
}
}
float maxd=-10000;
int maxd_i=0;
for(int kk=0;kk<decoy.size();kk++){
if(decoy.at(kk)>maxd){
maxd=decoy.at(kk);
maxd_i=kk;
}
}
if(maxt>=maxd) maxd=-1000000000;
if(maxt<maxd) maxt=-1000000000;
scores_file2<<"t\t"<<output[ii].sid<<"\t"<<(output[ii].peptide_names).at(maxt_i)<<"\t"<<maxt<<endl;
scores_file2<<"d\t"<<output[ii].sid<<"\t"<<(output[ii].peptide_names).at(maxd_i+lens/2)<<"\t"<<maxd<<endl;

if(maxt>=maxd) {

//para_file2<<output[ii].parameter[maxt_i]<<"\t"<<-10<<"\t"<<maxt<<endl;
para_file2<<ii<<"\t"<<maxt_i<<"\t"<<maxt<<endl;
if(maxt>162 && maxt<232){

para_file3<<"t\t"<<output[ii].edges[maxt_i].size()<<"\t"<<maxt<<endl;
for(int kk=0;kk<output[ii].edges[maxt_i].size();kk++){
para_file3<<output[ii].edges[maxt_i][kk]<<"\t";
}
para_file3<<endl;
}
}
if(maxt<maxd){

// para_file2<<-10<<"\t"<<output[ii].parameter[maxd_i+lens/2]<<"\t"<<maxd<<endl;
if(maxd>162){
para_file3<<"d\t"<<output[ii].edges[maxd_i+lens/2].size()<<"\t"<<maxd<<endl;
for(int kk=0;kk<output[ii].edges[maxd_i+lens/2].size();kk++){
para_file3<<output[ii].edges[maxd_i+lens/2][kk]<<"\t";
}
para_file3<<endl;
}

 }
float_int_s temp_float_int1;
temp_float_int1.s=maxt;
temp_float_int1.a=1;
all_scores.push_back(temp_float_int1);
float_int_s temp_float_int2;
temp_float_int2.s=maxd;
temp_float_int2.a=2;
all_scores.push_back(temp_float_int2);


//float maxd1=-10000;
//int nd=decoy.size();
//for(int kk=0;kk<nd/2;kk++){
//if(decoy.at(kk)>maxd1) maxd1=decoy.at(kk);
//p
//}
//float maxd2=-10000;
//for(int kk=nd/2;kk<nd;kk++){
//if(decoy.at(kk)>maxd2) maxd2=decoy.at(kk);
//}
//twodecoy1<<maxd1<<endl;
//twodecoy2<<maxd2<<endl;

//****************2014-10-26*****************************//

ii++;
}


sort(all_scores.begin(),all_scores.end(),sortBys_float_int);
//for(int kk=0;kk<100;kk++) cout<<all_scores.at(kk).s<<" "<<all_scores.at(kk).a<<endl;
float q[11];
for(int kk=0;kk<11;kk++) q[kk]=float(kk)/100;
float area=0;
float high_t=0;
float high_d=0;
int j=0;
for(int kk=0;kk<all_scores.size();kk++){
float s1=all_scores.at(kk).s;
float a1=all_scores.at(kk).a;
if(a1==1) high_t+=1;
else high_d+=1;
if(high_d/(float(kk))>q[j]) {j=j+1;area+=high_t;if(j>10) break;}



}
// Close any remaining files
//scores_file.close();
scores_file2.close();
infile.close();
//edges_file.close();
return area;
}
