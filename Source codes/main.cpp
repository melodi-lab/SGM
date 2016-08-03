#include<iostream>
#include<fstream>
#include<vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h> 
#include <cstdlib>
using namespace std;
float calculate_abs_ranking_ch2(vector<float> weight_b,vector<float> weight_y,vector<float>other_p,int k2_input,int kind,int index, int charge,vector<string>filenames,int n, int nn);
int main(int argc, char *argv[]){
	srand (time(NULL));
	vector<float> weight_b;
	vector<float> weight_b1;
	vector<float> weight_y;
	vector<float> weight_y1;
	float factor_t=0.5;
	float factor_t2=factor_t*factor_t;
	float factor_t3=factor_t2*factor_t;
	//int iso_mz_b[]={-28,-27,-19,-18,-17,-16,-15,-12,-1,0,1,2};
	float iso_b[]={0.1101*factor_t,0.0225*factor_t2,0.0121*factor_t2,0.3128*factor_t,0.2364*factor_t,0.0784*factor_t2,0.0112*factor_t3,0.0107*factor_t2,0.0481*factor_t,0.6122,0.2514*factor_t,0.0511*factor_t2};
//	float iso_b[]={0.1,0.02,0.01,0.1,0.1,0.01,0.01,0.01,0.01,0.5,0.1,0.01};
//	int iso_mz_y[]={-18,-17,-16,0,1,2};
	float iso_y[]={0.1364*factor_t,0.1179*factor_t,0.0345*factor_t2,1,0.4253*factor_t,0.0741*factor_t2};
//	float iso_y[]={0.2,0.2,0.01,1.0,0.2,0.01};
	for (int i=0;i<12;i++) weight_b.push_back(0.9*iso_b[i]);
	for (int i=0;i<6;i++) weight_y.push_back(iso_y[i]);
	ifstream input_file;
//	input_file.open("/s1/wrbai/codes/output/weight_learnt.txt");
//	for (int i=0;i<12;i++) input_file>>weight_b.at(i);
//	for (int i=0;i<6;i++) input_file>>weight_y.at(i);
	vector<float> other_p;
	other_p.push_back(1.0);//background
	other_p.push_back(1.0);//two
	other_p.push_back(1.0);//log
	other_p.push_back(1.0);//lambda1
	other_p.push_back(0.0);//test index; default 0 if not doing test
	other_p.push_back(0.0);//part index; default 0 if no partition 
	other_p.push_back(1.0);//lambda2
	other_p.push_back(1.0);//lambda3

	float lf=atof(argv[4]);

	other_p.at(3)=lf;

	lf=atof(argv[7]);
	other_p.at(6)=lf;

	lf=atof(argv[14]);
	other_p.at(7)=lf;


	lf=(float)atoi(argv[5]);
	other_p[4]=lf;
	cout<<other_p[4]<<"asdf\n";
	int part_arg=atoi(argv[6]);
	other_p.at(5)=(float)part_arg;
//	for (int i=0;i<3;i++) input_file>>other_p.at(i);
	weight_b1=weight_b;
	weight_y1=weight_y;
	vector<float> other_p1=other_p;
//	int k2_input=2;
	int k2_input=(int)(other_p[4]);
	k2_input=2;
	cout<<k2_input<<endl;
	float q=0;
	int kind=1;
	int index=1;
	int charge=2;
//	kind=1;
//	index=1;
//	q=calculate_abs_ranking_ch2(weight_b,weight_y,other_p,k2_input,kind,index);
//	kind=1;
//	index=2;
//	q=calculate_abs_ranking_ch2(weight_b,weight_y,other_p,k2_input,kind,index);
//	kind=1;
//	index=3;
//	q=calculate_abs_ranking_ch2(weight_b,weight_y,other_p,k2_input,kind,index);
//	kind=1;
//	index=4;
//	q=calculate_abs_ranking_ch2(weight_b,weight_y,other_p,k2_input,kind,index);
 	
 	kind=atoi(argv[1]);
	index=atoi(argv[2]);
	charge=atoi(argv[3]);
	vector<string> filenames;
	filenames.push_back(std::string(argv[8]));
	filenames.push_back(std::string(argv[9]));
	filenames.push_back(std::string(argv[10]));
	filenames.push_back(std::string(argv[11]));
	int n = atoi(argv[12]);
	n=(n-1)/2;
	int nn = atoi(argv[13]);
	q=calculate_abs_ranking_ch2(weight_b,weight_y,other_p,k2_input,kind,index,charge,filenames,n,nn);

//	kind=0;
//	index=2;
//	q=calculate_abs_ranking_ch2(weight_b,weight_y,other_p,k2_input,kind,index);
//	kind=0;
//	index=3;
//	q=calculate_abs_ranking_ch2(weight_b,weight_y,other_p,k2_input,kind,index);
//  	kind=0;
//	index=4;
//	q=calculate_abs_ranking_ch2(weight_b,weight_y,other_p,k2_input,kind,index);
  



/*	
	for(int t=0;t<200;t++){
	for(int i=0;i<weight_b.size();i++){
		float r=((float)rand())/((float)RAND_MAX);
	 weight_b1.at(i)=weight_b.at(i)*(1+(r-0.5)/10);
	}
	for(int i=0;i<weight_y.size();i++){
		float r=((float)rand())/((float)RAND_MAX);
	 weight_y1.at(i)=weight_y.at(i)*(1+(r-0.5)/10);
	}
	for(int i=0;i<other_p.size();i++){

		float r=((float)rand())/((float)RAND_MAX);
	 other_p1.at(i)=other_p.at(i)*(1+(r-0.5)/10);


	}
	
		float r=((float)rand())/((float)RAND_MAX);
		k2_input1=k2_input;
		if(r<0.1) k2_input1=k2_input+1;
		if(r>0.9) k2_input1=k2_input-1;
		if(k2_input1<1) k2_input1=1;

	float	q0=calculate_abs_ranking_ch2(weight_b1,weight_y1,other_p1,k2_input1,kind,index);

	cout<<q<<" "<<q0<<endl;
	if(q0>q) {
	q=q0;
	weight_b=weight_b1;
	weight_y=weight_y1;
	other_p=other_p1;
	k2_input=k2_input1;
		for(int kk=0;kk<weight_b.size();kk++) cout<<weight_b.at(kk)<<" ";
	cout<<endl;
	for(int kk=0;kk<weight_y.size();kk++) cout<<weight_y.at(kk)<<" ";
	cout<<endl;
	for(int kk=0;kk<other_p.size();kk++) cout<<other_p.at(kk)<<" ";
	cout<<endl;
	cout<<k2_input<<endl;
	}
	}
	for(int kk=0;kk<weight_b.size();kk++) cout<<weight_b.at(kk)<<" ";
	cout<<endl;
	for(int kk=0;kk<weight_y.size();kk++) cout<<weight_y.at(kk)<<" ";
	cout<<endl;
	for(int kk=0;kk<other_p.size();kk++) cout<<other_p.at(kk)<<" ";
	cout<<endl;
	ofstream outfile;
	outfile.open("/s1/wrbai/codes/output/weight_learnt.txt",ios::trunc);
	for(int kk=0;kk<weight_b.size();kk++) outfile<<weight_b.at(kk)<<" ";
	outfile<<endl;
	
	for(int kk=0;kk<weight_y.size();kk++) outfile<<weight_y.at(kk)<<" ";

	outfile<<endl;

	for(int kk=0;kk<other_p.size();kk++) outfile<<other_p.at(kk)<<" ";
	outfile<<endl;
	outfile<<k2_input<<endl;
	*/
}
