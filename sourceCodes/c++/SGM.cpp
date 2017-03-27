#include<iostream>
#include<fstream>
#include<vector>
#include<string>

#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <cstdlib>
#include <unistd.h>
#include "input_para.h"
using namespace std;
float calculate_abs_ranking_ch2(input_para paras);
float random_float(float m1, float m2) {
  float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  return m1 + (m2 - m1) * r;
}
void showhelpinfo(char *s)
{
  cout << "Usage:   " << s << " [-option] [argument]" << endl;
  cout << "option:  " << "-h  show help information" << endl;
  cout << "         " << "-c  charge" << endl;
  cout << "         " << "-i  input file" << endl;
  cout << "         " << "-o  output file" << endl;
  cout << "         " << "-m  ms2file" << endl;
}
int main(int argc, char *argv[]) {
  srand (time(NULL) * (getpid() + 1));

  input_para paras;
  paras.k1 = 10;
  paras.k2 = 2;

  /////////////////
  char tmp;
  /*if the program is ran witout options ,it will show the usgage and exit*/
  if (argc == 1)
  {
    showhelpinfo(argv[0]);
    exit(1);
  }
  /*use function getopt to get the arguments with option."hu:p:s:v" indicate
  that option h,v are the options without arguments while u,p,s are the
  options with arguments*/
  char* input = "";
  char* output = "";
  char* ms2file = "";
  char* n_s ="";
  char* nn_s = "";
  int set_id = 0;
  while ((tmp = getopt(argc, argv, "hc:i:o:m:n:N:p:")) != -1)
  {
    switch (tmp)
    {

    case 'h':
      showhelpinfo(argv[0]);
      break;

    case 'c':
      cout << "charge: " << optarg << endl;
      paras.charge = atoi(optarg);
      break;

    case 'i':
      cout << "input file: " << optarg << endl;
      input = optarg;

      break;

    case 'o':
      cout << "output file: " << optarg << endl;
      output = optarg;
      break;

    case 'm':
      cout << "ms2 file: " << optarg << endl;
      ms2file = optarg;
      break;

	case 'n':
      cout << "number of spectral: " << optarg << endl;
      nn_s = optarg;
      break;

	case 'N':
      cout << "number of spectral with at least one target spectrum: " << optarg << endl;
      n_s = optarg;
      break;


    default:
      showhelpinfo(argv[0]);
      break;
    }
  }
  string str1(input);
  paras.input_file = str1;
  string str2(output);
  paras.output_file = str2;
  string str3(ms2file);
  paras.input_ms2_file = str3;
  ////////////////
  paras.n = atoi(n_s);
  paras.nn = atoi(nn_s);
 	cout<<paras.n<<endl;
 	cout<<paras.nn<<endl;

  paras.lambda_vector.assign(8, 0);
  //lambda_vector[0]=random_float(0.0,1.0);
  //lambda_vector[1]=random_float(0.0,1.0);
  //lambda_vector[2]=(float)(rand()%4);
  //lambda_vector[3]=random_float(0.4,1.0);
  //lambda_vector[4]=(float)(rand()%4);
  //lambda_vector[5]=random_float(0.4,1.0);
  //lambda_vector[6]=random_float(0,40);
  //lambda_vector[7]=random_float(0,80);
  //other_p[3]=lambda_vector[0];
  //other_p[6]=lambda_vector[2];
  //other_p[7]=lambda_vector[3];
  ifstream fin_file;
  fin_file.open("/n/trombone/s1/wrbai/codes/pipeline/parameter_sets.txt");

  cout << "set_id " << set_id << endl;
  float temp_int = 10;
  for (int ii = 0; ii < paras.lambda_vector.size()*set_id; ii++) {
    fin_file >> temp_int;
  }
  for (int ii = 0; ii < paras.lambda_vector.size(); ii++) {
    fin_file >> temp_int;
    paras.lambda_vector[ii] = temp_int;
    cout << paras.lambda_vector[ii] << "---";
  }

  cout << "\n";

  paras.lambda_vector[2] = 3;
  paras.lambda_vector[4] = 3;

  calculate_abs_ranking_ch2(paras);

}

