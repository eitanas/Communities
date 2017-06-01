#include <iostream>
#include "networks.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <string>
#include <stdexcept>
#include <queue>
#include <ctime>
#include <algorithm>
#include <random>
#include <sstream>
#include <fstream>
#include <ostream>
#include <omp.h>

#define  P(A) cout << #A << ": " <<(A) << endl;

using namespace std;
typedef std::mersenne_twister_engine< uint32_t, 32, 351,
        175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253 > mt11213b;


//global variables:
float k_total = 4;              //k_inter + k_intra = k_total
int const M = 1000;
int const COMM_SIZE = 10000;
long int N = COMM_SIZE*M;
double a_min = 100000;           // alpha parameter, alpha = k_intra/k_inter
double a_max = 12e6;
int da = 50000;
int num_of__sim = 50;

//functions:
vector<double> build_alpha_vector(double min, double max, double da);
vector<float> calculteKinterKintra(long double a);
template <class T>
void printVector(vector<T>& v);
template  <class T>
T getVectorAvg(vector<T>& v);


int main() {
    double alpha;
    vector<double> alpha_v;
    vector<int> avg__gcc_nodes;
    vector<int> avg_secondGCC_nodes;
    vector<int> avg_gcc_modules;


    /*****************************************/
    cout<<"-------starting--------"<<endl;

    P(k_total); P(M); P(N);
    alpha_v = build_alpha_vector(a_min,a_max, da);

    for (int alpha_idx=0; alpha_idx<alpha_v.size(); alpha_idx++) {
        vector<int> gcc_nodes_storage;
        vector<int> secGcc_nodes_storage;
        vector<int> gcc_modules_storage;
        P(alpha);
        omp_set_num_threads(4);
        #pragma omp parallel for
        for (int s = 0; s < num_of__sim; s++) {
            cout<<"running simualtion "<<s<<" out of "<< num_of__sim<<endl;
            alpha = alpha_v[alpha_idx];

            network myNetwork(ER, alpha);
            myNetwork.build_network();

            vector<int> results;
            results = myNetwork.testConnectivity();

            printVector(results);
            gcc_nodes_storage.push_back(results[0]);
            gcc_modules_storage.push_back(results[1]);
            secGcc_nodes_storage.push_back(results[2]);
        }
        avg__gcc_nodes.push_back(getVectorAvg(gcc_nodes_storage));
        avg_gcc_modules.push_back(getVectorAvg(gcc_modules_storage));
        avg_secondGCC_nodes.push_back((getVectorAvg(secGcc_nodes_storage)));

    }

    cout<<"saving data"<<endl;
    ostringstream oss;
    oss << M <<"_"<<COMM_SIZE<<"_ER_modules_GCC_SGCC_"<<a_min<<"_"<<a_max<<".txt";
    string str = oss.str();
    ofstream data;
    string path =  "/home/eitan/Dropbox/M.Sc/results/" + str;
    //mac:        "/Users/eitanasher/Desktop/" + str;
    //linux:       "/home/eitan/Dropbox/M.Sc/results/" + str;
    //windows:    "G:\\dropbox\\Dropbox\\M.Sc\\results" + str;
    cout<<"path = "<<path<<endl;
    const char* c_path = path.c_str();
    data.open(c_path, ofstream::out );
    data<< "alpha      k inter      nodes in giant       nodes in second       modules"<< endl;

    for (int j=0; j< alpha_v.size(); j++){
        float k = calculteKinterKintra(alpha_v[j])[1];
        data<<  alpha_v[j]<<"      " <<k<<"       "<<avg__gcc_nodes[j]<<"          "
            <<avg_secondGCC_nodes[j]<<"         "<<avg_gcc_modules[j]<<endl;
    }
    data.close();
    cout<<"data saved successfully"<<endl;



    return 0;
}

vector<double>  build_alpha_vector(double min, double max, double da){
    double step = da;
    vector<double> v;
    v.push_back(min);
    double n = (max-min)/da;
    for (int i=1; i<n; i++){
        v.push_back(v[0] + da);
        da+=step;
    }
    cout<<"finished building alpha vector"<<endl;
    return v;
}
vector<float> calculteKinterKintra(long double a){
    //gets alpha as input and calculates k_inter and k_intra, assuming k_total is 4
    vector<float > res;
    float intra = float(k_total * a / (a + M - 1));
    float inter = float(k_total - intra);
    res.push_back(intra);
    res.push_back(inter);
    return res;
}
template <class T>
void printVector(vector<T>& v){
    cout<<"results:";
    int n = v.size();
    for (int i=0; i<n; i++){
        if (i%100==0)
            cout<<endl;
        cout<<v[i]<<" ";

    }cout<<endl;
}
template  <class T>
T getVectorAvg(vector<T>& v){
    T average = accumulate(v.begin(), v.end(), 0)/v.size();
    return  average;
}