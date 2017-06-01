// Created by eitan and la on 27/05/2017.
#ifndef NETWORKSFIX_NETWORKS_H
#define NETWORKSFIX_NETWORKS_H
#define  P(A) cout << #A << ": " <<(A) << endl;

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
using namespace std;

typedef std::mersenne_twister_engine< uint32_t, 32, 351,
        175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253 > mt11213b;

enum networkType {ER, LATTICE};
extern const int M, COMM_SIZE;
extern long int N;
extern float k_total;
mt11213b rnd;

class network{
public:
    //variables:
    networkType type_;
    vector<vector<int> > adj_list;
    vector <long int> distances;
    vector <long int> components;                   //for rvery node in the network, holds its module number
    vector <int> connected;
    double k_inter, k_intra;

    //functions:m
    network(networkType type, double a);
    void build_network();
    void assign_inter_links();
    void assignKinterKintra(double a);
    void addEdge(int u, int v);
    void removeEdge(int u, int v);
    template <class T>
    T vector_sum(vector<T>& v);
    template <class T>
    T print_vector(vector<T>& v);
    bool isConnected(int u, int v);
    int getModule(int node_id);
    vector<int> bfs(int s);
    long double MPL(vector<vector<int> >& list);                       //calculate the networks' mean path length
    void setNetworkType(networkType t);
    networkType getNetworkType();
    void assign_inter_links(vector<vector<int> >& adj_list);
    int getDegree(vector<vector<int> >& adjlist,int u);
    vector<int> testConnectivity();
    vector<int> getSecondGiant(int max_clstr);
    vector<int> createRandomVector(int max);
    int getLongestPath();
    int getDegree(int u);
    void deleteAll();
    ~network();
};

network::network(networkType t, double alpha) {
    setNetworkType(t);
    assignKinterKintra(alpha);
    switch (t) {
        case (ER):
            adj_list.resize(N);
            distances.resize(N);
            components.resize(N);
            connected.resize(N);
            break;
        case (LATTICE):
            components.resize(N);
            connected.resize(N);
            break;
    }
}
void network::setNetworkType(networkType t){
    type_ = t;
    //cout<<"network typre is "<<network.type<<endl;
}
networkType network::getNetworkType(){
    return type_;
}
void network::build_network() {
    rnd.seed(time(0));
    for (int m = 0; m < M; m++) {
        int num_of_edges;
        num_of_edges = int(round(COMM_SIZE * k_intra / 2));
        for (int i = 0; i < num_of_edges; i++) {
            int r1 = (m * COMM_SIZE) + (rnd() % COMM_SIZE);
            int r2 = (m * COMM_SIZE) + (rnd() % COMM_SIZE);
            //cout<<"r1 = "<<r1<<" r2 = "<<r2<<endl;
            if ((r1 != r2) && !isConnected(r1, r2)) {
                adj_list[r1].push_back(r2);
                adj_list[r2].push_back(r1);
                //cout<<"edge assigned, i = "<<i<<endl;
                i++;
            }
            i--;
        }
    }
    assign_inter_links();
}
void network::assign_inter_links(){

    networkType t = getNetworkType();
    switch (t) {
        case (ER):
            cout<<"number of inter-links = "<<round(k_inter * N / 2)<<endl;
            //double inter = round(k_inter * N / 2);
            for (int l = 0; l < round(k_inter * N / 2); l++) {
                int r1 = rnd() % N;
                int r2 = rnd() % N;
                //cout<<"r1 , r2 = "<<r1<<" & "<<r2<<endl;
                if (getModule(r1) != getModule(r2) && !isConnected(r1, r2)) {
                    adj_list[r1].push_back(r2);
                    adj_list[r2].push_back(r1);
                    //here need to add network of networks list
                }
            }
            cout<<"finished assigning inter-links"<<endl;
            break;

        case (LATTICE):
            cout<<"i built LATTICE and its not good"<<endl;
            break;
    }
}
int network::getModule(int u){
    return u/COMM_SIZE;
}
vector<int> network::bfs(int s){
    vector<int> results;
    bool visited[N];
    int longest_distance, num_p_degree, num_n_degree,degree=1;
    long int sum=0;
    long int counter=0;
    //visited vector initialization:
    for (int k=0; k<N; k++) visited[k]=0;
    visited[s]=1;

    queue <int> bfsq;
    bfsq.push(s);
    num_p_degree = 1;
    num_n_degree = 0;

    while (!bfsq.empty()){
        int u = bfsq.front();
        bfsq.pop();
        ////ask bnaya whats n and p
        if (num_p_degree==0){
            degree++;
            num_p_degree = num_n_degree;
            num_n_degree = 0;
        }
        //looking at all neighbors of u
        for (int i=0; i < adj_list.size(); i++){
            int v = adj_list[u][v];

            if (!visited[v]){
                bfsq.push(v);
                num_n_degree++;
                sum+=degree;
                counter+=1;

                distances[degree]++;
            }
        }
        num_p_degree--;
    }
    results.push_back(sum);
    results.push_back(counter);
    results.push_back(degree);
    return results;
}
vector<int> network::testConnectivity(){
    //this function uses BFS and returns 3 components: 1.number of modules in gcc, 2.number of nodes in gcc
// 3.second giant component size
    vector<int> results;
    queue <int> Q;
    int size, expCount = 0, cluster_id = 0, giant_size = 0, giant_id;
    //int explored[N];         // 0-white, 1-grey, 2-black
    vector<int> explored(N,0);
    //int *explored = new int[N];// 0-white, 1-grey, 2-black
    vector<int> community_gcc(M,0);             //this vector holds the communities that have nodes that belong to the gcc

    for (long int i=0; i<N; i++){
        //        if (active[i]==0)
//        {
//            explored[i] = 2;
//            expCount++;
//            components[i]=0;
//        }
        //this part of code takes care of all the unconnected nodes, i.e nodes with degree = 0
        if (getDegree(i) == 0){
            explored[i] = 2;
            expCount++;
            components[i] = cluster_id;
            cluster_id++;
        }
        else
            explored[i] = 0;
    }
    int node_id = 0;
    while (expCount < N ){
        while (explored[node_id]!=0){     //find a source for BFS algorithm
            node_id++;
        }
        cluster_id++;
        Q.push(node_id);
        components[node_id] = cluster_id;
        size = 1;
        explored[node_id] = 1;
        expCount++;
        while (!Q.empty()){
            int u = Q.front();
            Q.pop();
            explored[u] = 2;
            for(int i=0; i<adj_list[u].size(); i++){
                if ( explored[adj_list[u][i]] == 0){
                    //cout<<"explored[list[u][i] = "<< explored[list[u][i]]<<endl;
                    size++;
                    Q.push(adj_list[u][i]);
                    components[adj_list[u][i]] = cluster_id;
                    explored[adj_list[u][i]] = 1;
                    expCount++;
                }
            }
        }
        if (size>giant_size){
            giant_size = size;
            giant_id = cluster_id;
        }
    }
    for (int k=0; k<N; k++) {
        //cout<<"components[k] = "<<components[k]<<endl;
        if (components[k] == giant_id) {
            connected[k] = 1;
            community_gcc[getModule(k)] = 1;
            //cout<<"getModule(k,M) = "<<getModule(k,C)<<endl;
        } else {
            connected[k] = 0;

        }
        //if (k%1000==0){cout<<"k = "<<k<<endl;}
    }

    vector<int> second_gcc_results = getSecondGiant(cluster_id);
    int second_id  = second_gcc_results[0];
    int second_count = second_gcc_results[1];

    results.push_back(giant_size);
    results.push_back(vector_sum(community_gcc));
    results.push_back(second_count);

    second_gcc_results.clear();
    explored.clear();
    community_gcc.clear();
    return  results;
}
bool network::isConnected(int u, int v) {
    //l stands for list
    for (int k = 0; k < adj_list[u].size(); k++) {
        if (adj_list[u][k] == v) return true;
    }
    return false;
}
void network::assignKinterKintra(double a){
    //gets alpha as input and sets k_inter and k_intra, assuming k_total is 4
    k_intra = k_total * a / (a + M - 1);
    k_inter = k_total - k_intra;
}
vector<int> network::getSecondGiant( int max_cluster) {
    //cout<<"calculating second giant"<<endl;
    vector<int> results;
    vector<int> count(max_cluster+1, 0);
    for (int i=0; i < N; i++){
        count[components[i]]++;
    }
    // Traverse through the count[] and find second highest element.
    int first=0, second_id=0;
    for (int i=0; i<max_cluster; i++){
        /* If current element is higher than first then update both
          first and second */
        if (count[i] > count[first]){
            second_id = first;
            first = i;
        }
            /* If count[i] is in between first and second then update second  */
        else if (count[i] > count[second_id] &&
                 count[i] != count[first])
            second_id = i;
    }
    results.push_back(second_id);    //second is the second gcc id
    results.push_back(count[second_id]);

    count.clear();
    //cout<<"finished calculating second"<<endl;
    return results;  //this is the number of nodes in second gcc

}
int network::getDegree(int u) {
    return adj_list[u].size();
}
template <class T>
T network::vector_sum(vector<T>& v){
    T sum = 0;
    int n = v.size();
    for (int i=0; i<n; i++){
        sum += v[i];
    }
    return sum;
}


network::~network() {
    adj_list.clear();
    distances.clear();
    components.clear();
    connected.clear();
}


int network::getLongestPath(){
    int longest = 0, dist;
    for (int i = 0; i < N ; ++i) {
        if (i%5000==0){cout<<"traversing node "<<i<<endl;}
        dist = network::bfs(i)[2];
        if (dist>longest){
            longest = dist;
            cout<<"longest path found is = "<<longest<<endl;
        }
    }
    return longest;
}













#endif //NETWORKSFIX_NETWORKS_H
