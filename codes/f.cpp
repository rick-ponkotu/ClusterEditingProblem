#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <cmath>


using namespace std;
const static int inf = 1000000;
int k = 0;
typedef tuple<int, int, int> Type;


Edge EdgeBranching(Graph G);
Graph DataReduction(Graph G, int N);
vector<Type> ConflictTriple(Graph G, int N);

int main() {
    //ifstream ifs(argv[1]);// 第一引数のファイルを受け取り
    string input_file = "exact001.gr";
    ifstream ifs(input_file);
    string line;
    int N, M;

    if(!ifs){
        cerr << "Error: file not opened." << endl;
        return 1;
    }
    int num = 0;
    while(getline(ifs, line, ' ')){
        if(num == 2) N = stoi(line);
        if(num == 3) M = stoi(line);
        num ++;
    }
    ifs.clear();
    ifs.seekg(0, ios::beg);
    int tmp[M][2];
    string buf;

    for(int i=0; i < M+2; i++){
        if(i == 0) getline(ifs, line);
        if(i > 1){
            getline(ifs, line);
            stringstream ss;
            int a, b;
            ss << line;
            ss >> a;
            ss.ignore();
            ss >> b;
            tmp[i-2][0] = a;
            tmp[i-2][1] = b;
        }
    }

    ifs.close();

    Graph G(N+1);

    for(int a=1; a < N+1; a++) {
        for(int b=1; b < N+1; b++) {
            if(a != b) {
                int from = a;
                int to = b;
                int weight;
                for(int i=0; i < M; i++) {
                    if(from == tmp[i][0] && to == tmp[i][1]) {
                        weight = 1;
                        break;
                    }else if(from == tmp[i][1] && to == tmp[i][0]) {
                        weight = 1;
                        break;
                    }else {
                        weight = -1;
                    }
                }
                G[from].push_back(Edge(to, weight));
                //G[to].push_back(Edge(from, weight));
            }
        }
    }

/*
    for(int i=0; i < M; i++){
        int from = tmp[i][0];
        int to = tmp[i][1];
        int weight = 1;

        G[from].push_back(Edge(to, weight));
        G[to].push_back(Edge(from, weight));
    }
*/

    //時間計測開始
    int start = clock();
    for(int a=1; a < N+1; a++){
        for (int i = 0; i < G[a].size(); i++) {
            Edge e = G[a][i];
            //if(e.weight > 0) cout << "G[" << a << "] = {" << e.to << ", weight: " << e.weight << "}" << endl;
        }
    }
    int end = clock();
    //時間計測終了

    cout << "Execution Time of " << input_file << " : " << 1.0 * (end - start) / CLOCKS_PER_SEC << endl;
}


Graph DataReduction(Graph G, int N) {

    //Rule 2
    int icf[N+1][N+1];
    int icp[N+1][N+1];

    for(int u=1; u < N+1; u++) {
        for(int v=1; v < N+1; v++) {
            if(u != v) {
                vector<int> N_uv;
                vector<int> delta_uv;
                vector<int> union_uv;
                vector<int> N_u;
                vector<int> N_v;
                int sum_icf = 0;
                int sum_icp = 0;

                for(int i=1; i < G[u].size(); i++) {
                    Edge e = G[u][i];
                    if(e.weight > 0) N_u.push_back(e.to);
                }

                for(int i=1; i < G[v].size(); i++) {
                    Edge f = G[v][i];
                    if(f.weight > 0) N_v.push_back(f.to);
                }
                sort(N_u.begin(), N_u.end());
                sort(N_v.begin(), N_v.end());
                set_union(begin(N_u), end(N_u), begin(N_v), end(N_v), inserter(union_uv, end(union_uv)));
                set_intersection(begin(N_u), end(N_u), begin(N_v), end(N_v), inserter(N_uv, end(N_uv)));
                set_difference(begin(N_uv), end(N_uv), begin(union_uv), end(union_uv), inserter(delta_uv, end(delta_uv)));

                for(int i= 0; i < N_uv.size(); i++){
                    int w = N_uv[i];
                    Edge E = G[u][w];
                    Edge F = G[v][w];
                    sum_icf += min(E.weight, F.weight);
                }

                for(int i= 0; i < delta_uv.size(); i++){
                    int w = delta_uv[i];
                    Edge E = G[u][w];
                    Edge F = G[v][w];
                    sum_icf += min(abs(E.weight), abs(F.weight));
                }

                Edge E = G[u][v];
                icf[u][v] = max(0, E.weight) + sum_icf;
                icp[u][v] = max(0, -E.weight) + sum_icp;


            }
        }
    }

    for(int u=1; u < N+1; u++) {
        for(int v=1; v < N+1; v++) {
            if(u != v) {
                if(icf[u][v] > k && icp[u][v] > k) {
                    return G;
                }else if(icf[u][v] > k) {
                    Edge E = G[u][v];
                    E.weight = inf;
                }else if(icp[u][v] > k) {
                    Edge E = G[u][v];
                    E.weight = -inf;
                }
            }
        }
    }

    //Rule 3
}

vector<Type> ConflictTriples(Graph G, int N) {
    vector<int> N_u;
    vector<Type> tri;

    for(int u=1; u < N+1; u++) {
        for(int i=1; i < G[u].size(); i++) {
            Edge e = G[u][i];
            if(e.weight > 0) N_u.push_back(e.to);
        }
    }

    for(int u=1; u < N+1; u++) {
        if(G[u].size() > 1) {
            for(int i=0; i < N_u.size(); i++) {
                for(int j=0; j < N_u.size(); j++) {
                    int v = N_u[i];
                    int w = N_u[j];
                    Edge e = G[v][w];
                    if(e.weight < 0) {
                        tuple<int, int, int> uvw = make_tuple(u, v, w);
                        tri.push_back(uvw);
                    }
                }
            }
        }
    }

    return tri;
}
