#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "main.h"
#include "graph.h"

using namespace std;

Graph::Graph(){
  this-> num_nodes = 0;
}

Graph::Graph(int n){
  this->num_nodes = n;
  this->weight = vector <vector <int>>(n, vector <int>(n,-1));
  this->flag = vector <vector <int>>(n, vector <int>(n,0));

  this->node_names = vector <int>(n);
  this->node_pointers = vector <int>(n);
  this->neighbors = vector <vector <int>>(n);
  this->name_to_ind = vector <int>(n);

  REP(i, n){
    this->node_names[i] = i;
    this->node_pointers[i] = i;
    this->name_to_ind[i] = i;
  }
}

Graph::Graph(const char *filename){

  ifstream ifs(filename);
  string line;
  int N;

  if(!ifs){
    cerr << "Error: file not opened." << endl;
    return;
  }

  int num = 0;

  while(getline(ifs, line)){
    if(line[0] == 'c') continue;

    if(line[0] == 'p'){
      stringstream s{line};
      int cnt = 0;
      while(getline(s, line, ' ')){
	      if(cnt == 2){
	        N = stoi(line);
	        this->num_nodes = N;
	        this->weight = vector <vector <int>>(N, vector <int>(N,-1));
	        this->flag = vector <vector <int>>(N, vector <int>(N,0));
          this->node_names = vector <int>(N);
          this->node_pointers = vector <int>(N);
          this->neighbors = vector <vector <int>>(N);
          this->name_to_ind = vector <int> (N);

          REP(i, N){
            this->node_names[i] = i;
            this->node_pointers[i] = i;
            this->name_to_ind[i] = i;
          }
	      }
	      cnt++;
      }
    }
    else{
      stringstream ss;
      int a, b;
      ss << line;
      ss >> a;
      ss.ignore();
      ss >> b;
      this -> weight[a-1][b-1] = 1;
      this -> weight[b-1][a-1] = 1;
      this -> neighbors[a-1].push_back(b-1);
      this -> neighbors[b-1].push_back(a-1);
    }

    num++;
  }

  ifs.close();
}

void Graph::delete_node(int i, vector <edge>& sol, const Graph &G_orig){
  for(int &u : this->node_pointers) if(u == this->node_names[i]) u = -1;
  REP(j, this->num_nodes){
    if(j == i) continue;
    if(this->Flag(i,j) == 0) this->forbid(i,j,sol, G_orig);
  }
  //this->weight.erase(this->weight.begin()+i);
  //this->flag.erase(this->flag.begin()+i);
  this->name_to_ind[this->node_names[i]] = -1;
  this->num_nodes--;
  this->node_names.erase(this->node_names.begin()+i);
  FOR(j, i, this->num_nodes) this->name_to_ind[this->node_names[j]]--;

  /*
  REP(j, this->num_nodes){
    this->weight[j].erase(this->weight[j].begin()+i);
    this->flag[j].erase(this->flag[j].begin()+i);
  }
  */

  for(int u: this->neighbors[this->node_names[i]]){
    auto it = remove(this->neighbors[u].begin(), this->neighbors[u].end(), this->node_names[i]);
    this->neighbors[u].erase(it, this->neighbors[u].end());
  }

  return;
}

void Graph::delete_nodes(vector <int> U, vector <edge>& sol, const Graph &G_orig){
  sort(U.begin(), U.end());

  REP(i, (int) U.size()){
    int j = i+1;
    REP(k, this->num_nodes){
      if(k == U[i]) continue;

      if(j < (int) U.size() && U[j]==k){
        if(this->Flag(U[i],k) == 0) this->permanent(U[i], k, sol, G_orig);
        j++;
      }
      else{
        if(this->Flag(U[i],k) == 0) this->forbid(U[i], k, sol, G_orig);
      }
    }
  }

  for(int i = U.size()-1; i >= 0; i--){
    //this->weight.erase(this->weight.begin()+U[i]);
    //this->flag.erase(this->flag.begin()+U[i]);
    for(int a : this->neighbors[this->node_names[U[i]]]){
      auto it = remove(this->neighbors[a].begin(), this->neighbors[a].end(), this->node_names[U[i]]);
      this->neighbors[a].erase(it, this->neighbors[a].end());
    }
    this->name_to_ind[this->node_names[U[i]]] = -1;
    this->num_nodes--;
    this->node_names.erase(this->node_names.begin()+U[i]);
    FOR(j, U[i], this->num_nodes) this->name_to_ind[this->node_names[j]]--;
    /*
    REP(j, this->num_nodes){
      this->weight[j].erase(this->weight[j].begin()+U[i]);
      this->flag[j].erase(this->flag[j].begin()+U[i]);
    }
    */
  }

  return;
}


int Graph::merge_nodes(int a, int b, vector <edge> &sol, MergeData& mg, const Graph &G_orig){
  int ans = 0;
  if(this->node_names[b] < this->node_names[a]) SWAP(int, a,b);

  REP(k, this->num_nodes){
    mg.flag[this->node_names[k]] = this->Flag(a,k);
    mg.weight[this->node_names[k]] = this->Weight(a,k);
    if(k != b && this->Weight(b,k) > 0){
      auto it = remove(this->neighbors[this->node_names[k]].begin(),this->neighbors[this->node_names[k]].end(), this->node_names[b]);
      this->neighbors[this->node_names[k]].erase(it, this->neighbors[this->node_names[k]].end());
    }

    if(k == a || k == b) continue;

    if(this->Flag(a,k) * this->Flag(b,k) == -1) {
      cerr << "merge error " << a << ", " << b << ", " << k << endl;
      return -1;
    }
    else if(this->Flag(a,k) == 0 && this->Flag(b,k) == 0){
      if(this->Weight(k,a)*this->Weight(k,b) < 0){
        ans += min(abs(this->Weight(k,a)), abs(this->Weight(k,b)));
      }
      if(this->Weight(k,a) <= 0 && this->Weight(k,a) + this->Weight(k,b) > 0){
        this->neighbors[this->node_names[a]].push_back(this->node_names[k]);
        this->neighbors[this->node_names[k]].push_back(this->node_names[a]);
      }
      if(this->Weight(k,a) > 0 && this->Weight(k,a) + this->Weight(k,b) <= 0){
        auto it = remove(this->neighbors[this->node_names[a]].begin(), this->neighbors[this->node_names[a]].end(), this->node_names[k]);
        this->neighbors[this->node_names[a]].erase(it, this->neighbors[this->node_names[a]].end());
        it = remove(this->neighbors[this->node_names[k]].begin(), this->neighbors[this->node_names[k]].end(), this->node_names[a]);
        this->neighbors[this->node_names[k]].erase(it, this->neighbors[this->node_names[k]].end());
      }
      this->weight[this->node_names[k]][this->node_names[a]] += this->weight[this->node_names[k]][this->node_names[b]];
      this->weight[this->node_names[a]][this->node_names[k]] = this->weight[this->node_names[k]][this->node_names[a]];
    }
    else if(this->Flag(b,k) == -1 && this->Flag(a,k) == 0){
      this->forbid(a,k,sol,G_orig);
      if(this->Weight(a,k) > 0){
        ans += this->Weight(a,k);
        this->flip_edge(a,k);
      }
      this->weight[this->node_names[a]][this->node_names[k]] += this->weight[this->node_names[b]][this->node_names[k]];
      this->weight[this->node_names[k]][this->node_names[a]] = this->weight[this->node_names[a]][this->node_names[k]];
    }
    else if(this->Flag(a,k) == -1 && this->Flag(b,k) == 0){
      this->forbid(b,k,sol,G_orig);
      this->flag[this->node_names[k]][this->node_names[b]] = 0;
      this->flag[this->node_names[b]][this->node_names[k]] = 0;
      if(this->Weight(b,k) > 0){
        ans += this->Weight(b,k);
        this->weight[this->node_names[a]][this->node_names[k]] -= this->weight[this->node_names[b]][this->node_names[k]];
        this->weight[this->node_names[k]][this->node_names[a]] = this->weight[this->node_names[a]][this->node_names[k]];
      }
      else{
        this->weight[this->node_names[a]][this->node_names[k]] += this->weight[this->node_names[b]][this->node_names[k]];
        this->weight[this->node_names[k]][this->node_names[a]] = this->weight[this->node_names[a]][this->node_names[k]];
      }
    }
    else if(this->Flag(a,k) == 1 && this->Flag(b,k) == 0){
      this->permanent(b,k,sol,G_orig);
      this->flag[this->node_names[k]][this->node_names[b]] = 0;
      this->flag[this->node_names[b]][this->node_names[k]] = 0;

      if(this->Weight(b,k) < 0){
        ans -= this->Weight(b,k);
        this->weight[this->node_names[a]][this->node_names[k]] -= this->weight[this->node_names[b]][this->node_names[k]];
        this->weight[this->node_names[k]][this->node_names[a]] = this->weight[this->node_names[a]][this->node_names[k]];
      }
      else{
        this->weight[this->node_names[a]][this->node_names[k]] += this->weight[this->node_names[b]][this->node_names[k]];
        this->weight[this->node_names[k]][this->node_names[a]] = this->weight[this->node_names[a]][this->node_names[k]];
      }
    }
    else if(this->Flag(b,k) == 1 && this->Flag(a,k) == 0){
      this->permanent(a,k,sol,G_orig);
      if(this->Weight(a,k) <= 0){
        ans -= this->Weight(a,k);
        this->flip_edge(a,k);
      }
      this->weight[this->node_names[a]][this->node_names[k]] += this->weight[this->node_names[b]][this->node_names[k]];
      this->weight[this->node_names[k]][this->node_names[a]] = this->weight[this->node_names[a]][this->node_names[k]];
    }
    else if(this->Flag(b,k) * this->Flag(a,k) == 1){
        this->weight[this->node_names[a]][this->node_names[k]] += this->weight[this->node_names[b]][this->node_names[k]];
        this->weight[this->node_names[k]][this->node_names[a]] = this->weight[this->node_names[a]][this->node_names[k]];
    }
  }

  REP(u, G_orig.num_nodes) if(this->node_pointers[u] == this->node_names[b]){
    this->node_pointers[u] = this->node_names[a];
    mg.pointers.push_back(u);
  }
  /*
  this->weight.erase(this->weight.begin()+b);
  this->flag.erase(this->flag.begin()+b);
  */
  this->name_to_ind[this->node_names[b]] = -1;
  this->node_names.erase(this->node_names.begin()+b);
  if(b < this->num_nodes -1) FOR(j, b, this->num_nodes-1) this->name_to_ind[this->node_names[j]]--;
  this->num_nodes--;
  /*
  REP(j, this->num_nodes){
    this->weight[j].erase(this->weight[j].begin()+b);
    this->flag[j].erase(this->flag[j].begin()+b);
  }
  */
  return ans;
}

void Graph::expand_nodes(int a, int b, const MergeData& mg){
  if(a > b) SWAP(int, a, b);

  REP(k, this->num_nodes){
    if(this->weight[this->node_names[k]][b] > 0) this->neighbors[this->node_names[k]].push_back(b);
    if(this->node_names[k] == a) continue;

    if(this->weight[a][this->node_names[k]] <= 0 && mg.weight[this->node_names[k]] > 0){
      this->neighbors[a].push_back(this->node_names[k]);
      this->neighbors[this->node_names[k]].push_back(a);
    }

    if(this->weight[a][this->node_names[k]] > 0 && mg.weight[this->node_names[k]] <= 0){
      auto it = remove(this->neighbors[a].begin(),this->neighbors[a].end(), this->node_names[k]);
      this->neighbors[a].erase(it, this->neighbors[a].end());
      it = remove(this->neighbors[this->node_names[k]].begin(), this->neighbors[this->node_names[k]].end(), a);
      this->neighbors[this->node_names[k]].erase(it, this->neighbors[this->node_names[k]].end());
    }

    this->flag[a][this->node_names[k]] = mg.flag[this->node_names[k]];
    this->flag[this->node_names[k]][a] = mg.flag[this->node_names[k]];
    this->weight[a][this->node_names[k]] = mg.weight[this->node_names[k]];
    this->weight[this->node_names[k]][a] = mg.weight[this->node_names[k]];
  }

  this->name_to_ind[b] = num_nodes;
  this->num_nodes++;
  this->node_names.push_back(b);

  for(int i: mg.pointers) this->node_pointers[i] = b;
  return;
}


void Graph::forbid(int a, int b, vector <edge>& sol, const Graph& G_orig){
  if(this->Flag(a,b) == 1) cerr << "err: this pair cannnot be forbidden" << endl;
  if(this->Flag(a,b) == -1) cerr << "err: this pair has been already forbidden" << endl;
  this->flag[this->node_names[a]][this->node_names[b]] = -1;
  this->flag[this->node_names[b]][this->node_names[a]] = -1;

  vector <int> A, B;
  REP(i, (int) this->node_pointers.size()){
    if(this->node_pointers[i] == this->node_names[a]) A.push_back(i);
    if(this->node_pointers[i] == this->node_names[b]) B.push_back(i);
  }

  for(auto i: A) for(auto j: B){
    if (G_orig.Weight(i,j) > 0){
      sol.push_back(make_pair(i,j));
    }
  }
  return;
}

void Graph::permanent(int a, int b, vector <edge> & sol, const Graph &G_orig){
  if(this->Flag(a,b) == 1) {
    cerr << "err: this pair has been already permanent" << endl;
    cout << a << " " << b << endl;
  }
  if(this->Flag(a,b) == -1) cerr << "err: this pair cannot be permanent" << endl;
  this->flag[this->node_names[a]][this->node_names[b]] = 1;
  this->flag[this->node_names[b]][this->node_names[a]] = 1;

  vector <int> A, B;
  REP(i, (int) this->node_pointers.size()){
    if(this->node_pointers[i] == this->node_names[a]) A.push_back(i);
    if(this->node_pointers[i] == this->node_names[b]) B.push_back(i);
  }

  for(auto i: A) for(auto j: B){
    if (G_orig.Weight(i,j) < 0) sol.push_back(make_pair(i,j));
  }
  return;
}


void Graph::delete_edge(int a, int b){
  this->weight[this->node_names[a]][this->node_names[b]] = -1;
  this->weight[this->node_names[b]][this->node_names[a]] = -1;
  auto it = remove(this->neighbors[this->node_names[a]].begin(), this->neighbors[this->node_names[a]].end(), this->node_names[b]);
  this->neighbors[this->node_names[a]].erase(it, this->neighbors[this->node_names[a]].end());
  it = remove(this->neighbors[this->node_names[b]].begin(), this->neighbors[this->node_names[b]].end(), this->node_names[a]);
  this->neighbors[this->node_names[b]].erase(it, this->neighbors[this->node_names[b]].end());
  return;
}


void Graph::add_edge(int a, int b){
  this->weight[this->node_names[a]][this->node_names[b]] = 1;
  this->weight[this->node_names[b]][this->node_names[a]] = 1;
  this->neighbors[this->node_names[a]].push_back(this->node_names[b]);
  this->neighbors[this->node_names[b]].push_back(this->node_names[a]);
  return;
}

void Graph::flip_edge(int a, int b){
  if(this->Weight(a,b) > 0){
    auto it = remove(this->neighbors[this->node_names[a]].begin(), this->neighbors[this->node_names[a]].end(), this->node_names[b]);
    this->neighbors[this->node_names[a]].erase(it, this->neighbors[this->node_names[a]].end());
    it = remove(this->neighbors[this->node_names[b]].begin(), this->neighbors[this->node_names[b]].end(), this->node_names[a]);
    this->neighbors[this->node_names[b]].erase(it, this->neighbors[this->node_names[b]].end());
  }
  else{
    this->neighbors[this->node_names[a]].push_back(this->node_names[b]);
    this->neighbors[this->node_names[b]].push_back(this->node_names[a]);
  }
  this->weight[this->node_names[a]][this->node_names[b]] *= -1;
  this->weight[this->node_names[b]][this->node_names[a]] *= -1;
  return;
}

int Graph::Weight(int a, int b) const{
  return this->weight[this->node_names[a]][this->node_names[b]];
}

int Graph::Flag(int a, int b) const{
  return this->flag[this->node_names[a]][this->node_names[b]];
}

void Graph::reset_flag(int a, int b) {
  this->flag[this->node_names[a]][this->node_names[b]] = 0;
  this->flag[this->node_names[b]][this->node_names[a]] = 0;
}

/*
bool Graph::conflict_triple (vector <int>& triple) const {
  FOR(u, 0, this->num_nodes) FOR(v, 0, this->num_nodes-1){
    if (u == v || this->Weight(u,v) <= 0 || this->Flag(u,v) == -1) continue;

    FOR(w, v+1, this->num_nodes){
      if(u == w) continue;

      if(this->Weight(u,w) > 0 && this->Flag(u,w) >= 0 && (this->Weight(v,w) <= 0 || this->Flag(v,w) == -1)){
	      triple.clear();
	      triple.push_back(u);
	      triple.push_back(v);
	      triple.push_back(w);
	      return true;
      }
    }
  }

  return false;
}
*/

bool Graph::conflict_triple (vector <int>& triple) const {
  for(int u: this->node_names){
    if(this->neighbors[u].size() < 2) continue;
    for(auto v = this->neighbors[u].begin(); next(v) != this->neighbors[u].end(); v++) {
      for(auto w=next(v); w != this->neighbors[u].end(); w++){
        if(this->weight[(*v)][(*w)] <= 0 || this->flag[(*v)][(*w)] == -1){
          triple.clear();
          triple.push_back(this->name_to_ind[u]);
	        triple.push_back(this->name_to_ind[*v]);
	        triple.push_back(this->name_to_ind[*w]);
          /*
          triple[0] = u;
          triple[1] = *v;
          triple[2] = *w;
          */
  	      return true;
        }
      }
    }
  }

  return false;
}

void Graph::show () const{
  int n = this->num_nodes;

  FOR(i, 0, n){
    cout << i << ": ";
    FOR(j, 0, n){
      cout << j << "," << this->Weight(i,j) << "; ";
    }
    cout << endl;
  }

  return;
}

MergeData::MergeData(int n){
  this->flag = vector <int>(n);
  this->weight = vector <int>(n);
  this->pointers.clear();
}