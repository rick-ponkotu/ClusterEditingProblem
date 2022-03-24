#pragma once

#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

/*namespace {
	struct custom_hash { // could make it a singleton but whatever
		static uint64_t splitmix64(uint64_t x) {
			x += 0x9e3779b97f4a7c15;
			x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
			x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
			return x ^ (x >> 31);
		}

		std::size_t oprator()(uint64_t x) const {
			static const uint64_t fixedRand = std::chrono::steady_clock::now().time_since_epoch().count();
			return splitmix64(x + fixedRand);
		}
	} hasher;

	inline void hashCombine(std::size_t& seed, uint32_t v)
	{
		seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	}
}*/

namespace akt {
	class Graph
	{
	public:
		class Edge
		{
		public:
			// Ensures that there is a unique tuple corresponsding to each edge of the graph
			Edge(int u, int v) : u(std::min(u, v)), v(std::max(u, v)) { }
			//friend bool operator<=>(const Edge& e1, const Edge& e2) = default;
			friend bool operator==(const Edge& e1, const Edge& e2) { return std::tie(e1.u, e1.v) == std::tie(e2.u, e2.v); }
			friend bool operator<(const Edge& e1, const Edge& e2) { return std::tie(e1.u, e1.v) < std::tie(e2.u, e2.v); }
			//friend std::size_t hashEdge(const Edge& e1) { return hashCombine(hasher(e1.u), hasher(e1.v)); }

			int u, v;
		};

		Graph() = default;

		// Reads a graph in the DIMACS format
		Graph(std::istream& is)
		{
			std::string line;
			while (std::getline(is, line) && (line.front() == 'c'))
				;
			if (line.empty() || (line.front() != 'p'))
				throw std::invalid_argument("First non-comment line is not the problem line");
			auto iss = std::istringstream(line);
			std::string s;
			int n, m;
			if ((!(iss >> line >> line >> n >> m)) || (n < 0) || (m < 0))
				throw std::invalid_argument("Incomplete or invalid problem line");
			adj = std::vector<std::unordered_set<int>>(n);
			for (int i = 0; i < m; ++i) {
				while (std::getline(is, line) && (line.front() == 'c'))
					;
				iss.clear();
				iss.str(line);
				int u, v;
				if (!(iss >> u >> v))
					throw std::invalid_argument("An invalid line or not enough lines in the input. Offending line:\n" + line + "\n");
				--u; --v;
				connect(u, v);
			}
		}

		using EdgeSet = std::set<Edge>;

		EdgeSet removeAllEdgesHeuristic() const
		{
			EdgeSet res;
			for (int i = 0; i < adj.size(); ++i) {
				for (auto j : adj[i])
					res.emplace(i, j);
			}
			return res;
		}
		EdgeSet makeCliqueHeuristic() const
		{
			EdgeSet res;
			for (int i = 0; i < adj.size(); ++i) {
				for (int j = i + 1; j < adj.size(); ++j) {
					if (!connected(i, j))
						res.emplace(i, j);
				}
			}
			return res;
		}
		EdgeSet removeMinHeuristic() const { return neighborhoodHeuristic(0, adj.size(), 1); }
		EdgeSet removeMaxHeuristic() const { return neighborhoodHeuristic(int(adj.size()) - 1, -1, -1); }

		bool connected(int u, int v) const { if (!validIndices(u, v)) return false; return adj[u].count(v) > 0; }

		void connect(int u, int v) { if (!validIndices(u, v) || (u == v)) return; adj[u].insert(v); adj[v].insert(u); }
		void disconnect(int u, int v) { if (!validIndices(u, v)) return; adj[u].erase(v); adj[v].erase(u); }
		void symmetricDifference(const EdgeSet& es)
		{
			for (auto e : es) {
				if (connected(e.u, e.v))
					disconnect(e.u, e.v);
				else
					connect(e.u, e.v);
			}
		}

		void printGraph() const
		{
			for (int i = 0; i < adj.size(); ++i) {
				std::cout << i << ": ";
				for (auto j : adj[i])
					std::cout << j << ", ";
				std::cout << '\n';
			}
		}
		
		std::tuple<int,int,int> findP3() const 
		{ 
			auto [ccs, sizes] = connectedComponents();
			for (int i = 0; i < adj.size(); ++i) {
				if (adj[i].size() != (sizes[ccs[i]] - 1)) {
					for (auto v : adj[i]) {
						for (auto w : adj[v]) {
							if( w != i && adj[w].find(i) == adj[w].end())
								return std::make_tuple(i,v,w);
						}
					}

				}
			}
			return std::make_tuple(-1,-1,-1);
		}

		bool isCluster() const 
		{ 
			auto p3 = findP3();
			return std::get<0>(p3) == -1 && std::get<1>(p3) == -1 && std::get<2>(p3) == -1;
		}

	private:
		bool validIdx(int u) const { return (u >= 0) && (u < adj.size()); }
		bool validIndices() const { return true; }
		template<typename... Args>
		bool validIndices(int u, Args... rest) const { return validIdx(u) && validIndices(rest...); }

		EdgeSet neighborhoodHeuristic(int start, int stop, int step) const
		{
			EdgeSet res;
			std::unordered_set<int> visited;
			for (int i = start; validIdx(i) && (i != stop); i += step) {
				if (visited.count(i) > 0)
					continue;
				visited.insert(i);
				const auto& neigh = adj[i];
				// Make the neighbourhood of i a clique
				for (auto j : neigh) {
					if (visited.count(j))
						continue;
					visited.insert(j);
					// Add all edges needed to make the neighborhood of i a clique
					for (auto k : neigh) {
						if ((visited.count(k) > 0) || connected(j, k))
							continue;
						res.emplace(j, k);
						std::cout << j << " " << k << '\n';
					}
					// Remove all edges going outside of the neighborhood of i
					for (auto k : adj[j]) {
						if ((visited.count(k) > 0) || (neigh.count(k) > 0))
							continue;
						res.emplace(j, k);
						std::cout << j << " " << k << '\n';
					}
				}
			}

			return res;
		}


		void bfs(const int u, std::vector<int>& who, std::vector<int>& sizes) const 
		{
			std::queue<int> q;
			q.push(u);
			while (!q.empty()) {
				auto cur = q.front();
				q.pop();
				for (auto i : adj[cur]) {
					if (who[i] >= 0)
						continue;
					who[i] = who[u];
					++sizes[who[u]];
					q.push(i);
				}
			}
		}

		std::pair<std::vector<int>, std::vector<int>> connectedComponents() const 
		{
			auto who = std::vector<int>(adj.size(), -1);
			std::vector<int> sizes;
			for (int i = 0; i < adj.size(); ++i) {
				if (who[i] >= 0)
					continue;
				who[i] = sizes.size();
				sizes.push_back(1);
				bfs(i, who, sizes);
			}
			return { who, sizes };
		}

		std::vector<std::unordered_set<int>> adj;
	};
}
