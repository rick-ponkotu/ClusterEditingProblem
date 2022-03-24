#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include "Graph.h"

namespace {
	using akt::Graph;
}

int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cerr << "usage: verifer [graph file] [solution file]" << std::endl;
		return -1;
	}
	Graph orig;
	Graph::EdgeSet sol;
	if (std::ifstream ifsol{ argv[2] }; ifsol) {
		std::string line;
		while (std::getline(ifsol, line)) {
			auto iss = std::istringstream(line);
			int u,v;
			if (!(iss >> u >> v)) {
				std::cerr << "An invalid line or not enough lines in the input. Offending line:\n" + line + "\n" << std::endl;
			}
			--u; --v;
			sol.insert(Graph::Edge(u,v));
		}
		if (std::ifstream ifs{ argv[1] }; ifs) {
			try {
				orig = Graph(ifs);
				orig.symmetricDifference(sol);
				auto p3 = orig.findP3();
				if(std::get<0>(p3) == -1) std::cout << "OK" << std::endl;
				else std::cout << "P3:" << std::get<0>(p3) + 1 << " " << std::get<1>(p3) + 1 << " " << std::get<2>(p3) + 1 << std::endl;


			}
			catch (const std::invalid_argument& ia) {
				return -1;
			}
		}
		else { return -1; }
	}
	return 0;
}
