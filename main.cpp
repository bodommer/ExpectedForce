#include "stdafx.h"
using std::cin;
using std::cout;
using std::endl;

typedef std::vector<int> svi;
typedef std::vector<int>::iterator svii;
typedef TUNGraph::TNodeI nodeI;

/* HELPER FUNCTION
Converts a SNAP undirected graph to a sorted edgelist.
 @param[out] egos, alters: The edgelist.
 @param[in] G: The SNAP graph.
*/
int snap_to_edgelist(svi &egos, svi &alters,
	const PUNGraph G) {
	egos.clear(); alters.clear();
	int node_id, node_deg = 0;
	for (nodeI NI = G->BegNI(); NI < G->EndNI(); NI++) { //iterates over nodes
		node_id = NI.GetId();
		node_deg = NI.GetDeg();
		for (int out_cnt = 0; out_cnt < node_deg; out_cnt++) { //iterates over outgoing edges
			//augment the edgelist
			egos.push_back(node_id);
			alters.push_back(NI.GetOutNId(out_cnt));
			//std::cout << node_id << " - " << NI.GetOutNId(out_cnt) << "\n";
		}
	}
	return 0;
}



int main(int argc, char* argv[]) {

	cout << "This program determines the Expected Force of every node in a graph.\n Loads from the binary file 'test.graph'.\n Stores the results in 'results.txt'" << endl;
	system("pause");

	//load a SNAP graph from a SNAP binary file
	TFIn FIn("test.graph"); 
	PUNGraph G = TUNGraph::Load(FIn);
	
	/*
	// alternatively, generate a network using Forest Fire model
	PNGraph H = TSnap::GenForestFire(100, 0.35, 0.35);
	// convert to undirected graph
	PUNGraph G = TSnap::ConvertGraph<PUNGraph>(H);
	{ TFOut FOut("test.graph"); G->Save(FOut); } //save graph in binary form
	*/

	int nNodes;
	nNodes = G->GetNodes(); //number of nodes

	svi egosVect, altersVect;                
	snap_to_edgelist(egosVect, altersVect, G); //converts SNAP graph to sorted edgelist.

	std::ofstream outfile;
	outfile.open("results.txt");
	
	double EXF;
	for (int i = 0; i < nNodes; i++){
		EXF = exfcpp(egosVect, altersVect, i);
		outfile << std::to_string(i) << "  " << std::to_string(EXF) << std::endl;
	}
	outfile.close();
	cout << "Results saved as 'results.txt'" << endl;
	system("pause");

	return 0;
}
