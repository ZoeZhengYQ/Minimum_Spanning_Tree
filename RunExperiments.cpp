/*
CSE6140 HW2
 Name: Yingqiao Zheng
This is an example of how your experiments should look like.
Feel free to use and modify the code below, or write your own experimental code, as long as it produces the desired output.
*/
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <time.h>
#include <fstream>
#include <queue>
#include <unordered_set>
#include <limits>

using namespace std;

typedef pair<int, int> edgePair;

class Graph{
public:
    int vertexAmount;       // vertex amount in both original graph and MST
    int edgesAmount;        // = Graph.edges.size()
    int MSTweight;
    vector< pair<int, edgePair> > edges;      // structure to store graph
    vector< vector< pair<int, int> > > MSTadjList;
    
public:
    Graph() { }
    Graph(int vertex_amount) {
        MSTadjList = vector< vector< pair<int,int> > > (vertex_amount, vector< pair<int,int> >());
    }
    ~Graph() { }
    
    
    /* function to print original graph */
    void print_graph() {
        for (int i = 0; i < edges.size(); ++i) {
                cout << edges[i].second.first << " <--> " << edges[i].second.second
                    << " : [" << edges[i].first << "]" << endl;
        }
    }
    
    void print_MST(){
        for (int i = 0;i < MSTadjList.size(); ++i){
            cout<< i <<" -> ";
            for(auto edges:MSTadjList[i])
                cout << "[" << edges.first<<": " << edges.second << "]" << ", ";
            cout << endl;
        }
        cout << "Weight of MST: "<< MSTweight << endl;
    }
    
    
};

/*
 *  funtion to read graph from file and store
 */
Graph parseEdges(string &graph_file) {
    int vertexAmount = 0, edgeAmount = 0;       // number of vetex and edges in the graph
    int startV = -1, endV = -1;                 // start and end vertex of current edge
    int preStartV = -1, preEndV = -1;           // start and end vertex of last inreading edge
    int weight = 0, preWeight = 0;              // weight of current edge and last inreading edge
    
    
    ifstream infile;
    infile.open(graph_file);

    cout << "Read from file..." << endl;
    infile >> vertexAmount >> edgeAmount;
    cout << vertexAmount << " " << edgeAmount << endl;
    Graph graph(vertexAmount);
    graph.vertexAmount = vertexAmount;
    
    
    for (int i = 0; i < edgeAmount; ++i) {
        infile >> startV >> endV >> weight;
        /* only need to store the edge with lowest weight between same two vertices */
        if (startV == preStartV && endV == preEndV) {
            if (preWeight > weight) {
                preWeight = weight;
                graph.edges.pop_back();
            }
            else continue;
        }
        if (startV != preStartV || endV != preEndV) { preWeight = weight; }
        preStartV = startV;
        preEndV = endV;
        graph.edges.push_back(make_pair(weight, make_pair(startV, endV)));
    }
    graph.edgesAmount = (int)graph.edges.size();
    return graph;
}



int findFather(vector<int> father, int x) {
    int a = x;
    while (x != father[x])
        x = father[x];
    while (a != father[a]) {
        int z = a;
        a = father[a];
        father[z] = x;
        
    }
    return x;
}


int computeMST(Graph &g) {
    int MSTEdgeNum = 0;
    vector<int> father(g.vertexAmount);
    sort(g.edges.begin(), g.edges.end());
    
    for (int i = 0; i < father.size(); i++)
        father[i] = i;
    
    for (int i = 0; i < g.edgesAmount; ++i) {
        int faU = findFather(father, g.edges[i].second.first);
        int faV = findFather(father, g.edges[i].second.second);
        if (faU != faV) {
            father[faU] = faV;
            g.MSTweight += g.edges[i].first;
            MSTEdgeNum++;
            g.MSTadjList[g.edges[i].second.first].push_back(make_pair(g.edges[i].second.second, g.edges[i].first));
            g.MSTadjList[g.edges[i].second.second].push_back(make_pair(g.edges[i].second.first, g.edges[i].first));
            if (MSTEdgeNum == g.vertexAmount - 1)
                break;
        }
    }
    if (MSTEdgeNum != g.vertexAmount - 1)
        return -1;
    else
        return g.MSTweight;
}

bool MST_DFS(Graph &g, int currentNode, int endNode, vector<int> &circle, vector<bool>& visited) {
    visited[currentNode] = true;
    /* if reach the end node then the circle is found */
    if (currentNode == endNode) return true;
    for (auto edge:g.MSTadjList[currentNode]) {
        circle.push_back(edge.first);
        if (!visited[edge.first] && MST_DFS(g, edge.first, endNode, circle, visited)) {
            return true;
        }
        circle.pop_back();
    }
    return false;
}

int recomputeMST(int newStart, int newEnd, int newWeight, Graph &g){
    bool is_in_MST = false;
    for (auto& edges:g.MSTadjList[newStart]) {
        if (edges.first == newEnd) {
            is_in_MST = true;
            if (edges.second > newWeight) {
                g.MSTweight = g.MSTweight - edges.second + newWeight;
                edges.second = newWeight;
            }
            break;
        }
    }
    for (auto& edges:g.MSTadjList[newEnd]) {
        if (edges.first == newStart && edges.second > newWeight) {
            edges.second = newWeight;
            break;
        }
    }
    
    if (!is_in_MST) {
        vector<bool> visited(g.MSTadjList.size(), false);
        vector<int> circle;
        int maxWeight = -1;
        circle.push_back(newStart);
        pair<int, int> maxWeightEdge = make_pair(0, 0);
        if(MST_DFS(g, newStart, newEnd, circle, visited)) {
            for (int i = 0; i < circle.size() - 1; ++i) {
                for (auto edges:g.MSTadjList[circle[i]]) {
                    if (edges.first == circle[i + 1] && edges.second > maxWeight) {
                        maxWeight = edges.second;
                        maxWeightEdge.first = circle[i];
                        maxWeightEdge.second = edges.first;
                        break;
                    }
                }
            }
            if (maxWeight > newWeight) {
                int offset = 0;
                for (auto edges:g.MSTadjList[maxWeightEdge.first]) {
                    if(edges.first == maxWeightEdge.second) {
                        g.MSTadjList[maxWeightEdge.first].erase(g.MSTadjList[maxWeightEdge.first].begin() + offset);
                        break;
                    }
                    ++offset;
                }
                g.MSTadjList[newStart].push_back(make_pair(newEnd, newWeight));
    
                offset = 0;
                for (auto edges:g.MSTadjList[maxWeightEdge.second]) {
                    if(edges.first == maxWeightEdge.first) {
                        g.MSTadjList[maxWeightEdge.second].erase(g.MSTadjList[maxWeightEdge.second].begin() + offset);
                        break;
                    }
                    ++offset;
                }
                g.MSTadjList[newEnd].push_back(make_pair(newStart, newWeight));
                
                g.MSTweight = g.MSTweight - maxWeight + newWeight;
            }
        }
    }
    
    return g.MSTweight;
}





int main(int argc, char *argv[]) {

	/*
	1. inputs: graph file, change file, name of output file
	2. parseEdges to parse graph file
	3. calculate MST (returns integer, weight of MST); we print this integer to the output file
	4. loop through change file, call function pass in new edge and MST
	*/

	if (argc < 4) {
		cout << "Usage: " << argv[0] << " <graph_file> <change_file> <output_file>" << endl;
		return 1;
	}

	string graph_file = argv[1];
	string change_file = argv[2];
	string output_file = argv[3];

	ofstream output;
	output.open(output_file);

	//Write this function to parse edges from graph file to create your graph object
	Graph G = parseEdges(graph_file);
    // G.print_graph();

	//Run your MST function on graph G and collect as output the total weight of the MST
	clock_t startMST = clock();
	int MSTweight = computeMST(G);
	clock_t endMST = clock();

	//Subtract the start time from the finish time to get the actual algorithm running time
	// clock_t totalTime = 1000 * (endMST - startMST) / (float) CLOCKS_PER_SEC;
    
	//Write initial MST weight and time to output file
	output << MSTweight << "\t" << 1000 * double(endMST - startMST) / CLOCKS_PER_SEC << endl;
    cout << "compute: "<<1000 * double(endMST - startMST) / CLOCKS_PER_SEC << endl;
	//Iterate through changes file
	ifstream changes(change_file);

	int newMSTWeight = -1;

	if (changes.is_open()) {
		int numChanges;
		changes >> numChanges; //read number of changes

		int counter = 0;
        clock_t start = clock();
		while (counter < numChanges) {
			int u, v, weight;
			changes >> u; //read u
			changes >> v; //read v
			changes >> weight; //read w(u,v)

			//Run your recomputeMST function to recalculate the new weight of the MST given the addition of this new edge
			//Note: you are responsible for maintaining the MST in order to update the cost without recalculating the entire MST
			clock_t startNewMST = clock();
			newMSTWeight = recomputeMST(u, v, weight, G);
			clock_t endNewMST = clock();
            
			//clock_t totalNewMST = 1000 * double(endNewMST - startNewMST) / CLOCKS_PER_SEC;

			//Write new weight and time to output file
			output << newMSTWeight << "\t" << 1000 * double(endNewMST - startNewMST) / CLOCKS_PER_SEC << endl;

			counter++;
		}
        clock_t end = clock();
        // output << "recompute time:" << 1000 * double(end - start) / CLOCKS_PER_SEC << endl;
        cout << "recompute: "<< 1000 * double(end - start) / CLOCKS_PER_SEC << endl;
		changes.close();
	}
    // G.print_MST();

	output.close();
	return 0;

}
