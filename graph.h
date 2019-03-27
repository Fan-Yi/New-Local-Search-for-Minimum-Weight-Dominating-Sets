#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <string.h>
#include <unordered_set>
#include <unordered_map>
#include <math.h>

#include "config.h"
#include "myBijection.h"

//#define test_graph_construction_mode

using namespace std;

typedef long long LL;

class Vertex
{
private:
	LL *neighbors;
	//LL *adj_edges;
	LL degree;
	LL weight;

#if defined(two_level_cc_mode) || defined(scenario_hash_mode)
	#ifdef two_level_cc_mode
	LL *two_dist_neighbors;
	// LL *two_dist_edges;
	LL two_dist_neighbor_num;
	#endif
#endif

public:
	Vertex()
	{
		neighbors = NULL;
		//adj_edges = NULL;
		degree = 0;
		weight = 0;

#if defined(two_level_cc_mode) || defined(scenario_hash_mode)
	#ifdef two_level_cc_mode
		two_dist_neighbors = NULL;
		// two_dist_edges = NULL;
		two_dist_neighbor_num = 0;
	#endif
#endif
	}

	~Vertex()
	{
		free_neighborhood_space();

#if defined(two_level_cc_mode) || defined(scenario_hash_mode)
	#ifdef two_level_cc_mode
		free_two_dist_neighbor_space();
		// free_two_dist_edge_space();
	#endif
#endif
	}

	void allocate_neighborhood_space(LL dgr)
	{
		neighbors = new LL[dgr];
		//adj_edges = new LL[dgr];
	}

	void free_neighborhood_space()
	{
		delete[] neighbors;
		//delete[] adj_edges;
	}

#if defined(two_level_cc_mode) || defined(scenario_hash_mode)
	#ifdef two_level_cc_mode
	void allocate_two_dist_neighbor_space(LL two_dist_dgr)
	{
		two_dist_neighbors = new LL[two_dist_dgr];
	}

	void free_two_dist_neighbor_space()
	{
		delete[] two_dist_neighbors;
	}
/*
	void allocate_two_dist_edge_space(LL two_dist_dgr)
	{
		two_dist_edges = new LL[two_dist_dgr];
	}

	void free_two_dist_edge_space()
	{
		delete[] two_dist_edges;
	}
*/
	LL* get_two_dist_neighbors()
	{
		return two_dist_neighbors;
	}
/*
	LL* get_two_dist_edges()
	{
		return two_dist_edges;
	}
*/
	LL get_two_dist_neighbor_num()
	{
		return two_dist_neighbor_num;
	}
	
	void set_two_dist_neighbor_num(LL two_d_nb_n)
	{
		two_dist_neighbor_num = two_d_nb_n;
	}
	#endif
#endif


	void add_neighbor(LL name, LL index)
	{
		neighbors[index] = name;
	}
/*
	void add_adj_edge(LL name, LL index)
	{
		adj_edges[index] = name;
	}
*/

	LL *get_neighbors()
	{
		return neighbors;
	}

/*
	LL *get_adj_edges()
	{
		return adj_edges;
	}
*/
	void set_degree(LL d)
	{
		degree = d;
	}

	LL get_degree()
	{
		return degree;
	}

	LL get_weight()
	{
		return weight;
	}

	void set_weight(LL w)
	{
		weight = w;
	}

	void show_neighbors()
	{
		cout << "neighbors: ";
		for(LL i = 0; i < degree; i++)
		{
			cout << neighbors[i] << '\t';
		}
		cout << endl;
	}

};

class Edge
{
private:
	LL v1, v2;

public:
	Edge(){}

	void set_vertices(LL u1, LL u2)
	{
		v1 = u1;
		v2 = u2;
	}

	void get_vertices(LL &u1, LL &u2)
	{
		u1 = v1;
		u2 = v2;
	}

	~Edge(){}
};

class Graph
{
private:
	

protected:
	LL max_degree;
	LL max_weight;

protected:
	Vertex *vertices;
	Edge	*edges;
	LL v_num;
	LL e_num;
	LL within_two_dist_pair_num;
	unordered_set<LL> edge_hash_id_set;
	unordered_map<LL, LL> unordered_pair_hash_id_to_its_index;

#ifdef test_scenario_hash_mode
	Edge *within_two_dist_edges;
#endif

private:
	void encode_pairID(LL &pairID, LL n1, LL n2)
	{
		pairID = ((n1 + n2 + 1) * (n1 + n2) >> 1) + n2;
	}

	void decode_pairID(LL pairID, LL &n1, LL &n2)
	{
		LL w = LL((sqrt(double((pairID << 3) + 1)) - 1) / 2);
		LL t = (w * w + w) >> 1;
		n2 = pairID - t;
		n1 = w - n2; 
	}

	void encode_unordered_pairID(LL &pairID, LL n1, LL n2)
	{
		LL u, v;
		if(n1 < n2)
		{
			u = n1; v = n2;
		}
		else
		{
			u = n2; v = n1;
		}
		encode_pairID(pairID, u, v);		
	}

	void insertEdgeHashIDToSet(LL n1, LL n2)
	{
		LL edge_hash_id;
		encode_unordered_pairID(edge_hash_id, n1, n2);
		edge_hash_id_set.insert(edge_hash_id);
	}

	void insert_vertex_pair_to_edge_into_map(LL v1, LL v2, LL e)
	{
		LL pairID;
		encode_unordered_pairID(pairID, v1, v2);
		unordered_pair_hash_id_to_its_index[pairID] = e;
	}

#ifdef test_graph_construction_mode
void output_info_on_vertex(LL u)
{
cout << "the degree of " << u << " is " << vertices[u].get_degree() << endl;

cout << "its neighbor list: \n";
for(LL i = 0; i < vertices[u].get_degree(); i++)
{
	cout << vertices[u].get_neighbors()[i] << 't';
}
cout << endl;
/*
cout << "its incident edge list: \n";
for(LL i = 0; i < vertices[u].get_degree(); i++)
{
	cout << vertices[u].get_adj_edges()[i] << '\t';
}
*/
cout << endl;
}
#endif

public:
	Graph(char *filename)
	{
#ifdef debug_mode
cout << "begin to build graphs" << endl;
#endif
		ifstream infile(filename);
		if(!infile)
		{
			cout << "File " << filename << " cannot be opened" << endl;
			exit(1);
		}

		char line[1024];
		infile.getline(line, 1024);

		while(line[0] != 'p')
			infile.getline(line, 1024);

		char tempstr1[1024], tempstr2[1024];
		sscanf(line, "%s %s %lld %lld", tempstr1, tempstr2, &v_num, &e_num);

		if(strcmp(tempstr1, "p") != 0 || (strcmp(tempstr2, "edge") != 0 && strcmp(tempstr2, "edges") != 0))
		{
			cout << "format error occurs in reading p lines" << endl;
			exit(1);
		}

#ifdef debug_mode
cout << "have read the p line successfully" << endl;
#endif

		vertices = new Vertex[v_num + 2];
		edges = new Edge[e_num + 2];

		char ch_tmp;
		LL v;

		LL v1, v2;
		LL *v_degree_tmp = new LL[v_num + 2];
		memset(v_degree_tmp, 0, sizeof(LL) * (v_num + 2));

#ifdef debug_mode
cout << "having allocated memories for vertices and edges successfully" << endl;
#endif
		max_weight = 0;
		for(v = 1; v <= v_num; v++)
		{
			LL v_id;
			LL v_weight;
			infile >> ch_tmp >> v_id >> v_weight;
			if(v_id != v)
			{
				cout << "vertex id may have been read incorrectly" << endl;
				exit(1);
			}
			vertices[v].set_weight(v_weight);
			if(v_weight > max_weight)
			{
				max_weight = v_weight;
			}
		}

		LL e;
		for(e = 0; e < e_num; e++)
		{
			infile >> ch_tmp >> v1 >> v2;
			if((v1 < 0 || v1 > v_num) || (v2 < 0 || v2 > v_num))
			{
				cout << "have read an edge connecting " << v1 << " and " << v2 << endl;
				cout << "the vertex indexes are out of range" << endl;
				exit(1);
			}

			edges[e].set_vertices(v1, v2);
			v_degree_tmp[v1]++;
			v_degree_tmp[v2]++;
			insertEdgeHashIDToSet(v1, v2);
			insert_vertex_pair_to_edge_into_map(v1, v2, e);

#ifdef debug_mode
//cout << "having successfully dealt with line " << e << endl;
//cout << "v1: " << v1 << ", v2: " << v2 << endl;
//getchar();
#endif
		}

#ifdef debug_mode
cout << "having read e lines successfully" << endl;
#endif

		for(v = 1; v <= v_num; v++)
		{
			vertices[v].allocate_neighborhood_space(v_degree_tmp[v]);
			vertices[v].set_degree(v_degree_tmp[v]);
		}

#ifdef debug_mode
cout << "having allocated memories for vertex neighborhood successfully" << endl;
#endif

		// compute max_degree
		max_degree = v_degree_tmp[1];
		for(LL i = 2; i <= v_num; i++)
		{
			if(v_degree_tmp[i] > max_degree)
				max_degree = v_degree_tmp[i];
		}

		memset(v_degree_tmp, 0, sizeof(LL) * (v_num + 2));
		for(e = 0; e < e_num; e++)
		{
			edges[e].get_vertices(v1, v2);
			vertices[v1].add_neighbor(v2, v_degree_tmp[v1]);
			vertices[v2].add_neighbor(v1, v_degree_tmp[v2]);
			//vertices[v1].add_adj_edge(e, v_degree_tmp[v1]);
			//vertices[v2].add_adj_edge(e, v_degree_tmp[v2]);
			v_degree_tmp[v1]++;
			v_degree_tmp[v2]++;
		}

		delete[] v_degree_tmp;
		infile.close();

#ifdef debug_mode
cout << "the instance file has been successfully closed" << endl;
#endif

#if defined(two_level_cc_mode) || defined(scenario_hash_mode)
		Bijection* ptr_to_within_two_dist_neighbor_collection = new Bijection(v_num);
		within_two_dist_pair_num = e_num;

		for(v = 1; v <= v_num; v++)
		{

			for(LL i = 0; i < vertices[v].get_degree(); i++)
			{

				LL n1 = vertices[v].get_neighbors()[i];
				if(!ptr_to_within_two_dist_neighbor_collection->element_in(n1))
					ptr_to_within_two_dist_neighbor_collection->insert_element(n1); // add level-1 neighbors

				for(LL j = 0; j < vertices[n1].get_degree(); j++)
				{
					LL n2 = vertices[n1].get_neighbors()[j];
					if(!ptr_to_within_two_dist_neighbor_collection->element_in(n2) && n2 != v)
						ptr_to_within_two_dist_neighbor_collection->insert_element(n2); // add level-2 neighbors
				}

			} // for(LL i = 0; i < vertices[v].get_degree(); i++)
	#ifdef two_level_cc_mode
			vertices[v].set_two_dist_neighbor_num(ptr_to_within_two_dist_neighbor_collection->size());
			vertices[v].allocate_two_dist_neighbor_space(ptr_to_within_two_dist_neighbor_collection->size());
			vertices[v].allocate_two_dist_edge_space(ptr_to_within_two_dist_neighbor_collection->size());
	#endif
			for(LL i = 1; i <= ptr_to_within_two_dist_neighbor_collection->size(); i++)
			{
				LL two_dist_n = ptr_to_within_two_dist_neighbor_collection->at(i);
	#ifdef two_level_cc_mode
				vertices[v].get_two_dist_neighbors()[i] = two_dist_n;
	#endif
				if(!unordered_pair_exist(v, two_dist_n))
				{
					// insertEdgeHashIDToSet(v, two_dist_n);
	#ifdef two_level_cc_mode	
					// vertices[v].get_two_dist_edges()[i] = within_two_dist_pair_num;
	#endif
					insert_vertex_pair_to_edge_into_map(v, two_dist_n, within_two_dist_pair_num);
#ifdef debug_mode
cout << "having just inserted unordered pair (" << v << ", " << two_dist_n << ") into map, mapping to " << within_two_dist_pair_num << endl;
#endif
					within_two_dist_pair_num++;
				}
	#ifdef two_level_cc_mode
				else
				{
					// LL e_index = unordered_pair_index_of(v, two_dist_n);
					// vertices[v].get_two_dist_edges()[i] = e_index;
				}
	#endif			
			}
			ptr_to_within_two_dist_neighbor_collection->clear();

		} // for(v = 1; v <= v_num; v++)

#ifdef test_scenario_hash_mode
		within_two_dist_edges = new Edge[within_two_dist_pair_num];
		for (auto& x: unordered_pair_hash_id_to_its_index) 
		{
			LL v1, v2;
			decode_pairID(x.first, v1, v2);
    	within_two_dist_edges[x.second].set_vertices(v1, v2);
#ifdef debug_mode
cout << "edge: " << x.second << ", set vertices " << v1 << " and " << v2 << endl;
#endif
  	}

#ifdef debug_mode
		for(e = 0; e < within_two_dist_pair_num; e++)
		{
			LL v1, v2;
			within_two_dist_edges[e].get_vertices(v1, v2);
			cout << "Edge " << e << ": " << v1 << ", " << v2 << endl;
		}
#endif

#endif
		delete ptr_to_within_two_dist_neighbor_collection;
#endif

#ifdef test_graph_construction_mode
LL u;
u = 1;
output_info_on_vertex(u);
u = 2;
output_info_on_vertex(u);
u = v_num - 1;
output_info_on_vertex(u);
u = v_num;
output_info_on_vertex(u);
cout << "**************" << endl;
cout << "max degree: " << max_degree << endl;
float d_sum = 0;
float d_avg = 0;
for(u = 1; u <= v_num; u++)
{
	d_sum += float(vertices[u].get_degree());
}
d_avg = d_sum / float(v_num);
cout << "avg degree: " << d_avg << endl;

getchar();

#if defined(two_level_cc_mode) || defined(scenario_hash_mode)
	
	for(u = 1; u <= v_num; u++)
	{
	#ifdef two_level_cc_mode
		cout << "vertex " << u << ":\t" << endl;
		cout << "two dist neighbors:" << endl;
		for(LL i = 0; i < vertices[u].get_two_dist_neighbor_num(); i++)
		{
			cout << vertices[u].get_two_dist_neighbors()[i] << "\t";
		}
		cout << endl;
/*
		cout << "two dist edges:" << endl;
		for(LL i = 0; i < vertices[u].get_two_dist_neighbor_num(); i++)
		{
			cout << vertices[u].get_two_dist_edges()[i] << "\t";
		}
		cout << endl;
*/
	#endif
	}
	cout << "within_two_dist_pair_num: " << within_two_dist_pair_num << endl;
#endif

cout << "vertex weights:" << endl;
for(u = 1; u <= v_num; u++)
{
	cout << u << ": " << vertices[u].get_weight() << endl;
}

getchar();
#endif

	}

	Vertex* get_vertices()
	{
		return vertices;
	}

	LL get_max_degree()
	{
		return max_degree;
	}

	LL get_max_weight()
	{
		return max_weight;
	}

	bool isConnected(LL n1, LL n2)
	{
		LL edge_hash_id;
		encode_unordered_pairID(edge_hash_id, n1, n2);
		return edge_hash_id_set.count(edge_hash_id);
	}

	LL unordered_pair_index_of(LL n1, LL n2)
	{
		LL pairID;
		encode_unordered_pairID(pairID, n1, n2);
		return unordered_pair_hash_id_to_its_index[pairID];
	}

	bool unordered_pair_exist(LL n1, LL n2)
	{
		LL pairID;
		encode_unordered_pairID(pairID, n1, n2);
		return (unordered_pair_hash_id_to_its_index.find(pairID) != unordered_pair_hash_id_to_its_index.end());
	}

	~Graph()
	{
		delete[] vertices;
		delete[] edges;

#ifdef test_scenario_hash_mode
		delete[] within_two_dist_edges;
#endif
	}

};
#endif
