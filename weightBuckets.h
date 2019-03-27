#ifndef WEIGHT_BUCKET
#define WEIGHT_BUCKET
#include "config.h"

class WeightBucket{
private:
	LL *v_ascending_weight_array;
	LL *index_in_v_ascending_weight_array;
	LL *ptr_to_out_vertices_of_weight;
	LL *ptr_to_in_vertices_of_weight;
	LL max_weight;
	LL remain_vertex_num;

public:
	WeightBucket(LL v_num, LL max_w, Vertex* vertices)
	{
		max_weight = max_w;
		remain_vertex_num = 0;
		
		// vertices
		v_ascending_weight_array = new LL[v_num + 2];
		index_in_v_ascending_weight_array = new LL[v_num + 2];
		//cout << "index_in_v_ascending_weight_array: " << index_in_v_ascending_weight_array << endl;

		// pointers
		ptr_to_out_vertices_of_weight = new LL[max_weight + 2];
		ptr_to_in_vertices_of_weight = new LL[max_weight + 2];

		//
		LL* v_weight_num = new LL[max_weight + 1];
		LL* qtr_to_out_vertices_of_weight = new LL[max_weight + 2];

		//
		LL i;
		LL v;
		LL v_w;

		// make statistics
		memset(v_weight_num, 0, sizeof(LL) * (max_weight + 1));
		
		for(v = 1; v <= v_num; v++)
		{
			v_weight_num[vertices[v].get_weight()]++;
		}

		// partition
		ptr_to_out_vertices_of_weight[0] = ptr_to_in_vertices_of_weight[0] = 0;
		qtr_to_out_vertices_of_weight[0] = 0;

		for(v_w = 1; v_w <= max_weight; v_w++)
		{
			ptr_to_in_vertices_of_weight[v_w] = ptr_to_in_vertices_of_weight[v_w-1] + v_weight_num[v_w-1];
		}
		ptr_to_in_vertices_of_weight[max_weight + 1] = v_num;

		for(v_w = 1; v_w <= max_weight; v_w++)
		{
			ptr_to_out_vertices_of_weight[v_w] = ptr_to_in_vertices_of_weight[v_w+1]; // all vertices are in
		}

		memcpy(qtr_to_out_vertices_of_weight, ptr_to_in_vertices_of_weight, sizeof(LL) * (max_weight + 2));

		// place in
		for(v = 1; v <= v_num; v++)
		{
			v_w = vertices[v].get_weight();
			v_ascending_weight_array[qtr_to_out_vertices_of_weight[v_w]] = v;
			index_in_v_ascending_weight_array[v] = qtr_to_out_vertices_of_weight[v_w];
			qtr_to_out_vertices_of_weight[v_w]++;
		}

		//
		delete[] v_weight_num;
		delete[] qtr_to_out_vertices_of_weight;
	}

	~WeightBucket()
	{
#if 1
		//cout << "about to delete v_ascending..." << endl;
		delete[] v_ascending_weight_array;
		//cout << "succeed in deleting v_ascending_weight_array" << endl;
#endif
		//cout << "index_in_v_ascending_weight_array: " << index_in_v_ascending_weight_array << endl;
		delete[] index_in_v_ascending_weight_array;
		//cout << "having freed space" << endl;
#if 1
		delete[] ptr_to_out_vertices_of_weight;
		delete[] ptr_to_in_vertices_of_weight;
#endif
	}

public:
/*
	void clear()
	{
		memcpy(ptr_to_out_vertices_of_weight, ptr_to_in_vertices_of_weight, sizeof(LL) * (max_weight + 2));
	}
*/

	void placeInVertexToCover(const LL v, Vertex* vertices)
	{
		LL v_w = vertices[v].get_weight();
		LL v_loc = index_in_v_ascending_weight_array[v];

		LL boundary_loc = ptr_to_out_vertices_of_weight[v_w];
		LL boundary_vtx = v_ascending_weight_array[boundary_loc];

		v_ascending_weight_array[boundary_loc] = v;
		index_in_v_ascending_weight_array[v] = boundary_loc;

		v_ascending_weight_array[v_loc] = boundary_vtx;
		index_in_v_ascending_weight_array[boundary_vtx] = v_loc;

		ptr_to_out_vertices_of_weight[v_w]++; 
	}

	void placeOutVertexFromCover(const LL v, Vertex* vertices)
	{
		LL v_w = vertices[v].get_weight();
		LL v_loc = index_in_v_ascending_weight_array[v];

		LL boundary_loc = ptr_to_out_vertices_of_weight[v_w] - 1;
		LL boundary_vtx = v_ascending_weight_array[boundary_loc];

		v_ascending_weight_array[boundary_loc] = v;
		index_in_v_ascending_weight_array[v] = boundary_loc;

		v_ascending_weight_array[v_loc] = boundary_vtx;
		index_in_v_ascending_weight_array[boundary_vtx] = v_loc;

		ptr_to_out_vertices_of_weight[v_w]--;
	}

	LL vertex_outside_ds_of_smallest_weight_then_time_stamp(HugeInt* time_stamp)
	{
		LL v_w = 0;
		LL best_v = 0;

		do{
// cout << "to try weight " << v_w << endl;
			LL i = ptr_to_out_vertices_of_weight[v_w];

			if(i < ptr_to_in_vertices_of_weight[v_w + 1])
			{
				best_v = v_ascending_weight_array[i];
				i++;
			}

			for(; i < ptr_to_in_vertices_of_weight[v_w + 1]; i++)
			{
				LL v = v_ascending_weight_array[i];
				if(time_stamp[v] < time_stamp[best_v])
				{
					best_v = v;
				}
			}

// cout << "having obtained best_v = " << best_v << endl;

			if(v_w == max_weight) break;

			v_w++;

		}while(!best_v);

		return best_v;
	}

#ifdef init_reduction_mode
	void filter_out_dead_vertices(bool* is_dead, LL v_num)
	{
		LL i, j, l;

		l = 0;
		i = j = ptr_to_in_vertices_of_weight[0];

		while(i < v_num)
		{

			while(i < ptr_to_out_vertices_of_weight[l])
			{
				LL v = v_ascending_weight_array[i++];
				if(is_dead[v])
				{
#ifdef debug_mode
cout << "filter out " << v << endl;
#endif

					continue;
				}
				v_ascending_weight_array[j] = v;
				index_in_v_ascending_weight_array[v] = j;
				j++;
			}

			ptr_to_out_vertices_of_weight[l] = j;
			while(i < ptr_to_in_vertices_of_weight[l + 1])
			{
				LL v = v_ascending_weight_array[i++];
				if(is_dead[v])
				{
#ifdef debug_mode
cout << "filter out " << v << endl;
#endif
					//index_in_v_ascending_weight_array[v] = -1;
					continue;
				}
				v_ascending_weight_array[j] = v;
				index_in_v_ascending_weight_array[v] = j;
				j++;
			}

			ptr_to_in_vertices_of_weight[l + 1] = j;

			l++;

		}

		while(l < max_weight + 1)
		{
			ptr_to_out_vertices_of_weight[l] = j;
			ptr_to_in_vertices_of_weight[l+1] = j;
			l++;
		}

		remain_vertex_num = j;
	}
#endif

#if 1

	void show_weight_bucket()
	{
		LL i;
		for(LL v_w = 0; v_w <= max_weight; v_w++)
		{
			cout << "weight: " << v_w << endl;
			cout << "in vertices: " << endl;
			for(i = ptr_to_in_vertices_of_weight[v_w]; i < ptr_to_out_vertices_of_weight[v_w]; i++)
			{
				cout << v_ascending_weight_array[i] << "\t";
			}
			cout << endl;
			cout << "out vertices: " << endl;
			for(i = ptr_to_out_vertices_of_weight[v_w]; i < ptr_to_in_vertices_of_weight[v_w + 1]; i++)
			{
				cout << v_ascending_weight_array[i] << "\t";
			}
			cout << endl;
		}


	}

	bool is_in_vertex_in_bucket_of_weight(const LL v, const LL v_w)
	{
		if(ptr_to_in_vertices_of_weight[v_w] <= index_in_v_ascending_weight_array[v] && index_in_v_ascending_weight_array[v] < ptr_to_out_vertices_of_weight[v_w])
			return true;
		return false;
	}

	bool is_out_vertex_in_bucket_of_weight(const LL v, const LL v_w)
	{
		if(ptr_to_out_vertices_of_weight[v_w] <= index_in_v_ascending_weight_array[v] && index_in_v_ascending_weight_array[v] < ptr_to_in_vertices_of_weight[v_w+1])
			return true;
		return false;
	}

#endif


#if 1

	bool check_weight_bucket(Vertex* vertices, LL v_num, bool* is_dead)
	{
// cout << "checking weight buckets" << endl;
		LL i;
		for(i = 0; i < remain_vertex_num; i++)
		{
			if(index_in_v_ascending_weight_array[v_ascending_weight_array[i]] != i) // check bijection
			{
				cout << "bijection error occurs" << endl;
				cout << "at location " << i << ", the element is " << v_ascending_weight_array[i] << endl;
				cout << "but the reverse location is " << index_in_v_ascending_weight_array[v_ascending_weight_array[i]] << endl;
				return false;
			}

			if(is_dead[v_ascending_weight_array[i]]) // cannot be dead
			{
				cout << v_ascending_weight_array[i] << " is_dead but still in the bucket " << endl;
				return false;
			}
		}
		
		for(LL v_w = 0; v_w <= max_weight; v_w++)
		{
			for(i = ptr_to_in_vertices_of_weight[v_w]; i < ptr_to_out_vertices_of_weight[v_w]; i++)  // should be put in the respective bucket
			{

				if(vertices[v_ascending_weight_array[i]].get_weight() != v_w)
				{
					cout << "the weight of the vertex is " << vertices[v_ascending_weight_array[i]].get_weight() << endl;
					cout << "but it is in the bucket of weight " << v_w << endl;
					return false;
				}

			}

			for(i = ptr_to_out_vertices_of_weight[v_w]; i < ptr_to_in_vertices_of_weight[v_w + 1]; i++)  // should be put in the respective bucket
			{

				if(vertices[v_ascending_weight_array[i]].get_weight() != v_w)
				{
					cout << "the weight of the vertex is " << vertices[v_ascending_weight_array[i]].get_weight() << endl;
					cout << "but it is in the bucket of weight " << v_w << endl;
					return false;
				}

			}
		}

		return true;
/////////////////////////
	}

#endif


/*
	bool is_vertex_outside_bucket_of_weight(const LL v, const LL v_w, bool* is_dead)
	{
		
	}
*/
/*
	void clear()
	{
		
	}
*/
};
#endif
