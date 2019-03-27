#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

#include <cstdlib>
#include  <stdio.h>
#include <float.h>
#include <string>

#include <ctime>// include this header 

#include "graph.h"
// #include "myBijection.h"
#include "ansUpdate.h"
#include "hugeInt.h"
#include "operandSets.h"

#include "weightBuckets.h"
#include "dominationHash.h"
#include "scenarioHash.h"
#include "config.h"


class StochasticLocalSearch : private Graph{

private:
	string solver_name;
	string file_name;

	int seed;

	Bijection* ptr_to_dominating_set;
	Bijection* ptr_to_complementary_set;

#ifdef init_reduction_mode
	Bijection* ptr_to_live_vertices;
	bool* is_fixed;
#endif

	LL dominating_set_weight;
	LL* connect_dominating_set_degree;

	Bijection* ptr_to_uncov_vertices;
#ifdef ans_update_op_mode
	AnsChangeSet* ptr_to_moved_v;
#endif
	LL *score;
	HugeInt *time_stamp;

#ifdef edge_time_stamp_mode
	HugeInt* edge_time_stamp;
#endif

#if defined(nvcc_mode) || defined(score_changed_mode) || defined(two_level_cc_mode) || defined(two_level_cc_brute_force_mode)
	bool *confChange;
#endif

#ifdef cc2v3_mode
	char *conf;
#endif

#ifdef unlocking_relation_mode
	LL* unlocker;
#endif

	HugeInt step;
	LL time_limit;

	LL start, stop;
	double best_cmp_time;

	LL first_solution_weight;

	bool *best_in_dominating_set;
	LL best_dominating_set_size;
	LL best_dominating_set_weight;

	HugeInt best_solve_step;

#ifdef individual_analysis_on_init_sls_mode
	double init_time;
	double sls_time;
#endif

#ifdef init_reduction_mode
	Bijection* ptr_to_degree_0_vertices;
	Bijection* ptr_to_degree_1_vertices;
	Bijection* ptr_to_degree_2_vertices;
	bool having_confirmed_optimality;

#endif

#ifdef strategy_analysis_mode
	int ans_update_times;
	#ifdef restart_mode
		int restart_num;
		int restart_num_to_find_best_solution;
	#endif
#endif

	WeightBucket *ptr_to_weight_bucket;
	ConfchangedWeightAgePickSet* ptr_to_no_loss_vertices;

#ifdef domination_hash_mode
	DominationHash *ptr_to_hashed_domination;
#endif

#ifdef scenario_hash_mode
	ScenarioHash *ptr_to_hashed_scenario;
#endif

#ifdef detect_local_optimum_mode
	bool last_step_improved;
	bool penult_step_improved;
#endif

public:
	StochasticLocalSearch(char* solv_name, char *fname, int sd, LL cut_off) : Graph(fname)
	{
#ifdef debug_mode
cout << "the problem instance has been constructed" << endl;
#endif
		solver_name = string(solv_name);
		file_name = string(fname);

		seed = sd;

		step = 0;

		start = clock();
		best_cmp_time = DBL_MAX;

		time_limit = cut_off;

#ifdef individual_analysis_on_init_sls_mode
		init_time = 0;
#endif

		// in_dominating_set = new bool[v_num + 2];
		// memset(in_dominating_set, 0, sizeof(bool) * (v_num + 2));

		// dominating_set_size = v_num;
		dominating_set_weight = 0;

		connect_dominating_set_degree = new LL[v_num + 2];
		// memset(connect_dominating_set_degree, 0, sizeof(LL) * (v_num + 2));

		ptr_to_uncov_vertices = new Bijection(v_num);
		//num_of_uncov_vertices = v_num;

		ptr_to_dominating_set = new Bijection(v_num);
		ptr_to_complementary_set = new Bijection(v_num);

#ifdef init_reduction_mode
		ptr_to_live_vertices = new Bijection(v_num);
		is_fixed = new bool[v_num + 1];
		memset(is_fixed, 0, sizeof(bool) * (v_num + 1));

#endif

#ifdef debug_mode
// cout << 1 << endl;
#endif

		score = new LL[v_num + 2];
		for(LL v = 1; v <= v_num; v++)
		{
			// in_dominating_set[v] = 1;
			ptr_to_dominating_set->insert_element(v);
			dominating_set_weight += vertices[v].get_weight();
			connect_dominating_set_degree[v] = vertices[v].get_degree();
			//score[v] = vertices[v].get_degree() + 1;//
			score[v] = 0;

#ifdef init_reduction_mode
			ptr_to_live_vertices->insert_element(v);
#endif
		}

#ifdef debug_mode
// cout << 1.2 << endl;
#endif
		ptr_to_complementary_set = new Bijection(v_num);

		time_stamp = new HugeInt[v_num + 2];

#if defined(nvcc_mode) || defined(score_changed_mode) || defined(two_level_cc_mode) || defined(two_level_cc_brute_force_mode)
		confChange = new bool[v_num + 2];
#endif

#ifdef cc2v3_mode
		conf = new bool[v_num + 2];
#endif

		for(LL v = 1; v <= v_num; v++)
		{
#ifdef age_init_i_mode
			time_stamp[v] = -v_num - 1 + v;
#endif

#if defined(nvcc_mode) || defined(score_changed_mode) || defined(two_level_cc_mode) || defined(two_level_cc_brute_force_mode)
			confChange[v] = 1;
#endif

#ifdef cc2v3_mode
			conf[v] = 1;
#endif
		}

#ifdef debug_mode
// cout << 1.5 << endl;
#endif

#ifdef unlocking_relation_mode
		unlocker = new LL[v_num + 1];
		memset(unlocker, 0, sizeof(LL) * (v_num + 1));
#endif

		ptr_to_no_loss_vertices = new ConfchangedWeightAgePickSet(v_num);
		for(LL v = 1; v <= v_num; v++)
		{
			if(vertices[v].get_degree())
				ptr_to_no_loss_vertices->insert_element(v);
		}

#ifdef debug_mode
// cout << 2 << endl;
#endif

#ifdef domination_hash_mode
		// ptr_to_hashed_domination = new DominationHash(v_num);
#endif

#ifdef scenario_hash_mode
		ptr_to_hashed_scenario = new ScenarioHash(v_num, within_two_dist_pair_num);
		// last_unlock_edge = new LL[v_num + 1];
		// memset(last_unlock_edge, 0, sizeof(LL) * (v_num + 1));
#endif

#ifdef detect_local_optimum_mode
		last_step_improved = 0;
		penult_step_improved = 0;
#endif

#ifdef ans_update_op_mode
		ptr_to_moved_v = new AnsChangeSet(v_num);
#endif

		ptr_to_weight_bucket = new WeightBucket(v_num, max_weight, vertices);

#if 1

#ifdef restart_mode
	//first_init = true;
#endif

#ifdef strategy_analysis_mode
	ans_update_times = 0;
	#ifdef restart_mode
	restart_num = 0;
	restart_num_to_find_best_solution = 0;
	#endif
#endif

#ifdef debug_mode
//cout << 3 << endl;
#endif


		best_in_dominating_set = new bool[v_num + 2];
		memset(best_in_dominating_set, 1, sizeof(bool) * (v_num + 2));
		best_dominating_set_size = v_num;
		best_dominating_set_weight = dominating_set_weight;
		update_best_dominating_set();

#ifdef init_reduction_mode
		init_solution_with_reductions();
#else
		init_solution();
#endif

#ifdef debug_mode
cout << "having successfully initialized a solution*********************" << endl;
//getchar();
#endif




#endif
/*
cout << "having successfully initialized a solution*********************" << endl;
delete ptr_to_weight_bucket;
cout << "~~~~~~~~~~~~" << endl;
exit(1);
*/
	}

private:

	bool check_solution()
	{
		LL v;
		LL verified_best_dominating_set_size = 0;
		LL verified_best_dominating_set_weight = 0;

		for(v = 1; v <= v_num; v++)
		{
//cout << best_in_dominating_set[v] << endl;
			if(best_in_dominating_set[v])
			{
				verified_best_dominating_set_size++;
				verified_best_dominating_set_weight += vertices[v].get_weight();
			}
		}

		if(verified_best_dominating_set_size != best_dominating_set_size)
		{
			cout << "verified_best_dominating_set_size: " << verified_best_dominating_set_size << endl;
			cout << "best_dominating_set_size: " << best_dominating_set_size << endl;
			cout << "the dominating set size is computed incorrectly" << endl;
			return 0;
		}

		if(verified_best_dominating_set_weight != best_dominating_set_weight)
		{
			cout << "verified_best_dominating_set_weight: " << verified_best_dominating_set_weight << endl;
			cout << "best_dominating_set_weight: " << best_dominating_set_weight << endl;
			cout << "the dominating set weight is computed incorrectly" << endl;
			return 0;
		}

		for(v = 1; v <= v_num; v++)
		{
			if(best_in_dominating_set[v]) continue;
			bool flag = false;
			LL *nbs = vertices[v].get_neighbors();
			for(LL i = 0; i < vertices[v].get_degree(); i++)
			{
				if(best_in_dominating_set[nbs[i]])
				{
//cout << v << " is covered by " << nbs[i] << endl;
					flag = true;
					break;
				}
			}
			if(!flag)
			{
				cout << "vertex " << v << " is uncovered" << endl;
				return false;
			}
		}

		return true;
	}

	bool check_current_cand_solution()
	{
		LL v;
		LL verified_dominating_set_size = 0;
		LL verified_dominating_set_weight = 0;
		
		for(v = 1; v <= v_num; v++)
		{
			if(ptr_to_dominating_set->element_in(v))
			{
				verified_dominating_set_size++;
				verified_dominating_set_weight += vertices[v].get_weight();
			}
		}

		if(verified_dominating_set_size != ptr_to_dominating_set->size())
		{
			cout << "verified_dominating_set_size: " << verified_dominating_set_size << endl;
			cout << "dominating_set_size: " << ptr_to_dominating_set->size() << endl;
			cout << "the dominating set size is computed incorrectly" << endl;
			return 0;
		}

		if(verified_dominating_set_weight != dominating_set_weight)
		{
			cout << "verified_dominating_set_weight: " << verified_dominating_set_weight << endl;
			cout << "dominating_set_weight: " << dominating_set_weight << endl;
			cout << "the dominating set weight is computed incorrectly" << endl;
			return 0;			
		}
//cout << "before checking domination" << endl;
		LL verified_uncovered_vertex_num = 0;
		for(v = 1; v <= v_num; v++)
		{
			if(ptr_to_dominating_set->element_in(v)) continue;
			bool flag = false;
			LL *nbs = vertices[v].get_neighbors();
			for(LL i = 0; i < vertices[v].get_degree(); i++)
			{
				if(ptr_to_dominating_set->element_in(nbs[i]))
				{
					flag = 1;
					break;
				}
			}
			if(!flag)
			{
				verified_uncovered_vertex_num++;
			}
		}
//cout << "after checking domination" << endl;
		if(verified_uncovered_vertex_num != ptr_to_uncov_vertices->size())
		{
			cout << "verified_uncovered_vertex_num: " << verified_uncovered_vertex_num << endl;
			cout << "num_of_uncov_vertices: " << ptr_to_uncov_vertices->size() << endl;
			cout << "the number of uncovered vertices is computed incorrectly" << endl;
			return false;
		}

		return true;
	}

#ifdef init_reduction_mode
	void put_vertex_away_from_dominating_set_by_reductions(const LL v)
	{
#ifdef debug3_mode
cout << "having just entered the put_vertex_away function" << endl;
if(v == 583600)
{
	cout << "to put away " << 583600 << endl;
	getchar();
}
#endif

#if 0
		if(connect_dominating_set_degree[v] == 2)
			ptr_to_degree_2_vertices->delete_element(v);
		else if(connect_dominating_set_degree[v] == 1)
			ptr_to_degree_1_vertices->delete_element(v); 
#endif


#ifdef debug_mode
//cout << "dominating_set_weight: " << dominating_set_weight << endl;
//show_state();

if(!check_current_cand_solution())
{
	cout << "the current solution may be maintained incorrectly" << endl;
	exit(1);
}

#endif

#ifdef debug_mode
cout << "==================delete " << v << " loss: " << -score[v] << ", weight: " << vertices[v].get_weight() << endl;
if(score[v] != 0)
{
	cout << "this vertex will be put away incorrectly because it has a nonzero loss value" << endl;
	exit(1);
}
#endif

#ifdef scenario_hash_mode
		ptr_to_hashed_scenario->update_hash_wrt_remove(v); // no tabu, no unlocking relation, but vertex inclusion
#endif

#ifdef debug_mode
		LL dominating_set_weight_by_brute_force = 0;
		for(LL i = 1; i <= ptr_to_dominating_set->size(); i++)
		{
			dominating_set_weight_by_brute_force += vertices[ptr_to_dominating_set->at(i)].get_weight();
		}
//cout << "before the removal, dominating_set_weight_by_brute_force: " << dominating_set_weight_by_brute_force << endl;
//cout << "dominating_set_weight before being decreased by " << v << "\'s removal is " << dominating_set_weight << endl;
#endif

#ifdef debug_mode
cout << "to delete " << v << " from dominating sets" << endl;
#endif

		ptr_to_dominating_set->delete_element(v);
		// ptr_to_complementary_set->insert_element(v); // cannot enter complementary set, never add any more

		dominating_set_weight -= vertices[v].get_weight();

#ifdef debug_mode
		dominating_set_weight_by_brute_force = 0;
		for(LL i = 1; i <= ptr_to_dominating_set->size(); i++)
		{
			dominating_set_weight_by_brute_force += vertices[ptr_to_dominating_set->at(i)].get_weight();
		}
		//cout << "after removal, dominating_set_weight_by_brute_force: " << dominating_set_weight_by_brute_force << endl;
		//cout << "dominating_set_weight after being decreased by " << v << "\'s removal is " << dominating_set_weight << endl;
#endif


		ptr_to_weight_bucket->placeOutVertexFromCover(v, vertices);

#ifdef debug_mode
cout << 1 << endl;
cout << "to delete " << v << " from no-loss-vertex-set" << endl;
#endif
		ptr_to_no_loss_vertices->delete_element(v);
#ifdef debug_mode
cout << 2 << endl;
#endif

		LL* nbs = vertices[v].get_neighbors();
		LL dgr = vertices[v].get_degree();

		for(LL i = 0; i < dgr; i++)
		{
			LL n = nbs[i];

			connect_dominating_set_degree[n]--;
			
			if(!is_fixed[n])
			{			
				if(connect_dominating_set_degree[n] == 2)
				{
					ptr_to_degree_2_vertices->insert_element(n);
				}
				else if(connect_dominating_set_degree[n] == 1)
				{
#ifdef debug_mode
cout << 3 << endl;
#endif
					if(ptr_to_degree_2_vertices->element_in(n))					
						ptr_to_degree_2_vertices->delete_element(n);
#ifdef debug_mode
cout << 3.5 << endl;
#endif
					ptr_to_degree_1_vertices->insert_element(n);
				}
				else if(connect_dominating_set_degree[n] == 0)
				{
#ifdef debug_mode
cout << 4 << endl;
#endif
					ptr_to_degree_1_vertices->delete_element(n);
#ifdef debug_mode
cout << 4.5 << endl;
#endif
					ptr_to_degree_0_vertices->insert_element(n);
				}
			}

			if(ptr_to_dominating_set->element_in(n))
			{
				if(connect_dominating_set_degree[v] == 1)// v becomes critically covered
				{
					if(!score[n])
					{
#ifdef debug_mode
cout << 5 << endl;
#endif
						ptr_to_no_loss_vertices->delete_element(n);
#ifdef debug_mode
cout << 5.5 << endl;
#endif
					}
					score[n]--;  // inside
				}

				if(connect_dominating_set_degree[n] == 0)// n becomes critically covered
				{
					if(!score[n])
					{
#ifdef debug_mode
cout << 6 << endl;
#endif
						ptr_to_no_loss_vertices->delete_element(n);
#ifdef debug_mode
cout << 6.5 << endl;
#endif
					}
					score[n]--;  // inside
				}
			}
			else
			{
				if(connect_dominating_set_degree[v] == 0)// v becomes uncovered
				{
					score[n]++;  // outside
				}
				if(connect_dominating_set_degree[n] == 0)// n becomes uncovered
				{
					score[n]++;  // outside
					ptr_to_uncov_vertices->insert_element(n);
					LL* nbs_ = vertices[n].get_neighbors();
					LL dgr_ = vertices[n].get_degree();
					for(LL j = 0; j < dgr_; j++)
					{
						LL n_ = nbs_[j];
						if(!ptr_to_dominating_set->element_in(n_) && n_ != v)
						{
							score[n_]++;  // outside
						}
					}
				}
				else if(connect_dominating_set_degree[n] == 1)
				{
					LL* nbs_ = vertices[n].get_neighbors();
					LL dgr_ = vertices[n].get_degree();
					for(LL j = 0; j < dgr_; j++)
					{
						LL n_ = nbs_[j];
						if(ptr_to_dominating_set->element_in(n_))
						{
							if(!score[n_])
								ptr_to_no_loss_vertices->delete_element(n_);
							score[n_]--;  // inside
							break;
						}
					}					
				}
			}
		}
#ifdef ans_update_op_mode
		ptr_to_moved_v->ans_update(v);
#endif

#ifdef debug_mode
cout << "dominating_set_weight: " << dominating_set_weight << endl;
//show_state();

if(!check_current_cand_solution())
{
	cout << "the current solution may be maintained incorrectly" << endl;
	exit(1);
}

#endif
//show_state();
		step++;

#ifdef debug_mode
// cout << 2.9 << endl;
#endif		
	}


	void init_solution_with_reductions()
	{
		ptr_to_degree_0_vertices = new Bijection(v_num);
		ptr_to_degree_1_vertices = new Bijection(v_num);
		ptr_to_degree_2_vertices = new Bijection(v_num);

		LL* temp_array = new LL[max_degree];
		LL temp_array_size = 0;

		bool reduction_rule_applied = false;

		for(LL v = 1; v <= v_num; v++)
		{
			if(connect_dominating_set_degree[v] == 0)
			{
				is_fixed[v] = 1;
				ptr_to_live_vertices->delete_element(v);
#ifdef debug_mode
cout << "having put away " << v << " by Degree-0 Reduction Rule" << endl;
getchar();
#endif
			}
			else if(connect_dominating_set_degree[v] == 1)
			{
				ptr_to_degree_1_vertices->insert_element(v);
			}
			else if(connect_dominating_set_degree[v] == 2)
			{
				ptr_to_degree_2_vertices->insert_element(v);
			}
		}

		// analysis
#ifdef debug_mode
cout << "reductions" << endl;
#endif	
		do{

			reduction_rule_applied = false;

			while(ptr_to_degree_0_vertices->size()) // degree-0 rule
			{
				LL v = ptr_to_degree_0_vertices->at(1);
				ptr_to_degree_0_vertices->delete_element(v);

				is_fixed[v] = 1;
				ptr_to_no_loss_vertices->delete_element(v);

				ptr_to_live_vertices->delete_element(v);
				reduction_rule_applied = true;
#ifdef debug_mode
cout << "having put away " << v << " by Degree-0 Reduction Rule" << endl;
//getchar();
#endif

#ifdef debug3_mode
		for(LL i = 1; i <= ptr_to_no_loss_vertices->size(); i++)
			if(is_fixed[ptr_to_no_loss_vertices->at(i)])
			{
cout << "Degree-0 Rule being applied" << endl;
cout << ptr_to_no_loss_vertices->at(i) << " is dead but still in no_loss_set" << endl;
exit(1);
			}
#endif
			}
		
			while(ptr_to_degree_1_vertices->size()) // degree-1 rule
			{
				LL v;
				v = ptr_to_degree_1_vertices->at(1);

#ifdef debug_mode
cout << 3 << endl;
#endif
				ptr_to_degree_1_vertices->delete_element(v);

#ifdef debug_mode
cout << 4 << endl;
#endif

#ifdef debug_mode
cout << "considering " << v << " for Degree-1 reductions" << endl;
#endif


				if(is_fixed[v]) 
				{
#ifdef debug_mode
cout << "no considerations because " << v << " are fixed" << endl;
#endif
					continue;
				}

				LL u;

				for(LL j = 0; j < vertices[v].get_degree(); j++)
				{
					u = vertices[v].get_neighbors()[j];
					if(ptr_to_dominating_set->element_in(u))
					{
#ifdef debug_mode
cout << "having obtained its unique neighbor " << u << endl;
#endif
						break;
					} 
				}

				if(is_fixed[u])
				{
					continue;
				}
			
				LL w;
			
				for(LL j = 0; j < vertices[u].get_degree(); j++)
				{
					w = vertices[u].get_neighbors()[j];
					if(ptr_to_dominating_set->element_in(w) && connect_dominating_set_degree[w] == 1 && !is_fixed[w])
						temp_array[temp_array_size++] = w;
				}

				if(temp_array_size > 1)
				{
#ifdef debug_mode
cout << "considering Degree-1 Reduction Rule 2" << endl;
#endif
					LL neighbor_weight_sum = 0;

					for(LL i = 0; i < temp_array_size; i++)
					{
						neighbor_weight_sum += vertices[temp_array[i]].get_weight();
#ifdef debug_mode
cout << "considering a neighbor: " << temp_array[i] << ", weight: " << vertices[temp_array[i]].get_weight() << endl;
#endif
					}

					if(neighbor_weight_sum > vertices[u].get_weight())  // degree-1 rule 2
					{
						for(LL i = 0; i < temp_array_size; i++)
						{
							w = temp_array[i];
							put_vertex_away_from_dominating_set_by_reductions(w);
							is_fixed[w] = 1;
#ifdef debug_mode
cout << w << " is fixed" << endl;
#endif

#ifdef debug_mode
cout << 5 << endl;
#endif
							ptr_to_live_vertices->delete_element(w);
#ifdef debug_mode
cout << 6 << endl;
#endif

#if defined(debug_mode)
cout << "having put away " << w << " by Degree-1 Reduction Rule 2 " << endl;
//getchar();
#endif
						}
						is_fixed[u] = 1;
#ifdef debug_mode
cout << 7 << endl;
#endif
						//ptr_to_no_loss_vertices->delete_element(u);
#if defined(debug_mode)
cout << u << " is fixed" << endl;
#endif

#ifdef debug_mode
cout << 7 << endl;
#endif
						ptr_to_live_vertices->delete_element(u);
#ifdef debug_mode
cout << 8 << endl;
#endif
						reduction_rule_applied = true;

						
					}

#ifdef debug3_mode
		for(LL i = 1; i <= ptr_to_no_loss_vertices->size(); i++)
			if(is_fixed[ptr_to_no_loss_vertices->at(i)])
			{
cout << "Degree-1 Rule 2 being applied" << endl;
cout << ptr_to_no_loss_vertices->at(i) << " is dead but still in no_loss_set" << endl;
exit(1);
			}
#endif

					temp_array_size = 0;

					continue;
				}

				temp_array_size = 0;


				if(vertices[v].get_weight() > vertices[u].get_weight()) // degree-1 Rule 1
				{
#ifdef debug_mode
cout << "Degree-1 Reduction Rule 1 is being applied" << endl;
#endif
					put_vertex_away_from_dominating_set_by_reductions(v);
					is_fixed[u] = is_fixed[v] = true;
					ptr_to_live_vertices->delete_element(u);
					ptr_to_live_vertices->delete_element(v);
					//ptr_to_no_loss_vertices->delete_element(u);
					reduction_rule_applied = true;
#ifdef debug_mode
cout << "having put away " << v << " by Degree-1 Reduction Rule 1 " << endl;
cout << u << " and " << v << " have been fixed " << endl;
//getchar();
#endif

#ifdef debug3_mode
		for(LL i = 1; i <= ptr_to_no_loss_vertices->size(); i++)
			if(is_fixed[ptr_to_no_loss_vertices->at(i)])
			{
cout << "Degree-1 Rule 1 being applied" << endl;
cout << ptr_to_no_loss_vertices->at(i) << " is dead but still in no_loss_set" << endl;
exit(1);
			}
#endif

				}

			}

			while(ptr_to_degree_2_vertices->size()) // degree-2 rule
			{
				LL v1, v2;
				LL u;
				v1 = ptr_to_degree_2_vertices->at(1);

#ifdef debug_mode
cout << "considering " << v1 << " for Degree-2 reductions" << endl;
#endif
				ptr_to_degree_2_vertices->delete_element(v1);
#if defined(debug_mode)
cout << "having just deleted " << v1 << " from *ptr_to_degree_2_vertices" << endl;

#endif
				if(is_fixed[v1])
				{

#ifdef debug_mode
cout << "no considerations because " << v1 << " are fixed" << endl;
#endif				
					continue;
				}


				LL n1, n2;			
				bool flag = 0;
				for(LL i = 0; i < vertices[v1].get_degree(); i++)
				{
					LL v = vertices[v1].get_neighbors()[i];
#ifdef debug_mode
// cout << "considering whether " << v << " is a companion" << endl;
#endif
					
					if(ptr_to_dominating_set->element_in(v))
					{
						if(flag == 0)
						{
							n1 = v;
							flag = 1;
						}
						else
						{
							n2 = v;
							break;
						}
					}
				}
#ifdef debug_mode
cout << "n1: " << n1 << endl;
cout << "n2: " << n2 << endl;
#endif
				if(is_fixed[n1] || is_fixed[n2]) continue;

				if(connect_dominating_set_degree[n1] == 2)
				{
					v2 = n1;
					u = n2;
				}
				else if(connect_dominating_set_degree[n2] == 2)
				{
					v2 = n2;
					u = n1;
				}
				else
				{
					continue;
				}

				if(isConnected(u, v2))
				{
#if defined(debug_mode) || defined(debug3_mode)
// cout << "having obtained " << v2 << " and " << u << " as companions" << endl;
#endif					
					if((vertices[v1].get_weight() > vertices[u].get_weight()) && (vertices[v2].get_weight() > vertices[u].get_weight()))
					{
						put_vertex_away_from_dominating_set_by_reductions(v1);
#ifdef debug3_mode
if(ptr_to_no_loss_vertices->element_in(v1))
{
	cout << v1 << " should not be in no_loss_set after being put away " << v1 << endl;
	exit(1);
}
#endif
						put_vertex_away_from_dominating_set_by_reductions(v2);
#ifdef debug3_mode
if(ptr_to_no_loss_vertices->element_in(v1))
{
	cout << v1 << " should not be in no_loss_set after being put away " << v2 << endl;
	exit(1);
}
#endif
					//}

					is_fixed[u] = is_fixed[v1] = is_fixed[v2] = true;
#if defined(debug_mode) || defined(debug3_mode)
//cout << u << ", " << v1 << ", " << v2 << " are fixed" << endl;
#endif

					ptr_to_live_vertices->delete_element(u);
					ptr_to_live_vertices->delete_element(v1);
					ptr_to_live_vertices->delete_element(v2);
					//ptr_to_no_loss_vertices->delete_element(u);
					reduction_rule_applied = true;

#ifdef debug3_mode
		for(LL i = 1; i <= ptr_to_no_loss_vertices->size(); i++)
			if(is_fixed[ptr_to_no_loss_vertices->at(i)])
			{
cout << "Degree-2 Rule being applied" << endl;
cout << ptr_to_no_loss_vertices->at(i) << " is dead but still in no_loss_set" << endl;
exit(1);
			}
#endif

#ifdef debug_mode
cout << "having put away " << v1 << " and " << v2 << " away by Degree-2 Reduction Rule " << endl;
//getchar();
#endif
					}
				// continue;
				}
/*
				if(v1 has a neighbor in dominating set)
				{
					v2 <- v1's other neighbor;
					if(v2 is fixed) continue;
					if()
				}
*/
			}

#ifdef debug3_mode
		for(LL i = 1; i <= ptr_to_no_loss_vertices->size(); i++)
			if(is_fixed[ptr_to_no_loss_vertices->at(i)])
			{
cout << "reduction rules are being applied" << endl;
cout << ptr_to_no_loss_vertices->at(i) << " is dead but still in no_loss_set" << endl;
exit(1);
			}
#endif

		}while(reduction_rule_applied);

#ifdef debug_mode
cout << "reductions ended" << endl;
//getchar();
#endif

		ptr_to_weight_bucket->filter_out_dead_vertices(is_fixed, v_num);

#ifdef debug3_mode		
		if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
		{
			exit(1);
		}

		cout << "right after filtering, having checked weight buckets right after initialization with reductions. Success." << endl;
		for(LL i = 1; i <= ptr_to_no_loss_vertices->size(); i++)
			if(is_fixed[ptr_to_no_loss_vertices->at(i)])
			{
cout << ptr_to_no_loss_vertices->at(i) << " is dead but still in no_loss_set" << endl;
exit(1);
			}
		getchar();
#endif

#ifdef debug_mode
cout << "now pure rand construction" << endl;
#endif

	while(ptr_to_no_loss_vertices->size())
	{
		LL remove_v = ptr_to_no_loss_vertices->rand_element();
		remove(remove_v);
	}

#ifdef debug3_mode		
		if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
		{
			exit(1);
		}

		cout << "right after removing random redundante vertices, having checked weight buckets. Success." << endl;
		getchar();
#endif

#ifdef detect_local_optimum_mode
	last_step_improved = true;
#endif

	having_confirmed_optimality = true;
	for(LL i = 1; i <= ptr_to_dominating_set->size(); i++)
	{
#ifdef debug_mode
cout << "checking whether " << ptr_to_dominating_set->at(i) << " is fixed" << endl;
#endif
		if(!is_fixed[ptr_to_dominating_set->at(i)])
		// can I simply remove all fixed vertices from the neighbor lists
		{
			having_confirmed_optimality = false;
			break;
		}
	}

#ifdef debug_mode
	cout << "having_confirmed_optimality: " << having_confirmed_optimality << endl;
#endif

#ifdef strategy_analysis_mode
ans_update_times++;
#endif
		
	best_dominating_set_size = ptr_to_dominating_set->size();
	first_solution_weight = best_dominating_set_weight = dominating_set_weight;

#ifdef ans_update_op_mode

	#ifdef debug_mode
	cout << "answer changed vertices: " << endl;
	ptr_to_moved_v->show_elements();
	#endif

	ptr_to_moved_v->efficiently_update_best_dominating_set_vertices();
#else
	update_best_dominating_set();
#endif

	best_solve_step = step;


	stop = clock();
	best_cmp_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;


#ifndef only_test_init_mode
		if(having_confirmed_optimality)
#endif
		{

#ifdef ans_update_op_mode
			bool* in_dominating_set = new bool[v_num + 1];
			for(LL v = 1; v <= v_num; v++)
			{
				in_dominating_set[v] = ptr_to_dominating_set->element_in(v);
			}
			ptr_to_moved_v->compute_final_answer(v_num, in_dominating_set, best_in_dominating_set);
			delete[] in_dominating_set;
#endif

			if(check_solution())
				{
					cout << "o " << best_dominating_set_weight << endl;
					cout << "c size " << best_dominating_set_size << endl;
					cout << "c searchSteps " << best_solve_step <<endl;
					cout << "c solveTime(ms) " << best_cmp_time << endl;

#ifdef init_reduction_mode
					if(having_confirmed_optimality)
						cout << "c optimality confirmed" << endl;
					else
						cout << "c optimality non-confirmed" << endl;
#endif

				}
				else
				{
					cout << "the solution is wrong." << endl;
				}
				exit(0);

		}

#ifdef individual_analysis_on_init_sls_mode
	init_time = round(best_cmp_time * 100) / 100.0;
#endif	

#ifdef debug_mode
cout << "after initialization, " << endl;
show_state();
#endif

		delete[] temp_array;
		
#ifdef debug3_mode
		
		if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
		{
			exit(1);
		}

		cout << "at the end of init with reductions, having checked weight buckets. Success." << endl;
		getchar();
#endif

#ifdef debug_mode
cout << "end of initilization" << endl;
#endif
	}
#endif

	void clear()
	{
		while(ptr_to_complementary_set->size())
		{
			LL v = ptr_to_complementary_set->at(ptr_to_complementary_set->begin());
			add(v);
		}
	
		while(ptr_to_no_loss_vertices->size())
		{
			LL v = ptr_to_no_loss_vertices->rand_element();
			remove(v);
		}

#ifdef detect_local_optimum_mode
		last_step_improved = true;
#endif

#ifdef scenario_hash_mode
		ptr_to_hashed_scenario->forget_all_visits();
#endif
	}


void init_solution() //random construction mode
{
#ifdef debug_mode
cout << "pure rand construction" << endl;
#endif

	while(ptr_to_no_loss_vertices->size())
	{
		LL remove_v = ptr_to_no_loss_vertices->rand_element();
		remove(remove_v);
	}

#ifdef detect_local_optimum_mode
	last_step_improved = 1;
#endif

#ifdef strategy_analysis_mode
ans_update_times++;
#endif

#ifdef debug_mode
cout << "ptr_to_dominating_set->size(): " << ptr_to_dominating_set->size() << endl;
#endif
		
	best_dominating_set_size = ptr_to_dominating_set->size();
	first_solution_weight = best_dominating_set_weight = dominating_set_weight;

#ifdef ans_update_op_mode
	ptr_to_moved_v->efficiently_update_best_dominating_set_vertices();
#else
	update_best_dominating_set();
#endif

	best_solve_step = step;


	stop = clock();
	best_cmp_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;


#ifdef only_test_init_mode
		{
	#ifdef ans_update_op_mode
			bool* in_dominating_set = new bool[v_num + 1];
			for(LL v = 1; v <= v_num; v++)
			{
				in_dominating_set[v] = ptr_to_dominating_set->element_in(v);
			}
	
			ptr_to_moved_v->compute_final_answer(v_num, in_dominating_set, best_in_dominating_set);
	
			delete[] in_dominating_set;
	#endif

			if(check_solution())
				{
					cout << "o " << best_dominating_set_weight << endl;
					cout << "c size " << best_dominating_set_size << endl;
					cout << "c searchSteps " << best_solve_step <<endl;
					cout << "c solveTime(ms) " << best_cmp_time << endl;
				}
				else
				{
					cout << "the solution is wrong." << endl;
				}
				exit(0);
		}
#endif

#ifdef individual_analysis_on_init_sls_mode
	init_time = round(best_cmp_time * 100) / 100.0;
#endif	

#ifdef debug_mode
cout << "after initilization, " << endl;
show_state();
#endif
}



	void update_best_dominating_set()
	{
		for(LL v = 1; v <= v_num; v++)
		{
			best_in_dominating_set[v] = ptr_to_dominating_set->element_in(v);
		}
	}

#ifdef two_level_cc_mode
	void set_two_level_neighbors_free(const LL v)
	{
		for(LL i = 0; i < vertices[v].get_two_dist_neighbor_num(); i++)
		{
			LL n = vertices[v].get_two_dist_neighbors()[i];
			confChange[n] = 1;
		}
	}
#endif

#ifdef two_level_cc_brute_force_mode
	void set_two_level_neighbors_free(const LL v)
	{
		for(LL i = 0; i < vertices[v].get_degree(); i++)
		{
			LL n1 = vertices[v].get_neighbors()[i];
			confChange[n1] = true;
			for(LL j = 0; j < vertices[n1].get_degree(); j++)
			{
				LL n2 = vertices[n1].get_neighbors()[j];
				if(n2 == v) continue;
				confChange[n2] = true;
			}
		}
	}
#endif

#ifdef cc2v3_mode
	void cc2v3rule2(const LL v)
	{
		for(LL i = 0; i < vertices[v].get_degree(); i++)
		{
			LL n1 = vertices[v].get_neighbors()[i];
			conf[n1] = 2;
			for(LL j = 0; j < vertices[n1].get_degree(); j++)
			{
				LL n2 = vertices[n1].get_neighbors()[j];
				if(n2 == v) continue;
				conf[n2] = 2;
			}
		}

		for(LL i = 0; i < vertices[v].get_degree(); i++)
		{
			LL n = vertices[v].get_neighbors()[i];
			conf[n] = 1;
		}
	}

	void cc2v3rule3(const LL v)
	{
		conf[v] = 0;
		for(LL i = 0; i < vertices[v].get_degree(); i++)
		{
			LL n1 = vertices[v].get_neighbors()[i];
			conf[n1] = 2;
			for(LL j = 0; j < vertices[n1].get_degree(); j++)
			{
				LL n2 = vertices[n1].get_neighbors()[j];
				if(n2 == v) continue;
				conf[n2] = 2;
			}
		}
	}
#endif

	void unlock(const LL v, const LL n)
	{

#ifdef score_changed_mode
		if(confChange[n])
			return;

	#ifdef forbid_unlocking_twice_in_a_row_mode
		if(unlocker[n] != v)
	#endif
	#ifdef forbid_reverse_unlocking
		if(unlocker[v] != n)
	#endif
		{
#ifdef scenario_hash_mode
			LL last_unlock_e;
				
			if(unlocker[n] != 0) // not the first time to be unlocked, should delete older binaries
			{
				last_unlock_e = unordered_pair_index_of(n, unlocker[n]);
				if(unlocker[n] > n)
				{
					ptr_to_hashed_scenario->update_hash_wrt_delete_big_id_to_small_unlock_edge(last_unlock_e);
				}
				else // unlocker[n] < n
				{
					ptr_to_hashed_scenario->update_hash_wrt_delete_small_id_to_big_unlock_edge(last_unlock_e);
				}
			}
#endif
				
#ifdef unlocking_relation_mode
			unlocker[n] = v;
#endif

#ifdef scenario_hash_mode
			last_unlock_e = unordered_pair_index_of(v, n);
#endif

#ifdef debug_mode
// cout << v << " unlocks " << n << ", via edge: " << last_unlock_e << endl;
#endif

#ifdef scenario_hash_mode
			if(unlocker[n] > n)
			{
				ptr_to_hashed_scenario->update_hash_wrt_insert_big_id_to_small_unlock_edge(last_unlock_e);
			}
			else // unlocker[n] < n
			{
				ptr_to_hashed_scenario->update_hash_wrt_insert_small_id_to_big_unlock_edge(last_unlock_e);
			}
#endif

			confChange[n] = true;
#ifdef scenario_hash_mode
			ptr_to_hashed_scenario->update_hash_wrt_unlock(n);
#endif
				
		}
#endif		
	}

	void add(const LL v)
	{
#if defined(debug_mode) || defined(debug3_mode)
cout << "==================add " << v << " gain: " << score[v] << endl;
cout << "in dominating set of " << v << " is " << ptr_to_dominating_set->element_in(v) << endl;
if(score[v] != 0)
{
	cout << "in our algorithm, gain should be 0, but now is " << score[v] << endl;
	exit(1);
}
if(ptr_to_dominating_set->element_in(v))
{
	cout << "cannot add " << v << ", since it has already been in the candidate solution" << endl;
	exit(1);
}

#endif

		ptr_to_dominating_set->insert_element(v);
		ptr_to_complementary_set->delete_element(v);

#ifdef scenario_hash_mode
		ptr_to_hashed_scenario->update_hash_wrt_add(v);
#endif

#if defined(nvcc_mode) || defined(score_changed_mode) || defined(two_level_cc_mode) || defined(two_level_cc_brute_force_mode)
		if(confChange[v])
		{
			confChange[v] = false;
	#ifdef scenario_hash_mode
			ptr_to_hashed_scenario->update_hash_wrt_lock(v);
	#endif
		}
#endif

		// in_dominating_set[v] = 1;
		if(!connect_dominating_set_degree[v])//may add an uncovered vertex, so the connect degree may be 0
			ptr_to_uncov_vertices->delete_element(v);
			// num_of_uncov_vertices--;

		// dominating_set_size++;
		dominating_set_weight += vertices[v].get_weight();

		// ptr_to_partitioned_dominating_set->placeInVertexToCover(v, score);
		ptr_to_weight_bucket->placeInVertexToCover(v, vertices);

		//score[v] = -score[v]; in our algorithm, the score of the exchanged vertices are always 0
		//if(!score[v])
		ptr_to_no_loss_vertices->insert_element(v);

		LL* nbs = vertices[v].get_neighbors();
		LL dgr = vertices[v].get_degree();
		for(LL i = 0; i < dgr; i++)
		{
			LL n = nbs[i];
			connect_dominating_set_degree[n]++;
			if(ptr_to_dominating_set->element_in(n))
			{
				if(connect_dominating_set_degree[v] == 1)// v becomes strongly covered
				{
					// ptr_to_partitioned_dominating_set->lossMinusMinus(n, score);
					score[n]++;  // inside
#ifdef score_changed_mode
					unlock(v, n);
#endif
					if(!score[n])
#ifdef init_reduction_mode
						if(!is_fixed[n]) 
#endif
							ptr_to_no_loss_vertices->insert_element(n);
				}

				if(connect_dominating_set_degree[n] == 1)// n becomes strongly covered
				{
					// ptr_to_partitioned_dominating_set->lossMinusMinus(n, score);// one vertex may cover two vertices
					score[n]++;  // inside
#ifdef score_changed_mode
					unlock(v, n);
#endif
					if(!score[n])
#ifdef init_reduction_mode
						if(!is_fixed[n]) 
#endif
							ptr_to_no_loss_vertices->insert_element(n);
				}
			}
			else
			{
				if(connect_dominating_set_degree[v] == 0)// v becomes covered
				{
					// ptr_to_partitioned_dominating_set->gainMinusMinus(n, score);
					score[n]--;  // outside
#ifdef score_changed_mode
					unlock(v, n);
#endif
				}
				if(connect_dominating_set_degree[n] == 2)// n becomes strongly covered
				{
					LL* nbs_ = vertices[n].get_neighbors();
					LL dgr_ = vertices[n].get_degree();
					for(LL j = 0; j < dgr_; j++)
					{
						LL n_ = nbs_[j];
						if(n_ != v && ptr_to_dominating_set->element_in(n_))
						{
							// ptr_to_partitioned_dominating_set->lossMinusMinus(n_, score);// no longer causes loss of n
							score[n_]++;  // inside
#ifdef score_changed_mode
							unlock(v, n_);
#endif
							if(!score[n_])
#ifdef init_reduction_mode
								if(!is_fixed[n_]) 
#endif
								ptr_to_no_loss_vertices->insert_element(n_);
							break;
						}
					}							
				}
				else if(connect_dominating_set_degree[n] == 1)// n becomes critically covered
				{
					// ptr_to_partitioned_dominating_set->gainMinusMinus(n, score);//n itself no longer causes gain of n
					score[n]--;  // outside
#ifdef score_changed_mode
					unlock(v, n);
#endif
					ptr_to_uncov_vertices->delete_element(n);
					//num_of_uncov_vertices--;
					LL* nbs_ = vertices[n].get_neighbors();
					LL dgr_ = vertices[n].get_degree();
					for(LL j = 0; j < dgr_; j++)
					{
						LL n_ = nbs_[j];
						if(!ptr_to_dominating_set->element_in(n_))//n_ != v
						{
							// ptr_to_partitioned_dominating_set->gainMinusMinus(n_, score);//n_ no longer causes gain of n
							score[n_]--;  // outside
#ifdef score_changed_mode
							unlock(v, n_);
#endif
						}
					}
				}
#ifdef nvcc_mode
				confChange[n] = 1;
#endif
			}
		}
		time_stamp[v] = step;
#ifdef domination_hash_mode
		ptr_to_hashed_domination->update_hash_wrt_add(v);
#endif
#ifdef ans_update_op_mode
		ptr_to_moved_v->ans_update(v);
#endif

#if defined(two_level_cc_mode) || defined(two_level_cc_brute_force_mode)
		set_two_level_neighbors_free(v);
#endif

#ifdef cc2v3_mode
		cc2v3rule2(v);
#endif
//show_state();
		step++;

#ifdef debug3_mode
		if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
		{
			exit(1);
		}
#endif
	}

	void remove(const LL v)
	{

#ifdef debug3_mode
if(v == 583600)
{
cout << "to remove " << 583600 << endl;
getchar();
}
#endif

#if defined(debug_mode) || defined(debug3_mode)
//cout << "==================remove " << v << " loss: " << -score[v] << endl;

if(!ptr_to_dominating_set->element_in(v))
{
	cout << "cannot remove " << v << ", since it is not in the candidate solution" << endl;
	exit(1);
}

if(score[v])
{
	cout << "cannot remove " << v << ", since its score value is not 0" << endl;
	exit(1);
}

if(!vertices[v].get_degree())
{
	cout << "cannot remove " << v << ", its degree is 0. " << endl;
	exit(1);
}
#endif

		ptr_to_dominating_set->delete_element(v);
		ptr_to_complementary_set->insert_element(v);

#ifdef scenario_hash_mode
		ptr_to_hashed_scenario->update_hash_wrt_remove(v);
#endif

#if defined(nvcc_mode) || defined(score_changed_mode) || defined(two_level_cc_mode) || defined(two_level_cc_brute_force_mode)
		if(!confChange[v])
		{
			confChange[v] = true;
	#ifdef scenario_hash_mode
			ptr_to_hashed_scenario->update_hash_wrt_unlock(v);
	#endif
		}
#endif

		// in_dominating_set[v] = 0;
		if(!connect_dominating_set_degree[v])
			ptr_to_uncov_vertices->insert_element(v);
			//num_of_uncov_vertices++;

		// dominating_set_size--;
		dominating_set_weight -= vertices[v].get_weight();

		// ptr_to_partitioned_dominating_set->placeOutVertexFromCover(v, score);
		ptr_to_weight_bucket->placeOutVertexFromCover(v, vertices);

		//score[v] = -score[v];
		//if(!score[v])
		ptr_to_no_loss_vertices->delete_element(v);
// cout << 1 << endl;
		LL* nbs = vertices[v].get_neighbors();
		LL dgr = vertices[v].get_degree();

		for(LL i = 0; i < dgr; i++)
		{
// cout << 1.1 << endl;
			LL n = nbs[i];

			connect_dominating_set_degree[n]--;
			if(ptr_to_dominating_set->element_in(n))
			{
				if(connect_dominating_set_degree[v] == 1)// v becomes critically covered
				{
// cout << 1.4 << endl;
					// ptr_to_partitioned_dominating_set->lossPlusPlus(n, score);
					if(!score[n])
						ptr_to_no_loss_vertices->delete_element(n);
					score[n]--;  // inside
#ifdef score_changed_mode
					unlock(v, n);
#endif
// cout << 1.5 << endl;					
				}

				if(connect_dominating_set_degree[n] == 0)// n becomes critically covered
				{
// cout << 1.6 << endl;
					// ptr_to_partitioned_dominating_set->lossPlusPlus(n, score);
					if(!score[n])
						ptr_to_no_loss_vertices->delete_element(n);
					score[n]--;  // inside
#ifdef score_changed_mode
					unlock(v, n);
#endif
// cout << 2 << endl;
				}
			}
			else // if n is outside the current dominating set
			{
// cout << 2.1 << endl;
				if(connect_dominating_set_degree[v] == 0)// v becomes uncovered
				{
					// ptr_to_partitioned_dominating_set->gainPlusPlus(n, score);
					score[n]++;  // outside
#ifdef score_changed_mode
					unlock(v, n);
#endif
				}
// cout << 2.2 << endl;
				if(connect_dominating_set_degree[n] == 0)// n becomes uncovered
				{
// cout << 2.3 << endl;
					// ptr_to_partitioned_dominating_set->gainPlusPlus(n, score);
					score[n]++;  // outside
#ifdef score_changed_mode
					unlock(v, n);
#endif
					ptr_to_uncov_vertices->insert_element(n);
					// num_of_uncov_vertices++;
					LL* nbs_ = vertices[n].get_neighbors();
					LL dgr_ = vertices[n].get_degree();
					for(LL j = 0; j < dgr_; j++)
					{
						LL n_ = nbs_[j];
						if(!ptr_to_dominating_set->element_in(n_) && n_ != v)
						{
							// ptr_to_partitioned_dominating_set->gainPlusPlus(n_, score);// undirect effects
							score[n_]++;  // outside
#ifdef score_changed_mode
							unlock(v, n_);
#endif
						}
					}
				}
				else if(connect_dominating_set_degree[n] == 1)
				{
// cout << 2.4 << endl;
					LL* nbs_ = vertices[n].get_neighbors();
					LL dgr_ = vertices[n].get_degree();
					for(LL j = 0; j < dgr_; j++)
					{
						LL n_ = nbs_[j];
						if(ptr_to_dominating_set->element_in(n_))
						{
// cout << 2.9 << endl;
							// ptr_to_partitioned_dominating_set->lossPlusPlus(n_, score);
							if(!score[n_])
								ptr_to_no_loss_vertices->delete_element(n_);
							score[n_]--;  // inside
#ifdef score_changed_mode
							unlock(v, n_);
#endif
// cout << 3 << endl;
							break;
						}
					}					
				}
#ifdef nvcc_mode
				confChange[n] = 1;
#endif
			}
		}
// cout << 3.1 << endl;
		time_stamp[v] = step;
#ifdef domination_hash_mode
		ptr_to_hashed_domination->update_hash_wrt_remove(v);
#endif
#ifdef ans_update_op_mode
		ptr_to_moved_v->ans_update(v);
#endif

#if defined(two_level_cc_mode) || defined(two_level_cc_brute_force_mode)
		set_two_level_neighbors_free(v);
#endif

#ifdef cc2v3_mode
		cc2v3rule3(v);
#endif
//show_state();
// cout << 4 << endl;
		step++;

#ifdef debug3_mode
		if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
		{
			exit(1);
		}
#endif
	}

	void local_move()
	{

#ifdef debug3_mode
cout << "step&&&&&&: " << step << endl;
/*
if(!check_solution())
{
	cout << "the current best solution is incorrect" << endl;
	exit(1);
}
*/
if(!check_current_cand_solution())
		{
			cout << "current cand solution has not been maintained correctly" << endl;
			exit(1);
		}
cout << "dominating_set_size: " << ptr_to_dominating_set->size() << endl;
cout << "last_step_improved: " << last_step_improved << endl;
#endif
		
		LL best_remove_v;
		// exploitation
		do{

// cout << "no loss vertices: " << endl;
// ptr_to_no_loss_vertices->show_elements();
			best_remove_v = ptr_to_no_loss_vertices->best_element(confChange, vertices, time_stamp);
#ifdef debug3_mode
cout << "to remove " << best_remove_v << endl;
#endif
			if(best_remove_v == 0) break;
			remove(best_remove_v);
#ifdef detect_local_optimum_mode
			last_step_improved = 1;
#endif

#ifdef debug3_mode
cout << "having just removed a vertex" << endl;
#endif

#ifdef debug_mode
show_state();
#endif
		}while(best_remove_v != 0);
		

// cout << "dominating_set_weight: " << dominating_set_weight << endl;
// cout << "best_dominating_set_weight: " << best_dominating_set_weight << endl;

// cout << 1 << endl;
		if(last_step_improved) // 1-->0, increase-->decrease, local optimum
		{
#ifdef debug3_mode
cout << "local optimum detected" << endl;
#endif


#ifdef debug_mode
cout << "local optimum detected" << endl;
#endif
			if(dominating_set_weight < best_dominating_set_weight)
			{
#ifdef ans_update_op_mode
					ptr_to_moved_v->efficiently_update_best_dominating_set_vertices();
#else
					update_best_dominating_set();
#endif
					best_dominating_set_size = ptr_to_dominating_set->size();
					best_dominating_set_weight = dominating_set_weight;

					if(best_dominating_set_size == 1) return;

					stop = clock();
					best_cmp_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;

					best_solve_step = step;

#ifdef restart_mode
					restart_num_to_find_best_solution = restart_num;
#endif

#ifdef strategy_analysis_mode
					ans_update_times++;
#endif

#ifdef scenario_hash_mode
					ptr_to_hashed_scenario->forget_all_visits();
#endif
					return;			
			}

#ifdef scenario_hash_mode
			if(ptr_to_hashed_scenario->curr_hash_entry_visited())
			{
#ifdef debug3_mode
cout << "in the scenario revisiting scope" << endl;
#endif

#ifdef debug_mode
	#ifdef scenario_hash_mode
cout << "in a local optimum, scenario reoccurred" << endl;
getchar();
	#endif
#endif

#ifdef restart_mode
		
	#ifdef strategy_analysis_mode
		restart_num++;
	#endif

				clear();
				return;
#else 
	#ifdef random_walk_mode

				LL v = ptr_to_complementary_set->rand_element();
				add(v);
				return;
	#endif
#endif	
			
			}

			else
			{
				ptr_to_hashed_scenario->mark_hash_entry();
			}

#endif
		}// if(last_step_improved)

		

// cout << 2 << endl;
// ptr_to_weight_bucket->show_weight_bucket(max_weight);

		LL best_add_v;

#ifdef debug3_mode

#endif

		best_add_v = ptr_to_weight_bucket->vertex_outside_ds_of_smallest_weight_then_time_stamp(time_stamp);
// cout << 2.5 << endl;
		if(best_add_v != 0)
		{
#ifdef debug3_mode
cout << "in local_move(), right before adding one vertex, dominating_set_size: " << ptr_to_dominating_set->size() << endl;
LL element_in_count = 0;
for(LL i = 1; i <= v_num; i++)
if(ptr_to_dominating_set->element_in(i))
element_in_count++;
cout << "element_in_count: " << element_in_count << endl;
cout << "to insert " << best_add_v << endl;
if(ptr_to_dominating_set->element_in(best_add_v))
{
	cout << best_add_v << " has already been in the current dominating set" << endl;
	exit(1);
}
#endif
			add(best_add_v);
#ifdef debug3_mode
cout << "having just added " << best_add_v << endl;
#endif
#ifdef debug3_mode
cout << "in local_move(), right after adding one vertex, dominating_set_size: " << ptr_to_dominating_set->size() << endl;
element_in_count = 0;
for(LL i = 1; i <= v_num; i++)
if(ptr_to_dominating_set->element_in(i))
element_in_count++;
cout << "element_in_count: " << element_in_count << endl;
#endif
		}
		else
		{
#ifdef debug_mode
cout << "step: " << step << endl;
cout << "no vertices can change its state" << endl;
#endif
			clear();
			return;
		}
#ifdef detect_local_optimum_mode
		last_step_improved = false;
#endif
#if 0
		if(best_add_v != 0)
		{
			add(best_add_v);
		}
		else
		{
cout << 3 << endl;

ptr_to_weight_bucket->show_weight_bucket(max_weight);

			best_add_v = ptr_to_weight_bucket->vertex_outside_ds_of_smallest_weight_then_time_stamp(vertices, time_stamp); // may still fail to obtain any vertices, having to develop a suitable tabu strategy based on score change

			add(best_add_v);
		}
cout << 4 << endl;
#endif

#ifdef debug_mode
//show_state();
if(!check_solution())
{
	cout << "the current best solution is incorrect" << endl;
	exit(1);
}
if(!check_current_cand_solution())
		{
			cout << "current cand solution has not been maintained correctly" << endl;
			exit(1);
		}
#endif

#ifdef debug3_mode
cout << "step******: " << step << endl;
/*
if(!check_solution())
{
	cout << "the current best solution is incorrect" << endl;
	exit(1);
}
*/
if(!check_current_cand_solution())
		{
			cout << "current cand solution has not been maintained correctly" << endl;
			exit(1);
		}
getchar();
#endif
	}// void local_move()

	

public:

	~StochasticLocalSearch()
	{


		//delete[] in_dominating_set;
		delete[] connect_dominating_set_degree;

		delete ptr_to_dominating_set;
		delete ptr_to_complementary_set;
		delete ptr_to_uncov_vertices;

#ifdef init_reduction_mode
		delete ptr_to_live_vertices;
		delete[] is_fixed;
#endif

		delete[] time_stamp;

		delete[] score;

#if defined(nvcc_mode) || defined(score_changed_mode) || defined(two_level_cc_mode) || defined(two_level_cc_brute_force_mode)
		delete[] confChange;
#endif

#ifdef cc2v3_mode
		delete[] conf;
#endif

#ifdef unlocking_relation_mode
		delete[] unlocker;
#endif



		delete ptr_to_no_loss_vertices;

#ifdef domination_hash_mode
		delete ptr_to_hashed_domination;
#endif

#ifdef scenario_hash_mode
		delete ptr_to_hashed_scenario;
#endif



		delete[] best_in_dominating_set;

#ifdef ans_update_op_mode
		delete ptr_to_moved_v;
#endif


#if 1

	#ifdef init_reduction_mode
/*
		if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
		{
			cout << "something is wrong with the weight bucket" << endl;
			exit(1);
		}
*/
	#endif
		//ptr_to_weight_bucket->show_weight_bucket();

		delete ptr_to_weight_bucket;
#endif

#ifdef init_reduction_mode
		delete ptr_to_degree_0_vertices;
		delete ptr_to_degree_1_vertices;
		delete ptr_to_degree_2_vertices;
#endif


	}

	void cover_sls()
	{
		while(1)
		{
			LL step_bound = 1;
			if(LL(step / 1000) > step_bound)
			{
					step_bound++;

					stop = clock();
					double elap_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;

					if(LL(elap_time) >= time_limit)
					{

#ifdef ans_update_op_mode
						bool* in_dominating_set = new bool[v_num + 1];
						for(LL v = 1; v <= v_num; v++)
						{
							in_dominating_set[v] = ptr_to_dominating_set->element_in(v);
						}
						ptr_to_moved_v->compute_final_answer(v_num, in_dominating_set, best_in_dominating_set);
						delete[] in_dominating_set;
#endif

						return;
					}
			}

			local_move();

			if(best_dominating_set_size == 1)
			{
#ifdef ans_update_op_mode
						bool* in_dominating_set = new bool[v_num + 1];
						for(LL v = 1; v <= v_num; v++)
						{
							in_dominating_set[v] = ptr_to_dominating_set->element_in(v);
						}
						ptr_to_moved_v->compute_final_answer(v_num, in_dominating_set, best_in_dominating_set);
						delete[] in_dominating_set;
#endif
				return;
			}
		}
	}

	void show_results()
	{

		stop = clock();
		double elap_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;
#ifdef individual_analysis_on_init_sls_mode
		sls_time = elap_time - init_time;
#endif
		if(check_solution())
		{
			cout << "o " << best_dominating_set_weight << endl;
			cout << "c size " << best_dominating_set_size << endl;
			cout << "c solveTime(ms) " << best_cmp_time << endl;
			cout << "c searchSteps " << best_solve_step << endl;

#ifdef individual_analysis_on_init_sls_mode
			cout << "c init_time " << init_time << endl;
			cout << "c sls_time " << sls_time << endl;
			cout << "c stepSpeed(/ms) " << step / sls_time << endl;
#endif
			cout << "c seed " << seed << endl;
			cout << "c first_solution_weight " << first_solution_weight << endl;
#ifdef init_reduction_mode
			cout << "c non_fixed_vertex_num " << ptr_to_live_vertices->size() << endl;

			cout << "c fixed_rate " << 1 - round((float(ptr_to_live_vertices->size()) / float(v_num)) * 100) / 100.0 << endl;
#endif
			cout << "c restart_num_to_find_best_solution " << restart_num_to_find_best_solution << endl;
			cout << "c user cutoff(ms) " << time_limit << endl;
			cout << "c init_time(ms) " << init_time << endl;
			cout << "c actual_cutoff " << elap_time << endl;
			cout << "c actual_step " << step << endl;
			cout << "c v_num " << v_num << endl;
			cout << "c e_num " << e_num << endl;
			cout << "c avg_degree " << round(2 * float(e_num) * 100 / float(v_num)) / 100.0 << endl;
			cout << "c solver_name " << solver_name << endl;
			cout << "c file_path " << file_name << endl;
#ifdef restart_mode
			cout << "c restart_num " << restart_num << endl;
			cout << "c restart_period " << step / restart_num << endl;
#endif
#if 0
			cout << "c dominating set: ";
			for(int v = 1; v <= v_num; v++) if(best_in_dominating_set[v]) cout << v << '\t'; 
			cout << endl;
#endif

#ifdef strategy_analysis_mode
cout << "c the total number of steps: " << step << endl;
cout << "c answer update times: " << ans_update_times << endl;
#endif
		}
		else
		{
			cout << "sorry, something is wrong" << endl;
		}
	}	

	void show_state()
	{ 
		if(!check_solution())
		{
			cout << "the current best solution is incorrect" << endl;
			exit(1);
		}

		if(!check_current_cand_solution())
		{
			cout << "current cand solution has not been maintained correctly" << endl;
			exit(1);
		}


		LL v, e;
/*
		for(e = 0; e < e_num; e++)
		{
			cout << edges[e].get_v1() << " " << edges[e].get_v2() << endl;
		}

		for(v = 1; v <= v_num; v++)
		{
			cout << "vertex " << v << ": " << endl;
			vertices[v].show_neighbors();
		}
*/
		LL i, j;
		//cout << "step: " << step << endl;
#if 0
		cout << "the dominating set: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			if(in_dominating_set[v])
				cout << v << '\t';
		}
		cout << endl;
#endif

		cout << step << endl;
#if 0
		cout << "current dominating set: " << endl;
		for(v = 1; v <= v_num; v++)
			if(ptr_to_dominating_set->element_in(v))
				cout << v << "\t";
		cout << endl;
		cout << "loss: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			if(ptr_to_dominating_set->element_in(v))
				cout << -score[v] << '\t';
		}
		cout << endl;
		cout << "best dominating set: " << endl;
		for(v = 1; v <= v_num; v++)
			if(best_in_dominating_set[v])
				cout << v << "\t";
		cout << endl;

		cout << "dominating set weight: " << dominating_set_weight << endl;
		cout << "best dominating set weight: " << best_dominating_set_weight << endl;
#endif


#if 0
		cout << "connect dominating set degrees: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << v << '\t';
		}
		cout << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << connect_dominating_set_degree[v] << '\t';
		}
		cout << endl;
#endif

		cout << "vertices: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << v << '\t';
		}
		cout << endl;

		cout << "whether in dominating set: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << ptr_to_dominating_set->element_in(v) << '\t';
		}
		cout << endl;

		cout << "score: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << score[v] << '\t';
		}
		cout << endl;

		cout << "confChange: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << confChange[v] << '\t';
		}
		cout << endl;

#ifdef unlocking_relation_mode
		cout << "unlocker: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << unlocker[v] << '\t';
		}
		cout << endl;
#endif

#if defined(scenario_hash_mode)		
		cout << "last_unlock_edge: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			if(unlocker[v] != 0)
				cout << unordered_pair_index_of(v, unlocker[v]) << '\t';
			else
				cout << -1 << '\t';
		}
#endif

		cout << endl;
		cout << "weight: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << vertices[v].get_weight() << '\t';
		}
		cout << endl;

#ifdef scenario_hash_mode
	#ifdef test_scenario_hash_mode
		LL curr_hash_entry = 0;
		for(v = 1; v <= v_num; v++)
		{
			if(ptr_to_dominating_set->element_in(v))
				curr_hash_entry = (curr_hash_entry + ptr_to_hashed_scenario->get_vertex_inclusion_mod_array()[v]) % PRIME_NUM;
		}

		for(v = 1; v <= v_num; v++)
		{
			if(confChange[v] == 1)
				curr_hash_entry = (curr_hash_entry + ptr_to_hashed_scenario->get_vertex_freedom_mod_array()[v]) % PRIME_NUM;
		}

		for(e = 0; e < within_two_dist_pair_num; e++)
		{
	//cout << "update hash wrt. unlocking relation update" << endl;
			LL v1, v2;
			within_two_dist_edges[e].get_vertices(v1, v2);

			if(unlocker[v1] == v2)
			{
				if(v1 > v2)
					curr_hash_entry = (curr_hash_entry + ptr_to_hashed_scenario->get_small_id_to_big_unlock_relation_mod_array()[e]) % PRIME_NUM;
				else
					curr_hash_entry = (curr_hash_entry + ptr_to_hashed_scenario->get_big_id_to_small_unlock_relation_mod_array()[e]) % PRIME_NUM;
			}

			if(unlocker[v2] == v1)
			{
				if(v1 > v2)
					curr_hash_entry = (curr_hash_entry + ptr_to_hashed_scenario->get_big_id_to_small_unlock_relation_mod_array()[e]) % PRIME_NUM;
				else
					curr_hash_entry = (curr_hash_entry + ptr_to_hashed_scenario->get_small_id_to_big_unlock_relation_mod_array()[e]) % PRIME_NUM;
			}

		}
		cout << "expected curr_hash_value: \t";
		cout << curr_hash_entry << endl;
		cout << "actual curr_hash_value: \t";
		cout << ptr_to_hashed_scenario->curr_hash_entry_val() << endl;
		if(curr_hash_entry != ptr_to_hashed_scenario->curr_hash_entry_val())
		{
			cout << "the hash value is computed incorrectly" << endl;
			exit(1);
		}
	#endif
#endif
		// cout << "max_gain: " << ptr_to_partitioned_dominating_set->get_max_gain() << endl;

#if 0
		cout << "current best dominating set: " << endl;
		for(v = 1; v <= v_num; v++)
			if(best_in_dominating_set[v])
				cout << v << "\t";
		cout << endl;
		cout << "best_dominating_set_weight: " << best_dominating_set_weight << endl;
#endif

		//ptr_to_weight_bucket->show_weight_bucket(max_weight);
#ifdef init_reduction_mode
		if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
		{
			exit(1);
		}
#endif
#if 0
		cout << "no loss vertices:" << endl;
		ptr_to_no_loss_vertices->show_elements_together_with_confChange(confChange);

		cout << "answer changed vertices: " << endl;
		ptr_to_moved_v->show_elements();
#endif

/*
		cout << "confChange: " << endl;
		for(i = ptr_to_uncov_vertices->begin(); i < ptr_to_uncov_vertices->end(); i++)
		{
			int v = ptr_to_uncov_vertices->at(i);
			cout << confChange[v] << '\t';
		}
		cout << endl;
*/
#ifdef strategy_analysis_mode
	#ifdef restart_mode
		// cout << "restart_num: " << restart_num << endl;
	#endif
#endif
/*
		cout << "uncov edges: " << endl;
		for(i = ptr_to_uncov_edges->begin(); i < ptr_to_uncov_edges->end(); i++)
		{
			int e = ptr_to_uncov_edges->at(i);
			cout << e;
			int v1, v2;
			edges[e].get_vertices(v1, v2);
			cout << ", " << "endpoints: " << v1 << " and " << v2 << endl;
		}
		cout << endl;
*/

cout << "step: " << step << endl;
	
cout << "**************************************" << endl;

//getchar();
	}
};
#endif
