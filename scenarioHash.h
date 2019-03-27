#ifndef STATE_HASH_H
#define STATE_HASH_H

#include <bitset>
#include "constants.h"

//#define debug_mode

typedef long long LL;

class ScenarioHash{
private:
	LL* vertex_inclusion_mod_array;
	LL* vertex_freedom_mod_array;
	LL* big_id_to_small_unlock_relation_mod_array;
	LL* small_id_to_big_unlock_relation_mod_array;
	bitset<PRIME_NUM> hash_entries;
	LL curr_hash_entry;
	LL init_hash_entry;

public:
	ScenarioHash(LL v_num, LL e_num)
	{
		// the candidate solution
		vertex_inclusion_mod_array = new LL[v_num + 2];
		vertex_inclusion_mod_array[0] = 1;
		for(LL v = 1; v <= v_num; v++)
		{
			vertex_inclusion_mod_array[v] = (vertex_inclusion_mod_array[v - 1] * 2) % PRIME_NUM;
		}
		// the freedom status
		vertex_freedom_mod_array = new LL[v_num + 2];
		vertex_freedom_mod_array[0] = vertex_inclusion_mod_array[v_num];
		for(LL v = 1; v <= v_num; v++)
		{
			vertex_freedom_mod_array[v] = (vertex_freedom_mod_array[v - 1] * 2) % PRIME_NUM;
		}
		// the unlock relations (bigger unlocks smaller)
		big_id_to_small_unlock_relation_mod_array = new LL[e_num + 2];
		big_id_to_small_unlock_relation_mod_array[0] = (vertex_freedom_mod_array[v_num] * 2) % PRIME_NUM;
		for(LL e = 1; e < e_num; e++)
		{
			big_id_to_small_unlock_relation_mod_array[e] = (big_id_to_small_unlock_relation_mod_array[e - 1] * 2) % PRIME_NUM;
		}
		// the unlock relations (smaller unlocks bigger)
		small_id_to_big_unlock_relation_mod_array = new LL[e_num + 2];
		small_id_to_big_unlock_relation_mod_array[0] = (big_id_to_small_unlock_relation_mod_array[e_num - 1] * 2) % PRIME_NUM;
		for(LL e = 1; e < e_num; e++)
		{
			small_id_to_big_unlock_relation_mod_array[e] = (small_id_to_big_unlock_relation_mod_array[e - 1] * 2) % PRIME_NUM;
		}


//////////////////////////////////////////////////////////////////////////////////////
		curr_hash_entry = 0;
		// the candidate solution
		for(LL v = 1; v <= v_num; v++)
		{
			curr_hash_entry = (curr_hash_entry + vertex_inclusion_mod_array[v]) % PRIME_NUM; // all vertices are in
		}

		// the freedom status
		for(LL v = 1; v <= v_num; v++)
		{
			curr_hash_entry = (curr_hash_entry + vertex_freedom_mod_array[v]) % PRIME_NUM;
		}
#ifdef debug_mode
cout << "initial hash value: " << curr_hash_entry << endl;
#endif
		// the unlocking relation
		// nothing needs to be done
		init_hash_entry = curr_hash_entry;
	}
/*
	void clear()
	{
		hash_entries.reset();
		curr_hash_entry = init_hash_entry;
		// the candidate solution
		// nothing needs to be done
	}
*/
	void forget_all_visits()
	{
		hash_entries.reset();
	}

	~ScenarioHash()
	{
		delete[] vertex_inclusion_mod_array;
		delete[] vertex_freedom_mod_array;
		delete[] big_id_to_small_unlock_relation_mod_array;
		delete[] small_id_to_big_unlock_relation_mod_array;
	}

	void update_hash_wrt_add(const LL v)
	{
		curr_hash_entry = (curr_hash_entry + vertex_inclusion_mod_array[v]) % PRIME_NUM;
#ifdef debug_mode
//cout << "hash value (add): add into " << vertex_inclusion_mod_array[v] << endl;
//cout << "\tbecomes " << curr_hash_entry << endl;
#endif
	}

	void update_hash_wrt_remove(const LL v)
	{
		curr_hash_entry = (curr_hash_entry + PRIME_NUM - vertex_inclusion_mod_array[v]) % PRIME_NUM;
#ifdef debug_mode
//cout << "hash value (drop): del from " << vertex_inclusion_mod_array[v] << endl;
//cout << "\tbecomes " << curr_hash_entry << endl;
#endif
	}
///////////////////////////////////////////////////
	void update_hash_wrt_unlock(const LL v)
	{
		curr_hash_entry = (curr_hash_entry + vertex_freedom_mod_array[v]) % PRIME_NUM;
#ifdef debug_mode
//cout << "hash value (unlock): add into " << vertex_freedom_mod_array[v] << endl;
//cout << "\tbecomes " << curr_hash_entry << endl;
#endif
	}

	void update_hash_wrt_lock(const LL v)
	{
		curr_hash_entry = (curr_hash_entry + PRIME_NUM - vertex_freedom_mod_array[v]) % PRIME_NUM;
#ifdef debug_mode
//cout << "hash value (lock): del from " << vertex_freedom_mod_array[v] << endl;
//cout << "\tbecomes " << curr_hash_entry << endl;
#endif
	}	
//////////////////////////////////////////////////
	void update_hash_wrt_insert_big_id_to_small_unlock_edge(const LL e)
	{
		curr_hash_entry = (curr_hash_entry + big_id_to_small_unlock_relation_mod_array[e]) % PRIME_NUM;
#ifdef debug_mode
//cout << "hash value (insert): add into " << big_id_to_small_unlock_relation_mod_array[e] << endl;
//cout << "\tbecomes " << curr_hash_entry << endl;
#endif
	}

	void update_hash_wrt_delete_big_id_to_small_unlock_edge(const LL e)
	{
		curr_hash_entry = (curr_hash_entry + PRIME_NUM - big_id_to_small_unlock_relation_mod_array[e]) % PRIME_NUM;
#ifdef debug_mode
//cout << "hash value (delete): del from " << big_id_to_small_unlock_relation_mod_array[e] << endl;
//cout << "\tbecomes " << curr_hash_entry << endl;
#endif
	}
///////////////////////////////////////////////////
	void update_hash_wrt_insert_small_id_to_big_unlock_edge(const LL e)
	{
		curr_hash_entry = (curr_hash_entry + small_id_to_big_unlock_relation_mod_array[e]) % PRIME_NUM;
#ifdef debug_mode
//cout << "hash value (insert): add into " << small_id_to_big_unlock_relation_mod_array[e] << endl;
//cout << "\tbecomes " << curr_hash_entry << endl;
#endif
	}

	void update_hash_wrt_delete_small_id_to_big_unlock_edge(const LL e)
	{
		curr_hash_entry = (curr_hash_entry + PRIME_NUM - small_id_to_big_unlock_relation_mod_array[e]) % PRIME_NUM;
#ifdef debug_mode
//cout << "hash value (delete): del from " << small_id_to_big_unlock_relation_mod_array[e] << endl;
//cout << "\tbecomes " << curr_hash_entry << endl;
#endif
	}
//////////////////////////////////////////////////
	void mark_hash_entry()
	{
		hash_entries.set(curr_hash_entry);
	}
/////////////////////////////////////////////////

	bool curr_hash_entry_visited()
	{
		return hash_entries.test(curr_hash_entry);
	}

	LL curr_hash_entry_val()
	{
		return curr_hash_entry;
	}


////////////////////////////////////////////////

	LL *get_vertex_inclusion_mod_array()
	{
		return vertex_inclusion_mod_array;
	}

	LL *get_vertex_freedom_mod_array()
	{
		return vertex_freedom_mod_array;
	}

	LL *get_big_id_to_small_unlock_relation_mod_array()
	{
		return big_id_to_small_unlock_relation_mod_array;
	}

	LL *get_small_id_to_big_unlock_relation_mod_array()
	{
		return small_id_to_big_unlock_relation_mod_array;
	}

};

#endif
