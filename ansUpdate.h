#ifndef ANS_UPDATE_H
#define ANS_UPDATE_H

#include "myBijection.h"

class AnsChangeSet : public Bijection
{
public:
	AnsChangeSet(LL v_num) : Bijection(v_num)
	{
	}

	~AnsChangeSet()
	{
	}

	void ans_update(LL v)
	{
		if(element_in(v))
		{
			delete_element(v);
		}
		else
		{
			insert_element(v);
		}
	}

	void efficiently_update_best_dominating_set_vertices()
	{
		clear();
	}

	void compute_final_answer(LL v_num, bool* in_cover, bool* best_in_cover)
	{
		for(LL v = 1; v <= v_num; v++)
		{
			if(element_in(v))
			{
				best_in_cover[v] = 1 - in_cover[v];
			}
			else
			{
				best_in_cover[v] = in_cover[v];
			}
		}		
	}
};
#endif
