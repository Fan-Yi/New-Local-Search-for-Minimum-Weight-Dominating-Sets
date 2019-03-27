#ifndef OPERAND_SETS_H
#define OPERAND_SETS_H

#include "hugeInt.h"

class ConfchangedWeightAgePickSet : public Bijection
{
public:	
	ConfchangedWeightAgePickSet(LL sz) : Bijection(sz)
	{
	}

	~ConfchangedWeightAgePickSet()
	{
	}

	LL best_element(bool* confChange, Vertex* vertices, HugeInt* time_stamp)
	{
		if(!array_size) return 0;
		LL best_v = 0;
		LL i;
		for(i = begin(); i != end(); i++)
		{
			LL v = array[i]; 

			if(confChange[v])
			{
				best_v = v;

				i++;
				break;
			}
		}
		for(; i != end(); i++)
		{
			LL v = array[i];
			if(!confChange[v]) continue;
			if(vertices[v].get_weight() > vertices[best_v].get_weight() || (vertices[v].get_weight() == vertices[best_v].get_weight() && time_stamp[v] < time_stamp[best_v]))
				best_v = v;
		}
		return best_v;
	}

};
#endif
