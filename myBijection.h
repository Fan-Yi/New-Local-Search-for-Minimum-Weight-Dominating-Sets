#ifndef _MY_BIJECTION_H
#define _MY_BIJECTION_H

#include <cstring>

class Bijection
{
protected:
	LL *array;
	LL *index_in_array;
	LL array_size;
	LL array_capacity;

public:
	Bijection()
	{
		array = NULL;
		index_in_array = NULL;
		array_size = 0;
		array_capacity = 0;
	}

	Bijection(LL sz)
	{
		array_capacity = sz + 2;
		array = new LL[array_capacity];
		index_in_array = new LL[array_capacity];
		memset(array, 0, sizeof(LL) * (array_capacity));
		memset(index_in_array, 0, sizeof(LL) * (array_capacity));
		array_size = 0;
	}

	~Bijection()
	{
		delete[] array;
		array = NULL;
		delete[] index_in_array;
		index_in_array = NULL;
		array_size = 0;
		array_capacity = 0;
	}

	void clear()
	{
/*
		for(LL i = 1; i <= array_size; i++)
		{
			index_in_array[array[i]] = 0;
		}
*/
		memset(array, 0, sizeof(LL) * (array_capacity));
		memset(index_in_array, 0, sizeof(LL) * (array_capacity));
		array_size = 0;	
	}

	void insert_element(const LL e)
	{
#ifdef debug3_mode
if(element_in(e))
{
	cout << "Element " << e << " does exist and should not be inserted" << endl;
	exit(1);
}
#endif
//cout << 2.1 << endl;
//cout << "edge: " << e << endl;
//cout << "array_size: " << array_size << endl;
		++array_size;
//cout << "array_size_2: " << array_size << endl;
//cout << "*****" << endl;
		array[array_size] = e;
//cout << 2.2 << endl;

		index_in_array[e] = array_size;
//cout << 2.3 << endl;
	}

	void delete_element(const LL e)
	{
#ifdef debug_mode
if(!element_in(e))
{
	cout << "Element " << e << " does not exist and should not be deleted" << endl;
	exit(1);
}
#endif
		LL last_e = array[array_size--];
		LL e_pos = index_in_array[e];
		array[e_pos] = last_e;
		index_in_array[last_e] = e_pos;
		index_in_array[e] = 0;
	}

	LL index_of(const LL e)
	{
		return index_in_array[e];
	}

	bool element_in(const LL e)
	{
		return index_in_array[e];
	}

	LL size()
	{
		return array_size;
	}

	LL at(LL i)
	{
		return array[i];
	}

	LL begin()
	{
		return 1;
	}

	LL end()
	{
		return array_size + 1;
	}

	LL rand_element()
	{
		return array[rand() % array_size + 1];
	}

	void show_elements()
	{
		cout << "show elements: " << endl;
		for(LL i = 1; i <= array_size; i++)
		{
			cout << array[i] << "\t";
		}
		cout << endl;
	}

	void show_elements_together_with_confChange(bool* free)
	{
		cout << "show elements: " << endl;
		for(LL i = 1; i <= array_size; i++)
		{
			cout << array[i] << "\t";
		}
		cout << endl;
		cout << "confChange: " << endl;
		for(LL i = 1; i <= array_size; i++)
		{
			cout << free[array[i]] << "\t";
		}
		cout << endl;
/*
		cout << "weight: " << endl;
		for(LL i = 1; i <= array_size; i++)
		{
			cout << vtces[array[i]].get_weight() << "\t";
		}
		cout << endl;
*/
	}
};
#endif
