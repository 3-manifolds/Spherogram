/**************************************************
                  util.h

This header file is intended to hold common utilities not
related to any particular template class

**************************************************/

template <class T> T gcd (T i, T j)
{
    if (i == T(0))
		return (j);
    else
		if (j == T(0))
		    return (i);
		else
	    	return ( gcd<T>(j, i%j));
}

template <class T> inline void get_number (T& num, char* cptr)
{
	istringstream ss(cptr);
	ss >> num;
}

template <class T> inline void get_number (T& num, string s, string::size_type pos)
{
	istringstream ss(s.substr(pos));
	ss >> num;
}

/* Other utility functions defined in util.cpp */
int num_len (long n);
char* c_string(const string& s);
void print_state(istream& is);

// This function is used only for debugging .h files
void breakpoint();
