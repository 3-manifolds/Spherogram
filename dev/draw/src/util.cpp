/**************************************************
                  util.cpp

This file is intended to hold common utilities 

**************************************************/
using namespace std;

#include <iostream>
#include <ios>
#include <string>

int num_len (long n)
{
    int i = 0;

    if (n<0)
    {
        i++;
        n *= -1;
    }

    do
    {
        i++;
        n /= 10;
    } while(n != 0);

    return i;
}

char* c_string(const string& s)
{
	char* p = new char [s.length()+1];
	p[s.copy(p,string::npos)] = 0; // add terminator
	return p;
}

void print_state(istream& is)
{
cout << "\nprint_state:" << endl;
	ios_base::iostate s = is.rdstate();
cout << "\niostate: " << s << endl;
	if (s&ios_base::badbit) cout << "\nbad bit set";
	if (s&ios_base::failbit) cout << "\nfail bit set";
	if (s&ios_base::eofbit) cout << "\neof bit set";
	if (s&ios_base::eofbit) cout << "\ngood bit set";
	cout << endl;
}

//This function is used only for debugging .h files
void breakpoint() {}
