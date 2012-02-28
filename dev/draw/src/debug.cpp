/****************************************************************************
 debug.cpp defines the debug ofstream and provides the functions 
 to set debug booleans from the command line.
 ****************************************************************************/
using namespace std;


#include <iostream>
#include <fstream>
#include <cstring>

ofstream debug;

/* debug.h should be provided by the main programme to describe what debug
   capabilities the debug subsystem should provide.  These are specified in
   terms of the following preprocessor ditectives, which should be #define
   if required by the main programme
   
   DEBUG_BIGINT
   DEBUG_BIGREAL
   DEBUG_MATRIX
   DEBUG_POLYNOMIAL
   DEBUG_QUATERNION
   DEBUG_RATIONAL
   DEBUG_SURD

   
*/
#include <debug.h>


#include <util.h>
#include <algorithm>

//#ifdef DEBUG_SCALAR
//	#include <scalar.h>
//#endif


#ifdef DEBUG_BIGINT
	#include <bigint.h>
#endif

#ifdef DEBUG_BIGREAL
	#include <bigreal_control.h>
#endif


#ifdef DEBUG_POLYNOMIAL
	#include <polynomial.h>
#endif

#ifdef DEBUG_MATRIX	
	#include <matrix.h>
#endif

#ifdef DEBUG_QUATERNION
	#include <quaternion.h>
#endif

#ifdef DEBUG_RATIONAL
	#include <rational.h>
#endif

#ifdef DEBUG_SURD
	#include <surd_control.h>
#endif

void set_debug_option(char* start, char* end);
void check_debug_option_parameters(char* pptr, string option);
void set_debug_option_parameter(char* pptr, string option);

/* These functions have to be provided by the calling programme */
void set_main_debug_option(char* start, char* end);
void set_main_debug_option_parameter(char* pptr, string option);
void main_display_default_options();

/* The parameter to debug_setup is the argv[] pointer that contains the debug option
   The format of the debug option is #[option:option{parameter:parameter}:option] etc.
   and this function removes the [...] from the argument to allow the calling programme 
   to process other options independently of the debug options.
   
   The function returns true if [...] is present, and false otherwise, this allows the default 
   debug options to be set by the calling programme
*/
bool debug_setup(char* argv_ptr)
{
	char* cptr = strchr(argv_ptr,'#');
	
	if (*++cptr == '[')
	{
		char* debug_options_start = cptr; // used to remove the [...] later

		/* make sure there's an end to the option definitions*/
		if (!strchr(cptr,']'))
		{
			cout << "\nError in command line: debug options must be enclosed within []\n";
			exit(0);
		}

		do
		{
			cptr++; // step over '[', or ','

			if (!isalpha(*cptr))
			{
				cout << "\nError in command line: debug options must start with a keyword\n";
				exit(0);
			}

			char* end = cptr;

			while (*end != ':' && *end != ']')
			{						
				if (*end == '{')
				{					
					end++;
					/* move to the end of the option parameters */
					while (*end != '}')
					{
						 if (*end == ']')
						 {
							 cout << "\nError in command line: debug option parameters must be enclosed in {}\n";
							 exit(0);
						 }
						 else
							end++;									 
					}
					end++; // step over '}'
				}
				else
					end++;
			}

			if (end>cptr) // deals with the degenerate case -#[]
				set_debug_option(cptr, end-1);
			cptr = end;

		} while (*cptr != ']');

		/* Remove the debug options from the options string for further processing */
		strcpy(debug_options_start,cptr+1);
		
		return true;
	}
	else
		return false;
}

void set_debug_option(char* start, char* end)
{
	char  loc_buf[end-start+2];
	char* c1 = start;
	char* c2 = loc_buf;

	/* if both start and end are zero, display debug help information.
	   this has been included here in this manner so that each time a debug option
	   is added, the help will be updated (hopefully!)
	*/
	if (start == 0 && end == 0)
	{
#ifdef DEBUG_BIGINT
		cout << "\t\tbigint{san:+:-:*:/:%:rs:ls:bc:==:gt:out:in:sum:diff:num_len:gcd:all}, bitmap: no default" << endl;
#endif
#ifdef DEBUG_BIGREAL
		cout << "\t\tbigreal{+:-:*:/:==:gt:out:in:sum:diff:num_len:carry:all}, bitmap: no default" << endl;
#endif
#ifdef DEBUG_MATRIX
		cout << "\t\tmatrix{det:inv:*:all}, bitmap: no default" << endl;
#endif
#ifdef DEBUG_POLYNOMIAL
		cout << "\t\tpoly{no-gen:san:+:*:/:in:gcd:all}, bitmap: general debug included unless no-gen specified" << endl;
#endif
#ifdef DEBUG_QUATERNION
		cout << "\t\tquaternion, boolean" << endl;
#endif
#ifdef DEBUG_RATIONAL
		cout << "\t\trational, boolean" << endl;
#endif
#ifdef DEBUG_SURD
		cout << "\t\tsurd{+:-:*:==:gt:in:newton:all}, bitmap: no default" << endl;
#endif
		return;
	}

	do
	{
		*c2++ = *c1++;
	} while (c1 <= end);
	
	*c2 = '\0';

	char* pptr = strchr(loc_buf,'{');
	if (pptr)
		*pptr++ = '\0';

	/* now, even if there are parameters, the first part of loc_buf 
	   is a C-string that identifies the option */
	
	if (!strcmp(loc_buf,"bigint"))
	{
#ifdef DEBUG_BIGINT
		if (pptr)
		{
			/* check for any parameters */
			check_debug_option_parameters(pptr, loc_buf);			
		}
#endif
	}
	else if (!strcmp(loc_buf,"bigreal"))
	{
#ifdef DEBUG_BIGREAL
		if (pptr)
		{
			/* check for any parameters */
			check_debug_option_parameters(pptr, loc_buf);			
		}
#endif
	}
	else if (!strcmp(loc_buf,"matrix"))
	{
#ifdef DEBUG_MATRIX
		if (pptr)
		{
			/* check for any parameters */
			check_debug_option_parameters(pptr, loc_buf);			
		}
#endif
	}
	else if (!strcmp(loc_buf,"poly"))
	{
#ifdef DEBUG_POLYNOMIAL
		polynomial_control::DEBUG |= polynomial_control::general;
		debug << "debug: setting debug option polynomial_control::general\n";		

		if (pptr)
		{
			/* check for any parameters */
			check_debug_option_parameters(pptr, loc_buf);			
		}
#endif
	}
	else if (!strcmp(loc_buf,"rational"))
	{
#ifdef DEBUG_RATIONAL
		rational_control::DEBUG = true;
		debug << "debug: setting debug option rational_control::DEBUG\n";		

#endif
	}
	else if (!strcmp(loc_buf,"quaternion"))
	{
#ifdef DEBUG_QUATERNION
		quaternion_control::DEBUG = true;
		debug << "debug: setting debug option quaternion_control::DEBUG\n";		
#endif
	}
	else if (!strcmp(loc_buf,"surd"))
	{
#ifdef DEBUG_SURD
		if (pptr)
		{
			/* check for any parameters */
			check_debug_option_parameters(pptr, loc_buf);			
		}
#endif
	}	else
		set_main_debug_option(start,end); // function has to be provided by calling programme


	debug.flush();
}

/* work through the parameters from pptr and call set_debug_parameter_option for each one */
void check_debug_option_parameters(char* pptr, string option)
{
	char* cptr = pptr+1;
	bool end_of_parameters = false;

	do
	{
		if (*cptr == ':' || *cptr == '}')
		{
		
			if (*cptr == '}')
				end_of_parameters = true;
				
			*cptr = '\0';
			set_debug_option_parameter(pptr, option);

			if (end_of_parameters)
				break;
			else
				pptr = ++cptr;
		}
		else
			cptr++;

	} while(true); // we break out of this loop
}

void set_debug_option_parameter(char *pptr, string option)
{
	if (option == "bigint")
	{
#ifdef DEBUG_BIGINT
		if (!strcmp(pptr,"san"))
		{
			bigint_control::DEBUG |= bigint_control::sanitize;
			debug << "debug: setting debug option bigint_control::sanitize\n";		
		}
		else if (!strcmp(pptr,"+"))
		{
			bigint_control::DEBUG |= bigint_control::add;
			debug << "debug: setting debug option bigint_control::add\n";		
		}
		else if (!strcmp(pptr,"-"))
		{
			bigint_control::DEBUG |= bigint_control::subtract;
			debug << "debug: setting debug option bigint_control::subtract\n";		
		}
		else if (!strcmp(pptr,"*"))
		{
			bigint_control::DEBUG |= bigint_control::multiply;
			debug << "debug: setting debug option bigint_control::multiply\n";		
		}
		else if (!strcmp(pptr,"/"))
		{
			bigint_control::DEBUG |= bigint_control::divide;
			debug << "debug: setting debug option bigint_control::divide\n";		
		}
		else if (!strcmp(pptr,"%"))
		{
			bigint_control::DEBUG |= bigint_control::remainder;
			debug << "debug: setting debug option bigint_control::remainder\n";		
		}
		else if (!strcmp(pptr,"rs"))
		{
			bigint_control::DEBUG |= bigint_control::r_shift;
			debug << "debug: setting debug option bigint_control::r_shift\n";		
		}
		else if (!strcmp(pptr,"ls"))
		{
			bigint_control::DEBUG |= bigint_control::l_shift;
			debug << "debug: setting debug option bigint_control::l_shift\n";		
		}
		else if (!strcmp(pptr,"bc"))
		{
			bigint_control::DEBUG |= bigint_control::bool_conv;
			debug << "debug: setting debug option bigint_control::bool_conv\n";		
		}
		else if (!strcmp(pptr,"=="))
		{
			bigint_control::DEBUG |= bigint_control::equal;
			debug << "debug: setting debug option bigint_control::equal\n";		
		}
		else if (!strcmp(pptr,"gt"))
		{
			bigint_control::DEBUG |= bigint_control::greater;
			debug << "debug: setting debug option bigint_control::greater\n";		
		}
		else if (!strcmp(pptr,"out"))
		{
			bigint_control::DEBUG |= bigint_control::output;
			debug << "debug: setting debug option bigint_control::output\n";		
		}
		else if (!strcmp(pptr,"in"))
		{
			bigint_control::DEBUG |= bigint_control::input;
			debug << "debug: setting debug option bigint_control::input\n";		
		}
		else if (!strcmp(pptr,"sum"))
		{
			bigint_control::DEBUG |= bigint_control::sum;
			debug << "debug: setting debug option bigint_control::sum\n";		
		}
		else if (!strcmp(pptr,"diff"))
		{
			bigint_control::DEBUG |= bigint_control::diff;
			debug << "debug: setting debug option bigint_control::diff\n";		
		}
		else if (!strcmp(pptr,"num_len"))
		{
			bigint_control::DEBUG |= bigint_control::num_len;
			debug << "debug: setting debug option bigint_control::num_len\n";		
		}
		else if (!strcmp(pptr,"gcd"))
		{
			bigint_control::DEBUG |= bigint_control::gcd;
			debug << "debug: setting debug option bigint_control::gcd\n";		
		}
		else if (!strcmp(pptr,"all"))
		{
			bigint_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option bigint_control::all\n";		
		}
#endif
	}
	else if (option == "bigreal")
	{
#ifdef DEBUG_BIGREAL
//		if (!strcmp(pptr,"san"))
//		{
//			bigreal_control::DEBUG |= bigreal_control::sanitize;
//			debug << "debug: setting debug option bigreal_control::sanitize\n";		
//		}
		if (!strcmp(pptr,"+"))
		{
			bigreal_control::DEBUG |= bigreal_control::add;
			debug << "debug: setting debug option bigreal_control::add\n";		
		}
		else if (!strcmp(pptr,"-"))
		{
			bigreal_control::DEBUG |= bigreal_control::subtract;
			debug << "debug: setting debug option bigreal_control::subtract\n";		
		}
		else if (!strcmp(pptr,"*"))
		{
			bigreal_control::DEBUG |= bigreal_control::multiply;
			debug << "debug: setting debug option bigreal_control::multiply\n";		
		}
		else if (!strcmp(pptr,"/"))
		{
			bigreal_control::DEBUG |= bigreal_control::divide;
			debug << "debug: setting debug option bigreal_control::divide\n";		
		}
		else if (!strcmp(pptr,"=="))
		{
			bigreal_control::DEBUG |= bigreal_control::equal;
			debug << "debug: setting debug option bigreal_control::equal\n";		
		}
		else if (!strcmp(pptr,"gt"))
		{
			bigreal_control::DEBUG |= bigreal_control::greater;
			debug << "debug: setting debug option bigreal_control::greater\n";		
		}
		else if (!strcmp(pptr,"out"))
		{
			bigreal_control::DEBUG |= bigreal_control::output;
			debug << "debug: setting debug option bigreal_control::output\n";		
		}
		else if (!strcmp(pptr,"in"))
		{
			bigreal_control::DEBUG |= bigreal_control::input;
			debug << "debug: setting debug option bigreal_control::input\n";		
		}
		else if (!strcmp(pptr,"sum"))
		{
			bigreal_control::DEBUG |= bigreal_control::sum;
			debug << "debug: setting debug option bigreal_control::sum\n";		
		}
		else if (!strcmp(pptr,"diff"))
		{
			bigreal_control::DEBUG |= bigreal_control::diff;
			debug << "debug: setting debug option bigreal_control::diff\n";		
		}
		else if (!strcmp(pptr,"num_len"))
		{
			bigreal_control::DEBUG |= bigreal_control::num_len;
			debug << "debug: setting debug option bigreal_control::num_len\n";		
		}
		else if (!strcmp(pptr,"carry"))
		{
			bigreal_control::DEBUG |= bigreal_control::carry;
			debug << "debug: setting debug option bigreal_control::carry\n";		
		}
		else if (!strcmp(pptr,"all"))
		{
			bigreal_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option bigreal_control::all\n";		
		}
#endif
	}
	else if (option == "matrix")
	{
#ifdef DEBUG_MATRIX	
		if (!strcmp(pptr,"det"))
		{
			matrix_control::DEBUG |= matrix_control::determinant;
			debug << "debug: setting debug option matrix_control::determinant\n";		
		}
		if (!strcmp(pptr,"*"))
		{
			matrix_control::DEBUG |= matrix_control::multiply;
			debug << "debug: setting debug option matrix_control::multiply\n";		
		}
		if (!strcmp(pptr,"inv"))
		{
			matrix_control::DEBUG |= matrix_control::inverse;
			debug << "debug: setting debug option matrix_control::inverse\n";		
		}
		else if (!strcmp(pptr,"all"))
		{
			matrix_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option matrix_control::all\n";		
		}
#endif
	}
	else if (option == "poly")
	{
#ifdef DEBUG_POLYNOMIAL
		if (!strcmp(pptr,"no-gen"))
		{
			polynomial_control::DEBUG >>= 1;
			polynomial_control::DEBUG <<= 1;
			debug << "debug: clearing debug option polynomial_control::general\n";		
		}
		else if (!strcmp(pptr,"san"))
		{
			polynomial_control::DEBUG |= polynomial_control::sanitize;
			debug << "debug: setting debug option polynomial_control::sanitize\n";		
		}
		else if (!strcmp(pptr,"+"))
		{
			polynomial_control::DEBUG |= polynomial_control::add;
			debug << "debug: setting debug option polynomial_control::add\n";		
		}
		else if (!strcmp(pptr,"*"))
		{
			polynomial_control::DEBUG |= polynomial_control::multiply;
			debug << "debug: setting debug option polynomial_control::multiply\n";		
		}
		else if (!strcmp(pptr,"/"))
		{
			polynomial_control::DEBUG |= polynomial_control::divide;
			debug << "debug: setting debug option polynomial_control::divide\n";		
		}
		else if (!strcmp(pptr,"gcd"))
		{
			polynomial_control::DEBUG |= polynomial_control::gcd;
			debug << "debug: setting debug option polynomial_control::gcd\n";		
		}
		else if (!strcmp(pptr,"in"))
		{
			polynomial_control::DEBUG |= polynomial_control::input;
			debug << "debug: setting debug option polynomial_control::input\n";		
		}
		else if (!strcmp(pptr,"all"))
		{
			polynomial_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option polynomial_control::all\n";		
		}
#endif
	}	
	if (option == "surd")
	{
#ifdef DEBUG_SURD
		if (!strcmp(pptr,"+"))
		{
			surd_control::DEBUG |= surd_control::add;
			debug << "debug: setting debug option surd_control::add\n";		
		}
		else if (!strcmp(pptr,"-"))
		{
			surd_control::DEBUG |= surd_control::subtract;
			debug << "debug: setting debug option surd_control::subtract\n";		
		}
		else if (!strcmp(pptr,"*"))
		{
			surd_control::DEBUG |= surd_control::times;
			debug << "debug: setting debug option surd_control::times\n";		
		}
		else if (!strcmp(pptr,"=="))
		{
			surd_control::DEBUG |= surd_control::equal;
			debug << "debug: setting debug option surd_control::equal\n";		
		}
		else if (!strcmp(pptr,"gt"))
		{
			surd_control::DEBUG |= surd_control::greater;
			debug << "debug: setting debug option surd_control::greater\n";		
		}
		else if (!strcmp(pptr,"in"))
		{
			surd_control::DEBUG |= surd_control::input;
			debug << "debug: setting debug option surd_control::input\n";		
		}
		else if (!strcmp(pptr,"newton"))
		{
			surd_control::DEBUG |= surd_control::newton;
			debug << "debug: setting debug option surd_control::newton\n";		
		}
		else if (!strcmp(pptr,"all"))
		{
			surd_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option surd_control::all\n";		
		}
#endif
	}
	else
		set_main_debug_option_parameter(pptr, option); // has to be provided by the calling programme
}

void debug_help ()
{
	cout << "\n\tdebug options may be specified using the syntax -#[option:option{parameter:parameter}:option]" << endl;
	cout << "\tsupported options and parameters:" << endl;
	set_debug_option(0,0);
	set_main_debug_option(0,0);
	cout << "\tdefault options set by omitting [...]" << endl;
	cout << "\tdefault options:" << endl;
	main_display_default_options();
	cout << endl;
	exit(0);
}
