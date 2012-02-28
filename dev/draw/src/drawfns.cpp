/************************************************************************
                  Support functions for draw

                  
void print(generic_code_data& code_data, ostream& s, string prefix="")
string invstr(string& str)
void get_word (ifstream& input, string& buffer, string& title)
void read_immersion_code (generic_code_data& code_data, string input_string)
void read_peer_code (generic_code_data& code_data, string input_string)
void write_immersion_code(ostream& s, generic_code_data& code_data)
void write_peer_code(ostream& s, generic_code_data& code_data, bool labelled=true);
void write_code_data(ostream& s, generic_code_data& code_data)
void read_code_data (generic_code_data& code_data, string input_string)
bool valid_knotoid_input(generic_code_data& code_data, vector<int>& shortcut_crossing)
double badness(string vertex_file, generic_code_data& code_data)
bool first_occurrence(matrix<int>& code_table, matrix<int>& infinite_cycle_first_visit, int crossing, int position)
void write_metapost(ofstream& os, matrix<int>& code_table, string title, metapost_control& mp_control)
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
void print (metapost_control& mp_control, ostream& os, string prefix)

**************************************************************************/
using namespace std;
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <complex>

/********************* External variables ***********************/
extern string		title;

extern ifstream     input;
extern ofstream 	output;
extern ofstream     debug;

extern bool DRAW_IMMERSION_ONLY;
extern bool ONE_METAPOST_PATH;
extern double badness_threshold; 
extern int metapost_coordinate_scale_factor;
extern int metapost_hyperbolic_scale_factor;
extern int DRAW_IN_HYPERBOLIC_DISC;

#include <util.h>
#include <matrix.h>
#include <draw.h>

bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles);
void write_peer_code(ostream& s, generic_code_data& code_data, bool labelled=true);
void hyperbolic_representation(matrix<double>& centre, vector<double>& radius);
void mobius_transform_to_origin(complex<double>& a, complex<double>& b, complex<double> z);
complex<double> mobius_transformation(complex<double> a, complex<double> b, complex<double> z);
complex<double> inv_mobius_transformation(complex<double> a, complex<double> b, complex<double> z);

void print(generic_code_data& code_data, ostream& s, string prefix="")
{
	int num_crossings = code_data.num_crossings;
	s << prefix << "code_type = " << (code_data.type == generic_code_data::immersion_code? "immersion_code" : "peer_code") << endl;
	s << prefix << "head = " << code_data.head << endl;
	s << prefix << "num_crossings = " << num_crossings << endl;
	s << prefix << "num_components = " << code_data.num_components << endl;
	s << prefix << "type:      ";	
	matrix<int>& code_table = code_data.code_table;
	
	for (int i=0; i< num_crossings; i++)
		s << code_table[TYPE][i] << ' ';
	s << endl;
	
	s << prefix << "odd peer: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[OPEER][i] << ' ';
	s << endl;
	s << prefix << "even peer:  ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[EPEER][i] << ' ';
	s << endl;	
	s << prefix << "even term: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[EVEN_TERMINATING][i] << ' ';
	s << endl;
	s << prefix << "odd term:  ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[ODD_TERMINATING][i] << ' ';
	s << endl;
	s << prefix << "even orig: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[EVEN_ORIGINATING][i] << ' ';
	s << endl;
	s << prefix << "odd orig:  ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[ODD_ORIGINATING][i] << ' ';
	s << endl;
	s << prefix << "label:     ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[LABEL][i] << ' ';
	s << endl;	
	s << prefix << "component: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[COMPONENT][i] << ' ';
	s << endl;	

	s << prefix << "num_component_edges: ";
	for (int i=0; i< code_data.num_components; i++)
		s << code_data.num_component_edges[i] << ' ';
	s << endl;
	s << prefix << "first_edge_on_comonent: ";
	for (int i=0; i< code_data.num_components; i++)
		s << code_data.first_edge_on_component[i] << ' ';
	s << endl;
	s << prefix << "term_crossing: ";
	for (int i=0; i< 2*num_crossings; i++)
		s << code_data.term_crossing[i] << ' ';
	s << endl;
	s << prefix << "orig_crossing: ";
	for (int i=0; i< 2*num_crossings; i++)
		s << code_data.orig_crossing[i] << ' ';
	s << endl;	
}


/* invstr creates the inverse of src and returns it to the call.
   The string src must be composed of si, -si, and ti components only.
   The inverse of ti is again ti.
*/
string invstr(string& str)
{
    char*	src = c_string(str);
    char*	loc_buf = new char[2*str.size() + 1]; // a slight overestimate, but always enough
    char*	mark = src;
    bool	not_finished = true;

    /* start at the end of the string and work back */
    char* cptr = strchr(src,0);
    if (strlen(src))
	    cptr--;  /* if src is the null string, cptr = src */

    char* lptr = loc_buf;
    do
    {
		if (isdigit(*cptr))
	    	cptr--;
		else
		{
		    if (*cptr == 's')
	    	{
				if (cptr !=src && *(cptr-1) == '-')
			    	mark = cptr-1;
				else
				{
		    		*lptr++ = '-';
		    		mark = cptr;
				}
				*lptr++ = 's';
	    	}
	    	else if (*cptr == 't')
	    	{
				mark = cptr;
				*lptr++ = 't';
		    }

		    cptr++;
		    while (isdigit(*cptr))
				*lptr++ = *cptr++;
	    	if (mark == src)
				not_finished = false;
	    	else
				cptr = mark-1;
		}
    } while (not_finished);

    *lptr = '\0';
    return string(loc_buf);
}

/* get_word retrieves the next braid word, Gauss code, peer code or labelled 
   immersion code from the input file, carrying out any manipulation required, 
   places the resultant word in buffer and any title string in title.  The 
   function does not remove the '--' leader from the title and sets the title 
   to the empty string if none is found.  Any qualifiers are appended to the 
   input string, the function does not remove the braces {} from around qualifiers.
*/
void get_word (ifstream& input, string& buffer, string& title)
{
    string 	next_line;
	string qualifiers;
    string 	loc_buf;
	string	w1,w2,w3,w4,w5,w6,w7,w8;
	string* target;
	char*	lptr;
    bool   	word_found = false;
	bool	immersion_code = false;
	bool	peer_code = false;
    bool	not_finished;

    /* start by setting the buffer and title to the empty string */
    buffer.clear();
    title.clear();

	while (!word_found && getline(input, next_line))
    {

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "draw::get_word: read next line: " <<  next_line << endl;

		/* kill the <cr> at the end of the line, if one exists */
		string::size_type pos = next_line.find("\r");
		if (pos != string::npos)
		     next_line[pos] = ' ';
		
		/* remove any comments from the end of the line */
		pos = next_line.find(';');
		if (pos != string::npos)
			next_line = next_line.substr(0,pos);

 		
		if (next_line.find("exit") != string::npos)
			break;
			
		/* move any qualifiers from the end of the line into the qualifiers string 
		   by assigning the qualifier sting here any misplaced qualifiers in the input file 
		   are ignored.  Only when qualifiers are added to the end of a braid statement or the
		   last line of an immersion code or Gauss code will they be acted upon.
		   
		   Also, removing the qualifiers from the current line at this stage avoids any 
		   confusion between qualifiers and braid statement
		*/
		pos = next_line.find('{');
		if (pos != string::npos && next_line.find("--") == string::npos) // titles may contain braces
		{
			qualifiers = next_line.substr(pos);

if (draw_control::DEBUG >= draw_control::INTERMEDIATE)
    debug << "draw::get_word: read qualifiers " << qualifiers << " from line " << next_line << endl;

			next_line = next_line.substr(0,pos);
			
		}

		char* line_buf = c_string(next_line);
	    /* if this line contains a switch, or programme options ignore it */
	    if ( (strchr(line_buf,'[') &&
	          !strchr(line_buf,'\\') &&
	          !strchr(line_buf,'/')
	         )|| 
		     (  (strchr(line_buf,'s') || strchr(line_buf,'S')) && 
			     strchr(line_buf,'=') && 
				!strchr(line_buf,'w') //not an assignment statement
			 )
		   )
		{
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "draw::get_word: line is a switch or contains programme options" << endl;
			goto done_with_line;
		}
		
	    /* is there any whitespace at the beginning of the line? */
	    lptr = line_buf;
	    while (isspace(*lptr))
			lptr++;

	    if (strlen(lptr) && *lptr != '\n')
	    {
			
if (draw_control::DEBUG >= draw_control::INTERMEDIATE)
    debug << "draw::get_word: next word line = " << lptr << endl;

			/* check for a title line */
			char* cptr = lptr;
			if (*cptr == '-' && *(cptr+1) == '-')
			{
		    	title = lptr;

				goto done_with_line;
			}

			/* the target of the line parsing is either one of the word buffers 
			   w1-w8, or is 'buffer' itself, if this line is a braid statement, an 
			   immersion code, or a Gauss code, .  An immersion code is 
			   indicated by the presence of a '(' character, and we look for
			   that first. A Gauss code is indicated by the presence
			   of a '/' character but no '('.			
			
			   We check that the line is not part of a Gauss or immersion
			   code when looking for braid words, or is not part of a braid
			   word when looking for Gauss or immersion codes
			*/
			
			if (strchr(lptr,'('))
			{
				immersion_code = true;
			}
			else if (strchr(line_buf,'['))
			{
				peer_code = true;
			}
			
			
			/* Here we have the kind of input we're looking for on the next_line */
			if ((immersion_code || peer_code) && (strchr(lptr,'/') || strchr (lptr,'\\')))
			{
			    /* Take out line escapes and build up either the peer code,
				   labelled immersion code or the Gauss code in buffer
			    */
			    cptr = strchr(lptr,'\\');
			    if (cptr)
					*cptr = ' ';

			    /* copy the line into buffer and decide if there's more to come */
				buffer += string(lptr);

			    /* This test means we cannot break lines after
				   the '/' character in the input file
				*/
				if (strchr(lptr,'/'))
				{
					word_found = true;					
					
					/* Add the qualifiers to the buffer */

if (draw_control::DEBUG >= draw_control::INTERMEDIATE)
    debug << "draw::get_word: adding qualifiers " << qualifiers << " to buffer " << buffer << endl;;

					buffer += qualifiers;
				}
			    
				goto done_with_line;
			}
			
			/* Here we're looking for a braid word, first check whether this is 
			   an assignment statement or a braid statement, if it's a braid
			   statement, we're done once we've parsed this line. */
			cptr = strchr(lptr, '=');
			
			if (cptr) // assignment statement
			{
			    cptr = strchr(lptr,'w');
		    	switch (*++cptr)
		    	{
					case '1': target = &w1; break;
					case '2': target = &w2; break;
					case '3': target = &w3; break;
					case '4': target = &w4; break;
					case '5': target = &w5; break;
					case '6': target = &w6; break;
					case '7': target = &w7; break;
					case '8': target = &w8; break;
					default: target = &w1;
		    	}				
			}
			else // braid statement
			{
			    target = &buffer;
		    	word_found = true;
			}

			/* now parse the line into target, first move cptr to the
			   start of the line contents, i.e. past any = sign.
			*/

			cptr = strchr(lptr, '=');
			if (cptr)
			    cptr++;
			else
			    cptr = lptr;

			/* move through whitespace */
			while (isspace(*cptr))
			    cptr++;

			(*target).clear();
			not_finished = true;
			do
			{
			    if (isspace(*cptr) || *cptr == '\0')
					not_finished = false;
		    	else if (*cptr == 'w')
		    	{
					switch (*++cptr)
					{
			    		case '1': *target += w1; break;
			    		case '2': *target += w2; break;
			    		case '3': *target += w3; break;
			    		case '4': *target += w4; break;
			    		case '5': *target += w5; break;
			    		case '6': *target += w6; break;
			    		case '7': *target += w7; break;
			    		case '8': *target += w8; break;
			    		default: *target += w1;
					}
					cptr++;
		    	}
		    	else if (*cptr == '-' && *(cptr+1) == 'w')
		    	{
					cptr++; /* moves cptr to the w character */
					switch (*++cptr)
					{
					    case '1': *target += invstr(w1); break;
					    case '2': *target += invstr(w2); break;
				    	case '3': *target += invstr(w3); break;
				    	case '4': *target += invstr(w4); break;
			    		case '5': *target += invstr(w5); break;
				    	case '6': *target += invstr(w6); break;
				    	case '7': *target += invstr(w7); break;
				    	case '8': *target += invstr(w8); break;
				    	default: *target += invstr(w1);
					}
					cptr++;
		    	}
		    	else
		    	{
					/* copy the characters up to the next whitespace or 'w' into 
					   a local buffer and appent to target.
				    */
					char* copy_buf = new char[strlen(cptr)+1];
					char* sptr = copy_buf;
				
					while (!isspace(*cptr) && *cptr != '\0' && *cptr != 'w')
			    		*sptr++ = *cptr++;
					if (*cptr == 'w' && *(cptr-1) == '-')
					{
			    		/* move back one */
			    		cptr--;
			    		sptr--;
					}
					*sptr = '\0';
					*target += string(copy_buf);

					delete[] copy_buf;
		    	}
			} while (not_finished);
			
			if (word_found)
			{
				/* we have just parsed a braid statement into the buffer so we
				   append any qualifiers provided with that braid statement
				*/
if (draw_control::DEBUG >= draw_control::INTERMEDIATE)
    debug << "draw::get_word: adding qualifiers " << qualifiers << " to buffer " << buffer << endl;;

				buffer += qualifiers;
			}
	    }
done_with_line:
		delete[] line_buf;
    } //end of while (getline(input, next_line) && !word_found)

    if (!word_found)
    	buffer = "exit";
		
if (draw_control::DEBUG >= draw_control::INTERMEDIATE)
    debug << "draw::get_word: returning buffer " << buffer << endl;;
		
}


/* This function reads a peer code from the input_string into a generic_code_data object.  Within that
   object it creates a code_table that is a matrix<int>(9,num_crossings).  The function records in this
   matrix the crossing type, the peers of the even and odd edges 
   that terminate at the crossing, the component to which the naming edge belongs and the crossing label. 
   For each crossing it also records the odd and even terminating and originating vertices.
   It uses the #define TYPE, OPEER, EPEER, COMPONENT, LABEL, EVEN_TERMINATING, ODD_TERMINATING, 
   EVEN_ORIGINATING, ODD_ORIGINATING to index the rows of the matrix.
   If non-zero, the int* head is used to identify the first shortcut crossing in a knotoid (identified by 
   a '^' character after the crossing number).  This in turn identifies where in the peer code the head of 
   the knotoid is located. 
*/
void read_peer_code (generic_code_data& code_data, string input_string)
{
	code_data.type = generic_code_data::peer_code;
	code_data.head = -1;
	
	char* inbuf = c_string(input_string);
	int num_crossings = 0;
	char* cptr = strchr(inbuf,'/');
	cptr++;
	while (*cptr != '\0')
	{
		if (*cptr == '+' || *cptr == '-' || *cptr == '*' || *cptr == '#' )
			num_crossings++;

		cptr++;
	}	

if (draw_control::DEBUG >= draw_control::BASIC)
	debug << "read_peer_code: num_crossings determined from labelled peer code = " << num_crossings << endl;
	
	code_data.num_crossings = num_crossings;
	
	matrix<int> code_table(CODE_TABLE_SIZE,num_crossings);

	int num_edges = 2*num_crossings;
	int component = 0;

	/* assume all crossings are TYPE2 and set those 
	   indicated as TYPE1 acordingly.
	*/
	for (int i=0; i<num_crossings; i++)
		code_table[TYPE][i] = TYPE2;
		
	cptr = strchr(inbuf,'[');
	cptr++;
    for (int i=0; i< num_crossings; i++)
	{
		/* We write the odd-numbered peer of edge 2i in code_table[OPEER][i] and the 
		   even-numbered peer of edge 2i+1 in code_table[EPEER][i] for i=0,...,num_crossings-1
		*/
		bool crossing_complete = false;
		do
		{
			if (*cptr == '-')
			{
				code_table[TYPE][i] = TYPE1;
				cptr++;
			}
			else if (isdigit(*cptr))
			{
				char* mark = cptr;
				while (isdigit(*cptr)) cptr++;

				get_number (code_table[OPEER][i],mark);
				code_table[EPEER][(code_table[OPEER][i]-1)/2] = 2*i;
				code_table[COMPONENT][i] = component;
				
				if (*cptr == '^')
				{
					/* we've found the knotoid head, note the 
					corresponding crossing number */
					code_data.head = i;
					cptr++;
				}			
				crossing_complete = true;
			}
			else if (*cptr == ',')
			{
				/* we've come to the end of a component and are about to 
				   read the first crossing of the new component
				*/
				component++;
				cptr++;
			}
			else
				cptr++;		
		} while (!crossing_complete);
	}
	
	/* now do the labels */
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
			code_table[LABEL][count++] = POSITIVE;
		else if (*cptr == '-')
			code_table[LABEL][count++] = NEGATIVE;
		else if (*cptr == '*')
			code_table[LABEL][count++] = VIRTUAL;
		else if (*cptr == '#')
			code_table[LABEL][count++] = FLAT;
		cptr++;	
	}

	/* write the component data into code_data */
	int num_components = component+1;
	code_data.num_components = num_components;

	vector<int> num_component_edges(num_components);
	vector<int> first_edge_on_component(num_components);
	first_edge_on_component[0]=0;
	component=0;
	
	for (int i=1; i< num_crossings; i++)
	{
		if (code_table[COMPONENT][i] != code_table[COMPONENT][i-1])
		{
			num_component_edges[component] = 2*i-first_edge_on_component[component];
			component++;
			first_edge_on_component[component] = 2*i;
		}
	}
	
	num_component_edges[component] = num_edges - first_edge_on_component[component];
	
	/* now we can write the originating and terminating vertices and edges */
	vector<int> term_crossing(num_edges);
	vector<int> orig_crossing(num_edges);
	
	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[code_table[OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = code_table[COMPONENT][(code_table[OPEER][i]-1)/2];
		int peer_successor = (code_table[OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

		orig_crossing[peer_successor] = i;
		
		code_table[EVEN_TERMINATING][i] = 2*i;
		code_table[ODD_ORIGINATING][i] = 2*i+1;
		code_table[ODD_TERMINATING][i] = code_table[OPEER][i];
		code_table[EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = code_table;
	code_data.num_component_edges = num_component_edges;
	code_data.first_edge_on_component = first_edge_on_component;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;
	
if (draw_control::DEBUG >= draw_control::INTERMEDIATE)
{
	debug << "read_peer_code: code data produced from peer code:" << endl;
	print(code_data,debug,"read_peer_code: ");	
}

	delete[] inbuf;
}

/* This function uses head, num_crossings and the OPEER, TYPE, LABEL 
   and COMPONENT rows of the code_table in the generic_code_data structure.
*/
void write_peer_code(ostream& s, generic_code_data& code_data, bool labelled)
{
	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	
	s << "[";
	if (code_table[TYPE][0] == TYPE1)
		s << '-';
	s << code_table[OPEER][0] << ' ';
	for (int i=1; i<num_crossings; i++)
	{
		if (code_table[COMPONENT][i] != code_table[COMPONENT][i-1])
			s << ", ";
			
		if (code_table[TYPE][i] == TYPE1)
			s << '-';
		s << code_table[OPEER][i];
		if (i == code_data.head)
			s << "^ ";
		else if ( i < num_crossings-1 && code_table[COMPONENT][i] == code_table[COMPONENT][i+1])
			s << ' ';
	}
	s << "]";
	if (labelled)
	{
		s << "/";
		for (int i=0; i< num_crossings; i++)
		{
			if (code_table[LABEL][i] == POSITIVE)
				s << "+ ";
			else if (code_table[LABEL][i] == NEGATIVE)
				s << "- ";
			else if (code_table[LABEL][i] == VIRTUAL)
				s << "* ";
			else // (code_table[LABEL][i] == FLAT)
				s << "# ";
		}	
	}
}

void read_immersion_code (generic_code_data& code_data, string input_string)
{
	code_data.type = generic_code_data::immersion_code;
	code_data.head = -1;

	char* inbuf = c_string(input_string);
	int num_crossings = 0;
	char* cptr = strchr(inbuf,'/');
	cptr++;
	while (*cptr != '\0')
	{
		if (*cptr == '+' || *cptr == '-' || *cptr == '*')
			num_crossings++;

		cptr++;
	}	

if (draw_control::DEBUG >= draw_control::BASIC)
	debug << "read_immersion_code: num_crossings determined from labelled immersion code = " << num_crossings << endl;

	int num_edges = 2*num_crossings;
	code_data.num_crossings = num_crossings;
	code_data.num_components = 1;
	vector<int> num_component_edges(1);
	vector<int> first_edge_on_component(1);
	first_edge_on_component[0]=0;
	num_component_edges[0] = num_edges;
	
	matrix<int> code_table(CODE_TABLE_SIZE,num_crossings);
	
	for (int i=0; i< num_crossings; i++)
		code_table[COMPONENT][i] = 0;

	vector<int> perm(num_crossings);
	vector<int> type(num_crossings);
	
	cptr = strchr(inbuf,'(');
    while(*cptr != '/')
	{
		int count = 0;

		/* We write the crossing numbers in the order they appear in the immersion code 
		   into perm together with their corresponding type.  Then we unravel 
		   the permuatation writing the PERM and TYPE rows from INVERSE and LABEL.
		*/
		
		for (int i=0; i<num_crossings; i++)
			type[i] = TYPE2;

		while (*cptr != ')')
		{		
			if (*cptr == '-')
			{
				type[count] = TYPE1;
				cptr++;
			}
			else if (isdigit(*cptr))
			{
				char* mark = cptr;
				while (isdigit(*cptr)) cptr++;

				get_number (perm[count],mark);
				
				if (*cptr == '^')
				{
					/* we've found the knotoid head, note the 
					   corresponding crossing number */
					code_data.head = perm[count];
					cptr++;
				}
				
				count++;
			}
			else
				cptr++;
		}

		
if (draw_control::DEBUG >= draw_control::DETAIL)
{
	debug << "read_immersion_code: number of terms in this part of perm = " << count << endl;
	debug << "read_immersion_code: read absolute terms: ";
	for (int i=0; i< count; i++)
		debug << perm[i] << ' ';
	debug << endl;		
	debug << "read_immersion_code: signs: ";
	for (int i=0; i< count; i++)
		debug << type[i] << ' ';
	debug << endl;
}
		
		/* work out this part of the permutation */
		for (int i=0; i<count-1; i++)
		{
			code_table[OPEER][perm[i]] = (2*perm[i+1]-1+num_edges)%num_edges;
			code_table[EPEER][(code_table[OPEER][perm[i]]-1)/2] = 2*perm[i];
			code_table[TYPE][perm[i]] = type[i];
		}
		code_table[OPEER][perm[count-1]] = (2*perm[0]-1+num_edges)%num_edges;
		code_table[EPEER][(code_table[OPEER][perm[count-1]]-1)/2] = 2*perm[count-1];
		code_table[TYPE][perm[count-1]] = type[count-1];

		cptr++; // move over ')'
		while (*cptr == ' ') cptr++;

if (draw_control::DEBUG >= draw_control::DETAIL)
{
	debug << "read_immersion_code: interim code table produced from immersion code:" << endl;
	debug << "read_immersion_code:   type:    ";	
	for (int i=0; i< num_crossings; i++)
		debug << code_table[TYPE][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   OPEER: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[OPEER][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   EPEER: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[EPEER][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   label: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[LABEL][i] << ' ';
	debug << endl;
}

	}

	/* now do the labels */
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
			code_table[LABEL][count++] = POSITIVE;
		else if (*cptr == '-')
			code_table[LABEL][count++] = NEGATIVE;
		else if (*cptr == '*')
			code_table[LABEL][count++] = VIRTUAL;
		cptr++;	
	}
	
	/* now we can write the originating and terminating vertices and edges */
	vector<int> term_crossing(num_edges);
	vector<int> orig_crossing(num_edges);

	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[code_table[OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = code_table[COMPONENT][(code_table[OPEER][i]-1)/2];
		int peer_successor = (code_table[OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

		orig_crossing[peer_successor] = i;
		
		code_table[EVEN_TERMINATING][i] = 2*i;
		code_table[ODD_ORIGINATING][i] = 2*i+1;
		code_table[ODD_TERMINATING][i] = code_table[OPEER][i];
		code_table[EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = code_table;
	code_data.num_component_edges = num_component_edges;
	code_data.first_edge_on_component = first_edge_on_component;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;
	
if (draw_control::DEBUG >= draw_control::INTERMEDIATE)
{
	debug << "read_immersion_code: code data produced from immersion code:" << endl;
	print(code_data,debug,"read_immersion_code: ");	
}

	delete[] inbuf;
}

/* This function only uses head, num_crossings and the OPEER, TYPE, and LABEL rows 
   of the code_table in the generic_code_data structure.
*/
void write_immersion_code(ostream& s, generic_code_data& code_data)
{
	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	
	vector<int> flag(num_crossings); // used to track which crossings have been written
	for (int i=0; i<num_crossings; i++)
		flag[i] = 1;
		
	/* work out the permutation from the generic code data */
	vector<int> perm(num_crossings);
	for (int i=0; i< num_crossings; i++)
		perm[i] = ((code_table[OPEER][i]+1)/2)%num_crossings;
	
	bool found;
	do
	{
		found = false;
		int crossing;
		int start;
		
		/* look for another starting place in rptr */
		for (int i=0; i<num_crossings; i++)
		{
			if (flag[i])
			{
				crossing = start = i;
				flag[i] = 0;
				found = true;
				break;
			}
		}

		if (found)
		{
			s << "(";
			if (code_table[TYPE][crossing] == TYPE1)
				s << '-';
			s << crossing;

			int next_crossing;
			do
			{
				next_crossing = perm[crossing];
				s << " ";
				if (code_table[TYPE][next_crossing] == TYPE1)
					s << '-';
				s << next_crossing;

				if (next_crossing == code_data.head)
					s << '^';
					
				flag[next_crossing] = 0;
				crossing = next_crossing;
			} while (perm[crossing] != start);
			s << ")";
		}
	} while (found);

	s << " / ";
							
	for (int i=0; i< num_crossings; i++)
	{
		if (code_table[LABEL][i] == POSITIVE)
			s << "+ ";
		if (code_table[LABEL][i] == NEGATIVE)
			s << "- ";
		if (code_table[LABEL][i] == VIRTUAL)
			s << "* ";
	}	
}


void write_code_data(ostream& s, generic_code_data& code_data)
{
	if (code_data.type == generic_code_data::peer_code)
		write_peer_code(s, code_data);
	else
		write_immersion_code(s,code_data);
}

void read_code_data (generic_code_data& code_data, string input_string)
{
	if (input_string.find('[') != string::npos)
		read_peer_code(code_data, input_string);
	else
		read_immersion_code(code_data, input_string);
}

/* The generic code data code table and crossing number identify a valid knotoid if the
   specified crossing is the start of a shortcut.  A shortcut from the head
   to the leg of a knotoid must pass under all the strands of the knotoid; thus
   if the label associated with the head crossing is '+' then the head lies on the
   odd-numbered semi-arc of the diagram, and if the label is '-' it lies
   on the even-numbered semi-arc.  From this arc onwards all the crossings we 
   encounter must be approached on the under-arc if this is a valid knotoid 
   combination of code_table and head.
   
   As we check the validity of the input we note the shortcut crossings and also
   the algebraic number of intersections of the knotoid and the shortcut.  This is 
   the number of times the knotoid crosses the shortcut from right to left minus 
   the number of times it crosses from left to right.  However, for the definition
   the shortcut is oriented from the leg to the head and we will be traversing it in 
   the opposite direction.  
   
   We set shortcut_crossing[i] to be non-zero if the ith crossing is a shortcut crossing.
   We set it to 1 if the knotoid crosses the shortcut from right to left and -1 if it crosses from 
   left to right, orienting the shortcut from leg to head.
   
*/
bool valid_knotoid_input(generic_code_data& code_data, vector<int>& shortcut_crossing)
{
	matrix<int>& code_table = code_data.code_table;
	int head = code_data.head;
	
	/* we must have been given a crossing indicating the position of the head */
	if (head == -1)
	{
		
if (draw_control::DEBUG >= draw_control::BASIC)
	debug << "braid::valid_knotoid_input: no indication of knotoid head provided, returning false" << endl;
	
		return false;
	}
		
	bool valid = true;
	
	int semi_arc;
	int num_crossings = code_table.numcols();

if (draw_control::DEBUG >= draw_control::BASIC)
{
	debug << "braid::valid_knotoid_input: code_table " << endl;
	print(code_table,debug,3,"valid_knotoid_input: ");	
	debug << "braid::valid_knotoid_input: head =" << head << endl;
}

	if (code_table[LABEL][head] == POSITIVE)
		semi_arc = code_table[OPEER][head];
	else if (code_table[LABEL][head] == NEGATIVE)
		semi_arc = 2*head;
	else
	{
		/* the first shortcut crossing cannot be virtual */
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input: head " << 0 << " indicates first shortcut crosing is virtual" << endl;
		return false;
	}

if (draw_control::DEBUG >= draw_control::BASIC)
	debug << "braid::valid_knotoid_input: head " << head << " lies on semi-arc " << semi_arc << endl;

    /* set the shortcut crossing flag for the head crossing */
    if (semi_arc % 2)
    {
		if (code_table[TYPE][head] == TYPE1)				
		{
			shortcut_crossing[head] = -1; 
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the left" << endl;
		}
		else
		{
			shortcut_crossing[head] = 1;
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the right" << endl;
		}
	}
	else
	{
		if (code_table[TYPE][head] == TYPE1)				
		{
			shortcut_crossing[head] = 1; 
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the right" << endl;
		}
		else
		{
			shortcut_crossing[head] = -1;
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the left" << endl;
		}
	}
		

	for (int i = semi_arc+1; i< 2*num_crossings; i++)
	{

if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input: semi-arc " << i << endl;
	
		/* check that semi-arc i is the under-arc at the next crossing */
		if (i%2)
		{				
			int crossing = code_data.term_crossing[i];

if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   crossing " << crossing << endl;
			
			/* label for this crossing must be '+' if odd arc is under-arc */
			if (code_table[LABEL][crossing] != POSITIVE)
			{
				valid = false;
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   not positive, as required for valid knotoid" << endl;
				break;
			}
			else
			{
				if (code_table[TYPE][crossing] == TYPE1)				
				{
					shortcut_crossing[crossing] = -1; 
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   OK, knotoid crosses from the left" << endl;
				}
				else
				{
					shortcut_crossing[crossing] = 1;
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   OK, knotoid crosses from the right" << endl;
				}
					
			}
		}
		else
		{
			int crossing = i/2 ;
			
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   crossing " << crossing << endl;

			/* label for this crossing must be '-' if even arc is under-arc */
			if (code_table[LABEL][crossing] != NEGATIVE)
			{
				valid = false;
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   not negative, as required for valid knotoid" << endl;
				break;
			}
			else
			{
				if (code_table[TYPE][crossing] == TYPE1)				
				{
					shortcut_crossing[crossing] = 1; 
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   OK, knotoid crosses from the right" << endl;
				}
				else
				{
					shortcut_crossing[crossing] = -1;
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input:   OK, knotoid crosses from the left" << endl;
				}
			}
		}
	}
	
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "braid::valid_knotoid_input: returning " << (valid? "true": "false") << endl;
	return valid;
}


/* the badness function looks at consecutive vertices z_{n+1} and z_n in the vertex sequence
   and evaluates |z_{n+1}-z_n| and returns min{100, |z_{i+1} - z_i|}.  For links this is not
   perfect, since at component boundaries the vertex sequence does not reflect the path
   drawn by metapost but the distances between the last and first vertex at a component
   boundary is typically large, so it should not affect the value returned by the function.
*/
double badness(string vertex_file, generic_code_data& code_data)
{
	
	/* create the vertex sequence for writing the metapost path using the terminating 
	   crossing information in code_data for each edge.  The sequence starts with the 
	   type 2 vertex preceeding crossing zero, which is the first type 2 vertex in the 
	   numbering scheme.	       
	*/
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;
	int vertex_sequence_length = 2*num_edges;
	vector<int>& term_crossing = code_data.term_crossing;
	vector<int> vertex_sequence(vertex_sequence_length);
	for (int i=0; i< num_edges; i++)
	{
		vertex_sequence[2*i] = num_crossings+(i%num_edges); //the type_2_vertex at the midpoint of edge i
		vertex_sequence[2*i+1] = term_crossing[i];
	}

	double x1, x2, y1, y2, z, w;
	double worst=badness_threshold;
	
	ifstream input;
	input.open(vertex_file.c_str());
	
	if (input)
	{

		int num_type123_vertices;
		input >> num_type123_vertices;
		
		matrix<double> drawcoords(vertex_sequence_length,2);
		for (int i=0; i< vertex_sequence_length; i++)
		{
			input >> drawcoords[i][0];
			input >> drawcoords[i][1];
		}

		for (int i=1; i<vertex_sequence_length; i++)
		{
			x1 = drawcoords[i-1][0] * metapost_coordinate_scale_factor;
			x2 = drawcoords[i][0] * metapost_coordinate_scale_factor;
			y1 = drawcoords[i-1][1] * metapost_coordinate_scale_factor;
			y2 = drawcoords[i][1] * metapost_coordinate_scale_factor;
			w = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			z = sqrt(w);
			
			if (z < worst) 
				worst = z;
		}
	}
	else
	{
		cout << "\nError draw::badness failed to open " << vertex_file << endl;
		exit(0);
	}
	return worst;
}

/* first_occurrence detemines whether the given position at the given crossing is the first occurrence
   of the infinte region at the crossing or not.  It does this by identifying an edge incident with the
   given position and then checks to see if that edge is either of those recorded in 
   infinite_cycle_first_visit for that crossing.
*/
bool first_occurrence(matrix<int>& code_table, matrix<int>& infinite_cycle_first_visit, int crossing, int position)
{
	int num_edges = 2* code_table.numcols();
	bool first;

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "first_occurrence: crossing = " << crossing << " position = " << position << endl;
	
	int search_edge;
	
	switch (position)
	{
		case 0: search_edge = 2*crossing; // naming edge
				break;
		case 1: search_edge = (code_table[TYPE][crossing] == TYPE2? (2*crossing+1)%num_edges : 2*crossing);
				break;
		case 2: search_edge = (2*crossing+1)%num_edges; // naming_edge sucessor
				break;
		case 3: search_edge = (code_table[TYPE][crossing] == TYPE1? (2*crossing+1)%num_edges : 2*crossing);
				break;
	}
	
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "first_occurrence: search_edge = " << search_edge << endl;

	if (infinite_cycle_first_visit[crossing][0] == search_edge || infinite_cycle_first_visit[crossing][1] == search_edge)
		first = true;
	else
		first = false;
		
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "first_occurrence: first_occurrence = " << (first? "true" : "false") << endl;

	return first;
}

void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control)
{
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.num_components;
	int head = code_data.head;
	int num_edges = 2*num_crossings;
	int num_type12_vertices = num_edges+num_crossings;	

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "write_metapost: mp_control: " << endl;
    print (mp_control, debug, "write_metapost:   ");
}

if (draw_control::DEBUG >= draw_control::DETAIL)
{
    debug << "write_metapost: num_crossings = " << num_crossings << endl;
    debug << "write_metapost: num_edges = " << num_edges << endl;
    debug << "write_metapost: num_type12_vertices = " << num_type12_vertices << endl;
}

	matrix<int>& code_table	 = code_data.code_table;	
	vector<int> orig_crossing = code_data.orig_crossing;
	vector<int> term_crossing = code_data.term_crossing;
	vector<int> num_component_edges = code_data.num_component_edges;
	vector<int> first_edge_on_component = code_data.first_edge_on_component;

	/* record the first edge that enters each crossing */
	vector<int> first_edge(num_crossings,-1);
    for (int i=0; i< num_edges; i++)
    {
		if (first_edge[term_crossing[i]] == -1)
			first_edge[term_crossing[i]] = i;
	}

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "write_metapost: first_edge: ";
    for (int i=0; i< num_crossings; i++)
		debug << first_edge[i] << ' ';
	debug << endl;
}

	int semi_arc;
	if (head != -1)
	{
		if (code_table[LABEL][head] == POSITIVE)
			semi_arc = code_table[OPEER][head];
		else if (code_table[LABEL][head] == NEGATIVE)
			semi_arc = 2*head;

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "write_metapost: head " << head << " lies on semi-arc " << semi_arc << endl;
	
	}
	else
	{
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "write_metapost: no head semi-arc" << endl;
	}

	
	ifstream input;
	
	input.open("__circlepack.out");
	if (!input)
	{
		cout << "\nError opening output file __circlepack.out\n";
		exit(0);
    }
    
//	int num_vertices; //number of type 1,2,3 & 4 vertices
//   input >> num_vertices;
	int num_type123_vertices;
	input >> num_type123_vertices; // number of type 1,2 & 3 vertices
	
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
	debug << "write_metapost: number of type 1,2 & 3 vertices = " << num_type123_vertices << endl;
}
    
    matrix<double> vcoords(num_type123_vertices,2);
    vector<double> vertex_radius(num_type123_vertices);

	for (int i=0; i< num_type123_vertices; i++)
	{
		for (int j=0; j< 2; j++)
			input >> vcoords[i][j];
	}
	for (int i=0; i< num_type123_vertices; i++)
		input >> vertex_radius[i];
		
	/* If we're drawing in the hyperbolic plane, evaluate the Euclidean centre
	   and radius of each circle
	*/
	if (DRAW_IN_HYPERBOLIC_DISC)
		hyperbolic_representation(vcoords, vertex_radius);

	/* scale vertices and radii here */
	
	for (int i=0; i< num_type123_vertices; i++)
	{
		for (int j=0; j< 2; j++)
			vcoords[i][j] *= metapost_coordinate_scale_factor;
	}
	for (int i=0; i< num_type123_vertices; i++)
		vertex_radius[i] *= metapost_coordinate_scale_factor;
	
if (draw_control::DEBUG >= draw_control::DETAIL)
{
    debug << "write_metapost: scaled vertex_coordinates: " << endl;
    for (int i=0; i< num_type123_vertices; i++)
    {
		debug << "write_metapost:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << vcoords[i][j] << ' ';
		debug << endl;
	}	
    debug << "write_metapost: scaled vertex_radius: " << endl;
    for (int i=0; i< num_type123_vertices; i++)
		debug << "write_metapost:   vertex " << i << ": " << vertex_radius[i] << endl;
}
	
	
	double minx = vcoords[0][0];
	double maxx = vcoords[0][0];
	double miny = vcoords[0][1];
	double maxy = vcoords[0][1];
	
	for (int i=1; i< num_type12_vertices; i++)
	{
		for (int j=0; j< 2; j++)
		{
			if (j==0 && vcoords[i][j] < minx)
				minx = vcoords[i][j];
			
			if (j==1 && vcoords[i][j] < miny)
				miny = vcoords[i][j];
			
			if (j==0 && vcoords[i][j] > maxx)
				maxx = vcoords[i][j];
			
			if (j==1 && vcoords[i][j] > maxy)
				maxy = vcoords[i][j];
		}
	}
	
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "write_metapost: minx = " << minx << ", maxx = " << maxx << ", miny = " << miny << ", maxy = " << maxy << endl;

	if (mp_control.rotate)
	{
		if (mp_control.implicit_rotation_centre)
		{
			mp_control.rotation_centre_x = vcoords[mp_control.rotation_centre_z][0];
			mp_control.rotation_centre_y = vcoords[mp_control.rotation_centre_z][1];
		}
		else if (!mp_control.explicit_rotation_centre)
		{
			mp_control.rotation_centre_x = (minx+maxx)/2;
			mp_control.rotation_centre_y = (miny+maxy)/2;
		}
		
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "write_metapost: rotating by " << mp_control.rotation_degrees << " around ";
    if (mp_control.explicit_rotation_centre)
		debug << "explicit centre (" << mp_control.rotation_centre_x << "," << mp_control.rotation_centre_y << ")" << endl;
	else if (mp_control.implicit_rotation_centre)
	{
		debug << "implicit centre z" << mp_control.rotation_centre_z << " = (";
		debug << vcoords[mp_control.rotation_centre_z][0] << "," << vcoords[mp_control.rotation_centre_z][1] << ")" << endl;
	}
	else
		debug << "default centre (" << mp_control.rotation_centre_x << "," << mp_control.rotation_centre_y << ")" << endl;
}

		/* we convert to radians by calculating rotation_degrees/360 * 2 pi and use 355/113 as an approximation to pi */
	
		double radians = mp_control.rotation_degrees * 355 / 180 / 113;
		double M00 = cos(radians);
		double M10 = sin(radians);
		double M01 = M10 * -1;
		double M11 = M00;

if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "write_metapost: cos(rotation_degrees = " << radians << ") = " << M00 << ", sin(rotation_degrees) = " << M10 << endl;

		for (int i=0; i< num_type123_vertices; i++)
		{
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "write_metapost: rotating z" << i << " = (" << vcoords[i][0] << "," << vcoords[i][1] << ")" << endl;
    
			double shift_x = vcoords[i][0] - mp_control.rotation_centre_x;
			double shift_y = vcoords[i][1] - mp_control.rotation_centre_y;

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "write_metapost:   shifted to (" << shift_x << "," << shift_y << ")" << endl;

			vcoords[i][0] = M00*shift_x + M01*shift_y;
			vcoords[i][1] = M10*shift_x + M11*shift_y;

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "write_metapost:   rotated to (" << vcoords[i][0] << "," << vcoords[i][1] << ")" << endl;
    
			vcoords[i][0] += mp_control.rotation_centre_x;
			vcoords[i][1] += mp_control.rotation_centre_y;

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "write_metapost:   shifted back to (" << vcoords[i][0] << "," << vcoords[i][1] << ")" << endl;
			
		}

		/* re-evaluate the min/max values after the rotation */
		minx = vcoords[0][0];
		maxx = vcoords[0][0];
		miny = vcoords[0][1];
		maxy = vcoords[0][1];
		
		for (int i=1; i< num_type12_vertices; i++)
		{
			for (int j=0; j< 2; j++)
			{
				if (j==0 && vcoords[i][j] < minx)
					minx = vcoords[i][j];
				
				if (j==1 && vcoords[i][j] < miny)
					miny = vcoords[i][j];
				
				if (j==0 && vcoords[i][j] > maxx)
					maxx = vcoords[i][j];
				
				if (j==1 && vcoords[i][j] > maxy)
					maxy = vcoords[i][j];
			}
		}

	}

if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "write_metapost: after rotation minx = " << minx << ", maxx = " << maxx << ", miny = " << miny << ", maxy = " << maxy << endl;

	input.close();
	
	os << "\n\n";
	
	if (title.length())
		os << "% " << title << endl;

	matrix<int>* Cycle = new matrix<int>(num_crossings+2, num_edges+1);
	int num_left_cycles;
	int num_cycles;
	
	os << "% ";
	write_code_data(os, code_data);
	os << endl;
	calculate_turning_cycles(code_data, *Cycle, num_left_cycles, num_cycles);

	matrix<int>& cycle = *Cycle;
	
    os << "% left turning cycles:";   
	for (int i=0; i<num_left_cycles; i++)
	{
		os << "\n% cycle " << i << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			os << cycle[i][j] << " ";
	}
	
    os << "\n% right turning cycles:";
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
		os << "\n% cycle " << i << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			os << cycle[i][j] << " "	;
	}
	
	os << endl;
	delete Cycle;
	
	/* controls for coordinate output */
	int num_decimal_points = 3;
	int output_field_width = 12;
	output.setf(ios::fixed,ios::floatfield);
	output.precision(num_decimal_points);

	os << "\nbeginfig(fignum);" << endl;
	os << "fignum:=fignum+1;" << endl;
	os << "numeric u,d;" << endl;
	float unit_points = mp_control.unit_size;
	os << "u=" << unit_points/100 << "pt;" << endl;
	os << "d=" << mp_control.disc_size << "u;" << endl;
	os << "path p[];" << endl;
	os << "pickup pencircle scaled " << mp_control.pen_size*0.5 << "pt;" << endl;

	/* vertex coordinates are numbered z<vertex> */
	for (int i=0; i< num_type123_vertices; i++)
		os << "z" << i << "=(" << setw(output_field_width) << vcoords[i][0] << "u," << setw(output_field_width) << vcoords[i][1] << "u);" << endl;
	
	os << endl;
	

	/* create the vertex sequence for writing the metapost path using the terminating 
	   crossing information in code_data for each edge.  The sequence starts with the 
	   type 2 vertex preceeding crossing zero, which is the first type 2 vertex in the 
	   numbering scheme.	       
	*/
	vector<int> vertex_sequence(2*num_edges);
	for (int i=0; i< num_edges; i++)
	{
		vertex_sequence[2*i] = num_crossings+(i%num_edges); //the type_2_vertex at the midpoint of edge i
		vertex_sequence[2*i+1] = term_crossing[i];
	}
	

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "write_metapost: vertex_sequence: ";
    for (int i=0; i< 2*num_edges; i++)
		debug << vertex_sequence[i] << ' ';
	debug << endl;
}

	if (ONE_METAPOST_PATH)
	{
		/*  This approach can produce inflexion points that are difficult to remove */
		for (int i=0; i< num_components; i++)
		{
			os << "p" << i+1 << " = ";
			for (int j=first_edge_on_component[i]; j< first_edge_on_component[i]+num_component_edges[i]; j++)
				os << "z" << vertex_sequence[2*j] << "..z" << vertex_sequence[2*j+1] << "..";
			os << "z" << vertex_sequence[2*first_edge_on_component[i]] << "..cycle;" << endl;
		}
	}
	else
	{
		for (int i=0; i< num_components; i++)
		{
			os << "p" << 2*i+1 << " = ";
			for (int j=first_edge_on_component[i]; j< first_edge_on_component[i]+num_component_edges[i]; j++)
				os << "z" << vertex_sequence[2*j] << "..z" << vertex_sequence[2*j+1] << "..";
			os << "z" << vertex_sequence[2*first_edge_on_component[i]] << ";" << endl;
		}
		for (int i=0; i< num_components; i++)
		{
			os << "p" << 2*i+2 << " = ";
			for (int j=first_edge_on_component[i]; j< first_edge_on_component[i]+num_component_edges[i]; j++)
				os << "z" << vertex_sequence[2*j] << "..z" << vertex_sequence[2*j+1] << "..";
			os << "{direction 0 of p" << 2*i+1 << "}z" << vertex_sequence[2*first_edge_on_component[i]] << ";" << endl;
		}
	}
	
	if (head == -1)
	{
		for (int i=0; i< num_components; i++)
		{
			if (mp_control.draw_oriented)
			{
				if (ONE_METAPOST_PATH)
					os << "drawarrow subpath(0,0.5) of p" << i+1 << ";" << endl;
				else
					os << "drawarrow subpath(0,0.5) of p" << 2*i+2 << ";" << endl;
			}
			
			if (ONE_METAPOST_PATH)
				os << "draw p" << i+1 << ";" << endl;
			else
				os << "draw p" << 2*i+2 << ";" << endl;
		}
	}
	else
	{	
		if (ONE_METAPOST_PATH)
		{
			os << "drawarrow subpath(0,0.75) of p1;" << endl;
			os << "draw subpath(0.75," << 2*semi_arc-1 << ".5) of p1;" << endl;
		}
		else
		{
			os << "drawarrow subpath(0,0.75) of p2;" << endl;
			os << "draw subpath(0.75," << 2*semi_arc-1 << ".5) of p2;" << endl;
		}

		
		if (mp_control.draw_shortcut)
		{
			if (ONE_METAPOST_PATH)
			{
				os << "draw subpath(0,0.5) of p1 dashed " << (mp_control.dash_with_dots? "withdots" : "evenly") << ";" << endl;
				os << "draw subpath(" << 2*semi_arc-1 << ".5," << 2*num_component_edges[0]+1 << ") of p1 dashed " 
				   << (mp_control.dash_with_dots? "withdots" : "evenly") << ";" << endl;
			}
			else
			{
				os << "draw subpath(0,0.5) of p2 dashed " << (mp_control.dash_with_dots? "withdots" : "evenly") << ";" << endl;
				os << "draw subpath(" << 2*semi_arc-1 << ".5," << 2*num_component_edges[0]+1 << ") of p2 dashed " 
				   << (mp_control.dash_with_dots? "withdots" : "evenly") << ";" << endl;
			}
		}

		for (int i=1; i< num_components; i++)
		{
			if (mp_control.draw_oriented)
			{
				if (ONE_METAPOST_PATH)				
					os << "drawarrow subpath(0,0.5) of p" << i+1 << ";" << endl;
				else
					os << "drawarrow subpath(0,0.5) of p" << 2*i+2 << ";" << endl;
			}
			
			if (ONE_METAPOST_PATH)
				os << "draw p" << i+1 << ";" << endl;
			else
				os << "draw p" << 2*i+2 << ";" << endl;
		}
	}
	
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "write_metapost: drawing crossing features" << endl;

	if (!DRAW_IMMERSION_ONLY)
	{
		for (int i=0; i< num_crossings; i++)
		{
	
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "write_metapost:   crossing " << i;
			
			if (code_table[LABEL][i] == VIRTUAL)
			{
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << " is virtual" << endl;
		
				os << "draw fullcircle scaled d shifted z" << vertex_sequence[2*first_edge[i]+1] << ";" << endl;
			}
			else if (code_table[LABEL][i] != FLAT)
			{
				os << "fill fullcircle scaled d shifted z" << vertex_sequence[2*first_edge[i]+1] << " withcolor 1white;" << endl;
				
				/* draw the over-arc back in, if the label is positive we want the naming edge
				   otherwise the other edge
				*/
				int edge;
				if (code_table[LABEL][i] == POSITIVE)	
				{
	
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << " has positive label, ";
				
					edge = 2*i;
				}
				else
				{
	
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << " has negative label, ";
				
					edge = code_table[OPEER][i];
				}
					
				int component;
				if ( edge%2 )
					component = code_table[COMPONENT][(edge-1)/2];
				else
					component = code_table[COMPONENT][edge/2];
					
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "over-arc starts on edge " << edge << " on component " << component << endl;
	
				/* we are at the midpoint vertex of edge on path p at time 2*(edge-first_edge_on_component[component])
				   and at the terminating crossing at time 2*(edge-first_edge_on_component[component])+1
				*/
				os << "draw subpath(" << 2*(edge-first_edge_on_component[component]) << ".5,";
				if (ONE_METAPOST_PATH)
					os << 2*(edge-first_edge_on_component[component])+1 << ".5) of p" << component+1 << ";" << endl;
				else
					os << 2*(edge-first_edge_on_component[component])+1 << ".5) of p" << 2*component+2 << ";" << endl;
			}
			else
			{
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << " is flat " << endl;
			}
		}
	}
	
	if (mp_control.draw_labels)
	{
		for (int i=0; i< num_edges; i++)
			os << "label(btex $" << (mp_control.label_edges_from_one? i+1: i) << "$ etex, z" << vertex_sequence[2*i] << ");" << endl;
	}

	if (mp_control.label_vertices)
	{
		for (int i=0; i< num_type123_vertices; i++)
			os << "label(btex $z" << i << "$ etex, z" << i << ");" << endl;
			
		os << "draw (0," << miny << "u)--(0," << maxy << "u);" << endl;
		os << "draw (" << minx << "u,0)--(" << maxx << "u,0);" << endl;
		os << "fill fullcircle scaled 0.5d shifted (0," << miny << "u);" << endl;
		os << "fill fullcircle scaled 0.5d shifted (0," << maxy << "u);" << endl;
		os << "fill fullcircle scaled 0.5d shifted (" << minx << "u,0);" << endl;
		os << "fill fullcircle scaled 0.5d shifted (" << maxx << "u,0);" << endl;

		os << "label.lft(btex (0," << int(miny) << "u) etex,(0," << miny << "u));" << endl;
		os << "label.lft(btex (0," << int(maxy) << "u) etex,(0," << maxy << "u));" << endl;
		os << "label.bot(btex (" << int(minx) << "u,0) etex,(" << minx << "u,0));" << endl;
		os << "label.bot(btex (" << int(maxx) << "u,0) etex,(" << maxx << "u,0));" << endl;
	}

	if (mp_control.circle_packing)
	{
		for (int i=0; i< num_type123_vertices; i++)
			os << "draw fullcircle scaled " << 2*vertex_radius[i] << "u shifted z" << i << ";" << endl;
	}

	os << "endfig;" << endl;
}

/* calculate_turning_cycles looks for left and right turning cycles in the given code_data and returns
   false if it is unable to calculate a set of cycles that may be realized.  From Eulers formula a
   realizable diagram cannot have more than num_crossings+2 turning cycles.  If the search therefore 
   exceeds this number of cycles the code_data cannot be realizable and the function returns false.
   If the search completes without reaching this limit the function returns true but this does not indicate
   that the code_data is necessarily realizable.  Note that the num_crossings+2 limit can only be breached
   when looking for right turning cycles, given that we look for left turning cycles first.
*/
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
{

if (draw_control::DEBUG >= draw_control::DETAIL)
{
	debug << "calculate_turning_cycles: code data: ";
	write_peer_code (debug, code_data);
	debug << endl;
	print(code_data,debug,"calculate_turning_cycles: ");	
}

	matrix<int>& code_table = code_data.code_table;
	vector<int> term_crossing = code_data.term_crossing;
	vector<int> orig_crossing = code_data.orig_crossing;
	
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;

	num_cycles = 0;
	
	/* First look for left turning cycles */
	for (int i=0; i<2*num_crossings; i++)
	{
		
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "calculate_turning_cycles: edge = " << i;

    	/* does edge i already appear in a cycle ? */
		bool found = false;
		for (int j=0; j<num_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
				if (abs(cycle[j][k]) == i)
//				if (cycle[j][k] == i)
				{

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << " found in left turning cycle " << j << endl;

					found = true;
					break;
				}
			}
		}

		if (!found)
		{

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << " not found in current left turning cycles "<< endl;

			/* start a new cycle */
			int column = 1;

			/* we always traverse odd edges backwards */
			int edge = (i%2? -i: i);
//			int edge = i;
			cycle [num_cycles][column++] = edge;
			bool complete = false;

			/* a cycle cannot be longer than num_edges so we check that we do indeed
			   cycle within that limit.  This is used to check that the component map
			   within the code_table is valid, since an unrealizable component map can
			   result in infinite loops being calculated in the terminating and 
			   originating edges of the crossings where we never return to the start
			   of a cycle.
			*/
			for (int j=0; !complete && j< num_edges; j++)
			{
				/* determine which vertex the current edge takes us to and 
	   			the next edge we turn onto */
				if (edge % 2) // edge is odd (and negative)
				{
					int vertex = orig_crossing[-edge];
//					int vertex = orig_crossing[edge];
					if (code_table[TYPE][vertex] == TYPE1)
		    			edge = code_table[EVEN_ORIGINATING][vertex];
					else
						edge = -code_table[ODD_TERMINATING][vertex];
//						edge = code_table[ODD_TERMINATING][vertex];

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}
				else // edge is even (and positive)
				{
					int vertex = term_crossing[edge];
					if (code_table[TYPE][vertex] == TYPE1)
						edge = -code_table[ODD_TERMINATING][vertex];
//						edge = code_table[ODD_TERMINATING][vertex];
					else
		    			edge = code_table[EVEN_ORIGINATING][vertex];

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;
				}

    
				if (edge == cycle[num_cycles][1])
				{
					complete = true;
					cycle[num_cycles][0] = column-1;
					num_cycles++;
				}
				else
				{
					cycle[num_cycles][column++] = edge;
				}				
			}
			
			if (!complete)
			{
				/* we've encounterd an infinte loop */
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "calculate_turning_cycles: exceeded the maximum possible length of a turning cycle in a realizable code" << endl;
				return false;
			}
		}
	}	

	/* record the number of left cycles */
	num_left_cycles = num_cycles;
		
if (draw_control::DEBUG >= draw_control::DETAIL)
{
    debug << "calculate_turning_cycles: number of left turning cycles = " << num_left_cycles;
	for (int i=0; i<num_left_cycles; i++)
	{
		debug << "\ncalculate_turning_cycles: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " ";
	}
	debug << endl;
}
		
	/* Now look for right turning cycles */

	for (int i=0; i<2*num_crossings; i++)
	{

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "calculate_turning_cycles: edge = " << i;

		/* does edge i already appear in a right cycle ? */
		bool found = false;
		for (int j=num_left_cycles; j<num_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
				if (abs(cycle[j][k]) == i)
//				if (cycle[j][k] == i)
				{
					
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << " found in right turning cycle " << j << endl;
    
    					found = true;
					break;
				}
			}
		}

		if (!found)
		{

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << " not found in current right turning cycles "<< endl;

			/* check we've not exceeded the maximum number of turning cycles */
			if (num_cycles == num_crossings+2)
			{
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "calculate_turning_cycles: exceeded the maximum possible number of turning cycles in a realizable code" << endl;
				return false;
			}
			
			/* start a new cycle */
			int column = 1;
			
			/* we always traverse odd edges backwards */
			int edge = (i%2? -i: i);
//			int edge = i;
			cycle [num_cycles][column++] = edge;
			bool complete = false;

			do
			{
				/* determine which vertex the current edge takes us to and 
	   			the next edge we turn onto */
				if (edge % 2) // edge is odd (and negative)
				{
					int vertex = orig_crossing[-edge];
//					int vertex = orig_crossing[edge];
					if (code_table[TYPE][vertex] == TYPE1)
						edge = -code_table[ODD_TERMINATING][vertex];
//						edge = code_table[ODD_TERMINATING][vertex];
					else
		    			edge = code_table[EVEN_ORIGINATING][vertex];

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}
				else // edge is even (and positive)
				{
					int vertex = term_crossing[edge];
					if (code_table[TYPE][vertex] == TYPE1)
		    			edge = code_table[EVEN_ORIGINATING][vertex];
					else
						edge = -code_table[ODD_TERMINATING][vertex];
//						edge = code_table[ODD_TERMINATING][vertex];

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}

				if (edge == cycle[num_cycles][1])
				{
					complete = true;
					cycle[num_cycles][0] = column-1;
					num_cycles++;
				}
				else
				{
					cycle[num_cycles][column++] = edge;
				}				
			} while(!complete);			
		}
	}

if (draw_control::DEBUG >= draw_control::DETAIL)
{
    debug << "calculate_turning_cycles: number of right turning cycles = " << num_cycles - num_left_cycles;
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
		debug << "\ncalculate_turning_cycles: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " "	;
	}
	debug << endl;
}
	return true;
}

void print (metapost_control& mp_control, ostream& os, string prefix)
{
	os << prefix << "rotate = " << (mp_control.rotate? "true": "false") << endl;
	os << prefix << "explicit_rotation_centre = " << (mp_control.explicit_rotation_centre? "true": "false") << endl;
	os << prefix << "implicit_rotation_centre = " << (mp_control.implicit_rotation_centre? "true": "false") << endl;
	os << prefix << "rotation_degrees = " << mp_control.rotation_degrees << endl;
	os << prefix << "rotation_centre_x = " << mp_control.rotation_centre_x << endl;
	os << prefix << "rotation_centre_y = " << mp_control.rotation_centre_y << endl;
	os << prefix << "rotation_centre_z = " << mp_control.rotation_centre_z << endl;
	os << prefix << "unit_size = " << mp_control.unit_size << endl;
	os << prefix << "disc_size = " << mp_control.disc_size << endl;
	os << prefix << "pen_size = " << mp_control.pen_size << endl;
	os << prefix << "infinite_cycle = " << mp_control.infinite_cycle << endl;
	os << prefix << "dash_with_dots = " << (mp_control.dash_with_dots? "true": "false") << endl;
	os << prefix << "knotoid = " << (mp_control.knotoid? "true": "false") << endl;
	os << prefix << "draw_immersion_only = " << (mp_control.draw_immersion_only? "true": "false") << endl;
	os << prefix << "draw_labels = " << (mp_control.draw_labels? "true": "false") << endl;
	os << prefix << "draw_oriented = " << (mp_control.draw_oriented? "true": "false") << endl;
	os << prefix << "draw_shortcut = " << (mp_control.draw_shortcut? "true": "false") << endl;
	os << prefix << "label_vertices = " << (mp_control.label_vertices? "true": "false") << endl;
	os << prefix << "circle_packing = " << (mp_control.circle_packing? "true": "false") << endl;
}

/* The function hyperbolic_representation produces the Euclidean centres and radii
   for a set of circles to be drawn in a representation of the disc model of the
   hyperbolic plane.

   For each circle, we determine the Mobius transformation that takes the hyperbolic 
   centre to zero.  There are two points of intersection, X and -X, of the image of
   circle under this trandformation and the x-axis, and the hyperbolic distance between
   0 and X is the hyperbolic radius.
   
   From the formula log (1+t_2)(1-t_1)/(1-t_2)(1+t_1) for hyperbolic distance we can
   evaluate the Euclidean distance of X from zero as |X| = |(e^{r_h} - 1)/(e^{r_h}+1)|.
   
   The inverse of the Mobius transformation then takes X and -X to points z_1 and z_2.
   The Euclidean centre of the circle is the point (z_1+z_2)/2 and the Euclidean radius
   is |z_1-z_2|/2
   
   Note that the hyperbolic radii produced by the circle packing are recorded in the
   form e^{-2h} for the hyperbolic radius h.  The factor of 2 is based on the choice of 
   definition for the hyperbolic metric but note that 
   
   (e^{-d} - 1)/(e^{-d}+1) = -(e^{d} - 1)/(e^{d}+1)
   
   so we may calculate the value of |X| directly from radius[i].
   
*/
void hyperbolic_representation(matrix<double>& centre, vector<double>& radius)
{
	int num_circles = radius.size();
	
	for (int i=0; i< num_circles; i++)
	{
		/* set e = e^{r_h} */
//		double e = exp(radius[i]);

		/* set t to be the modulus of the hyperbolic centre */
//		double t = sqrt(centre[i][0]*centre[i][0]+centre[i][1]*centre[i][1]);
		
		complex<double> a;
		complex<double> b;                     
		complex<double> z(centre[i][0],centre[i][1]);
		
		mobius_transform_to_origin(a,b,z);
		
//		double x = (e - 1)/(e + 1)/2;
		double x = abs((radius[i] - 1)/(radius[i] + 1));
		complex<double> z1 = inv_mobius_transformation(a,b, complex<double>(x,0));
		complex<double> z2 = inv_mobius_transformation(a,b, complex<double>(-x,0));

if (draw_control::DEBUG >= draw_control::DETAIL)
{
	debug << "hyperbolic_representation: circle " << i << ", centre " << z << ", radius " << radius[i] << endl;
	debug << "hyperbolic_representation: Mobius transformation parameters a = " << a << ", b = " << b << endl;
	debug << "hyperbolic_representation: Euclidean value of radius (X) = " << x << endl;
	debug << "hyperbolic_representation: pre-image of X = " << z1 << ", pre-image of -X = " << z2 << endl;
}    
		
		radius[i] = abs(z1-z2)/2*metapost_hyperbolic_scale_factor;

		z1 += z2;
		z1 /= 2;
		
		centre[i][0] = z1.real()*metapost_hyperbolic_scale_factor;
		centre[i][1] = z1.imag()*metapost_hyperbolic_scale_factor;
		

		/* we always evaluate tQ, even if the hyperbolic centre is the origin
		  because we need it for the Euclidean radius
		
		double tQ = e + e*t - 1 + t;
		double d  = e + e*t + 1 - t;
		tQ /= d;

		if (t != 0)
		{
			
			double tP = 1 + t - e + e*t; 
			d =         1 + t + e - e*t;
			tP /= d;
			
			double te = (tP + tQ)/2;
			
//			scale the hyperbolic centre to get the Euclidean centre
			double f = te/t;
			
			centre[i][0] *= f;
			centre[i][1] *= f;
		}
		*/
		
	}

if (draw_control::DEBUG >= draw_control::DETAIL)
{
    debug << "hyperbolic_representation: hyperbolic vertex_coordinates: " << endl;
    for (int i=0; i< num_circles; i++)
    {
		debug << "hyperbolic_representation:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << centre[i][j] << ' ';
		debug << endl;
	}	
    debug << "hyperbolic_representation: hyperbolic vertex_radius: " << endl;
    for (int i=0; i< num_circles; i++)
		debug << "hyperbolic_representation:   vertex " << i << ": " << radius[i] << endl;
}

}

/* The function mobius_transform_to_origin evaluates in a and b the parameters
   for the Mobius transformation that takes z to the origin
*/
void mobius_transform_to_origin(complex<double>& a, complex<double>& b, complex<double> z)
{
	double r_squared = 1/(1-norm(z)); // norm is the square of abs
	double r = sqrt(r_squared);
	double s = sqrt(r_squared - 1);
	b = complex<double>(s,0);
	double theta = 355/113 - arg(z);
	a = complex<double>(r*cos(theta), r*sin(theta));
}

/* the function mobius_transformation evaluates the image of z under the Mobius transformation
   that fixing the unit disc determined by a and b
*/
complex<double> mobius_transformation(complex<double> a, complex<double> b, complex<double> z)
{
	return (a*z+b)/(conj(b)*z + conj(a));
}

/* the function mobius_transformation evaluates the image of z under the inverse of the
   Mobius transformation fixing the unit disc determined by a and b
*/
complex<double> inv_mobius_transformation(complex<double> a, complex<double> b, complex<double> z)
{
	return (b-conj(a)*z)/(conj(b)*z - a);
}

