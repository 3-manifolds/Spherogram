/**************************************************************************
draw generates metapost code to draw diagrams from labelled immersion codes.  
It is currently able to draw classical and virtual knots and knotoids.  It 
was motivated by the requirement to draw knotoids.

The algorithm is taken from the drawing routine in Knotscape, by Jim Hoste and Morwen 
Thistlethwaite, and in particular Ken Stephenson's implementation of Thurston's circle 
packing algrithm.

My code provides the environment to read immersion codes from an input stream, then 
applies the algorithms of Knotscape to produce a sequence of coordinates at key 
points as we traverse the immersion.  Armed with these coordinates it is a simple 
matter to produce metapost output that can draw the immersion with the appropriate 
over, under and virtual crossings, knotoid shortcuts, and so on.

Specifically, the following key aspects have been taken from Knotscape.  I am very 
grateful to the authors of Knotscape for making this code available in the public
domain and fully acknowledge all rights they may have in that code and it's methods.

    triangulate.c - I have written code to produce the same triangulation as described in 
    this Knotscape module starting from a labelled immersion code.  The output of my code 
    is a file with the same file format as that produced by the Knotscape version.
    
    ken.c - this module has been used as is from Knotscape.  I create a programme called
    circle_pack by modifying the beginning of the Knotscape module draw_main.c as follows

		//triang(argv[1], argv[2]);
		//ken(argv[2], argv[3]);
		ken(argv[1], argv[2]);
		exit(0);

	then from the Knotscape distribution I run 'make draw', and rename the resultant 
	binary circle-pack and move to directory containing my draw executable.  The above 
	procedure isolates the circle packing algorithm from the rest of Knotscape.
	
	nodeseq.c - this module has been re-written to take into account the fact that the 
	sequence of vertices I am interested in is determined by a labelled immersion code;
	also, the use of metapost means I do not need all of the coordinates used by Knotscape.
	Thus, my nodeseq is a form of the corresponding Knotscape model optimized for immersion codes, 
	however, as with the other modules above, the algorithms (in particular the coordinate 
	scaling algorithm) are all taken from Knotscape.
	
	badness - I have incorporated the same assessment of a set of coordinates as used by
	Knotscape via the badness function.
	
	In addition to the above, wherever the Knotscape original modules call auxiliary functions
	I have written corresponding functions in support of my versions.  These are also 
	acknowledged effctively to have been taken directly from Knotscape.

                  Coding started 21st April 2010
 
   Version 1: support for Knotoids
   Version 2: added immersions and rotation
   Version 3: added labelled peer codes and Metapost control qualifiers
      
**************************************************************************/
using namespace std;
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <stdio.h>
#include <valarray>
#include <ctime>
#include <csignal>
#include <list>
#include <iomanip>
#include <vector>

/******************** Global Variables **********************/
string 		version  = "3.0";
   

extern ofstream    debug;
ofstream 	output;
ifstream    input;

int printf_bool=0; // debug for circle_pack.c
int DRAW_IN_HYPERBOLIC_DISC = 0;//used by the C module circle_pack
bool DRAW_IMMERSION_ONLY = false;
bool ONE_METAPOST_PATH = false;

#include <util.h>
#include <matrix.h>
#include <draw.h>

/********************* Function prototypes ***********************/
bool get_next_input_string(string input_file,string& input_string, string& title);
void help_info();
void debug_help();
void set_programme_long_option(char* cptr, string source, metapost_control& mp_control);
void set_programme_short_option(char* cptr, metapost_control& mp_control);
bool debug_setup(char* argv_ptr);
void check_debug_option_parameters(char* pptr, string option);
void get_word (ifstream& input, string& buffer, string& title);
void read_code_data (generic_code_data& code_data, string input_string);
void triangulate (generic_code_data& code_data, matrix<int>& cycle, int num_cycles, int num_left_cycles, int infinite_region);
double badness(string vertex_file, generic_code_data& code_data);
//double badness(string vertex_file, int vertex_sequence_length);
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control);
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles);
bool valid_knotoid_input(generic_code_data& code_data, vector<int>& shortcut_crossing);

extern "C" int circle_pack (char* inputfile, char* outputfile);

int draw_control::DEBUG = draw_control::OFF;

void sigfpe_handler (int sig) 
{
	cout << "Error, draw programme received SIGFPE!" << endl;
	exit(1);
}


/* scale factor applied to circle_pack output coordinates */
int metapost_coordinate_scale_factor = 3000;

/* additional scaling when drawing in a representation of the hyperbolic disc */
int metapost_hyperbolic_scale_factor = 1; 

/* badness controls */

/* If the minimum distance between successive vertices in the vertex sequence
   is less than badness_threshold, we look for an alternative infinite turning
   cycle, unless a specific infinite cycle has been provided.
*/
double badness_threshold = 100; 

/* candidate turning cycles must have minimum_infinite_cycle_length to 
   be considered as an infinite cycle
*/
int minimum_infinite_cycle_length = 3;

/******************* Main Function ************************/
int main (int argc, char* argv[])
{
    string 	input_file;
    string 	output_file;
	string 	input_string;
	string 	title;
    bool    input_file_provided = false;
    bool    output_file_provided = false;

	char* triangulation_output_file = c_string("__triangulate.out");
	char* circlepack_output_file =  c_string("__circlepack.out");

    
	/* default_metapost_control is initialized with the default settings */
    metapost_control default_metapost_control;

	signal(SIGFPE, sigfpe_handler); // register the SIGFPE handler

    /* Determine command line options and whether input is to
       come from the screen or a file */
    if (argc == 1)
	{
		help_info();
	}
	else
    {
		for (int i=1; i < argc; i++)
		{
			if (*(argv[i]) == '-')
			{
				if (*(argv[i]+1) == '-')
					set_programme_long_option(argv[i], "command line",default_metapost_control);
				else 
					set_programme_short_option(argv[i], default_metapost_control);
			}
			else if (!input_file.length()) 
			{
				input_file = argv[i];
	    		input_file_provided = true;
			}
			else if (!output_file.length()) 
			{
				output_file = argv[i];
		    	output_file_provided = true;
			}
			else
			{
		    	cout << "Usage draw --<task> [-<options>][<infile>[<outfile>]], type draw -h for help.\n";
		    	exit(0);
			}
		}	
    }

//if (draw_control::DEBUG >= draw_control::DETAIL)
//	print_prog_params(debug, "detail");

	if (input_file_provided)
	{
    	input.open(input_file.c_str());
    	if(!input)
    	{
			cout << "\nError opening input file\n";
			exit(0);
    	}	
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "\ndraw::main: input file: " << input_file << "\n";

		/* Read any options included in the input file */
		string next_line;
		while(getline(input,next_line))
		{
			char* line_buf = c_string(next_line);
			char* cptr = strchr(line_buf,';');
			if (cptr)
	    		*cptr = '\0';

			if (strlen(line_buf))
			{
				cptr = strchr(line_buf,'[');
				if (cptr && !strchr(line_buf,'\\') && !strchr(line_buf,'/'))
				{		
					/* this is an options line not a line containing a peer code
					   make sure there's an end to the option definitions on this line 
					*/
					if (!strchr(line_buf,']'))
					{
						cout << "\nError in " << input_file << ": programme options must be enclosed within []\n";
						exit(0);
					}
					do
					{
						cptr++;
						/* skip whitespace */
						while (isspace(*cptr))
							cptr++;
						set_programme_long_option(cptr++, "input file",default_metapost_control);
						while (*cptr != ',' && *cptr != ']')
							cptr++;
					} while (*cptr != ']');
				}
			}
			
			delete[] line_buf;
    	}
	}

    /* prepare output file */
	if (!output_file_provided)
		output_file = "draw.out";
	
    output.open (output_file.c_str());

    if (!output)
    {
		cout << "\nError opening output file " << output_file << endl;
		exit(0);
    }
    else
    {
		output << "% Output from draw v" << version << "\n";

		if (input_file_provided)
		{
			output << "% Input file: " << input_file << "\n";
		}
		
		output << "\ninput boxes\n" << endl;
		output << "\nprologues:=3;\n" << endl;		

		output << "numeric fignum; fignum=1;" << endl;
    }


    if (!input_file_provided)
    {
	
		if (argc > 1)
		{
			cout << "\nThis is A.Bartholomew's draw programme, v" << version << endl;
		}
		
		cout << "\n\nThe programme will produce a metapost file for drawing the knot or knotoid";
		cout << " specified by a labelled immersion code.\n";
		
		cout <<"\nType help at the input prompt to view the help screens, type q to exit.\n";

		cout << "\nEnter the labelled immersion code";		
    }
	
try {

    if (input_file_provided)
    {
		/* reset input ready for read */
		input.clear(); // state flags
		input.seekg(0);
	}

	while (get_next_input_string(input_file,input_string, title))
	{
		metapost_control mp_control(default_metapost_control);
		
		/* first remove any unwanted qualifiers from the input string */
		string::size_type pos = input_string.find('{');
		if (pos != string::npos)
		{
			char* qualifier_buf = c_string(input_string.substr(pos));
			char* cptr = qualifier_buf;

if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main: input_string qualifiers = " << qualifier_buf << endl;


			if (!strchr(qualifier_buf,'}'))
			{
				cout << "\nError in " << input_string << ": qualifiers must be enclosed within {}\n";
			}
			do
			{
				cptr++;
				/* skip whitespace */
				while (isspace(*cptr))
					cptr++;
				set_programme_long_option(cptr++, "qualifiers",mp_control);
				while (*cptr != ',' && *cptr != '}')
					cptr++;
			} while (*cptr != '}');

			delete qualifier_buf;
			input_string = input_string.substr(0,pos);
		}

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "draw::main: metapost_control for input_string:" << endl;
    print(mp_control, debug, "draw::main:   ");
}

		generic_code_data code_data;
		matrix<int>* Cycle;
		int num_cycles = 0;
		int num_left_cycles;
		int num_crossings;
		int num_edges;

		read_code_data(code_data,input_string);		
		
		/* There cannot be more than num_crossings + 2 cycles, and they cannot 
		   be more than num_edges in length, we add 1 to the length of a 
	   	   cycle as we store the actual length in column 0 */
		num_crossings = code_data.num_crossings;
		num_edges = 2*num_crossings;
		Cycle = new matrix<int>(num_crossings+2, num_edges+1);
		calculate_turning_cycles(code_data, *Cycle, num_left_cycles, num_cycles);

		if (code_data.head != -1)
		{
			vector<int> shortcut_crossing(num_crossings);
			if (!valid_knotoid_input(code_data,shortcut_crossing))
			{
				
				cout << "Input is not a valid knotoid, indicated shortcut does not pass under the rest of the diagram" << endl;
				continue;
				
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main: invalid knotoid input detected" << endl;
			}
			else
			{
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main: valid knotoid input detected" << endl;
			}
		}

		matrix<int>& cycle = *Cycle;

		int infinite_region = -1;
		if (mp_control.knotoid)
		{
			/* Select the initial choice of infinite region, this will be 
			   one bounded of maximal length that includes edge zero
			*/
			for (int i=1; i< num_cycles; i++)
			{
				int length = cycle[i][0];
				
				bool contains_zero = false;
				
				for (int j=1; j<= length; j++)
				{
					if (cycle[i][j] == 0)
					{
						contains_zero = true;
						break;
					}
				}
				
				if (contains_zero)
				{
					if (infinite_region != -1)
					{
						if (length > cycle[infinite_region][0])
							infinite_region = i;
					}
					else
						infinite_region = i;
				}
			}
		}
		else if (mp_control.infinite_cycle != -1)
		{
			infinite_region = mp_control.infinite_cycle;
		}
		else
		{
			infinite_region = 0;
			
			/* Select the initial choice of infinite region, this will be 
			   one bounded by a turning cycle of maximal length
			*/
			for (int i=1; i< num_cycles; i++)
			{
				if (cycle[i][0] > cycle[infinite_region][0])
					infinite_region = i;
			}
		}
			
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main: initial infinite_region = " << infinite_region << endl;
				

		triangulate(code_data, cycle, num_cycles, num_left_cycles,infinite_region);

//		circle_pack ("__triangulate.out", "__circlepack.out");
		circle_pack (triangulation_output_file, circlepack_output_file);
		
//		double badness_rating = badness("__vertex_coords.out", 2*num_edges);
		double badness_rating = badness("__circlepack.out", code_data);
			
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main: badness_rating when using cycle " << infinite_region << " as infinite_region = " << badness_rating << endl;
			
		if (!mp_control.knotoid && badness_rating < badness_threshold && mp_control.infinite_cycle == -1)
		{

if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main: badness_rating less than " << badness_threshold << ", looking for other candidates" << endl;

			int new_infinite_region = infinite_region;
			
			for (int i=0; i< num_cycles; i++)
			{
				if (i != infinite_region && cycle[i][0] >= minimum_infinite_cycle_length)
				{
					triangulate(code_data, cycle, num_cycles, num_left_cycles,i);
//					circle_pack ("__triangulate.out", "__circlepack.out");
					circle_pack (triangulation_output_file, circlepack_output_file);
				
					
//					double trial_badness = badness("__vertex_coords.out", 2*num_edges);
					double trial_badness = badness("__circlepack.out", code_data);

if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main:   badness_rating when using cycle " << i << " as infinite_region = " << badness_rating << endl;

					if (trial_badness > badness_rating)
					{
						badness_rating = trial_badness;
						new_infinite_region = i;
					}
				}
				else
				{
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main:   cycle " << i << " not a viable alternative candidate as infinite_region" << endl;
				}
			}

if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main: best infinite_region selection is cycle " << new_infinite_region << endl;

			triangulate(code_data, cycle, num_cycles, num_left_cycles,new_infinite_region);
//			circle_pack ("__triangulate.out", "__circlepack.out");
			circle_pack (triangulation_output_file, circlepack_output_file);
		}
		else
		{
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "draw::main: badness_rating acceptable" << endl;
		}
			
		write_metapost(output, code_data, title, mp_control);		
		
		delete Cycle;
	}
}
catch (bad_alloc)
{
	cerr << "ERROR! Out of memory, bad_alloc thown" << endl;
}
catch (...)
{
	cerr << "ERROR! Exception thrown that is not caught explicitly" << endl;
}


	output << "\nend" << endl;

	if (input_file_provided)
		input.close();
	
    output.flush();
    output.close();


if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug.close();

	delete[] triangulation_output_file;
	delete[] circlepack_output_file;

	return 0;
} /* End of main */

void set_programme_long_option(char* cptr, string source, metapost_control& mp_control)
{
	char  loc_buf[40];
	char* c1;
	char* c2;

	c1 = cptr;

	/* take out any leading -- */
	while (*c1 == '-')
		c1++;
	
	c2 = loc_buf;
	while (isalpha(*c1) || isdigit(*c1) || *c1 == '-' || *c1 == '_')
		*c2++ = *c1++;
    *c2 = '\0';
	
	if (!strcmp(loc_buf,"centre"))
	{

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: rotation centre read from " << source << endl;

		bool error = false;
		
		if (*c1 == '=')
		{
			c1++;
			
			if (*c1 == 'z')
			{				
				get_number(mp_control.rotation_centre_z,++c1);
				mp_control.implicit_rotation_centre = true;

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
	debug << "set_programme_long_option: rotation_centre_z read from " << source 
	      << ", rotation_centre_z = " << mp_control.rotation_centre_z << endl;
}	    
			}
			else
			{
				if (*c1 == '(')
				{
					get_number(mp_control.rotation_centre_x,++c1);
					while (isdigit(*c1))
						c1++;
						
					if (*c1 == ',')
					{
						get_number(mp_control.rotation_centre_y,++c1);
						
						while (isdigit(*c1))
							c1++;
						
						if (*c1 == ')')
						{
							mp_control.explicit_rotation_centre = true;
			
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
	debug << "set_programme_long_option: explicit rotation centre (" << mp_control.rotation_centre_x 
	      << "," << mp_control.rotation_centre_y << ")" << endl;
}	    
						}
						else
							error=true;
					}
					else
						error=true;
				}
				else
					error=true;				
			}
		}
		else
			error = true;
			
		if (error)
		{
			cout << "\nYou must specify the required rotation centre if you use the centre option." << endl;
			cout << "Specify the coordinates of the centre explicitly, e.g. centre=(123.45,67.89), or implicitly, e.g. centre=z21." << endl;
			exit(0);
		}
				
	}
	else if (!strcmp(loc_buf,"cycle"))
	{

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: setting infinite region turning cycle as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.infinite_cycle,++c1);
		}
		else
		{
			cout << "\nYou must specify a turning cycle if you use the cycle option, e.g. 'cycle=2'" << endl;
			exit(0);
		}
		
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
	debug << "set_programme_long_option: infinite_cycle read from " << source << ", infinite_cycle = " 
	      << mp_control.infinite_cycle << endl;
}	
	}
	else if (!strcmp(loc_buf,"dots"))
	{
    	mp_control.dash_with_dots = true;
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: dash_with_dots read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"disc-size"))
	{

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: setting disc_size as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.disc_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a unit size if you use the disc-size option, e.g. 'disc-size=20'" << endl;
			exit(0);
		}
		
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
	debug << "set_programme_long_option: disc_size read from " << source << ", disc_size = " 
	      << mp_control.disc_size << endl;
}	
	}
	else if (!strcmp(loc_buf,"hyperbolic"))
	{
		DRAW_IN_HYPERBOLIC_DISC = 1;
		
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: drawing diagram in the hyperbolic plane" << endl;
	
	}	
	else if (!strcmp(loc_buf,"knotoid"))
	{
    	mp_control.knotoid = true;
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: knotoid read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"labels"))
	{
    	mp_control.draw_labels = true;
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: draw_labels read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"oriented"))
	{
    	mp_control.draw_oriented = true;
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: draw_oriented read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"packing"))
	{
    	mp_control.circle_packing = true;
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: circle_packing read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"pen-size"))
	{

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: setting pen_size as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.pen_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a unit size if you use the pen-size option, e.g. 'pen-size=2'" << endl;
			exit(0);
		}
		
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
	debug << "set_programme_long_option: pen_size read from " << source << ", pen_size = " 
	      << mp_control.pen_size << endl;
}	
	}
	else if (!strcmp(loc_buf,"rotate"))
	{

		mp_control.rotate = true;
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: rotate read from " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.rotation_degrees,++c1);
		}
		else
		{
			cout << "\nYou must specify the required rotation in degrees if you use the rotate option, e.g. 'rotate=15'" << endl;
			exit(0);
		}
		
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
	debug << "set_programme_long_option: rotation_degrees read from " << source << ", rotation_degrees = " 
	      << mp_control.rotation_degrees << endl;
}	
	}
	else if (!strcmp(loc_buf,"shortcut"))
	{
    	mp_control.draw_shortcut = true;
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: draw_shortcut read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"unit"))
	{

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: setting unit as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.unit_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a unit size if you use the unit-size option, e.g. 'unit-size=30'" << endl;
			exit(0);
		}
		
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
	debug << "set_programme_long_option: unit_size read from " << source << ", unit_size = " 
	      << mp_control.unit_size << endl;
}	      
	}
	else if (!strcmp(loc_buf,"vertices"))
	{
    	mp_control.label_vertices = true;
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_long_option: label_vertices read from " << source << endl;
	}
}


void set_programme_short_option(char* cptr, metapost_control& mp_control)
{

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_short_option: provided with argument: " << cptr << endl;

	char* c1;

	if (strchr(cptr, 'H') && strchr(cptr, '!'))
	{
		if (strchr(cptr,'#'))
		{
			debug_help();
		}
		else
		{
			cout << "\nUsage draw [--<metapost-control>][-<option>][<infile>[<outfile>]]\n\n";

			cout << "<long-option> =\n";
			cout << "  centre=<centre>: specify the centre of rotation <centre> = (<x>,<y>)|z<n>\n";
			cout << "  cycle=<cycle>: specify the turning cycle to bound the infinite region\n";
			cout << "  dots: draw knotoid shortcuts dashed with dots rather than dashed evenly\n";
			cout << "  disc-size: set the unit multiplier to determine the crossing disc diameter\n";
			cout << "  knotoid: draw knotoids with the leg in the unbounded component of the immersion's complement\n";
			cout << "  labels: add edge labels to the diagram\n";
			cout << "  oriented: draw orientation for knots (orientation always shown for knotoids)\n";
			cout << "  packing: draw the underlying circle packing that determines the diagram\n";
			cout << "  pen-size: set the pen multiplier n to determine the pencircle scale n*0.5pt\n";
			cout << "  rotate=<degrees>: rotate the diagram anti-clockwise by the specified number of degrees\n";
			cout << "  shortcut: draw the shortcut if the immersion code is a knotoid\n";
			cout << "  unit: set the unit dimension nn to determine u=0.nn points\n";
			cout << "  vertices: label the vertices and draw coordinate axes\n";

			cout << "<short-option> =\n";
	    	cout << "  #: debug\n";
			cout << "  c=<centre>: same as centre\n";
			cout << "  C=<cycle>: same as cycle\n";
			cout << "  D: same as dots\n";
			cout << "  d: same as disc-size\n";
	    	cout << "  H!: this help screen\n";
	    	cout << "  #H!: display debug help screen\n";
			cout << "  i: draw immersion only\n";
			cout << "  k: same as knotoid\n";
			cout << "  l: same as labels\n";
			cout << "  L: start from one when labelling edges\n";
			cout << "  o: same as oriented\n";
			cout << "  O: create metapost with one cycled path for each component\n";
			cout << "  p: same as pen-size\n";
			cout << "  P: same as packing\n";
			cout << "  r=<degrees>: same as rotate\n";
			cout << "  s: same as shortcut\n";
			cout << "  u: same as unit\n";
			cout << "  v: same as vertices\n";
			cout << endl;
	    	exit(0);
		}
	}

	char* dptr = strchr(cptr, '#');
	if (dptr)
	{

    	/* establish a debug file */
    	debug.open ("draw.dbg"); 

    	if (!debug)
    	{
        	cout << "\nError opening debug file\n";
        	exit(0);
    	}
		else
			debug << "Debug information from draw version " << version << "\n\n";

if (!debug_setup(cptr))  // could probably be dptr, but the original code used cptr
{
	draw_control::DEBUG = draw_control::SUMMARY;
	debug << "set_programme_short_option: default debug options set" << endl;
}

	}

	c1=strchr(cptr, 'c');
	if (c1)
	{
		bool error = false;
		
		if (*++c1 == '=')
		{
			c1++;
			
			if (*c1 == 'z')
			{				
				get_number(mp_control.rotation_centre_z,++c1);
				mp_control.implicit_rotation_centre = true;

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_short_option: implicit rotation centre z" << mp_control.rotation_centre_z << endl;
	
			}
			else
			{
				if (*c1 == '(')
				{
					get_number(mp_control.rotation_centre_x,++c1);
					while (isdigit(*c1))
						c1++;
					if (*c1 == ',')
					{
						get_number(mp_control.rotation_centre_y,++c1);
						
						while (isdigit(*c1))
							c1++;
						
						if (*c1 == ')')
						{
							mp_control.explicit_rotation_centre = true;
			
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
	debug << "set_programme_short_option: explicit rotation centre (" << mp_control.rotation_centre_x 
	      << "," << mp_control.rotation_centre_y << ")" << endl;
}	    
						}
						else
							error = true;

				}
					else
						error=true;
				}
				else
					error=true;				
			}
		}
		else
			error = true;
			
		if (error)
		{
			cout << "\nYou must specify the required rotation centre if you use the c option." << endl;
			cout << "Specify the coordinates of the centre explicitly, e.g. -c=(123.45,67.89), or implicitly, e.g. -c=z21." << endl;
			exit(0);
		}
	}

	c1 = strchr(cptr, 'C');
	if (c1)
	{ 
		if (*++c1 == '=')
	    {
			get_number(mp_control.infinite_cycle,++c1);
	    }
		else
		{
			cout << "\nYou must specify a turning cycle if you use the I option, e.g. '-I=2'" << endl;
			exit(0);
		}
	
if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_short_option: set infinite_cycle = " << mp_control.infinite_cycle << endl;

	}

	if (strchr(cptr, 'D'))
	{
		mp_control.dash_with_dots = true;

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_short_option: set dash_with_dots" << endl;
	}

	c1 = strchr(cptr, 'd');
	if (c1)
	{ 
		if (*++c1 == '=')
	    {
			get_number(mp_control.disc_size,++c1);
	    }
		else
		{
			cout << "\nYou must specify a unit size if you use the d option, e.g. '-d=20'" << endl;
			exit(0);
		}

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_short_option: set disc_size = " << mp_control.disc_size << endl;

	}

	if (strchr(cptr, 'i'))
		DRAW_IMMERSION_ONLY = true;

	if (strchr(cptr, 'k'))
		mp_control.knotoid = true;

	if (strchr(cptr, 'l'))
		mp_control.draw_labels = true;

	if (strchr(cptr, 'L'))
		mp_control.label_edges_from_one = true;

	if (strchr(cptr, 'o'))
		mp_control.draw_oriented = true;

	if (strchr(cptr, 'O'))
		ONE_METAPOST_PATH = true;

	if (strchr(cptr, 's'))
		mp_control.draw_shortcut = true;

	if (strchr(cptr, 'h'))
		help_info();

	if (strchr(cptr, 'P'))
		mp_control.circle_packing = true;

    c1 = strchr(cptr, 'p');
	if (c1)
	{
	    if (*++c1 == '=')
	    {
			get_number(mp_control.pen_size,++c1);
	    }
		else
		{
			cout << "\nYou must specify a unit size if you use the p option, e.g. '-p=4'" << endl;
			exit(0);
		}

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_short_option: set pen_size = " << mp_control.pen_size << endl;
	}

	c1=strchr(cptr, 'r');
	if (c1)
	{
		mp_control.rotate = true;
	    if (*++c1 == '=')
	    {
			get_number(mp_control.rotation_degrees,++c1);

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_short_option: rotation of " << mp_control.rotation_degrees << " degrees" << endl;

	    }
		else
		{
			cout << "\nYou must specify the required rotation in degrees if you use the r option, e.g. '-r=15'" << endl;
			exit(0);
		}
	}

    c1 = strchr(cptr, 'u'); 
	if (c1)
	{
	    if (*++c1 == '=')
	    {
			get_number(mp_control.unit_size,++c1);
	    }
		else
		{
			cout << "\nYou must specify a unit size if you use the u option, e.g. '-u=12'" << endl;
			exit(0);
		}

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_short_option: set unit_size = " << mp_control.unit_size << endl;

	}

	if (strchr(cptr, 'v'))
		mp_control.label_vertices = true;

    c1 = strchr(cptr, 'x');
	if (c1)
	{
	    if (*++c1 == '=')
	    {
			get_number(printf_bool,++c1);
	    }
		else
		{
			cout << "\nYou must specify the circle packing debug level if you use the x option, e.g. '-x=3'" << endl;
			exit(0);
		}

if (draw_control::DEBUG >= draw_control::SUMMARY)
	debug << "set_programme_short_option: set printf_bool = " << printf_bool << endl;
	}

}

bool get_next_input_string(string input_file,string& input_string, string& title)
{
	bool success = true;

	if (input_file.length())
	{
		/* get next word from input */
		get_word(input, input_string, title);
	}
	else
	{
		do
		{
			cout << "\n\ninput: ";
			getline(cin, input_string); // has to be getline so there can be spaces in the string
		} while (input_string.length() == 0);
	}

	if (input_string == "exit" || input_string == "q" || input_string =="Q")
		success = false;
	else if (input_string == "help")
	{
	    help_info();
		exit(0);
	}
	else    	
	{

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "draw::get_next_input_string: got string " << input_string << endl;
}
	}	
	
	if (success)
	{ 
		cout << "\n\n";
		if (title.length() != 0)
			cout << title << endl;
		cout << input_string << endl;
	}
		
	return success;
}

void help_info()
{
	cout << "\nUsage draw [--<metapost-control>][-<option>][<infile>[<outfile>]]\n\n";

	cout << "<metapost-control> =\n";
	cout << "  centre=<centre>: specify the centre of rotation <centre> = (<x>,<y>)|z<n>\n";
	cout << "  cycle=<cycle>: specify the turning cycle to bound the infinite region\n";
	cout << "  dots: draw knotoid shortcuts dashed with dots rather than dashed evenly\n";
	cout << "  disc-size: set the unit multiplier to determine the crossing disc diameter, default 20\n";
	cout << "  knotoid: draw knotoids with the leg in the unbounded component of the immersion's complement\n";
	cout << "  labels: add edge labels to the diagram\n";
	cout << "  oriented: draw orientation for knots (orientation always shown for knotoids)\n";
	cout << "  pen-size: set the pen multiplier N to determine the pencircle scale n*0.5pt\n";
	cout << "  rotate=<degrees>: rotate the diagram anti-clockwise by the specified number of degrees\n";
	cout << "  shortcut: draw the shortcut if the immersion code is a knotoid\n";
	cout << "  unit: set the unit dimension N to determine u=N/100 points, default 20\n";
	cout << "  vertices: label the vertices and draw coordinate axes\n";


	cout << "<option> =\n";
	cout << "  c=<centre>: same as centre\n";
	cout << "  C=<cycle>: same as cycle\n";
	cout << "  D: same as dots\n";
	cout << "  d: same as disc-size\n";
	cout << "  h: this help screen\n";
	cout << "  i: draw immersion only\n";
	cout << "  k: same as knotoid\n";
	cout << "  l: same as labels\n";
	cout << "  o: same as oriented\n";
	cout << "  p: same as pen-size\n";
	cout << "  r=<degrees>: same as rotate\n";
	cout << "  s: same as shortcut\n";
	cout << "  u: same as unit\n";
	cout << "  v: same as vertices\n";
			cout << endl;
	exit(0);
}

/***************  Functions required for setting up draw specific debug **************/

void set_main_debug_option_parameter(char* pptr, string option);
void set_main_debug_option(char* start, char* end)
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
		set_main_debug_option_parameter(0,"draw");
//		cout << "\t\tvogel, boolean" << endl;
		return;
	}

	do
	{
		*c2++ = *c1++;
	} while (c1 <= end);
	
	*c2 = '\0';

	char* pptr = strchr(loc_buf,'{');
	if (pptr)
	{
		*pptr++ = '\0';
	}

	/* now, even if there are parameters, the first part of loc_buf 
	   is a C-string that identifies the option */
	
	if (!strcmp(loc_buf,"draw"))
	{
		draw_control::DEBUG = draw_control::SUMMARY;
		debug << "main::set_main_debug_option: setting debug option draw_control::DEBUG = draw_control::SUMMARY\n";		
		
		if (pptr)
		{
			/* check for any parameters */
			check_debug_option_parameters(pptr, "draw");
		}
	}

	debug.flush();

}

void set_main_debug_option_parameter(char* pptr, string option)
{
	if (option == "draw")
	{
		if (!pptr)
		{
			cout << "\t\tdraw{summary|1:basic|2:intermediate|3:detail|4}, integer: default 0=off, no parameters sets summary" << endl;
		}
		else
		{
			if (!strcmp(pptr,"summary") || !strcmp(pptr,"1") )
			{
				draw_control::DEBUG = draw_control::BASIC;
				debug << "main::set_main_debug_option_parameter: setting debug option draw_control::DEBUG = draw_control::SUMMARY\n";		
			}
			if (!strcmp(pptr,"basic") || !strcmp(pptr,"2") )
			{
				draw_control::DEBUG = draw_control::BASIC;
				debug << "main::set_main_debug_option_parameter: setting debug option draw_control::DEBUG = draw_control::BASIC\n";		
			}
			else if (!strcmp(pptr,"intermediate") || !strcmp(pptr,"3"))
			{
				draw_control::DEBUG = draw_control::INTERMEDIATE;
				debug << "main::set_main_debug_option_parameter: setting debug option draw_control::DEBUG = draw_control::INTERMEDIATE\n";		
			}
			else if (!strcmp(pptr,"detail") || !strcmp(pptr,"4"))
			{
				draw_control::DEBUG = draw_control::DETAIL;
				debug << "main::set_main_debug_option_parameter: setting debug option draw_control::DEBUG = draw_control::DETAIL\n";		
			}
		}
	}
}

void main_display_default_options()
{
	cout << "\t\tdraw{summary}" << endl;
}
