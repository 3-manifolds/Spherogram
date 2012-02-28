/************************************************************************
The function nodeseq is based on the corresponding Knotscape function and uses
the same scaling algorithm as Knotscape.  It reads the coordinates produced 
by the circle packing algorithm, re-scales the type 1 & 2 vertex
coordinates (immersion crossings and edge mid-points) and writes the scaled
type 1 and type 2 vertex coordinates to the local file __vertex_coords.out 
in the correct sequence as we trace around the immersion.  These coordinates
will be used used as the coordinate points for a metapost path by the function
write_metapost.
**************************************************************************/
using namespace std;
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>

/********************* External variables ***********************/

extern ofstream     debug;


#include <util.h>
#include <matrix.h>
#include <draw.h>

void nodeseq (generic_code_data& code_data)
{
	/* controls for coordinate output */
	int num_decimal_points = 3;
	int output_field_width = 12;

	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;
	int vertex_sequence_length = 2*num_edges; // this is the value used by knotscape
	
	ifstream triangulation;
	triangulation.open("__triangulate.out");
	int num_vertices;
	char ch;
	do
	{
		triangulation >> ch;
	} while (!isdigit(ch));
	triangulation.putback(ch);
	triangulation >> num_vertices;


if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "nodeseq: num_vertices = " << num_vertices << endl;

//	vector<int>& orig_crossing=code_data.orig_crossing;
	vector<int>& term_crossing=code_data.term_crossing;

	vector<int> vertex_sequence(vertex_sequence_length+1);
	for (int i=1; i<= num_edges; i++)
	{
		vertex_sequence[2*i-1] = term_crossing[i-1];
		vertex_sequence[2*i] = num_crossings+(i%num_edges); //the type_2_vertex at the midpoint of edge i
	}
	
	vertex_sequence[0] = vertex_sequence[vertex_sequence_length];	

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "nodeseq: vertex_sequence: ";
    for (int i=0; i< vertex_sequence_length+1; i++)
		debug << vertex_sequence[i] << ' ';
	debug << endl;
}
	
	/* type_12_vertex is an array of flags used to identify
	   those vertices that are type 1 or type 2.  It is used to 
	   determine the minimum and maximum x and y coordinates of
	   these vertex types so the coordinates may be scaled
	   appropriately.
	*/
	vector<int> type_12_vertex(num_vertices);
	for (int i=1; i<=vertex_sequence_length; ++i) 
		type_12_vertex[vertex_sequence[i]] = 1;

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "nodeseq: type_12_vertex: ";
    for (int i=0; i<= vertex_sequence_length; i++)
		debug << type_12_vertex[i] << ' ';
	debug << endl;
}

	float minx,maxx,miny,maxy;
	minx=maxx=miny=maxy=0;

	matrix<float> vertex_coordinates(num_vertices,2);
	ifstream coordinates;
	coordinates.open("__circlepack.out");
	
	if (!coordinates)
	{
		cout << "\nError opening __circlepack.out" << endl;
		exit(0);
	}
	
	/* determine the min/max X and Y coordinates of the type 1 and 2 vertices */
	for (int i=0; i< num_vertices; i++)
	{
		coordinates >> vertex_coordinates[i][0];
		coordinates >> vertex_coordinates[i][1];

		if (vertex_coordinates[i][0] < minx && type_12_vertex[i] == 1) 
			minx = vertex_coordinates[i][0];
		
		if (vertex_coordinates[i][0] > maxx && type_12_vertex[i] == 1) 
			maxx = vertex_coordinates[i][0];  
		
		if (vertex_coordinates[i][1] < miny && type_12_vertex[i] == 1) 
			miny = vertex_coordinates[i][1];  
		
		if (vertex_coordinates[i][1] > maxy && type_12_vertex[i] == 1) 
			maxy = vertex_coordinates[i][1];  
	}

if (draw_control::DEBUG >= draw_control::DETAIL)
{
    debug << "nodeseq: vertex_coordinates: " << endl;
    for (int i=0; i< num_vertices; i++)
    {
		debug << "nodeseq:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << vertex_coordinates[i][j] << ' ';
		debug << endl;
	}	
}
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "nodeseq: minx = " << minx << ", maxx = " << maxx << ", miny = " << miny << ", maxy = " << maxy << endl;

	/* scale so that picture fits nicely in canvas */
	float sc_scale = 600;
	float tr_scale = -300;
	
	float sc_x = sc_scale/(maxx - minx);
	float tr_x = tr_scale - sc_x*minx;
	float sc_y = sc_scale/(maxy - miny);
	float tr_y = tr_scale - sc_y*miny;
	float xt,yt;
		
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "nodeseq: sc_x = " << sc_x << ", tr_x = " << tr_x << ", sc_y = " << sc_y << ", tr_y = " << tr_y << endl;
		
	/* scale the vertex coordinates */
	if (sc_x < sc_y)
	{
		for (int i=0; i<num_vertices; i++)
		{
			xt = vertex_coordinates[i][0]*sc_x + tr_x + 320;
			yt = vertex_coordinates[i][1]*sc_x + (sc_x/sc_y)*tr_y + 320;
			vertex_coordinates[i][0] = xt;
			vertex_coordinates[i][1] = yt;
		}		

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "nodeseq: scaled vertex_coordinates (sc_x < sc_y): " << endl;
    for (int i=0; i< num_vertices; i++)
    {
		debug << "nodeseq:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << vertex_coordinates[i][j] << ' ';
		debug << endl;
	}	
}
	}
	else
	{
		for (int i=0; i< num_vertices; i++)
		{
			xt = vertex_coordinates[i][0]*sc_y + (sc_y/sc_x)*tr_x + 320;
			yt = vertex_coordinates[i][1]*sc_y + tr_y + 320;
			vertex_coordinates[i][0] = xt;
			vertex_coordinates[i][1] = yt;
		}

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "nodeseq: scaled vertex_coordinates (sc_x >= sc_y): " << endl;
    for (int i=0; i< num_vertices; i++)
    {
		debug << "nodeseq:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << vertex_coordinates[i][j] << ' ';
		debug << endl;
	}	
}
	}
	

	/* store polygon (type 1&2 vertices) vertex_coordinates in (vertex sequence) order in drawcoords */
	matrix<float> drawcoords(vertex_sequence_length+1,2);
	
	for (int i=1; i<=vertex_sequence_length; i++)
	{
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "nodeseq: drawcoords i = " << i << endl;

		for (int j=0; j<2; j++)
		{
			drawcoords[i][j]=vertex_coordinates[vertex_sequence[i]][j];
if (draw_control::DEBUG >= draw_control::DETAIL)
	debug << "nodeseq:   j = " << j << ", sequence vertex = " << vertex_sequence[i] << ", gives " << drawcoords[i][j] << endl;
		}
	}
	
	for (int j=0; j<2; ++j)
		drawcoords[0][j]=drawcoords[vertex_sequence_length][j];
		  
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "nodeseq: drawcoords: " << endl;
    for (int i=0; i< vertex_sequence_length+1; i++)
    {
		debug << "nodeseq:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << drawcoords[i][j] << ' ';
		debug << endl;
	}	
}

	/* write the scaled drawcoords to the output file.  They will appear
	   in the correct sequence ready to be turned into a metapost path */
	ofstream output;
	output.setf(ios::fixed,ios::floatfield);
	output.precision(num_decimal_points);
	
	output.open("__vertex_coords.out");
	if (!output)
	{
		cout << "\nError opening output file __vertex_coords.out\n";
		exit(0);
    }
    else
    {	
		for (int i=0; i< vertex_sequence_length+1; i++)
		{
			output << setw(output_field_width) << drawcoords[i][0] << ' ' << 
			          setw(output_field_width) << drawcoords[i][1] << endl;
		}
		
		output.close();  	
	}
}

