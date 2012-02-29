/************************************************************************
Triangulate calculates the triangulation of the disc determined by an 
immersion code.  This triangulation is the one used by Knotscape: a set of
triangulated sub-discs centred at the immersion crossings and barycentres of 
those components of the immersion's complement having more than two edes 
in their boundary.  The output of the triangulation is written in exactly 
the same form as that produced by the knotscape function triang.

The triangulation is constructed from a set of distinct left and right turning
cycles, that are in 1-1 correspondance with the components of the immersion's 
complement in S^2.  These turning cycles must be calculated by the calling code 
and are passed as a parameter to the function defined here.  A turning cycle is
selected to bound the infinte region of the complement when it is considered 
to be in R^2; that is, the choice of this cycle fixes the point at infinity.  
Either the infinite region may be selected by the calling code or left for
the triangulation to choose, in which case a turning cycle of maximal lngth will 
be selected.  In R2 the "infinite" turning cycle contains in the compact component 
of its complement the remainder of the immersion.  

Note that the order in which the turning cycles are enumerated determines 
a nummbering for the regions of the immersion's complement.

There are four types of vertex in the triangulation:

1. Each crossing of the immersion corresponds to a vertex, crossing i is 
   numbered as vertex i.
   
2. The midpoint of each edge of the immersion is a vertex.  The midpoint of
   edge i is vertex number (i + num_crossings)
   
3. The barycentry of each region in the immersion's complement having more 
   than two edges in it's boundary is a vertex.  These regions are identified 
   by turning cycles of length greater than two and their corresponding 
   vertices are numbered based on the ordering of the turning cycles.

4. There is a ring of vertices constructed in the infinte region in 1-1
   correspondance with the edges in the turning cycle bounding the infinte
   region.  These vertices are placed radially outwards from the mid-point 
   of the edges in this turning cycle.  These vertices are numbered sequentially
   according to the vertex that corresponds to the mid-point of each edge 
   as we proceed around the turning cycle.
   
Type 4 vertices are connected by edges of the triangulation to one type 1
vertex, two type 2 vertices and two other type 4 vertices, as follows.  

 a) Each type 4 vertex is connected radially to the corresponding mid-point
    in an edge of the infinite region turning cycle. 
 b) A type 4 vertex is connected to the type 1 vertex reached by following the
    infinite turning cycle from the type 2 vertex in a)
 c) A type 4 vertex is connected to the type 2 vertex in the next edge of the 
    infinite turning cycle following the type 1 vertex in b)
 d) The type 4 vertices are connected in a ring (the boundary of our triangulated disc)
 
Note that the choice of turning cycle to bound the infinite region may be a left 
or right turning cycle.  If it is a left turning cycle, as we follow the turning cycle
we shall be moving in a clockwise direction, and if it is a right turning cycle an 
anti-clockwise direction.  This is the case regardless of which region is selected
to contain the point at infinity.

The compact regions of the immersion's complement are triangluated as follows.  For a
region having two sides the midpoint vertices of the two edges are connected by an edge.
For a region having more than two sides each midpoint is connected by an edge to the 
midpoint of the adjacent edges in the turning cycle and to the barycentre type 3 vertex.
   
**************************************************************************/
using namespace std;
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <iomanip>

/********************* External variables ***********************/

extern ofstream     debug;


#include <util.h>
#include <matrix.h>
#include <draw.h>
pair<int,int> adjacent_edges(matrix<int>& cycle, int num_left_cycles, int region, int edge);
bool first_occurrence(matrix<int>& code_table, matrix<int>& infinite_cycle_first_visit, int crossing, int position);

void triangulate (generic_code_data& code_data, matrix<int>& cycle, int num_cycles, int num_left_cycles, int infinite_region = -1)
{
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;
	
	matrix<int>& code_table = code_data.code_table;	
	vector<int>& orig_crossing = code_data.orig_crossing;
	vector<int>& term_crossing = code_data.term_crossing;

	/* It is possible that the infinite region meets a crossing in two positions
	   (note that only the infinite region can do this since the components of
	   the immersion's complement are connected).  In this case we shall want to
	   know for each edge in the infinite region's turning cycle which edges 
	   meet these crossings first and which meet them second.  To determine this
	   we record the two edges of the cycle that first meet each crossing.  
	*/
	matrix<int> infinite_cycle_first_visit(num_crossings,2);
	for (int i=0; i< num_crossings; i++)
	{
		for (int j=0; j<2; j++)
			infinite_cycle_first_visit[i][j] = -1;
	}
	
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: calculating infinite_cycle_first_visit:" << endl;
    
	for (int i=1; i<= cycle[infinite_region][0]; i++)
	{
		int this_edge = cycle[infinite_region][i];
		int next_edge = (i<cycle[infinite_region][0]? cycle[infinite_region][i+1]: cycle[infinite_region][1]);

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   this edge  = " << this_edge << ", next_edge = " << next_edge << endl;

		int crossing = (this_edge < 0 ? orig_crossing[abs(this_edge)] : term_crossing[this_edge]);

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   intermediate crossing  = " << crossing << endl;
	
		if (infinite_cycle_first_visit[crossing][0] < 0)
		{
			infinite_cycle_first_visit[crossing][0] = abs(this_edge);
			infinite_cycle_first_visit[crossing][1] = abs(next_edge);

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   setting first visit edges" << endl;
		}
		else
		{
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   first visit edges already set" << endl;
		}
	}

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "triangulate: infinite_cycle_first_visit:" << endl;
    for (int i=0; i< num_crossings; i++)
    {
		debug << "triangulate:   ";
		debug << infinite_cycle_first_visit[i][0] << ' ' << infinite_cycle_first_visit[i][1] << endl;
	}
	debug << endl;
}

	/* Determine the regions that meet at each crossing.  We shall number the
	   regions according to their position in cycle, since each turning cycle
	   corresponds to a region.
	   
	   The regions meeting at a crossing will be stored as a row of crossing_region
	   as follows:
	             
	     \ 3 /        
          \ / 
        0  X  2
          / \ 
         / 1 \
              
         
     where the crosing is drawn int he usual manner, so that the region between 
     the two ingress edges is stored in position 0, and so on.
     
	 We complete crossing_region by considering the left and right turning cycles in turn.
	 Note that we may tell whether an edge in a turning cycle is odd or even by its sign, 
	 since we always traverse even edges following the orientation and odd edges against
	 the orientation.  Also note that a crossing has an odd and even ingress (egress) edge,
	 so if for a left(right) turning cycle we know the type of the crossing and whether
	 we're arriving on an odd or even edge, we know which region around the crossing we're
	 following.
   */
	
	matrix<int> crossing_region(num_crossings,4);
	
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: calculate crossing regions" << endl;
    
	/* First the right turning cycles */
	for (int i=0; i<num_left_cycles; i++)
	{
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   (left) turning cycle " << i << endl;

		for (int j=1; j<=cycle[i][0]; j++)
		{

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     edge " << cycle[i][j] << endl;

			int crossing;
			
			if (cycle[i][j] > 0)
				crossing = cycle[i][j]/2;
			else
				crossing = (abs(cycle[i][j])-1)/2; // we're going backwards along odd edges

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     takes us to crossing " << crossing;
    
			if (cycle[i][j] < 0)
			{
				if (code_table[TYPE][crossing] == TYPE1)
				{
					crossing_region[crossing][2] = i;
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << ", position 2" << endl;
				}
				else
				{
					crossing_region[crossing][1] = i;
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << ", position 1" << endl;
				}
			}
			else
			{
				if (code_table[TYPE][crossing] == TYPE1)
				{
					crossing_region[crossing][0] = i;
					
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << ", position 0" << endl;
				}
				else
				{
					crossing_region[crossing][3] = i;
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << ", position 3" << endl;
				}
			}
		}
	}

	/* Now the right turning cycles */
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   (right) turning cycle " << i << endl;

		for (int j=1; j<=cycle[i][0]; j++)
		{

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     edge " << cycle[i][j] << endl;

			int crossing;
			
			if (cycle[i][j] > 0)
				crossing = cycle[i][j]/2;
			else
				crossing = (abs(cycle[i][j])-1)/2; // we're going backwards along odd edges

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     takes us to crossing " << crossing;
    
			if (cycle[i][j] < 0)
			{
				if (code_table[TYPE][crossing] == TYPE1)
				{
					crossing_region[crossing][3] = i;
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << ", position 3" << endl;
				}
				else
				{
					crossing_region[crossing][2] = i;
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << ", position 2" << endl;
				}
			}
			else
			{
				if (code_table[TYPE][crossing] == TYPE1)
				{
					crossing_region[crossing][1] = i;
					
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << ", position 1" << endl;
				}
				else
				{
					crossing_region[crossing][0] = i;
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << ", position 0" << endl;
				}
			}
		}
	}

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "triangulate: crossing_region: " << endl;
    for (int i=0; i< num_crossings; i++)
    {
		debug << "triangulate:  ";
		for (int j=0; j<4; j++)
			debug << crossing_region[i][j] << ' ';
		debug << endl;
	}
}

	/* note the type 2 vertices corresponding to the mid-point of edges,
	   see notes at the top of the file
	*/
	vector<int> type_2_vertex(num_edges);
	for (int i=0; i< num_edges; i++)
		type_2_vertex[i] = i+num_crossings;

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "triangulate: type_2_vertex: ";
    for (int i=0; i< num_edges; i++)
		debug << type_2_vertex[i] << ' ';
	debug << endl;
}

	/* note the type 3 vertices corresponding to the barycentres
	   of regions having more than two edges.  Regions with two edges 
	   will be recorded as vertex zero but there is no ambiguity since
	   vertex zero always corresponds to crossing zero.
	*/
	vector<int> type_3_vertex(num_cycles);
	int temp = 3*num_crossings; // this is the next vertex number to allocate
	for (int i=0; i< num_cycles; i++)
	{
		if (cycle[i][0] > 2 && i != infinite_region )
			type_3_vertex[i] = temp++;
		else
			type_3_vertex[i] = 0;
	}

if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "triangulate: type_3_vertex: ";
    for (int i=0; i< num_cycles; i++)
		debug << type_3_vertex[i] << ' ';
	debug << endl;
}

	int num_type_123_vertices = temp;
if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "triangulate: num_type_123_vertices = " << num_type_123_vertices << endl;

	/* finally note the type 4 vertices corresponding to the mid-point of the 
	   edges in the turning cycle bounding the infinite region and the inverse
	   mapping of type 4 vertex to type 2 vertex
	*/
	int num_type_4_vertices = cycle[infinite_region][0];
	
	vector<int> type_4_vertex(num_type_4_vertices);

	for (int i=0; i< num_type_4_vertices; i++)
		type_4_vertex[i] = temp++;
	
	int num_vertices = temp;
	
		
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "triangulate: num_vertices = " << num_vertices << endl;
    debug << "triangulate: num_type_4_vertices = " << num_type_4_vertices << endl;
    debug << "triangulate: type_4_vertex: ";
    for (int i=0; i< num_type_4_vertices; i++)
		debug << type_4_vertex[i] << ' ';
	debug << endl;
}


	/* create a matrix to record the neighbours of each vertex, that is 
	   those other vertices to which it is connected by an edge of the 
	   triangulation.  Type 1 vertices are connected to at most five 
	   neighbours, type 2 to at most eight, type 3 vertices are connected
	   to at most the length of the longest turning cycle, and type 4 
	   vertices are connected to exactly five neighbours.  We therefore 
	   have the maximum number of neighbours being 
	   max(8,length of longest turning cycle).  We allow an additional 
	   column (0) to record the number of neighbours.
	*/
	int max_neighbours = 8;
	for ( int i=0; i< num_cycles; i++)
	{
		if (max_neighbours < cycle[i][0])
			max_neighbours = cycle[i][0];
	}

if (draw_control::DEBUG >= draw_control::SUMMARY)
    debug << "triangulate: max_neighbours = " << max_neighbours << endl;
    
    matrix<int> neighbour(num_vertices,max_neighbours+1);
    
    /* Start by noting the neighbours of the type 4 vertices, since this enables us to 
       note to which type 1 vertex each of these vertices is joined to, which we shall 
       need to complete the neighbours of the type 1 vertices.
       
    */
    matrix<int> type_4_vertex_corresponding_to_type_1_vertex(num_crossings,2);
   
    for (int i=num_type_123_vertices; i< num_vertices; i++)
    {
		neighbour[i][0] = 5; //type 4 vertices always have exactly five neighbours

		if (i==num_type_123_vertices)
			neighbour[i][1] = num_vertices-1;
		else
			neighbour[i][1] = i-1;

		int edge = cycle[infinite_region][i-num_type_123_vertices+1];
		neighbour[i][2] = type_2_vertex[abs(edge)];
		if (edge < 0)
		{
			neighbour[i][3] = (abs(edge)-1)/2;
			if (type_4_vertex_corresponding_to_type_1_vertex[(abs(edge)-1)/2][0])
				type_4_vertex_corresponding_to_type_1_vertex[(abs(edge)-1)/2][1] = i;
			else
				type_4_vertex_corresponding_to_type_1_vertex[(abs(edge)-1)/2][0] = i;
		}
		else
		{
			neighbour[i][3] = edge/2;
			if (type_4_vertex_corresponding_to_type_1_vertex[edge/2][0])
				type_4_vertex_corresponding_to_type_1_vertex[edge/2][1] = i;
			else
				type_4_vertex_corresponding_to_type_1_vertex[edge/2][0] = i;
		}
			
		if (i < num_vertices -1)
			edge = cycle[infinite_region][i-num_type_123_vertices+2];
		else
			edge = cycle[infinite_region][1];
			
		neighbour[i][4] = type_2_vertex[abs(edge)];
		
		if (i < num_vertices -1)
			neighbour[i][5] = i+1;
		else
			neighbour[i][5] = num_type_123_vertices;
	}
   
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "triangulate: type_4_vertex_corresponding_to_type_1_vertex: " << endl;
    for (int i=0; i< num_crossings; i++)
    {
		debug << "triangulate:  ";
		for (int j=0; j<2; j++)
			debug << type_4_vertex_corresponding_to_type_1_vertex[i][j] << ' ';
		debug << endl;
	}
	debug << endl;
}
    

	/* determine the neighbours of the type 1 vertices; i.e the immersion crossings.
	   We start from the region in potition 0 and work anti-clockwise around the crossing,
	   only the infinite region contributes a neighbour, and each edge contributes
	   its mid-point.
	*/
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: type 1 vertex neighbours" << endl;
	for (int i=0; i< num_crossings; i++)
	{

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   vertex " << i << endl;

		int nbr = 1; // index into neighbour
		
		/* region in position 0 */
		int region = crossing_region[i][0];

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     region" << region << endl;

		if (region == infinite_region)
		{
			if (type_4_vertex_corresponding_to_type_1_vertex[i][1])
			{
				
				if (first_occurrence(code_table, infinite_cycle_first_visit, i, 0))
					neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][0];				
				else
					neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][1];				
				
			}
			else
			{
				neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][0];
			}
		}
		
		/* following edge */
		int edge;
		if (code_table[TYPE][i] == TYPE1)
			edge = code_table[EVEN_TERMINATING][i];
		else
			edge = code_table[ODD_TERMINATING][i];
		
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     edge" << edge << endl;

		neighbour[i][nbr++] = type_2_vertex[edge];
		
		/* region in position 1 */
		
		region = crossing_region[i][1];

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     region " << region << endl;

		if (region == infinite_region)
		{
			if (type_4_vertex_corresponding_to_type_1_vertex[i][1])
			{
				
				if (first_occurrence(code_table, infinite_cycle_first_visit, i, 1))
					neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][0];				
				else
					neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][1];				
				
			}
			else
			{
				neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][0];
			}
		}
		
		/* following edge */
		if (code_table[TYPE][i] == TYPE1)
			edge = code_table[EVEN_ORIGINATING][i];
		else
			edge = code_table[ODD_ORIGINATING][i];
		
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     edge " << edge << endl;

		neighbour[i][nbr++] = type_2_vertex[edge];
		
		/* region in position 2 */
		region = crossing_region[i][2];

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     region " << region << endl;

		if (region == infinite_region)
		{
			if (type_4_vertex_corresponding_to_type_1_vertex[i][1])
			{
				
				if (first_occurrence(code_table, infinite_cycle_first_visit, i, 2))
					neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][0];				
				else
					neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][1];				
				
			}
			else
			{
				neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][0];
			}
		}
		
		/* following edge */
		if (code_table[TYPE][i] == TYPE1)
			edge = code_table[ODD_ORIGINATING][i];
		else
			edge = code_table[EVEN_ORIGINATING][i];
		
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     edge " << edge << endl;

		neighbour[i][nbr++] = type_2_vertex[edge];

		/* region in position 3 */
		region = crossing_region[i][3];

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     region " << region << endl;

		if (region == infinite_region)
		{
			if (type_4_vertex_corresponding_to_type_1_vertex[i][1])
			{
				
				if (first_occurrence(code_table, infinite_cycle_first_visit, i, 3))
					neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][0];				
				else
					neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][1];				
				
			}
			else
			{
				neighbour[i][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[i][0];
			}
		}
		
		/* following edge */
		if (code_table[TYPE][i] == TYPE1)
			edge = code_table[ODD_TERMINATING][i];
		else
			edge = code_table[EVEN_TERMINATING][i];
		
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     edge " << edge << endl;

		neighbour[i][nbr++] = type_2_vertex[edge];
		
		/* set the number of neighbours */
		neighbour[i][0] = nbr-1;
	}
	
	/* For type 2 vertices we start with the vertex reached by following the 
	   underlying orientation of the edge and then work anti-clockwise.  Note that
	   crossing_region tells us the two regions incident with the edge
	   underlying the type 2 vertex.
	*/
	
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: type 2 vertex neighbours" << endl;
	for (int edge=0; edge< num_edges; edge++)
	{
		int vertex = type_2_vertex[edge];
		
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: edge " << edge << ", vertex " << vertex << endl;

		int nbr = 1; // index into neighbour
			
		/* start with the terminating crossing */
		

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   terminating crossing " << term_crossing[edge] << endl;
    
		neighbour[vertex][nbr++] = term_crossing[edge];
    
		/* then the left hand region */
		int region;
		if (edge % 2)
		{			
			if (code_table[TYPE][term_crossing[edge]] == TYPE1)
				region = crossing_region[term_crossing[edge]][3];
			else
				region = crossing_region[term_crossing[edge]][0];
		}
		else
		{			
			if (code_table[TYPE][term_crossing[edge]] == TYPE1)
				region = crossing_region[term_crossing[edge]][0];
			else
				region = crossing_region[term_crossing[edge]][3];
		}

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   left hand region " << region << endl;
    
		if (region == infinite_region)
		{
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     infinite region" << endl;

			/* locate the edge in the turning cycle bounding the infinite region */
			int offset = 0;
			for (int i=1; i<= cycle[infinite_region][0]; i++)
			{
				if (abs(cycle[infinite_region][i]) == edge)
					break;
				else
					offset++;
			}
			
			/* We have two neighbours in the infinite region, the type 4 vertex associated with the
			   edge and the type 4 vertex associated with a crossing.  If the edge is odd it is the 
			   type 4 vertex corresponding to the crossing we reach along this edge .  If the edge 
			   is even we have as a neighbour the type 4 vertex corresponding  to the crossing 
			   we've just left.
			   
			   We are enumerating our neighbours in an anti-clockwise manner, so if the turning cycle 
			   bounding the infinite region is a left turning cycle we want the type 4 vertex 
			   corresponding to the edge first, and otherwise the type 4 vertex corresponding to
			   the crossing first.
			*/
			int crossing = (edge % 2? term_crossing[edge]: orig_crossing[edge]);
			int occurrence;
			
			if (infinite_cycle_first_visit[crossing][0] == edge || infinite_cycle_first_visit[crossing][1] == edge)
				occurrence = 0;
			else
				occurrence = 1;

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     crossing = " << crossing << ", occurrence = " << occurrence << endl;
			

			if (infinite_region < num_left_cycles)
			{
				neighbour[vertex][nbr++] = type_4_vertex[offset];			
				neighbour[vertex][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[crossing][occurrence];								
			}
			else
			{					
				neighbour[vertex][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[crossing][occurrence];								
				neighbour[vertex][nbr++] = type_4_vertex[offset];			
			}
		}
		else if (cycle[region][0] == 2)
		{
			/* identify the other edge bounding this region, since the neighbour
			   we want is the type 2 vertex on that edge
			*/
			int peer_edge;
			if (abs(cycle[region][1]) == edge)
				peer_edge = abs(cycle[region][2]);
			else
				peer_edge = abs(cycle[region][1]);

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     cycle_length 2, peer edge " << peer_edge << endl;

			neighbour[vertex][nbr++] = type_2_vertex[peer_edge];			

		}
		else
		{

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     non-infinite region of cycle_length > 2" << endl;

			/* identify neighbouring type 2 vertices in the turning cycle */
			pair<int,int> type_2_neighbours = adjacent_edges(cycle, num_left_cycles, region, edge);

if (draw_control::DEBUG >= draw_control::DETAIL)
{
    debug << "triangulate:     type 2 neighbours in turning cycle: " << 
          type_2_neighbours.first << ' ' << type_2_neighbours.second << endl;
}
			
			neighbour[vertex][nbr++] = type_2_vertex[type_2_neighbours.second];			

			neighbour[vertex][nbr++] = type_3_vertex[region];			

			neighbour[vertex][nbr++] = type_2_vertex[type_2_neighbours.first];			
		}

		/* next the originating crossing */
		
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   originating crossing " << orig_crossing[edge] << endl;
    
		neighbour[vertex][nbr++] = orig_crossing[edge];
    

		/* finally the right hand region */
		
		if (edge % 2)
		{
			if (code_table[TYPE][term_crossing[edge]] == TYPE1)
				region = crossing_region[term_crossing[edge]][0];
			else
				region = crossing_region[term_crossing[edge]][1];
		}
		else
		{			
			if (code_table[TYPE][term_crossing[edge]] == TYPE1)
				region = crossing_region[term_crossing[edge]][1];
			else
				region = crossing_region[term_crossing[edge]][0];
		}

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:   right hand region " << region << endl;
    
		if (region == infinite_region)
		{
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     infinite region" << endl;
    
			/* locate the edge in the turning cycle bounding the infinite region */
			int offset = 0;
			for (int i=1; i<= cycle[infinite_region][0]; i++)
			{
				if (abs(cycle[infinite_region][i]) == edge)
					break;
				else
					offset++;
			}
			
			/* We have two neighbours in the infinite region, the type 4 vertex associated with the
			   edge and the type 4 vertex associated with a crossing.  If the edge is odd it is the 
			   type 4 vertex corresponding to the crossing we reach along this edge .  If the edge 
			   is even we have as a neighbour the type 4 vertex corresponding  to the crossing 
			   we've just left.
			   
			   We are enumerating our neighbours in an anti-clockwise manner, so if the turning cycle 
			   bounding the infinite region is a left turning cycle we want the type 4 vertex 
			   corresponding to the edge first, and otherwise the type 4 vertex corresponding to
			   the crossing first.
			*/

			int crossing = (edge % 2? term_crossing[edge]: orig_crossing[edge]);
			int occurrence;
			
			if (infinite_cycle_first_visit[crossing][0] == edge || infinite_cycle_first_visit[crossing][1] == edge)
				occurrence = 0;
			else
				occurrence = 1;

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     crossing = " << crossing << ", occurrence = " << occurrence << endl;
			

			if (infinite_region < num_left_cycles)
			{
				neighbour[vertex][nbr++] = type_4_vertex[offset];			
				neighbour[vertex][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[crossing][occurrence];								
			}
			else
			{					
				neighbour[vertex][nbr++] = type_4_vertex_corresponding_to_type_1_vertex[crossing][occurrence];								
				neighbour[vertex][nbr++] = type_4_vertex[offset];			
			}
		}
		else if (cycle[region][0] == 2)
		{
			/* identify the other edge bounding this region, since the neighbour
			   we want is the type 2 vertex on that edge
			*/
			int peer_edge;
			if (abs(cycle[region][1]) == edge)
				peer_edge = abs(cycle[region][2]);
			else
				peer_edge = abs(cycle[region][1]);

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     cycle_length 2, peer edge " << peer_edge << endl;

			neighbour[vertex][nbr++] = type_2_vertex[peer_edge];			

		}
		else
		{

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate:     non-infinite region of cycle_length > 2" << endl;

			/* identify neighbouring type 2 vertices in the turning cycle */
			pair<int,int> type_2_neighbours = adjacent_edges(cycle, num_left_cycles, region, edge);

if (draw_control::DEBUG >= draw_control::DETAIL)
{
    debug << "triangulate:     type 2 neighbours in turning cycle: " << 
          type_2_neighbours.first << ' ' << type_2_neighbours.second << endl;
}
			neighbour[vertex][nbr++] = type_2_vertex[type_2_neighbours.second];			

			neighbour[vertex][nbr++] = type_3_vertex[region];			

			neighbour[vertex][nbr++] = type_2_vertex[type_2_neighbours.first];			

		}

		/* set the number of neighbours */
		neighbour[vertex][0] = nbr-1;
	}
	

	/* Finally the type 3 vertices, whose neighbours we enumerate anti-clockwise around
	   the vertex noting the mid-point vertex for each edge in the corresponding 
	   region's turning cycle.  Enumerating anti-clockwise means that for a right turning
	   cycles we have to work backwards along the turning cycle
	*/

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: type 3 vertex neighbours" << endl;

	for (int i=0; i< num_left_cycles; i++)
	{
		int vertex = type_3_vertex[i];
		
		if (vertex)
		{
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: vertex " << vertex << endl;
    
			int nbr = 1; // index into neighbour

			for (int j=1; j<= cycle[i][0]; j++)
			{
				int edge = abs(cycle[i][j]);

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: edge " << edge << endl;

				neighbour[vertex][nbr++] = type_2_vertex[edge];			
			}			

			/* set the number of neighbours */
			neighbour[vertex][0] = nbr-1;
		}
	}

	for (int i=num_left_cycles; i< num_cycles; i++)
	{
		int vertex = type_3_vertex[i];
		
		if (vertex)
		{


if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: vertex " << vertex << endl;
   
			int nbr = 1; // index into neighbour

			for (int j=cycle[i][0]; j>= 1; j--)
			{
				int edge = abs(cycle[i][j]);

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "triangulate: edge " << edge << endl;

				neighbour[vertex][nbr++] = type_2_vertex[edge];						

			}			
			
			/* set the number of neighbours */
			neighbour[vertex][0] = nbr-1;

		}
	}
	
if (draw_control::DEBUG >= draw_control::SUMMARY)
{
    debug << "triangulate: neighbour: " << endl;
    for (int i=0; i< num_vertices; i++)
    {
		debug << "triangulate:   vertex " << i << " num neighbours = " << neighbour[i][0] << ": ";
		for (int j=1; j<= neighbour[i][0]; j++)
			debug << neighbour[i][j] << ' ';
		debug << endl;
	}
}

	/* Write out the neighbours to outputfile.  The output file is required to number
	   the vertices from 1 not from zero as we have done internally.  
	   
	   The output file should record the neighbours in an anti-clockwise manner around the 
	   vertex.   We have done this in our calculation of neighbour except in the case that the 
	   turning cycle bounding the infinite region is a right turning cycle.  In this case
	   we have recorded them in a clockwise manner.
	   
	   Finally note that the output file requires the first neighbour to appear at the 
	   end of the list of neighbours as well as at the beginning.  
	   
	   The format of the output is taken from knotscape triang.c.  It includes the value
	   GAMMA that in triang.c is vertex[a[2]], that is, the terminating vertex (numbered from 1) 
	   of the edge paired in the Dowker code with edge 2 (edges numbered from 1).
	   Here this the the terminating vertex of edge 1.
	*/
	ofstream output;
	output.open("__triangulate.out");
	if (!output)
	{
		cout << "\nError opening output file __triangulate.out\n";
		exit(0);
    }
    else
    {
	
		output << "NODECOUNT: " << num_vertices << endl;
		output << "ALPHA/BETA/GAMMA: 1 " << num_type_123_vertices+1 << ' ' << term_crossing[1]+1 << endl;
        output << "GEOMETRY: hyperbolic" << endl;
        output << "FLOWERS:\n" << endl;
        
//		for (int i=0; i< (infinite_region < num_left_cycles? num_vertices : num_type_123_vertices); ++i)
		for (int i=0; i< num_type_123_vertices; ++i)
		{
			output << i+1 << " " << neighbour[i][0] << ' ';
            for (int j=1; j<=neighbour[i][0]; j++) 
				output << neighbour[i][j]+1 << ' ';
			output << neighbour[i][1]+1;
			output << endl; 
		}
		
		if (infinite_region < num_left_cycles)
		{
			for (int i=num_type_123_vertices; i < num_vertices; ++i)
			{
				output << i+1 << " " << neighbour[i][0]-1 << ' ';
//				output << neighbour[i][1]+1 << ' ';
				for (int j=1; j<=neighbour[i][0]; j++) 
					output << neighbour[i][j]+1 << ' ';
				output << endl;
			}
		}
		else
		{
			for (int i=num_type_123_vertices; i < num_vertices; ++i)
			{
				output << i+1 << " " << neighbour[i][0]-1 << ' ';
//				output << neighbour[i][1]+1 << ' ';
				for (int j=neighbour[i][0]; j>=1; j--) 
					output << neighbour[i][j]+1 << ' ';
				output << endl;
			}
		}
		
		output << "END" << endl;
		output.close();  	
	}
}

/* adjacent_edges identifies the pair of edges adjacent to edge in the turning cycle 
   corresponding to region.  The adjacent edges are listed so that the sequence
   first, edge, second appears anti-clockwise (with respect to the infinite region).
*/
pair<int,int> adjacent_edges(matrix<int>& cycle, int num_left_cycles, int region, int edge)
{
	pair<int,int> neighbours;
	
	int offset=1;
	int cycle_length = cycle[region][0];

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "adjacent_edges: cycle_length = " << cycle_length << endl;

	for (int i=1; i<= cycle_length; i++)
	{
		if (abs(cycle[region][i]) == edge)
			break;
		else
			offset++;
	}

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "adjacent_edges: offset = " << offset << endl;

    if (region < num_left_cycles)
    {
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "adjacent_edges: left turning cycle" << endl;

		if (offset > 1)
			neighbours.first = abs(cycle[region][offset-1]);
		else
			neighbours.first = abs(cycle[region][cycle_length]);

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "adjacent_edges: first = " << neighbours.first << endl;
		
		if (offset < cycle_length)
			neighbours.second = abs(cycle[region][offset+1]);
		else
			neighbours.second = abs(cycle[region][1]);		

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "adjacent_edges: second = " << neighbours.second << endl;
	}
    else
    {
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "adjacent_edges: right turning cycle" << endl;

		if (offset > 1)
			neighbours.second = abs(cycle[region][offset-1]);
		else
			neighbours.second = abs(cycle[region][cycle_length]);
		
if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "adjacent_edges: second = " << neighbours.second << endl;

		if (offset < cycle_length)
			neighbours.first = abs(cycle[region][offset+1]);
		else
			neighbours.first = abs(cycle[region][1]);		

if (draw_control::DEBUG >= draw_control::DETAIL)
    debug << "adjacent_edges: first = " << neighbours.first << endl;
	}
	
	return neighbours;
}
