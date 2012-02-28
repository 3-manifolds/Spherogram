/**********************************************************************

This header file defines local programme types used to control objects
used by the programme 
**********************************************************************/

#define TYPE1   -1
#define TYPE2   1

#define FLAT 2
#define POSITIVE 1
#define NEGATIVE -1
#define VIRTUAL 0

#define EPEER 0
#define OPEER 1
#define TYPE    2
#define LABEL	3
#define EVEN_TERMINATING 4
#define ODD_TERMINATING 5
#define EVEN_ORIGINATING 6
#define ODD_ORIGINATING 7
#define COMPONENT 8
#define CODE_TABLE_SIZE 9


#define IMMERSION_CODE true
#define PEER_CODE false

struct draw_control 
{
	enum level {OFF, SUMMARY, BASIC, INTERMEDIATE, DETAIL};
	static int DEBUG;
};

class metapost_control
{

public:

	bool rotate;
	bool explicit_rotation_centre;
	bool implicit_rotation_centre;
	double rotation_degrees;
	double rotation_centre_x;
	double rotation_centre_y;
	int rotation_centre_z;	
	int unit_size;
	int disc_size;
	int pen_size;
	int infinite_cycle;
	bool dash_with_dots;
	bool knotoid;
	bool draw_immersion_only;
	bool draw_labels;
	bool label_edges_from_one;
	bool draw_oriented;
	bool draw_shortcut;
	bool label_vertices;
	bool circle_packing;
	
	metapost_control(): rotate(false), explicit_rotation_centre(false), implicit_rotation_centre(false),
	                    rotation_degrees(0), rotation_centre_x(0) , rotation_centre_y(0), rotation_centre_z(0),
	                    unit_size(20), disc_size(20), pen_size(1), infinite_cycle(-1), dash_with_dots(false), 
	                    knotoid(false), draw_immersion_only(false), draw_labels(false), label_edges_from_one(false), 
	                    draw_oriented(false), draw_shortcut(false), label_vertices(false), circle_packing(false) {}	

};

void print (metapost_control& mp_control, ostream& os, string prefix);

class generic_code_data
{

public:

	enum code_type
	{
		immersion_code,
		peer_code,
	};
	
	int type;
	int head; // used for knotoids
	int num_crossings;
	int num_components;

	matrix<int> code_table;
	vector<int> num_component_edges;
	vector<int> first_edge_on_component;
	vector<int> term_crossing;
	vector<int> orig_crossing;

	generic_code_data(): code_table(matrix<int>(0,0)) {}
};

