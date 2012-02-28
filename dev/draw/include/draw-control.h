/**********************************************************************

This header file defines the #defines required by the programme
**********************************************************************/

struct draw_control 
{
	enum level {OFF, SUMMARY, BASIC, INTERMEDIATE, DETAIL};
	static int DEBUG;
};
