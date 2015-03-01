#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif

#include "PMdef.h"

/******************************
 this function displays the command line modifiers
******************************/
void pmHelp(void)
{
  printf("PlanarMap v1.2  Random sampling of planar maps and graphs\n");
  printf("    by Gilles Schaeffer, CNRS, Jan 2001\n");
  printf("    mailto:Gilles.Schaeffer@lix.polytechnique.fr\n\n");
  printf("usage : planarmap -CQMBNEVFRGIOSX<num> -lcspdhvw\n\n");
  printf("  map selection:\n");
  printf("   -C<2,3,4>: 2,3,4-edge-connected cubic maps\n");
  //  printf("        = singular, usual, irreducible triangulations\n");
  //  printf("   -M<1,(2),(3),4,(5),6>:  general, loopless, simple\n");
  //  printf("                     2-connected, simple 2-c, 3-connected\n");
  printf("   -Q<1,2,3>: 2,4,6-edge-connected quartic maps\n");
  printf("   -Q4      : bipartite quartic maps\n");
  printf("   -M<1,2,3>: 1,2,3-connected maps\n");
  printf("   -M4      : bipartite bicolor maps\n");
  printf("   -B<1,2>  : 2,3-edge-connected bipartite cubic\n");
  //  printf("        =  1,2-connected bipartite maps\n");
  //  printf("        =  singular, usual bicolored triangulations\n");
  //  printf("   -D<max> <d1> <d2> ... <dmax>: bipartite + face distrib\n");
  printf("  parameters:\n");
  printf("   -N<nb>: number of maps generated (default is one)\n");
  printf("   -E<nb>: number of edges\n");
  printf("   -V<nb>: number of vertices (=edges for -Mk)\n");
  printf("   -F<nb>: number of faces\n");
  printf("   -R<nb>: number of red faces in Q2 (=vertices for -M2)\n");
  printf("   -G<nb>: number of green faces in Q2 (=faces for -M2)\n");
  printf("   -I<nb>: approximate size [N+/-I]\n");
  printf("  methods:\n");
  printf("   -l extraction by largest component (default)\n");
  printf("   -c extraction by core \n");
  printf("   -s suppress pic optimisation in extraction\n");
  printf("  output:\n");
  printf("   -p: print map (human readable format is default)\n");
  printf("   -d: print dual map\n");
  printf("   -O<num>: output format\n");
  printf("       1= human, 2=edgelist, 3=maple, 4=graphcode,\n");
  printf("       5=leda dimacs, 6=adjacency matrix, 7=minor of laplacian\n");
  printf("       8=edgelist tulip\n");
  printf("  stats:\n");
  printf("   -S<num>: 1 = number of loops\n");
  printf("           2,3,4 = sizes (core, two largest, all),\n"); 
  printf("           5,6,7,8 = distances, 9 = face degrees  \n"); 
  printf("           10,11 = Gauss components (nb, max)  \n"); 
  printf("  system:\n");
  printf("   -h: help\n");
  printf("   -v: verbose (better put it first)\n");
  printf("   -w: super verbose (better put it first)\n");
  printf("   -X<seed>: seed of random generator\n"); 
  exit(0); 
}


/******************************/
/* This function parses the command line */
/******************************/
int pmParseArgs(int argc, char *argv[], 
		pmSize *Size, pmMethod *Meth,
		pmOutput *Outp, pmStats *Stat)
{
  long param;
  char modifier;
  long i,j;

  /* Size preinit */
  Size->m = 5; Size->m = 5;
  Size->e = 0; Size->v = 0; Size->f = 0;
  Size->r = 0; Size->g = 0; Size->d = 0;
  Size->t = -1; Size->dgArr = NULL;
  /* Meth preinit */
  Meth->core = 0; Meth->pic = 0;
  Meth->seed = getpid();
  Meth->verbose = 0; 
  /* Output preinit */
  Outp->format = 1;  Outp->transform = 0;
  Outp->map = 0;     Outp->dual = 0;
  /* Stats preinit */
  Stat->nb = 1;
  Stat->stats = 0; Stat->loop = 0;
  Stat->core = 0;  Stat->dist = 0;
  Stat->facedeg =0; Stat->gauss = 0;
  Stat->gaussmax = 0;

  /* main loop to parse args */

  for (i = 1; i < argc; i++){

    /* numerical arguments */

    if (sscanf(argv[i],"-%c%ld", &modifier, &param) == 2){
      switch(modifier){
      case 'C': //cubic maps
	if (param == 2){ // dually 2-connected
	  printvf("# 2-edge-connected cubic maps\n");
	  Size->m = 1; Size->b = 1;
	}else if (param == 3){ // dually 3-connected
	  printvf("# 3-edge-connected cubic maps\n");
	  Size->m = 2; Size->b = 1;
	}else if (param == 4){ // dually 4-connected
	  printvf("# 4-edge-connected cubic maps\n");
	  Size->m = 3; Size->b = 1;
	}else{
	  fprintf(stderr,"unknown kind of cubic\n");
	  exit(2);
	}break;
      case 'Q': //4-regular maps
	if (param == 1){ 
	  printvf("# 2-edge-connected quartic maps\n");
	  Size->m = 4; Size->b = 4;
	}else if (param == 2){
	  printvf("# 4-edge-connected quartic maps\n");
	  Size->m = 5; Size->b = 5;
	}else if (param == 3){
	  printvf("# 6-edge-connected quartic maps\n");
	  Size->m = 6; Size->b = 5;
	}else if (param == 4){
	  printvf("# bi-quartic maps");
	  Size->m = 9; Size->b = 9;
	}else{
	  fprintf(stderr,"unknown kind of quartic\n");
	  exit(2);
	}break;
      case 'M': //general maps
	Outp->transform = 1;
	if (param == 1){ 
	  printvf("# general map\n");
	  Size->m = 4; Size->b = 4;
	}else if (param == 2){
	  printvf("# nonseparable\n");
	  Size->m = 5; Size->b = 5;
	}else if (param == 3){
	  printvf("# 3-connected\n");
	  Size->m = 6; Size->b = 5;
	}else if (param == 4){
	  printvf("# bipartite bicolor\n");
	  Size->m = 9; Size->b = 9;
	}else if (param == 7){
	  printvf("# simple nonseparable\n");
	  Size->m = 13; Size->b = 5;
	}else if (param == 8){
	  printvf("# loopless maps\n");
	  Size->m = 11; Size->b = 4;
	}else if (param == 9){
	  printvf("# loopless simple maps\n");
	  Size->m = 12; Size->b = 4;
	}else{
	  fprintf(stderr,"unknown kind of map\n");
	  exit(2);
	}break;
      case 'B': //bipartie
	if (param == 1){ // 1-c
	  printvf("# 2-edge-connected bipartite cubic\n");
	  Size->m = 7; Size->b = 7;
	}else if (param == 2){ // 2-c
	  printvf("# 3-edge-connected bipartite cubic\n");
	  Size->m = 8; Size->b = 7;
	}else{
	  fprintf(stderr,"unknown kind of bicubic\n");
	  exit(2);
	}break;
      case 'D': // given degrees
	if (param > 0) { // max degree
	  printvf("# eulerian with reduced degree distribution:");
	  Size->dgArr = (long *)calloc(param, sizeof(long));
	  for (j=0; j < param; j++){
	    if (sscanf(argv[++i],"%ld",Size->dgArr+j)){
	      printvf("#  %ld ", Size->dgArr[j]);
	    } else {
	      fprintf(stderr,"\nincoherent degree distribution\n");
	      exit(2);
	    }	    
	  }
	  printvf("# \n");
	  Size->d = param;
	  Size->m = 10; Size->b = 10;
	}
      case 'N': // more than one map ?
	if (param > 0){
	  printvf("# Number of maps: %ld\n", param);
	  Stat->nb = param;
	}break;
      case 'E': // number of edges
	if (param > 0){
	  printvf("# Number of edges: %ld\n",param);
	  Size->e = param;
	}break;
      case 'V': // number of vertices
	if (param > 0){
	  printvf("# Number of vertices: %ld\n",param);
	  Size->v = param;
	}break;
      case 'F': // number of faces
	if (param > 0){
	  printvf("# Number of faces: %ld\n",param);
	  Size->f = param;
	}break;
      case 'R': // red ones
	if (param > 0){
	  printvf("# Number of red ones: %ld\n",param);
	  Size->r = param;
	}break;
      case 'G': // green ones
       	if (param > 0){
	  printvf("# Number of green ones: %ld\n",param);
	  Size->g = param;
	}break;
      case 'I': // error allowed on size
	if (param >=0){
	  printvf("# error allowed on size: %ld\n",param);
	  Size->t = param;
	}break;
      case 'O': // output format
	if (param == 1){
	  printvf("# human readable format\n");
	  Outp->format = 1;
	}if (param == 2){
	  printvf("# edge list output format\n");
	  Outp->format = 2;
	}else if (param == 3){
	  printvf("# maple output format\n");
	  Outp->format = 3;
	}else if (param == 4){
	  printvf("# graphcode output format\n");
	  Outp->format = 4;
	}else if (param == 5){
	  printvf("# LEDA dimacs format\n");
	  Outp->format = 5;
	}else if (param == 6){
	  printvf("# adjacency matrix\n");
	  Outp->format = 6;
	}else if (param == 7){
	  printvf("# principal minor of laplacian\n");
	  Outp->format = 7;
	}else if (param == 8){
	  printvf("# tulip format\n");
	  Outp->format = 8;
	}break;
      case 'S': // statistic
	if (param > 0 && param < 12){
	  printvf("# statistic output of type %ld enabled\n", param);
	  Stat->stats = param;
	  switch (param){
	  case 5:
	    Stat->dist = 1; break;
	  case 6:
	    Stat->dist = 2; break;
	  case 7:
	    Stat->dist = 3; break;
          case 8:
	    Stat->dist = 4; break;
	  case 9:
	    Stat->facedeg = TRUE; break;
	  case 10:
	    Stat->gauss = 1; break;
	  case 11:
	    Stat->gaussmax = 1; break;
	  default:
	      break;
	  }
	}else
	  fprintf(stderr,"unknown statistic\n");
	break;
      case 'X': // given seed of random generator
	Meth->seed = param;
	break;
      default :
	printvf("# unused arg -%c%ld\n", modifier, param);
      }
    }else if (sscanf(argv[i],"-%c", &modifier)){

      /* simple switches */

      switch(modifier){
      case 'v': // verbose
	Meth->verbose = TRUE;
	break;
      case 'w': // superverbose
	Meth->verbose = TRUE+TRUE;
	break;
      case 'h': // usage
	pmHelp();
	break;
      case 'p': // print map(s)
	Outp->map = 1;
	break;
      case 'd': // print dual map(s)
	Outp->dual = 1;
	break;
      case 'c': // core method
	Meth->core = 1;
	break;
      case 'l': // largest component method
	Meth->core = 2;
	break;
      case 's': // suppress pick correction
	Meth->pic = 2;
	printvf("# pic optimization disabled\n");
	break;
      default:
	printvf("# unused arg -%c\n", modifier);
      }
    }
  }
  return(TRUE);
}
