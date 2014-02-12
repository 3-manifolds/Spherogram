#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PMdef.h"
#include "PMconjugation.h"
#include "PMenlight.h"
#include "PMextract.h"
#include "PMplanmap.h"

/******************************/
/* This function checks and completes parameters with defaults */
/******************************/
int pmSetParameters(pmSize *Size, pmMethod *Meth)
{
  long i;

  /* Default parameters */

  if (Size->m == 0) {                 // type of maps generated 
    Size->m = 5; Size->b = 5;          // Q2 : 2-c quartic
    printvf("# 2-c quartic by default\n");
  }

  // compute face, vertice and edge numbers...

  if (Size->e + Size->v + Size->f + Size->r + Size->g == 0 && 
      Size->dgArr == NULL){
    fprintf(stderr,"Size must be given somehow\n");
    return(FALSE);
  }

  switch(Size->m){

  case 1:                          // all type of cubic maps
  case 2:
  case 3:
  case 7:
  case 8:
    if (Size->v){
      if (Size->v % 2 == 0){
	Size->e = 3 * Size->v / 2;
	Size->f  = Size->e + 2 - Size->v;
      }else{
	fprintf(stderr,"nb vertex must be even for cubic\n");
	return(FALSE);
      }
    }else if (Size->f){
      Size->e = 3 * Size->f - 6;
      Size->v = 2 * Size->e / 3;
    }else if (Size->e){
      if (Size->e % 3 == 0){
	Size->f = Size->e / 3 - 2;
	Size->v = 2*Size->e / 3;
      }else{
	fprintf(stderr,"nb edges must be multiple of three for bicubic\n");
	return(FALSE);
      }
    }else{
      fprintf(stderr,"degrees or colors not available for cubic\n");
      return(FALSE);
    }break;

  case 9:                     // all type of quartic
    if (Size->v % 2 != 0){
      fprintf(stderr,"vtx number must be even for bi-quartic\n");
      return(FALSE);
    }// no break here !
  case 4:                     
    if (Size->g || Size->r){
      fprintf(stderr,"color control only for 2-c or 3-c quartic\n");
      return(FALSE);
    }// no break here !
  case 5:                     
  case 6: 
    if (Size->g && Size->r){
      Size->f = Size->g + Size->r;
      Size->v = Size->f - 2;
      Size->e = 2 * Size->v;
    }else if (Size->v){
      Size->f = Size->v + 2;
      Size->e = 2 * Size->v;
    }else if (Size->e){
      if (Size->e % 2 == 0){
	Size->v = Size->e / 2;
	Size->f = Size->v + 2;
      }else{
	fprintf(stderr,"nb edges must be even for quartic\n");
	return(FALSE);
      }
    }else if (Size->f){
      Size->v = Size->f - 2;
      Size->e = 2 * Size->v;
    }
    if (Size->g) Size->r = Size->f - Size->g;
    if (Size->r) Size->g = Size->f - Size->r;
    if (Size->g < 0 || Size->r < 0){
      fprintf(stderr,"not enough faces for colors \n");
      return(FALSE);
    }break;

  case 10:                         // eulerian 
    if (Size->g || Size->r){
      fprintf(stderr,"color control not implemented for eulerian\n");
      return(FALSE);
    }else{
      Size->e = 0; 
      Size->v = 0;
      for (i=0; i<Size->d; i++){
	Size->v = Size->v + Size->dgArr[i];
	Size->e = Size->e + (i+1) * Size->dgArr[i];
      }
      Size->e = Size->e / 2;
      Size->f = Size->e + 2 - Size->v;
    }
  }
  printvf("# Edges: %ld ; Faces: %ld ; Vtx: %ld ; Greens: %ld ; Reds: %ld\n",
	  Size->e, Size->f, Size->v, Size->g, Size->r);
  
  if (Size->m == 2 || Size->m == 3 || // extraction methods for core maps 
      Size->m == 6 || Size->m == 8) {     
    if (!Meth->core) Meth->core = 2;  // largest component enabled
    if (!Meth->pic)  Meth->pic  = 1;  // pic optimization enabled
    if (Size->t == -1)                // default tolerence n^(2/3) 
      Size->t = (long) exp(2*log(Size->v)/3);
    printvf("# Size interval: %ld, %ld\n",
	    Size->v - Size->t, Size->v + Size->t);
  }
  

  return(TRUE);
}


/******************************/
/* this function pre-compute the basic memory requirements */
/******************************/
int pmMemoryInit(pmSize *S, pmMethod *Meth, pmMemory *M)
{
  
  switch(S->m){

  case 1:                            // cubic
  case 2:
  case 3:
    M->dTree = 3;                    // uses a ternary tree T
                                     // nb nodes in T
    if (S->m == 1)                   //    for dually 2-c 
      M->sTree = S->v / 2;              
    else if (S->m == 2){             //    for dually 3-c
      if (Meth->pic == 1) M->sTree = (long) (S->v)+2-0.77*exp(2*log(S->v)/3);
      else M->sTree = S->v;
    }else if (S->m == 3){             //    for dually 4-c
      if (Meth->pic == 1) M->sTree = (long) 2*(S->v)+2-0.77*exp(2*log(S->v)/3);
      else M->sTree = S->v * 2;
    }      
    M->sWrd  = 3 * M->sTree + 1;     // nb letters in word
    M->sEdge = 8 * M->sTree + 2;     // nb half edges used
    M->sVtx  = 4 * M->sTree + 2;     // nb vertices+faces
    M->sLeaf = 2 * M->sTree + 2;     // nb leaf in stack for balance
    break;

  case 4:                            // quartic
    M->dTree = 2;                    // uses a binary tree T
    M->sTree = S->v;                 // nb nodes in T
    M->sWrd  = 2* M->sTree + 1;      // nb letters in word
    M->sEdge = 4* M->sTree + 2;      // nb half edges used
    M->sVtx  = 2* M->sTree + 3;      // nb vtx+faces 
    M->sLeaf = 2* M->sTree + 2;      // nb leaf in stack for balance
    break;
    
  case 5:                                 // quartic 2-c and 3-c
  case 6:
    M->dTree = 3;                         // uses a ternary tree T
                                          // nb nodes in T
    if (S->m == 5){                       //   for 2-c
      if (S->r){                          //     with color control
	M->rTree = S->r - 1;
	M->gTree = S->g - 2;
      }else M->rTree = M->gTree = 0;
      M->sTree = S->v - 1;                //     without 
    }else if (S->m == 6){                 //   for 3-c 
      if (Meth->pic == 1) M->sTree =(long)3*(S->v)-1.22*exp(2*log(3*S->v)/3); 
      else M->sTree = 3 * S->v;
      if (S->r){                          //   with color control
	M->rTree = S->r * 3;                 //   (pic correction non ready)
	M->gTree = S->g * 3;
	M->sTree = S->v * 3;
      }else M->rTree = M->gTree = 0;
    }
    if (S->r) M->sWrd = 6 * M->sTree + 3; // nb letters in word
    else      M->sWrd = 3 * M->sTree + 1;  
    M->sEdge = 6 * M->sTree + 2;          // nb half edges used
    M->sVtx  = 2 * M->sTree + 4;          // nb vertices+faces used
    M->sLeaf = 4 * M->sTree;              // nb leaf in stack for balance
    break;

  case 7:                               // bicubic
  case 8:
    M->dTree = 2;                       // uses a binary tree T
                                        // nb nodes in T
    if (S->m == 7)                      //   for 1-c
      M->sTree = S->v / 2; 
    else if (S->m == 8){                //   for 2-c  
      if (Meth->pic == 1) 
	M->sTree = (long) 9*(S->v)/5+2-0.77*exp(2*log(9*S->v/5)/3);
      else M->sTree = 9*(S->v)/5+2;
    }
    M->sWrd  = 2 * M->sTree + 1;        // nb letters in word
    M->sEdge = 6 * M->sTree;            // nb half edges used
    M->sVtx  = 3 * M->sTree + 2;        // nb vertices+faces 
    M->sLeaf = 2 * M->sTree + 1;        // nb leaf in stack for balance
    break;
  case 9:                               //biquartic
    M->dTree = 3;                       // uses a ternary tree T
    M->sTree = S->v / 2;
    M->sWrd  = 3 * M->sTree + 1;
    M->sEdge = 8 * M->sTree;
    M->sVtx  = 4 * M->sTree + 2;
    M->sLeaf = 4 * M->sTree + 1;
    break;
  default:
    fprintf(stderr,"unknown type of map %d", (int) S->m);
    return(FALSE);
  }
  printvf("# Size of tree: %ld\n", M->sTree);
  printvf("# Memory       : %ld vtx, %ld edgs\n", M->sVtx, M->sEdge);
  
  return(TRUE);
}

/******************************/
/* this function extends the number of edges created in case of extraction */
/******************************/
int pmExtendMemory(pmSize *S, pmMethod *Meth, pmMemory *M, char OtherReason)
{
  char map = S->m;
  
  if (OtherReason == TRUE &&
      (map == 1 || map == 2 || map == 4 || map == 5 || map == 7)) map++;

  switch(map){                   // cubic 2c
  case 2:
    M->sVtx  = 2 * M->sVtx;
    break;
  case 3:                        // cubic 3c
    M->sEdge = 8 * M->sEdge;
    M->sVtx  = 8 * M->sVtx;
    break;
  case 5:                        // quartic 2-c 
    M->sVtx  = 2 * M->sVtx;
    break;
  case 6:                        // quartic 3-c
    M->sEdge = 8 * M->sEdge;
    M->sVtx  = 8 * M->sVtx;
    break;
  case 8:                        // bicubic
    M->sVtx  = 2 * M->sVtx;
    break;
  }
  printvf("# Memory (extd): %ld vtx, %ld edgs\n", 
	  M->sVtx, M->sEdge);

  return(TRUE);
}

/******************************/
/* This function uses conjugation of tree to sample a map */
/******************************/
int pmTreeConjugation(pmSize *S, pmMemory *M, 
		      pmMap *Map)
{
  long pos;
  char *Wrd;
  pmStck Stack;
  pm_edge *Root;
    

  pmCreateWrd(M->sWrd, &Wrd);     // allocation of the word
                                  
  if (S->b == 5 && M->rTree)      // generation of a random word
    pos = pmLuka3(M->rTree, M->gTree, Wrd);
  else                          
    pos = pmLuka1(M->sTree, M->dTree, Wrd);


  pmCreateEdge(M->sEdge);         // allocation of edges 
  pmCreateVtx(M->sVtx);           // and vertices in global variables

  if (S->b == 5 && M->rTree)      // transformation word -> tree
  Root = pmChottin2tree(pos, Wrd);
  else
  Root = pmLuka2tree(pos, Wrd);
  
  pmFreeWrd(Wrd);                   // finished with the word

  
  switch(S->b){                   // adjunction of buds
  case 1: pmSpring4(Root); break;
  case 4: pmSpring1(Root); break;
  case 5: pmSpring3(Root); break;
  case 7: pmSpring2(Root); break;
  case 9: pmSpring5(Root); break;
  }

  Root = pmBalance(Root);           // the conjugation

  pmCreateStck(M->sLeaf, &Stack);   // the closure, using a stack
  Root = pmClosure(Root, &Stack);
  Root = pmSuppress(Root);            
  pmFreeStck(Stack);      

  
  Map->e = pmLabelCanon(Root);      // set labels and create faces
  Map->v = pmChainVtx(Root);
  Map->f = pmAddFaces(Root);
  Map->root = Root;


  return(TRUE);
}

/******************************/
/* this function extract the core or largest component */
/******************************/
void pmExtract(pmSize *S, pmMethod *Meth, pmMemory *M, pmMap *Map)
{
  pm_edge *Root, *maxEdge;
  long maxVnb, maxEnb;

  Root = Map->root;

  pmCreatePost(M->sEdge); pmCreateSeed(M->sEdge);   // compoment lists
  pmCreateBloc(M->sEdge); pmCreateComp(M->sEdge);

  if (Meth->core == 1){                     // core extraction method
    switch (S->m){
    case 2:
    case 8:
      Root = pmGet3tri(Root); break;
    case 3:
      Root = pmGet4tri(Root); break;
    case 6:
      Root = pmGet3c(Root); break;
    }
    if (Root != NULL){
      Map->e = pmLabelCanon(Root);
      Map->v = pmChainVtx(Root);
      Map->f = pmChainFaces(Root);
      pmLabelFaces(Map->root->face);
    }else  
      Map->v = 0;
    Map->root = Root;
  }

  if (Meth->core == 2){                      // largest component method
    switch (S->m){
    case 2:
    case 8:
      pmFull2to3tri(Root); break;            // all components are extracted
    case 3:
      pmFull2to4tri(Root); break;
    case 6:
      pmFull2to3c(Root); break;
    }

    maxVnb = 0; maxEnb = 0; maxEdge = Root;  // the largest is selected
    pmFirstComp();
    while(pmIsComp()){
      Root = pmNextComp();
      Map->e = pmLabelCanon(Root);
      Map->v = pmChainVtx(Root);
      if (Map->v >= maxVnb){
	maxVnb = Map->v; 
	maxEnb = Map->e;
	maxEdge = Root;
      }
    }
    Map->v = maxVnb;
    Map->e = maxEnb;
    Map->root = maxEdge;
    Map->f = pmChainFaces(Map->root);
    pmLabelFaces(Map->root->face);
  }
  pmFreeBloc();pmFreePost();pmFreeSeed();pmFreeComp();
}


/******************************/
/* this function is the main sampling function */
/******************************/
int pmPlanMap(pmSize *S, pmMethod *Meth, pmMemory *M, pmMap *Map)
{
  long numTry;
  
  if ((S->m == 1 || S->m == 4 ||                    // basic families
       S->m == 5 || S->m == 7 || S->m == 9))
    pmTreeConjugation(S, M, Map);

  else if (S->m == 2 || S->m == 3 ||                // extracted families
	   S->m == 6 || S->m == 8){
    numTry = 0;
    do{
      pmTreeConjugation(S, M, Map);
      pmExtract(S, Meth, M, Map);
      if (Map->v < S->v - S->t || Map->v > S->v + S->t){
	pmFreeEdge(); pmFreeVtx();
      }
      numTry++;
    }while (Map->v < S->v - S->t || Map->v > S->v + S->t);
    printwf("# NbTry%ld = %ld; Final Size = %ld;\n",Map->i, numTry, Map->v);
  }
  
  return(TRUE);
}

int pmFreeMap(pmMap *Map)
{
  Map->root=NULL;
  pmFreeEdge();
  pmFreeVtx();
  return (TRUE);
}






