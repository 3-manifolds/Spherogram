#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "PMdef.h"
#include "PMenlight.h"
#include "PMdisplay.h"

/* fonctions d'affichage des cartes */

/******************************/
/* this function MAPLE-display the map in depth first search order */
/******************************/
void pmPrintMap(pm_edge *Root, long num)
{
  pm_edge *Cur1, *Cur2;
  long mark;
  
  mark = pmNewMark();

  Root->from->mark = mark;
  Root->mark = mark;
  printf("Map%ld := [\n", num);
  printf(" [ %ld, [ %ld", Root->from->label, Root->label);
  for (Cur1 = Root->next; Cur1 != Root; Cur1 = Cur1->next)
    printf(", %ld",Cur1->label);
  printf("]]");

  if (Root->oppo->from != Root->from)
    Cur1 = Root->oppo->next;
  else 
    Cur1 = Root->next;
  while (Cur1 != Root){
    if (Cur1->from->mark != mark){
      Cur1->from->mark  = mark;
      printf(",\n [ %ld, [%ld",Cur1->from->label, Cur1->prev->label);
      for (Cur2 = Cur1; Cur2 != Cur1->prev; Cur2 = Cur2->next)
	printf(", %ld",Cur2->label);
      printf("]]");    
    }
    if ((Cur1->oppo->mark == mark && Cur1->oppo->from != Cur1->from) ||
	Cur1->oppo->from->mark != mark){
      Cur1->mark = mark;
      Cur1 = Cur1->oppo;
    }
    Cur1 = Cur1->next;
  }
  printf("];\n");
}
/******************************/
/* this function MAPLE-display the dual in depth first search order */
/******************************/
void pmPrintDualMap(pm_edge *Root, long num)
{
  pm_edge *Cur1, *Cur2;
  long mark;
  
  mark = pmNewMark();

  Root->face->mark = mark;
  Root->mark = mark;
  printf("Dual%ld := [\n", num);
  printf(" [ %ld, [ %ld", Root->face->label, Root->label);
  for (Cur1 = Root->Next; Cur1 != Root; Cur1 = Cur1->Next)
    printf(", %ld",Cur1->label);
  printf("]]");

  if (Root->oppo->face != Root->face)
    Cur1 = Root->oppo->Next;
  else 
    Cur1 = Root->Next;
  while (Cur1 != Root){
    if (Cur1->face->mark != mark){
      Cur1->face->mark  = mark;
      printf(",\n [ %ld, [%ld",Cur1->face->label, Cur1->Prev->label);
      for (Cur2 = Cur1; Cur2 != Cur1->Prev; Cur2 = Cur2->Next)
	printf(", %ld",Cur2->label);
      printf("]]");    
    }
    if ((Cur1->oppo->mark == mark && Cur1->oppo->face != Cur1->face) ||
	Cur1->oppo->face->mark != mark){
      Cur1->mark = mark;
      Cur1 = Cur1->oppo;
    }
    Cur1 = Cur1->Next;
  }
  printf("];\n");
}

/******************************/
/* this function USER-display the map in vertex-chain order */
/******************************/
void pmPrintChndVtx(pm_vertex *Vtx)
{
  pm_edge *Cur1;
  
  printf("Vertex %ld: ", Vtx->label);
  Cur1 = Vtx->root;
  while(Cur1 != Vtx->root->prev){
    printf(" %ld", Cur1->label);
    Cur1 = Cur1->next;
  }
  printf(" %ld\n", Cur1->label);
  while(Vtx->next != NULL){
    Vtx = Vtx->next;
    printf("Vertex %ld: ", Vtx->label);
    Cur1 = Vtx->root;
    while(Cur1 != Vtx->root->prev){
      printf(" %ld", Cur1->label);
      Cur1 = Cur1->next;
    }
    printf(" %ld\n", Cur1->label);
  }
}


/******************************/
/* this function USER-display the dual in face-chain order */
/******************************/
void pmPrintChndFaces(pm_vertex *Face)
{
  pm_edge *Cur1;
  
  printf("Face %ld - %d: ", Face->label,Face->type);
  Cur1 = Face->root;
  while(Cur1 != Face->root->prev->oppo){
    printf(" %ld", Cur1->label);
    Cur1 = Cur1->oppo->next;
  }
  printf(" %ld\n", Cur1->label);
  while(Face->next != NULL){
    Face = Face->next;
    printf("Face %ld - %d: ", Face->label, Face->type);
    Cur1 = Face->root;
    while(Cur1 != Face->root->prev->oppo){
      printf(" %ld", Cur1->label);
      Cur1 = Cur1->oppo->next;
    }
    printf(" %ld\n", Cur1->label);
  }
}

/******************************/
/* this function display the map as edge list (graph structure only) */
/******************************/
void pmPrintEdgeList(pm_vertex *Vtx) 
{
  pm_edge *Cur1; 

  printf("\n");  
  Cur1 = Vtx->root; 
  while(Cur1 != Vtx->root->prev){ 
    if (Cur1->label > 0) 
      printf(" %ld - %ld \n", Cur1->from->label, Cur1->oppo->from->label); 
    Cur1 = Cur1->next; 
  }
  if (Cur1->label > 0) 
    printf(" %ld - %ld \n", Cur1->from->label, Cur1->oppo->from->label); 
  while(Vtx->next != NULL){ 
    Vtx = Vtx->next; 
    Cur1 = Vtx->root; 
    while(Cur1 != Vtx->root->prev){ 
      if (Cur1->label > 0) 
	printf(" %ld - %ld \n", Cur1->from->label, Cur1->oppo->from->label); 
      Cur1 = Cur1->next; 
    }
    if (Cur1->label > 0) 
      printf(" %ld - %ld \n", Cur1->from->label, Cur1->oppo->from->label); 
  }
  printf("\n");
}

/******************************/
/* this function display the dual as edge list (graph structure only) */
/******************************/
void pmPrintDualEdgeList(pm_vertex *Face) 
{
  pm_edge *Cur1; 
  
  //  printf("Face %d: ", Face->label); 
  printf("\n");
  Cur1 = Face->root; 
  while(Cur1 != Face->root->prev->oppo){ 
    if (Cur1->label > 0) 
      printf(" %ld - %ld \n", Cur1->face->label, Cur1->oppo->face->label); 
    Cur1 = Cur1->oppo->next; 
  }
  if (Cur1->label > 0) 
    printf(" %ld - %ld \n", Cur1->face->label, Cur1->oppo->face->label); 
  while(Face->next != NULL){ 
    Face = Face->next; 
    //    printf("Face %ld: ", Face->label);
    Cur1 = Face->root; 
    while(Cur1 != Face->root->prev->oppo){ 
      if (Cur1->label > 0) 
	printf(" %ld - %ld \n", Cur1->face->label, Cur1->oppo->face->label); 
      Cur1 = Cur1->oppo->next; 
    }
    if (Cur1->label > 0) 
      printf(" %ld - %ld \n", Cur1->face->label, Cur1->oppo->face->label); 
  }
  printf("\n");
}


/******************************/
/* this function display the map in GraphCode */
/******************************/
void pmPrintGraphCode(pm_edge *Root) 
{
  pm_edge *Cur1; 
  pm_vertex *Vtx  = Root->from;
  pm_vertex *Face = Root->face;
  long mark = pmNewMark();
  long i,d;
  
  printf("v 0 0 0\n");
  while(Vtx->next != NULL){
    printf("v 0 0 0\n");
    Vtx = Vtx->next;
  }
  printf("\n");  
  for (Cur1 = Root->Next, d=1; Cur1 != Root; Cur1 = Cur1->Next, d++);
  Root->from->label = 1; Root->from->mark = mark;
  printf("vt %f %f 0\n", cos(2*3.1415/d),sin(2*3.1415/d));
  for (Cur1 = Root->Next, i=2; Cur1 != Root; Cur1 = Cur1->Next){
    if (Cur1->from->mark != mark){
      Cur1->from->label = i++;
      Cur1->from->mark = mark;
      printf("vt %f %f 0\n", cos(2*3.1415/d*i),sin(2*3.1415/d*i));
    }
  }  
  Vtx = Root->from;
  while(Vtx->next != NULL){
    Vtx = Vtx->next;
    if (Vtx->mark != mark){
      Vtx->label = i++;
      printf("vt 0 0 0\n");
    }
  }
  printf("\n");

  while(Face->next != NULL){ 
    Face = Face->next; 
    printf("f %ld/%ld ", Face->root->from->label, Face->root->from->label);
    for(Cur1 = Face->root->Next; Cur1 != Face->root; Cur1 = Cur1->Next)
      printf("%ld/%ld ", Cur1->from->label, Cur1->from->label);
    printf("\n");
  }
  printf("\n");
}


/******************************/
/* this function display the dual in GraphCode */
/******************************/
void pmPrintDualGraphCode(pm_edge *Root) 
{
  pm_edge *Cur1; 
  pm_vertex *Vtx  = Root->from;
  pm_vertex *Face = Root->face;
  long mark = pmNewMark();
  long i,d;
  
  printf("v 0 0 0\n");
  while(Face->next != NULL){
    printf("v 0 0 0\n");
    Face = Face->next;
  }
  printf("\n");  
  for (Cur1 = Root->next, d=1; Cur1 != Root; Cur1 = Cur1->next, d++);
  Root->face->label = 1; Root->face->mark = mark;
  printf("vt %f %f 0\n", cos(2*3.1415/d),sin(2*3.1415/d));
  for (Cur1 = Root->next, i=2; Cur1 != Root; Cur1 = Cur1->next){
    if (Cur1->face->mark != mark){
      printf("vt %f %f 0\n", cos(2*3.1415/d*i),sin(2*3.1415/d*i));
      Cur1->face->label = i++;
      Cur1->face->mark = mark;
    }
  }  
  Face = Root->face;
  while(Face->next != NULL){
    Face = Face->next;
    if (Face->mark != mark){
      Face->label = i++;
      printf("vt 0 0 0\n");
    }
  }
  printf("\n");

  while(Vtx->next != NULL){ 
    Vtx = Vtx->next; 
    printf("f %ld/%ld ", Vtx->root->face->label, Vtx->root->face->label);
    for(Cur1 = Vtx->root->next; Cur1 != Vtx->root; Cur1 = Cur1->next)
      printf("%ld/%ld ", Cur1->face->label, Cur1->face->label);
    printf("\n");
  }
  printf("\n");
}


/******************************/
/* this function display the map for LEDA xlman dimacs  */
/******************************/
void pmPrintEdgeListLEDA(pm_vertex *Vtx) 
{
  pm_edge *Cur1; 
  pm_vertex *Vtx1 = Vtx;
  long n=0, i=0;

  printf("c LEDA maxflow problem\n");
  while (Vtx1 != NULL){
    i++;
    for (Cur1 = Vtx1->root->next; Cur1 != Vtx1->root; Cur1 = Cur1->next)
      i++;
    Vtx1 = Vtx1->next;
    n++;
  }
  printf("p max %ld %ld\n", n,i/2);
  printf("n 1 s\n");
  printf("n %ld t\n",n);
  
  Cur1 = Vtx->root; 
  while(Cur1 != Vtx->root->prev){ 
    if (Cur1->label > 0) 
      printf("a %ld %ld 0\n", Cur1->from->label, Cur1->oppo->from->label); 
    Cur1 = Cur1->next; 
  }
  if (Cur1->label > 0) 
    printf("a %ld %ld 0\n", Cur1->from->label, Cur1->oppo->from->label); 
  while(Vtx->next != NULL){ 
    Vtx = Vtx->next; 
    Cur1 = Vtx->root; 
    while(Cur1 != Vtx->root->prev){ 
      if (Cur1->label > 0) 
	printf("a %ld %ld 0\n", Cur1->from->label, Cur1->oppo->from->label); 
      Cur1 = Cur1->next; 
    }
    if (Cur1->label > 0) 
      printf("a %ld %ld 0\n", Cur1->from->label, Cur1->oppo->from->label); 
  }
}


/******************************/
/* this function display the dual for LEDA xlman dimacs  */
/******************************/
void pmPrintDualEdgeListLEDA(pm_vertex *Face) 
{
  pm_edge *Cur1; 
  pm_vertex *Vtx1 = Face;
  long n=0, i=0;

  printf("c LEDA maxflow problem\n");
  while (Vtx1 != NULL){
    i++;
    for (Cur1 = Vtx1->root->Next; Cur1 != Vtx1->root; Cur1 = Cur1->Next)
      i++;
    Vtx1 = Vtx1->next;
    n++;
  }
  printf("p max %ld %ld\n", n,i/2);
  printf("n 1 s\n");
  printf("n %ld t\n",n);

  Cur1 = Face->root; 
  while(Cur1 != Face->root->prev->oppo){ 
    if (Cur1->label > 0) 
      printf("a %ld %ld 0\n", Cur1->face->label, Cur1->oppo->face->label); 
    Cur1 = Cur1->oppo->next; 
  }
  if (Cur1->label > 0) 
    printf("a %ld %ld 0\n", Cur1->face->label, Cur1->oppo->face->label); 
  while(Face->next != NULL){ 
    Face = Face->next; 
    Cur1 = Face->root; 
    while(Cur1 != Face->root->prev->oppo){ 
      if (Cur1->label > 0) 
	printf("a %ld %ld 0\n", Cur1->face->label, Cur1->oppo->face->label); 
      Cur1 = Cur1->oppo->next; 
    }
    if (Cur1->label > 0) 
      printf("a %ld %ld 0\n", Cur1->face->label, Cur1->oppo->face->label); 
  }
}


/******************************/
/* this function display the map as adjacency matrix (graph structure) */
/* the syntax is that of octave. switch for minor of laplacian */
/******************************/
void pmPrintAdjacency(pmMap *Map, char minor)
{
  pm_edge *Root, *Cur1;
  pm_vertex *Vtx;
  long i,d,n,shift=0;
  long *Adj;

  Root = Map->root;
  n = Map->v;
  if (minor){
    shift=1;
    printf("# name: MinorMap%ld\n", Map->i + 1);
  }else
    printf("# name: map%ld\n", Map->i + 1);
  printf("# type: matrix\n# rows: %ld\n# columns: %ld\n",n-shift,n-shift);

  Adj = (long *)calloc(n+1,sizeof(long));
  
  Vtx = Root->from;

  if (! minor){
    Cur1 = Vtx->root;
    while(Cur1 != Vtx->root->prev){
      Adj[Cur1->oppo->from->label]++;
      Cur1 = Cur1->next;
    }
    Adj[Cur1->oppo->from->label]++;
    for (i=1 ; i<=n ; i++){
      printf(" %ld", Adj[i]);
      Adj[i]=0;
    }
  }
  while(Vtx->next != NULL){
    Vtx = Vtx->next;
    printf("\n");
    d=1;
    Cur1 = Vtx->root;
    while(Cur1 != Vtx->root->prev){
      Adj[Cur1->oppo->from->label]++;
      d++;
      Cur1 = Cur1->next;
    }
    Adj[Cur1->oppo->from->label]++;
    if (minor) Adj[Vtx->label]-=d;
    for (i=1+shift; i<=n ; i++){
      printf(" %ld", Adj[i]);
      Adj[i]=0;
    }
  }
  printf("\n");
  free(Adj);
}

/******************************/
/* this function display the dual as adjacency matrix (graph structure) */
/* the syntax is that of octave. switch for minor of laplacian */
/******************************/
void pmPrintDualAdjacency(pmMap *Map, char minor)
{
  pm_edge *Root, *Cur1;
  pm_vertex *Fce;
  long i,d,n,shift =0;
  long *Adj;

  Root = Map->root;
  n = Map->f;
  if (minor){
    shift = 1;
    printf("# name: MinorDual%ld\n", Map->i + 1);
  }else
    printf("# name: Dual%ld\n", Map->i + 1);
  printf("# type: matrix\n# rows: %ld\n# columns: %ld\n",n-shift,n-shift);
  
  Adj = (long *)calloc(n+1,sizeof(long));
  
  Fce = Root->face;

  if (! minor){
    Cur1 = Fce->root;
    while(Cur1 != Fce->root->Prev){
      Adj[Cur1->oppo->face->label]++;
      Cur1 = Cur1->Next;
    }
    Adj[Cur1->oppo->face->label]++;
    for (i=1 ; i<=n ; i++){
      printf(" %ld", Adj[i]);
      Adj[i]=0;
    }
  } 
  while(Fce->next != NULL){
    Fce = Fce->next;
    printf("\n");
    d=1;
    Cur1 = Fce->root;
    while(Cur1 != Fce->root->Prev){
      Adj[Cur1->oppo->face->label]++;
      d++;
      Cur1 = Cur1->Next;
    }
    Adj[Cur1->oppo->face->label]++;
    if (minor) Adj[Fce->label]-=d;
    for (i=1+shift; i<=n ; i++){
      printf(" %ld", Adj[i]);
      Adj[i]=0;
    }
  }
  printf("\n");
  free(Adj);
}



/******************************/
/* this function display the map as edge list (graph structure only) */
/******************************/
void pmPrintEdgeListTulip(pm_vertex *Vtx) 
{
  pm_edge *Cur1; 
  pm_vertex *Vtx1;
  
  printf("\n(nodes %ld", Vtx->label);  
  
  if (Vtx1 -> next !=NULL){
    Vtx1 = Vtx->next;
    while (Vtx1->next !=NULL){
      printf(" %ld", Vtx1->label);
      Vtx1 = Vtx1->next;
    }
  }
  printf(" %ld )\n", Vtx1->label);
  
  Cur1 = Vtx->root; 
  while(Cur1 != Vtx->root->prev){ 
    if (Cur1->label > 0) 
      printf("(edge %ld %ld %ld) \n", Cur1->label,
	     Cur1->from->label, Cur1->oppo->from->label); 
    Cur1 = Cur1->next; 
  }
  if (Cur1->label > 0) 
    printf("(edge %ld %ld %ld) \n", Cur1->label,
	   Cur1->from->label, Cur1->oppo->from->label); 
  while(Vtx->next != NULL){ 
    Vtx = Vtx->next; 
    Cur1 = Vtx->root; 
    while(Cur1 != Vtx->root->prev){ 
      if (Cur1->label > 0) 
	printf("(edge %ld %ld %ld) \n", Cur1->label,
	       Cur1->from->label, Cur1->oppo->from->label); 
      Cur1 = Cur1->next; 
    }
    if (Cur1->label > 0) 
      printf("(edge %ld %ld %ld) \n", Cur1->label,
	     Cur1->from->label, Cur1->oppo->from->label); 
  }
  printf("\n");
}

/******************************/
/* this function display the dual as edge list (graph structure only) */
/******************************/
void pmPrintDualEdgeListTulip(pm_vertex *Face) 
{
  pm_edge *Cur1; 
  
  pm_vertex *Fce1;
  
  printf("\n(nodes %ld", Face->label);  
  
  if (Face -> next !=NULL){
    Fce1 = Face->next;
    while (Fce1->next !=NULL){
      printf(" %ld", Fce1->label);
      Fce1 = Fce1->next;
    }
  }
  printf(" %ld )\n", Fce1->label);
  //  printf("Face %d: ", Face->label); 

  Cur1 = Face->root; 
  while(Cur1 != Face->root->prev->oppo){ 
    if (Cur1->label > 0) 
      printf("(edge %ld %ld %ld) \n", Cur1->label,
	     Cur1->face->label, Cur1->oppo->face->label); 
    Cur1 = Cur1->oppo->next; 
  }
  if (Cur1->label > 0) 
    printf("(edge %ld %ld %ld) \n", Cur1->label,
	   Cur1->face->label, Cur1->oppo->face->label); 
  while(Face->next != NULL){ 
    Face = Face->next; 
    //    printf("Face %ld: ", Face->label);
    Cur1 = Face->root; 
    while(Cur1 != Face->root->prev->oppo){ 
      if (Cur1->label > 0) 
	printf("(edge %ld %ld %ld) \n", Cur1->label,
	       Cur1->face->label, Cur1->oppo->face->label); 
      Cur1 = Cur1->oppo->next; 
    }
    if (Cur1->label > 0) 
      printf("(edge %ld %ld %ld) \n", Cur1->label,
	     Cur1->face->label, Cur1->oppo->face->label); 
  }
  printf("\n");
}


/******************************/
/* this function choose the right display */
/******************************/
void pmPrint(pmMap *Map, pmOutput *Outp)
{
  pm_edge *Root;
  Root = Map->root;
  
  if (Outp->transform){
    pmBicolorFaces(Map->root);
    pmEdgeMap(Map);
  }
  switch (Outp->format){
  case 1:
    if (Outp->map) pmPrintChndVtx(Root->from);
    if (Outp->dual) pmPrintChndFaces(Root->face);
  break;
  case 2:
    if (Outp->map) pmPrintEdgeList(Root->from);
    if (Outp->dual) pmPrintDualEdgeList(Root->face);
    break;
  case 3:
    if (Outp->map) pmPrintMap(Root, Map->i);
    if (Outp->dual) pmPrintDualMap(Root, Map->i);
    break;
  case 4:
    if (Outp->map) pmPrintGraphCode(Root);
    if (Outp->dual) pmPrintDualGraphCode(Root);
    break;
  case 5:
    if (Outp->map) pmPrintEdgeListLEDA(Root->from);
    if (Outp->dual) pmPrintDualEdgeListLEDA(Root->face);
    break;
  case 6:
    if (Outp->map) pmPrintAdjacency(Map, FALSE);
    if (Outp->dual) pmPrintDualAdjacency(Map, FALSE);
    break;
  case 7:
    if (Outp->map) pmPrintAdjacency(Map, TRUE);
    if (Outp->dual) pmPrintDualAdjacency(Map, TRUE);
    break;
  case 8:
    if (Outp->map) pmPrintEdgeListTulip(Root->from);
    if (Outp->dual) pmPrintDualEdgeListTulip(Root->face);
    break;
  }
}



