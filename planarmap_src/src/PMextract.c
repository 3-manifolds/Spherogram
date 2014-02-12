#include <stdlib.h>
#include <stdio.h>

#include "PMdef.h"
#include "PMenlight.h"
#include "PMextract.h"

/* ici les procedures d'extraction :
   en fait c'est toujours la meme, on elimine les cycles
   maximaux de longueur k, 1<=k<=4.
*/


/**********************************************************************/
// elimination des 4-cocycles dans une 4-reguliere sans 2-cocycle 


pm_edge *pmVide4cocycle(pm_edge *Root, pm_edge *Cot1, pm_edge *Cot2, pm_edge *Cot3)
{
  pm_edge *Curr, *Cur1, *Cur2, *Cur3, *Cur4,
       *Inn1, *Inn2, *Inn3, *Inn4;
    
  pm_vertex *Vtx = pmNewVtx(NULL);
    
  Cur1 = pmNewEdge(Vtx,NULL,NULL,Root->oppo,INNER); //remplace Root
  Cur2 = pmNewEdge(Vtx,NULL,Cur1,Cot1->oppo,INNER); //remplace Cot1
  Cur3 = pmNewEdge(Vtx,NULL,Cur2,Cot2->oppo,INNER);
  Cur4 = pmNewEdge(Vtx,Cur1,Cur3,Cot3->oppo,INNER);
  Vtx->root = Cur1; Vtx->label = Root->from->label;
  Cur1->next = Cur4;
  Cur1->prev = Cur2;
  Cur2->prev = Cur3;
  Cur3->prev = Cur4;
  Cur1->label = Root->label;
  Cur2->label = Cot1->label;
  Cur3->label = Cot2->label;
  Cur4->label = Cot3->label;
  // maintenant il faut cauteriser...
  // d'abord l'exterieur
  Vtx = pmNewVtx(NULL);
  Inn1 = pmNewEdge(Vtx,NULL,NULL,Root,INNER); // face Root
  Inn2 = pmNewEdge(Vtx,Inn1,NULL,Cot1,INNER);  
  Inn3 = pmNewEdge(Vtx,Inn2,NULL,Cot2,INNER);
  Inn4 = pmNewEdge(Vtx,Inn3,Inn1,Cot3,INNER);
  Vtx->root = Inn1; Vtx->label = Root->from->label+1;
  Inn1->prev = Inn4; Inn1->next = Inn2;
  Inn2->next = Inn3; Inn3->next = Inn4;
  Inn1->label = Root->oppo->label;
  Inn2->label = Cot1->oppo->label;
  Inn3->label = Cot2->oppo->label;
  Inn4->label = Cot3->oppo->label;
  Inn1->face = Root->oppo->face;
  Inn2->face = Cot1->oppo->face;
  Inn3->face = Cot2->oppo->face;
  Inn4->face = Cot3->oppo->face;
  Root->oppo = Inn1; Cot1->oppo = Inn2;
  Cot2->oppo = Inn3; Cot3->oppo = Inn4;
  Root->face->root = Root; Cot1->face->root = Cot1;
  Cot2->face->root = Cot2; Cot3->face->root = Cot3;
  // puis l'interieur
  Cur1->oppo->oppo = Cur1;
  Cur2->oppo->oppo = Cur2;
  Cur3->oppo->oppo = Cur3;
  Cur4->oppo->oppo = Cur4;
  Cur1->face = pmNewFace(Cur1); Cur1->face->label = Root->face->label;
  for(Curr = Cur1->Next; Curr != Cur1; Curr = Curr->Next)
    Curr->face = Cur1->face;
  Cur2->face = pmNewFace(Cur2); Cur2->face->label = Cot1->face->label;
  for(Curr = Cur2->Next; Curr != Cur2; Curr = Curr->Next)
    Curr->face = Cur2->face;
  Cur3->face = pmNewFace(Cur3); Cur3->face->label = Cot2->face->label;
  for(Curr = Cur3->Next; Curr != Cur3; Curr = Curr->Next)
    Curr->face = Cur3->face;
  Cur4->face = pmNewFace(Cur4); Cur4->face->label = Cot3->face->label;
  for(Curr = Cur4->Next; Curr != Cur4; Curr = Curr->Next)
    Curr->face = Cur4->face;
  
  return(Cur1);  
}


// on commence par une procedure pour decomposer en somme horiz.
// si la composante est indecomposable on renvoie TRUE, sinon FALSE
// et les composantes resultantes sont ajoutees dans Fut

int pmInSum(pm_edge *Root){
  pm_edge *Edge, *Edg1, *Edg2;
  long mark = pmNewMark();
  short indec = TRUE;

  Edg1 = Root;
  Edg2 = Root->next;
  
  // marquage des faces etiquetees 2 et detection des diagonales
  for(Edge = Root->Next;
      Edge != Root->Prev;
      Edge = Edge->Next){
    Edge->oppo->face->mark = mark;
    Edge->oppo->face->root = Edge->Oppo;
  }
      
  for(Edge = Root->next->oppo->Prev;
      Edge != Root->next->next;
      Edge = Edge->Prev){
    if (Edge->Oppo->face->mark == mark){
      // diagonal detected
      pmNewBloc(pmVide4cocycle(Edg1, Edge->oppo->face->root, Edge, Edg2));
      Edg1 = Edge->oppo->face->root->oppo;
      Edg2 = Edge->oppo;
      indec = FALSE;
    } else
      Edge->oppo->face->mark = mark;
  }
  if (!indec) 
    pmNewBloc(pmVide4cocycle(Edg1, Root->prev, Root->next->next, Edg2));
  return(indec);  
}

// formerly local function

int pmCheck1(pm_edge *Edge){
    long label = pmNewLabel();
    pm_edge *Edg1, *Edg2;

    for (Edg1 = Edge->Next; Edg1->from->label == 0;
	 Edg1 = Edg1->Next){
      Edg1->oppo->face->label = label;
      Edg1->oppo->face->root  = Edg1->oppo;
    }  
    for (Edg1 = Edge->oppo->Prev; Edg1->from->label == 0;
	 Edg1 = Edg1->Prev);
    for (; Edg1->oppo->from->label == 0;
	 Edg1 = Edg1->Next){
      for (Edg2 = Edg1->oppo->Prev; 
	   Edg2->from->label == 0 && Edg2 != Edg1->oppo->Next;
	   Edg2 = Edg2->Prev);
      for (; Edg2 != Edg1->oppo; Edg2 = Edg2->Next)
	if (Edg2->oppo->face->label == label &&
	    Edge->oppo->from != Edg2->oppo->from &&
	    Edge->from != Edg2->from){
	  pmNewBloc(pmVide4cocycle(Edge,Edg2->oppo->face->root,Edg2,Edg1));
	  return(1);
	}
    }
    return(0);  printf("kes tu fous la dans check1 ??\n");
}

void pmC3kernel(pm_edge *Root){
  
  pm_edge *Edge, *Edg1;
  
  long level;
  short i;
  long mark = pmNewMark();


  pmResetPost();
  level = 1;
  Root->from->label = level;
  for (i=0, Edge = Root; i < 4; i++, Edge = Edge->next) 
    pmNewPost(Edge);
  while(pmIsPost()){
    pmCopyPostSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
      pmCheck1(Edge);
    }
    level++;
    pmFirstSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
#if 0 // FM's comment: doesn't seem to be useful and is C-false
      if (Edge->from->label == level - 1 &&
	  Edge->oppo->from->label == 0);
#endif
      Edge->oppo->from->label = level;
    }
    pmFirstSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
      for (Edg1 = Edge->oppo->next; Edg1 != Edge->oppo; Edg1 = Edg1->next)
	  if (Edg1->from->label == level &&
	      Edg1->oppo->from->label == 0 &&
	      Edg1->mark != mark){
	    Edg1->mark = mark;
	    pmNewPost(Edg1);
	  }
    }
    //     printf("%ld -> %ld\n",SeedEnd,PostEnd);
  }
}

void pmFull2to3c(pm_edge *Root){
  pm_edge *Edge;
  
  pmClearLblFace(Root->face);
  pmClearLblVtx(Root->from);
  
  pmNewBloc(Root);
  while(pmIsBloc()){
    Edge = pmNextBloc();
    if (Edge->oppo->from != Edge->next->next->oppo->from){
      // non trivial
      if (pmInSum(Edge) && pmInSum(Edge->next)){
	// vraiment indecomposable
	pmNewComp(Edge);    
	pmC3kernel(Edge);
      }
    }
  }
}

pm_edge *pmGet3c(pm_edge *Root){
  
  pmClearLblFace(Root->face);
  pmClearLblVtx(Root->from);
  
  if (pmInSum(Root) && pmInSum(Root->next)){
    pmC3kernel(Root);
    return(Root);
  }else
    return(NULL);
}

/**********************************************************************/
/*************************************************************/
/**********************************************************************/

// Elimination des 2cocycles dans les cubiques

pm_edge *pmVide2cocycle(pm_edge *Root, pm_edge *Cot1)
{
  pm_edge *Curr, *Inn1, *Inn2;
    
  Inn1 = Root->oppo;
  Inn2 = Cot1->oppo;

  Inn1->oppo = Inn2;
  Inn2->oppo = Inn1;
  Root->oppo = Cot1;
  Cot1->oppo = Root;

  Inn1->face = pmNewFace(Inn1); 
  Inn2->face = pmNewFace(Inn2);

  Inn1->face->label = Cot1->face->label;
  for(Curr = Inn1->Next; Curr != Inn1; Curr = Curr->Next)
    Curr->face = Inn1->face;

  Inn2->face->label = Root->face->label;
  for(Curr = Inn2->Next; Curr != Inn2; Curr = Curr->Next)
    Curr->face = Inn2->face;
  
  return(Inn1);  
}

// formerly local fnct
int pmCheck2(pm_edge *Edge){

  pm_edge *Edg1;
  
  Edge->face->root = Edge;
  for (Edg1 = Edge->oppo->Next; Edg1 != Edge->oppo; Edg1 = Edg1->Next)
    if (Edg1->oppo->face->root == Edge)
      pmNewBloc(pmVide2cocycle(Edge,Edg1));
  return(1);
}

void pmTri3kernel(pm_edge *Root){
  
  pm_edge *Edge, *Edg1;
  
  long level;
  short i;
  long mark = pmNewMark();


  pmResetPost();
  level = 1;
  Root->from->label = level;
  for (i=0, Edge = Root; i < 3; i++, Edge = Edge->next) 
    pmNewPost(Edge);
  while(pmIsPost()){
    pmCopyPostSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
      pmCheck2(Edge);
    }
    level++;
    pmFirstSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
#if 0 // FM's comment: idem
      if (Edge->from->label == level - 1 &&
	  Edge->oppo->from->label == 0);
#endif
      Edge->oppo->from->label = level;
    }
    pmFirstSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
      for (Edg1 = Edge->oppo->next; Edg1 != Edge->oppo; Edg1 = Edg1->next)
	  if (Edg1->from->label == level &&
	      Edg1->oppo->from->label == 0 &&
	      Edg1->mark != mark){
	    Edg1->mark = mark;
	    pmNewPost(Edg1);
	  }
    }
    //     printf("%ld -> %ld\n",SeedEnd,PostEnd);
  }
}



pm_edge *pmGet3tri(pm_edge *Root){
  
  pmClearLblFace(Root->face);
  pmClearLblVtx(Root->from);
  
  pmTri3kernel(Root);
  return(Root);
}

void pmFull2to3tri(pm_edge *Root){
  pm_edge *Edge;
  
  pmClearLblFace(Root->face);
  pmClearLblVtx(Root->from);
  
  pmNewBloc(Root);
  while(pmIsBloc()){
    Edge = pmNextBloc();
    pmNewComp(Edge);    
    if (Edge->oppo->from != Edge->next->oppo->from ||
        Edge->oppo->from != Edge->prev->oppo->from){
      // non reduit a un triangle
      pmTri3kernel(Edge);
    }
  }
}



/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

// Elimination des 3cocycles dans les cubiques

pm_edge *pmVide3cocycle(pm_edge *Root, pm_edge *Cot1, pm_edge *Cot2)
{
  pm_edge *Curr, *Cur1, *Cur2, *Cur3,
       *Inn1, *Inn2, *Inn3;
    
  pm_vertex *Vtx = pmNewVtx(NULL);
    
  Cur1 = pmNewEdge(Vtx,NULL,NULL,Root->oppo,INNER); //remplace Root
  Cur2 = pmNewEdge(Vtx,NULL,Cur1,Cot1->oppo,INNER); //remplace Cot1
  Cur3 = pmNewEdge(Vtx,Cur1,Cur2,Cot2->oppo,INNER);
  Vtx->root = Cur1; Vtx->label = Root->from->label;
  Cur1->next = Cur3;
  Cur1->prev = Cur2;
  Cur2->prev = Cur3;
  Cur1->label = Root->label;
  Cur2->label = Cot1->label;
  Cur3->label = Cot2->label;
  // maintenant il faut cauteriser...
  // d'abord l'exterieur
  Vtx = pmNewVtx(NULL);
  Inn1 = pmNewEdge(Vtx,NULL,NULL,Root,INNER); // face Root
  Inn2 = pmNewEdge(Vtx,Inn1,NULL,Cot1,INNER);  
  Inn3 = pmNewEdge(Vtx,Inn2,Inn1,Cot2,INNER);
  Vtx->root = Inn1; Vtx->label = Root->from->label+1;
  Inn1->prev = Inn3; Inn1->next = Inn2; Inn2->next = Inn3; 
  Inn1->label = Root->oppo->label;
  Inn2->label = Cot1->oppo->label;
  Inn3->label = Cot2->oppo->label;
  Inn1->face = Root->oppo->face;
  Inn2->face = Cot1->oppo->face;
  Inn3->face = Cot2->oppo->face;
  Root->oppo = Inn1; Cot1->oppo = Inn2; Cot2->oppo = Inn3; 
  // puis l'interieur
  Cur1->oppo->oppo = Cur1;
  Cur2->oppo->oppo = Cur2;
  Cur3->oppo->oppo = Cur3;
  Cur1->face = pmNewFace(Cur1); Cur1->face->label = Root->face->label;
  for(Curr = Cur1->Next; Curr != Cur1; Curr = Curr->Next)
    Curr->face = Cur1->face;
  Cur2->face = pmNewFace(Cur2); Cur2->face->label = Cot1->face->label;
  for(Curr = Cur2->Next; Curr != Cur2; Curr = Curr->Next)
    Curr->face = Cur2->face;
  Cur3->face = pmNewFace(Cur3); Cur3->face->label = Cot2->face->label;
  for(Curr = Cur3->Next; Curr != Cur3; Curr = Curr->Next)
    Curr->face = Cur3->face;

  return(Cur1);  
}

// formerly local fct
int pmCheck3(pm_edge *Edge){
    long label = pmNewLabel();
    pm_edge *Edg1;

    Edge->face->root = Edge;
    for (Edg1 = Edge->oppo->Next; Edg1 != Edge->oppo; Edg1 = Edg1->Next)
      if (Edg1->oppo->face->root == Edge)
	pmNewBloc(pmVide2cocycle(Edge,Edg1));
    
    for (Edg1 = Edge->Next; Edg1->from->label == 0;
	 Edg1 = Edg1->Next){
      Edg1->oppo->face->label = label;
      Edg1->oppo->face->root  = Edg1->oppo;
    }  
    for (Edg1 = Edge->oppo->Prev; Edg1->from->label == 0;
	 Edg1 = Edg1->Prev);
    for (; Edg1 != Edge->oppo->Prev; Edg1 = Edg1->Next)
      if (Edg1->oppo->face->label == label &&
	  Edge->from != Edg1->from){
	pmNewBloc(pmVide3cocycle(Edge,Edg1->oppo->face->root,Edg1));
	return(1);
      }
    return(0);  printf("kes tu fous la dans check3 ??\n");
}

void pmTri4kernel(pm_edge *Root){
  
  pm_edge *Edge, *Edg1;
  
  long level;
  short i;
  long mark = pmNewMark();

  pmResetPost();
  level = 1;
  Root->from->label = level;
  for (i=0, Edge = Root; i < 3; i++, Edge = Edge->next) 
    pmNewPost(Edge);
  while(pmIsPost()){
    pmCopyPostSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
      pmCheck3(Edge);
    }
    level++;
    pmFirstSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
#if 0 // FM's comment: idem
      if (Edge->from->label == level - 1 &&
	  Edge->oppo->from->label == 0);
#endif
      Edge->oppo->from->label = level;
    }
    pmFirstSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
      for (Edg1 = Edge->oppo->next; Edg1 != Edge->oppo; Edg1 = Edg1->next)
	  if (Edg1->from->label == level &&
	      Edg1->oppo->from->label == 0 &&
	      Edg1->mark != mark){
	    Edg1->mark = mark;
	    pmNewPost(Edg1);
	  }
    }
    //     printf("%ld -> %ld\n",SeedEnd,PostEnd);
  }
}


void pmFull2to4tri(pm_edge *Root){
  pm_edge *Edge;
  
  pmClearLblFace(Root->face);
  pmClearLblVtx(Root->from);
  
  pmNewBloc(Root);
  while(pmIsBloc()){
    Edge = pmNextBloc();
    pmNewComp(Edge);    
    if (Edge->oppo->from != Edge->next->oppo->from ||
        Edge->oppo->from != Edge->prev->oppo->from){
      // non reduit a un triangle
      //      pmNewComp(Edge);    
      pmTri4kernel(Edge);
    }
  }
}



pm_edge *pmGet4tri(pm_edge *Root){
  
  pmClearLblFace(Root->face);
  pmClearLblVtx(Root->from);
  
  pmTri4kernel(Root);
  return(Root);
}





/**********************************************************************/
/**********************************************************************/
/**********************************************************************/



// Elimination des 2cocycles des 4regular

pm_edge *pmVide2cocycle4r(pm_edge *Root, pm_edge *Cot1)
{
  pm_edge *Curr, *Inn1, *Inn2;
    
  Inn1 = Root->oppo;
  Inn2 = Cot1->oppo;

  Inn1->oppo = Inn2;
  Inn2->oppo = Inn1;
  Root->oppo = Cot1;
  Cot1->oppo = Root;

  Inn1->face = pmNewFace(Inn1); 
  Inn2->face = pmNewFace(Inn2);

  Inn1->face->label = Cot1->face->label;
  for(Curr = Inn1->Next; Curr != Inn1; Curr = Curr->Next)
    Curr->face = Inn1->face;

  Inn2->face->label = Root->face->label;
  for(Curr = Inn2->Next; Curr != Inn2; Curr = Curr->Next)
    Curr->face = Inn2->face;
  
  return(Inn1);  
}




int pmCheck4(pm_edge *Edge){

  pm_edge *Edg1;
  
  Edge->face->root = Edge;
  for (Edg1 = Edge->oppo->Next; Edg1 != Edge->oppo; Edg1 = Edg1->Next)
    if (Edg1->oppo->face->root == Edge)
      pmNewBloc(pmVide2cocycle4r(Edge,Edg1));
  return(1);
}


void pmSimplekernel(pm_edge *Root){
  
  pm_edge *Edge, *Edg1;
  
  long level;
  short i;
  long mark = pmNewMark();

  pmResetPost();
  level = 1;
  Root->from->label = level;
  for (i=0, Edge = Root; i < 4; i++, Edge = Edge->next) 
    pmNewPost(Edge);
  while(pmIsPost()){
    pmCopyPostSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
      pmCheck4(Edge);
    }
    level++;
    pmFirstSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
#if 0 // FM's comment: idem
      if (Edge->from->label == level - 1 &&
	  Edge->oppo->from->label == 0);
#endif
      Edge->oppo->from->label = level;
    }
    pmFirstSeed();
    while(pmIsSeed()){
      Edge = pmNextSeed();
      for (Edg1 = Edge->oppo->next; Edg1 != Edge->oppo; Edg1 = Edg1->next)
	  if (Edg1->from->label == level &&
	      Edg1->oppo->from->label == 0 &&
	      Edg1->mark != mark){
	    Edg1->mark = mark;
	    pmNewPost(Edg1);
	  }
    }
    //     printf("%ld -> %ld\n",SeedEnd,PostEnd);
  }
}



pm_edge *pmGetsimple4r(pm_edge *Root){
  
  pmClearLblFace(Root->face);
  pmClearLblVtx(Root->from);
  
  pmSimplekernel(Root);
  return(Root);
}

void pmFull4rtosimple(pm_edge *Root){
  pm_edge *Edge;
  
  pmClearLblFace(Root->face);
  pmClearLblVtx(Root->from);
  
  pmNewBloc(Root);
  while(pmIsBloc()){
    Edge = pmNextBloc();
    pmNewComp(Edge);    
    if (Edge->from != Edge->oppo->from ||
        Edge->next->next->from != Edge->next->next->oppo->from){
      // non reduit a un sommet et 2 boucles
      //      pmNewComp(Edge);    
      pmSimplekernel(Edge);
    }
  }
}

