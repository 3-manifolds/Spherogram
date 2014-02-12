#include <stdlib.h>
#include <stdio.h>

#include "PMdef.h"
#include "PMenlight.h"

/* ici on passe au dual, aux triangulations ou autres variantes. */

long pmLabelCanon(pm_edge *Root)
{
  pm_edge *Cur1, *Cur2;
  long mark;
  long vtx=1, edge=1;
  
  mark = pmNewMark();

  Root->from->mark = mark;
  Root->from->label = vtx++;
  Cur1 = Root->prev;
  do {
    Cur1 = Cur1->next;
    if (Cur1->oppo->from->mark != mark){
      Cur1->label       = edge;
      Cur1->oppo->label = -edge;
      edge++; 
    }else if (Cur1->oppo->mark != mark){
      Cur1->label       = edge;
      Cur1->oppo->label = -edge;
      edge++;
      Cur1->mark = mark;
    }
  }while (Cur1 != Root->prev);
  Root->mark = mark;
  
  if (Root->oppo->from->mark != mark)
    Cur1 = Root->oppo->next;
  else 
    Cur1 = Root->next;
  while (Cur1 != Root){
    if (Cur1->from->mark != mark){
      Cur1->from->mark  = mark;
      Cur1->from->label = vtx++;
      Cur2 = Cur1->prev;
      do {
	Cur2 = Cur2->next;
	if (Cur2->oppo->from->mark != mark){
	  Cur2->label       = edge;
	  Cur2->oppo->label = -edge;
	  edge++;
	}else if (Cur2->oppo->from == Cur2->from &&
		  Cur2->oppo->mark != mark){
	  Cur2->label       = edge;
	  Cur2->oppo->label = -edge;
	  edge++;
	  Cur2->mark = mark;
	}
      }while (Cur2 != Cur1->prev);
    }
    if ((Cur1->oppo->mark == mark && Cur1->oppo->from != Cur1->from) ||
	Cur1->oppo->from->mark != mark){
      Cur1->mark = mark;
      Cur1       = Cur1->oppo;
    }
    Cur1 = Cur1->next;
  }
  return(edge);
}

// chain vertices

long pmChainVtx(pm_edge *Root)
{
  pm_edge *Cur1;
  long mark;
  pm_vertex *Vtx;
  long nb=1;
  
  mark = pmNewMark();
  
  Root->from->mark = mark;          
  Vtx = Root->from;
  Cur1 = Root;
  do{
    if (Cur1->oppo->mark == mark){
      Cur1 = Cur1->oppo;
    }else
    if (Cur1->oppo->from->mark != mark){
      Cur1->mark = mark;
      Cur1 = Cur1->oppo;
      Cur1->from->mark  = mark;
      Vtx->next = Cur1->from;
      Vtx = Cur1->from;
      nb++;
    }
    Cur1 = Cur1->next;
  }while (Cur1 != Root);
  Vtx->next = NULL;
  return(nb);
}

// chainage des faces lorsqu'elles existent

long pmChainFaces(pm_edge *Root)
{
  pm_edge *Cur1;
  long mark;
  pm_vertex *Face;
  long nb=1;
  
  mark = pmNewMark();
  
  Root->face->mark = mark;          
  Face = Root->face;
  Cur1 = Root;
  do{
    if (Cur1->oppo->mark == mark){
      Cur1 = Cur1->oppo;
    }else
    if (Cur1->oppo->face->mark != mark){
      Cur1->mark = mark;
      Cur1 = Cur1->oppo;
      Cur1->face->mark  = mark;
      Face->next = Cur1->face;
      Face = Cur1->face;
      nb++;
    }
    Cur1 = Cur1->Next;
  }while (Cur1 != Root);
  Face->next = NULL;
  return(nb);
}

// creation et chainage des faces d'une carte.

void pmMakeaFace(pm_edge *C1, pm_vertex *Fce, long mark, long *nbf){
  pm_edge *C2;
  Fce->mark = mark;
  Fce->label= ++(*nbf);
  C2 = C1->oppo->next;
  while(C2 != C1){
    C2->face = Fce;
    C2 = C2->oppo->next;
  }
}

long pmAddFaces(pm_edge *Root)
{
  pm_edge *Cur1;
  long mark;
  pm_vertex *Face;
  long nbf=0;
  
  mark = pmNewMark();
  
  Face = pmNewFace(Root);
  Cur1 = Root;
  pmMakeaFace(Cur1, Face, mark, &nbf);
  do{
    if (Cur1->oppo->mark == mark){
      Cur1 = Cur1->oppo;
    }else
    if (Cur1->oppo->face == NULL ||
	Cur1->oppo->face->mark != mark){
      Cur1->mark = mark;
      Cur1 = Cur1->oppo;
      Face->next = pmNewFace(Cur1);
      Face       = Face->next;
      pmMakeaFace(Cur1, Face, mark, &nbf);
    }
    Cur1 = Cur1->oppo->next;
  }while (Cur1 != Root);
  Face->next = NULL;
  return(nbf);
}



void pmLabelFaces(pm_vertex *Face){
  int i=1;
  for (; Face->next != NULL; Face = Face->next)
    Face->label = i++;
  Face->label =i;
}
  
    
void pmClearLblFace(pm_vertex *Face)
{
  do{
    Face->label = 0;
    Face = Face->next;
  }while(Face != NULL);
}
void pmClearLblVtx(pm_vertex *Vtx)
{
  do{
    Vtx->label = 0;
    Vtx = Vtx->next;
  }while(Vtx != NULL);
}





long pmBicolorFaces(pm_edge *Root)
{
  pm_edge *Cur1;
  long mark, nb=1;
  short t=1;
  
  mark = pmNewMark();
  
  Root->face->mark = mark;          
  Root->face->type = t;
  Cur1 = Root;
  do{
    if (Cur1->oppo->mark == mark){
      Cur1 = Cur1->oppo;
      t = (t%2) +1;
    }else
    if (Cur1->oppo->face->mark != mark){
      Cur1->mark = mark;
      Cur1 = Cur1->oppo;
      t = (t%2) +1;
      Cur1->face->mark  = mark;
      Cur1->face->type = t;
      if (t%2) nb++;
    }
    Cur1 = Cur1->Next;
  }while (Cur1 != Root);
  return(nb);
}



void pmEdgeMap(pmMap *Map)
{
  pm_edge *Cur1, *Root = Map->root;
  pm_vertex *Face = Root->face;
  pm_vertex *Vtx  = Root->from;
  short t = Face->type;

  //  bicolorFaces(Root);

  Root->prev = Root->Next;
  Root->from = Root->face;
  Cur1 = Root->Next;
  while(Cur1 != Root){
    Cur1->prev = Cur1->Next;
    Cur1->from->root = Cur1;
    Cur1->from = Cur1->face;
    Cur1 = Cur1->Next;
  }
  while(Face->next != NULL){
    Face = Face->next;
    if (Face->type == t){
      Face->root->prev = Face->root->Next;
      Face->root->from->root = Face->root;
      Face->root->from = Face->root->face;
      Cur1 = Face->root->Next;
      while(Cur1 != Face->root){
	Cur1->prev = Cur1->Next;
	Cur1->from->root = Cur1;
	Cur1->from = Cur1->face;
	Cur1 = Cur1->Next;
      }
    }
  }
  
  Vtx->type = DLTD;
  Root->face = Root->oppo->face;
  Root->face->root = Root;
  Root->oppo = Root->next->next;
  Root->next->next->face = Root->next->next->oppo->face;
  Root->next->next->face->root = Root->next->next;
  Root->next->next->oppo = Root;
  Root->next->type = DLTD; Root->next->next->next->type = DLTD;
  while(Vtx->next != NULL){
    Vtx = Vtx->next;
    Vtx->type = DLTD;
    Cur1 = Vtx->root;
    Cur1->face = Cur1->oppo->face;
    Cur1->face->root = Cur1;
    Cur1->oppo = Cur1->next->next;
    Cur1->next->next->face = Cur1->next->next->oppo->face;
    Cur1->next->next->face->root = Cur1->next->next;
    Cur1->next->next->oppo = Cur1;
    Cur1->next->type = DLTD; Cur1->next->next->next->type = DLTD;
  }
  
  Face = Root->from;
  Cur1 = Root->prev;
  Cur1->next = Root;
  while (Cur1 != Root){
    Cur1->prev->next = Cur1;
    Cur1 = Cur1->prev;
  }
  while (Face->next != NULL){
    Face = Face->next;
    if (Face->type == t){
      Cur1 = Face->root->prev;
      Cur1->next = Face->root;
      while (Cur1 != Face->root){
	Cur1->prev->next = Cur1;
	Cur1 = Cur1->prev;
      }
    }
  }
  Map->e = pmLabelCanon(Root);
  Map->v = pmChainVtx(Root);
  Map->f = pmChainFaces(Root);
  pmLabelFaces(Root->face);
}
    


