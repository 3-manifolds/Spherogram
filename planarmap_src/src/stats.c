#include <stdio.h>
#include <stdlib.h>

#include "PMdef.h"
#include "stats.h"

void pmStatPrint(long i, char *Name, long *Table)
{
  long d;
  printf("%s%ld:=[", Name, i);
  for (d=1; d<Table[0]; d++)
    if (Table[d]) printf("[%ld,%ld],",d,Table[d]);
  printf("[%ld,%ld]];\n",d,Table[d]);
}

long pmStatMaxGauss(pmMap *Map)
{
  pm_edge *Cur1;
  pmStck Stck;
  long maxSize, size ;
  
  pmNewMark();
  pmCreateStck(Map->e,&Stck);

  maxSize = 2; size = 0;
  for (Cur1 = Map->root; Cur1!=NULL; Cur1 = pmStckOut(&Stck)) {
    if (Cur1->mark != pmCurMark()){
    while(Cur1->mark != pmCurMark()){
	Cur1->mark = pmCurMark();
	Cur1->oppo->mark = pmCurMark();
	size++;
	if (Cur1->next->mark != pmCurMark()) 
	  pmStckIn(Cur1->next, &Stck);
	Cur1 = Cur1->next->next->oppo;
      }
    if (size > maxSize) maxSize = size;
    size=0;
    }
  }
  pmFreeStck(Stck);
  return(maxSize/2);
}

long pmStatGauss(pmMap *Map)
{
  pm_edge *Cur1;
  pmStck Stck;
  long nbComp;
  
  pmNewMark();
  pmCreateStck(Map->e,&Stck);

  nbComp = 0;
  for (Cur1 = Map->root; Cur1!=NULL; Cur1 = pmStckOut(&Stck)) {
    if (Cur1->mark != pmCurMark()){
      nbComp++;
      while(Cur1->mark != pmCurMark()){
	Cur1->mark = pmCurMark();
	Cur1->oppo->mark = pmCurMark();
	if (Cur1->next->mark != pmCurMark()) 
	  pmStckIn(Cur1->next, &Stck);
	Cur1 = Cur1->next->next->oppo;
      }
    }
  }
  pmFreeStck(Stck);
  return(nbComp);
}

void pmStatCumulGauss(long n, long **cumul){
  long *extend, i;

  if (*cumul == NULL){
    *cumul = (long *)calloc(n+1,sizeof(long));
    (*cumul)[0]=n;
  }else if (n > (*cumul)[0]){
    extend = (long *)calloc(n+1,sizeof(long));
    for (i=1; i <= (*cumul)[0]; i++) extend[i]=(*cumul)[i];
    extend[0]=n;
    free(*cumul);
    *cumul=extend;
  }
  (*cumul)[n]++;
}

void pmStatFaceDeg(pmMap *Map, long **deglist)
{
  pm_edge *Cur1;
  pm_vertex *Fce, *Face;
  long number=0, d, dmax=0;
  
  Face = Map->root->face;
  for (Fce = Face; Fce != NULL; Fce = Fce->next) {
    number ++;
    d = 1;
    for (Cur1 = Fce->root; Cur1 != Fce->root->prev->oppo;
	 Cur1 = Cur1->oppo->next)
      d++;
    if (d>dmax) dmax=d;
  }
  *deglist = (long *)calloc(dmax+1,sizeof(long));
  for (Fce = Face; Fce != NULL; Fce = Fce->next) {
    d = 1;
    for (Cur1 = Fce->root; Cur1 != Fce->root->prev->oppo;
	 Cur1 = Cur1->oppo->next)
      d++;
    (*deglist)[d]++;
  }
  (*deglist)[0]=dmax;
}  

void pmStatCumulDist(long *ddist, pmCumul *C){
  long i;
  long *extendmaxi;

  if (C->maxDist == NULL){
    C->maxDist = (long *)calloc(ddist[0]+1,sizeof(long));
    (C->maxDist)[ddist[0]]=1;
    (C->maxDist)[0]=ddist[0];
    C->allDist = ddist;
  }else{
    if (ddist[0]>(C->allDist)[0]){
      extendmaxi = (long *)calloc(ddist[0]+1,sizeof(long));
      for (i=1; i<=(C->maxDist)[0]; i++)
	extendmaxi[i]=(C->maxDist)[i];
      extendmaxi[0]=ddist[0];
      free(C->maxDist);
      C->maxDist=extendmaxi;
      (C->maxDist)[ddist[0]]++;
      for (i=1; i<=(C->allDist)[0]; i++)
	ddist[i]=ddist[i]+(C->allDist)[i];
      free(C->allDist);
      C->allDist=ddist;
    }else{
      (C->maxDist)[ddist[0]]++;
      for (i=1; i<=ddist[0]; i++)
	(C->allDist)[i]=(C->allDist)[i]+ddist[i];
      free(ddist);
    }      
  }
}

void pmStatDistVtx(pm_edge *Root, long **distances){
  long mark = pmNewMark();
  pm_vertex *Vtx;
  long number=0, current, last;
  pm_edge **later, *Cur1;
  long dmax;
  long *ddist;
  for (Vtx = Root->from; Vtx != NULL; Vtx = Vtx->next) 
    number ++;
  later = (pm_edge **)calloc(number+1,sizeof(pm_edge *));
  
  later[0]=Root; last=1;
  Root->from->mark=mark;
  Root->from->label=0;
  if (Root->oppo->from->mark != mark){
    Root->oppo->from->mark = mark;
    Root->oppo->from->label = 1;
    later[last++] = Root->oppo;
  }
  for(current = 0; current < number; current++){
    for (Cur1=later[current]->next; Cur1!=later[current];
	 Cur1=Cur1->next){
      if (Cur1->oppo->from->mark != mark){
	Cur1->oppo->from->mark = mark;
	Cur1->oppo->from->label = Cur1->from->label+1;
	later[last++] = Cur1->oppo;
      }
    }
  }
  free(later);
  dmax = Cur1->from->label;
  ddist = (long *)calloc(dmax+1,sizeof(long));
  ddist[0] = dmax;
  for (Vtx = Root->from->next; Vtx != NULL; Vtx = Vtx->next) 
    ddist[Vtx->label]++;
  *distances = ddist;
}

void pmStatDistDual(pm_edge *Root, long **distances){
  long mark = pmNewMark();
  pm_vertex *Fce;
  long number=0, current, last;
  pm_edge **later, *Cur1;
  long dmax;
  long *ddist;
  for (Fce = Root->face; Fce != NULL; Fce = Fce->next) 
    number ++;
  later = (pm_edge **)calloc(number+1,sizeof(pm_edge *));
  
  later[0]=Root; last=1;
  Root->face->mark=mark;
  Root->face->label=0;
  if (Root->oppo->face->mark != mark){
    Root->oppo->face->mark = mark;
    Root->oppo->face->label = 1;
    later[last++] = Root->oppo;
  }
  for(current = 0; current < number; current++){
    for (Cur1=later[current]->Next; Cur1!=later[current];
	 Cur1=Cur1->Next){
      if (Cur1->oppo->face->mark != mark){
	Cur1->oppo->face->mark = mark;
	Cur1->oppo->face->label = Cur1->face->label+1;
	later[last++] = Cur1->oppo;
      }
    }
  }
  free(later);
  dmax = Cur1->face->label;
  ddist = (long *)calloc(dmax+1,sizeof(long));
  ddist[0] = dmax;
  for (Fce = Root->face->next; Fce != NULL; Fce = Fce->next) 
    ddist[Fce->label]++;
  *distances = ddist;
}


void pmStatistic(pmMap *Map, pmStats *Stat, pmCumul *Cumul)
{
  long *distances, *degrees;
  
  if (Stat->facedeg){
    pmStatFaceDeg(Map, &degrees);
    pmStatPrint(Map->i,"statDegrees", degrees);
    free(degrees);
  }  
  if (Stat->dist == 1 || Stat->dist == 3){
    if (Stat->dist ==1) pmStatDistVtx(Map->root, &distances);
    else pmStatDistDual(Map->root, &distances);
    pmStatCumulDist(distances, Cumul);
    if (Map->i + 1 == Stat->nb){
      pmStatPrint(Map->i+1,"cumulDist", Cumul->allDist);
      pmStatPrint(Map->i+1,"cumulRadius", Cumul->maxDist);
      free(Cumul->maxDist); free(Cumul->allDist);
    }
  }else if (Stat->dist == 2 || Stat->dist == 4){ 
    if (Stat->dist ==2 ) pmStatDistVtx(Map->root, &distances);
    else pmStatDistDual(Map->root, &distances);
    pmStatPrint(Map->i,"statDist", distances);
    free(distances);
  }
  if (Stat->gauss){
    pmStatCumulGauss(pmStatGauss(Map), &(Cumul->gauss));
    if (Map->i + 1 == Stat->nb){
      pmStatPrint(Map->i+1,"cumulGauss", Cumul->gauss);
      free(Cumul->gauss);
    }
  }
  if (Stat->gaussmax){
    pmStatCumulGauss(pmStatMaxGauss(Map), &(Cumul->maxgauss));
    if (Map->i + 1 == Stat->nb){
      pmStatPrint(Map->i+1,"cumulSizeGauss", Cumul->maxgauss);
      free(Cumul->gauss);
    }
  }
}



