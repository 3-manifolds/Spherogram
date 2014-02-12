#include <stdlib.h>
#include <stdio.h>

#include "PMdef.h"
#include "PMconjugation.h"

/* ici on construit des mots de lukacievicz */


/* tirage d'un mot et calcul a la volee de sa factorisation
 *
 * long luka1 (long n,        nombre de noeuds
 *             long k,        l'arite des noeuds (si 0 => variable)
 *             char LkWrd[])  le mot resultat
 * renvoie la position de depart du mot conjugue dans LkWrd.
 * il faut que le tableau LkWrd soit assez grand !
 * Les lettres sont 'a'+k pour l'arite k.
 */

long pmLuka1 (long n, long k, char *LkWrd )
{
  /* ici on fait un unranking tout bete. */
  /* au fur et a mesure on calcule la hauteur */
  long s,l=k*n+1;
  long start=0, min=-1, h=-1;
  LkWrd[l]='\0';
  for (s=l-1; s>=0; s--){
    if (pmRandom(s+1)<=n) {
      LkWrd[s]=k+'a'; n--;
      h-=k-1;
    }else {
      LkWrd[s]=0+'a'; 
      h++;
    }
    if (h<=min) { min=h; start=s; }
  }    
  return start;
}

/*
 * long luka2 (long l,        la longueur du mot ; ie somme des DgArr[k]
 *	       long DgArr[],  table des degres dans le cas variable
 *             char LkWrd[])  le mot resultat
 * renvoie la position de lecture du debut du mot conjugue.
 *
 * il faut que le tableau LkWrd soit assez grand !
 * il faut aussi que la somme des degres des noeuds soit 
 * egale a l.
 */ 

long pmLuka2 (long l, long DgArr[], char LkWrd[] )
{
  /* cette fonction utilise le tirage d'une permutation
     aleatoire en temps lineaire pour trier un tableau
     contenant les lettres. La permutation est construite
     via sa table d'inversion, par transpositions successives. */

  char c;
  long i,k=0;
  long start, min=0, h=0;
  /* construction du mot trie */
  for (i=0 ; i<l; i++){
    while (!(DgArr[k]--)) k++;
    LkWrd[i]='a'+k;
  }    
  /* permutation */
  /* pour i de 1 a l on multiplie par (i,k) avec k<=i */
  for (i=0 ; i<l; i++){
    k=pmRandom(i+1)-1;
    c=LkWrd[k];
    LkWrd[k]=LkWrd[i];
    LkWrd[i]=c;
  }
  /* ici on ne peut pas facilement calculer le point de conjugaison 
     on the fly donc on le fait apres coup */
  for (i=0 ; i<l; i++){
    h+=LkWrd[i]-'a'-1;
    if (h<min) { min=h; start=i+1; }
  }
  return start%l;
}   

/* tirage d'un mot bicolore, par chottin
 *
 * long luka3 (long i,        nombre de noeuds blancs
 *             long j,        nombre de noeuds noirs
 *             char LkWrd[])  le mot resultat
 * renvoie la position de depart du mot conjugue dans LkWrd.
 * il faut que le tableau LkWrd soit assez grand !
 * Attention : chaque lettre est de la forme = (B(b*)a(b*)a)*A
 *    et son arite totale est le nombre de b moins 1.
 */

long pmLuka3 (long i, long j, char LkWrd[] )
{
  /* cette fonction utilise l'algo de chottin, tel que decrit dans
     ma these p151-154. */
  int v;
  long s=i+2*j,t=2*i+j-1,l=3*(i+j)+1,p;
  long start=0, min=0, h=0;
  for (p=0; p<l; p++, s--){
    if (pmRandom(s)<=i) {
      LkWrd[p]='B'; i--;
      for (v=0 ; v<2; v++){
	while (pmRandom(t--)<=j){
	  j--;    
	  LkWrd[++p]='b';
	  h+=2;
	} 
	p++;
	LkWrd[p]='a';
      }
    }else {
      LkWrd[p]='A'; 
      h--;
    }
    if (h<min) { min=h; start=p+1; }
  }    
  LkWrd[l+1]='\0';
  return start%l;
}






/* ici on prend les mots de lukacievicz et on construit les
   arbres associes. */

/*
 * Cette fonction construit l'arbre associe a un 
 * code de lukacievicz normal quelconque. 
 *
 *
 * pm_edge *luka2tree(                 renvoie le 1/2 brin racine
 *                 long st,         position du conjugue dans
 *                 char LkWrd[])    le mot de lukacievicz
 *
 *
 * Attention cet arbre est enracine sur une 1/2 arete supplementaire
 * non donnee par le code. 
 *
 * Attention encore, les feuilles ne sont pas crees :
 * il y a des 1/2 brins.  
 *
 */

pm_edge *pmLuka2tree(long st, char LkWrd[])
{
  pm_edge *Root, *Cur1, *Cur2;
  pm_vertex *Vtx;
  int arity;
  int i;
  
  Root = pmEmptyEdge();
  Root->type = ROOT;
  
  /* creation du sommet racine */

  arity = LkWrd[st]-MSK;
  
  Vtx = pmNewVtx(Root);
  Cur1 = Root;
  while (arity--){
    Cur1->next = pmNewEdge(Vtx,Cur1,NULL,NULL,BLACK);
    Cur1       = Cur1->next;
  }
  Cur1->next = Root;
  Root->prev = Cur1;
  Cur1       = Root->next;
    
  /* creation des autres sommets */
  i=st+1;
  if (!LkWrd[i]) i=0;
  while(i!=st){
    arity = LkWrd[i]-MSK;
    if (arity) {
      Cur1->type = INNER;
      Cur1->oppo = pmNewEdge(NULL,NULL,NULL,Cur1,INNER);
      Cur1       = Cur1->oppo;
      Vtx        = pmNewVtx(Cur1);
      Cur2       = Cur1;
      while (arity--){
	Cur2->next = pmNewEdge(Vtx,Cur2,NULL,NULL,BLACK);
	Cur2       = Cur2->next;	
      }
      Cur2->next = Cur1;
      Cur1->prev = Cur2;
      Cur1       = Cur1->next;
    }else{
      Cur1 = Cur1->next;
      while (Cur1->oppo != NULL)
	Cur1=Cur1->oppo->next;
    }
    i++;
    if (!LkWrd[i]) i=0;
  }
  return Root;
}

/* Cette fonction construit l'arbre bicolore associe a un 
   code de Chottin. */

pm_edge *pmChottin2tree(long st, char LkWrd[])
{
  pm_edge *Cur0, *Cur1, *Cur2;
  pm_vertex *Vtx;
  long mark = pmNewMark();
  int i=st;
  pm_edge Factice1,Factice2;
  /* pour amorcer et areter on utilise deux arete factice */
  Factice1.next = &Factice2;
  Factice1.prev = &Factice2;
  Factice2.next = &Factice1;
  Factice2.prev = &Factice1;
  Factice1.oppo = NULL;
  Cur0 = &Factice2;
  
  do{
    Cur1 = Cur0;
    /* Reconstitution d'une "lettre" */
    while(LkWrd[i++]!='A'){
      if (!LkWrd[i]) i=0;
      Cur1->type = INNER;
      Cur1->oppo = pmNewEdge(NULL,NULL,NULL,Cur1,INNER);
      Cur1       = Cur1->oppo;
      Vtx        = pmNewVtx(Cur1);
      Cur1->next = pmNewEdge(Vtx,Cur1,NULL,NULL,BLACK2);
      Cur1->prev = pmNewEdge(Vtx,NULL,Cur1,NULL,BLACK2);
      Cur1->next->next = pmNewEdge(Vtx,Cur1->next,Cur1->prev,NULL,BLACK);
      Cur1->prev->prev = Cur1->next->next;
      /* fils gauche */
      Cur2 = Cur1->next;
      while (LkWrd[i++]!='a'){
	if (!LkWrd[i]) i=0;
	Cur2->type = INNER;
	Cur2->oppo = pmNewEdge(NULL,NULL,NULL,Cur2,INNER);
	Cur2       = Cur2->oppo;
	Vtx        = pmNewVtx(Cur2);
	Cur2->next = pmNewEdge(Vtx,Cur2,NULL,NULL,BLACK);
	Cur2->prev = pmNewEdge(Vtx,NULL,Cur2,NULL,BLACK);
	Cur2->next->next = pmNewEdge(Vtx,Cur2->next,Cur2->prev,NULL,BLACK2);
	Cur2->prev->prev = Cur2->next->next;
	Cur2 = Cur2->next->next;
      }
      if (!LkWrd[i]) i=0;
      Cur2->mark = mark; //CLSD;
      /* fils droit */
      Cur2 = Cur1->prev;
      while (LkWrd[i++]!='a'){
	if (!LkWrd[i]) i=0;
	Cur2->type = INNER;
	Cur2->oppo = pmNewEdge(NULL,NULL,NULL,Cur2,INNER);
	Cur2       = Cur2->oppo;
	Vtx        = pmNewVtx(Cur2);
	Cur2->next = pmNewEdge(Vtx,Cur2,NULL,NULL,BLACK);
	Cur2->prev = pmNewEdge(Vtx,NULL,Cur2,NULL,BLACK);
	Cur2->next->next = pmNewEdge(Vtx,Cur2->next,Cur2->prev,NULL,BLACK2);
	Cur2->prev->prev = Cur2->next->next;
	Cur2 = Cur2->next->next;
      }
      if (!LkWrd[i]) i=0;
      Cur2->mark = mark; //CLSD;
      /* la suite de la lettre*/
      Cur1 = Cur1->next->next;
    }
    if (!LkWrd[i]) i=0;
    Cur1->mark = mark; //CLSD;
    
    /* on avance jusqu'a la prochaine feuille ouverte */
    Cur0 = Cur0->prev;
    do{
      Cur0 = Cur0->next;
      while(Cur0->oppo != NULL)
	Cur0 = Cur0->oppo->next;
    }while(Cur0->mark == mark);//CLSD);
  }while(Cur0 != &Factice1);   
  Factice2.oppo->oppo = NULL;
  Factice2.oppo->type = ROOT;
  return Factice2.oppo;
}






/* Ici on ajoute les bourgeons de toutes sortes de facons */


/* arbre binaire -> bourgeons pour cartes planaires normale
 *   1 bourgeon par sommets:
 *      3 possibilites au hasard,
 *        sauf si le sommet porte une indication, auquel cas
 *        apres sa racine locale.
 */

void pmSpring1(pm_edge *Root)
{
  pm_edge *Cur1, *Cur2, *Cur3;
  pm_vertex *Vtx;
  
  Cur1 = Root->next;
  while (Cur1 != Root){
    Vtx = Cur1->from;
    switch (Vtx->type){
    case DONE: break;
    case FORCED:
      Vtx->type  = DONE;
      Cur2       = Cur1->from->root;
      Cur3       = pmNewEdge(Vtx,Cur2,Cur2->next,NULL,WHITE);
      Cur2->next->prev = Cur3;
      Cur2->next = Cur3;
      break;
    default: 
      Vtx->type  = DONE;
      switch (pmRandom(3)){
      case 1: Cur2 = Cur1; break;
      case 2: Cur2 = Cur1->next; break;
      case 3: Cur2 = Cur1->prev; 
      }
      Cur3       = pmNewEdge(Vtx,Cur2,Cur2->next,NULL,WHITE);
      Cur2->next->prev = Cur3;
      Cur2->next = Cur3;
    }
    if (Cur1->oppo != NULL)
      Cur1 = Cur1->oppo;
    Cur1 = Cur1->next;
  }
}

/* arbre binaire -> bourgeons pour cartes bicubiques
 *   1 bourgeons par aretes internes
 *      2 possibilites au hasard
 *   je ne sais pas ce que voudrait dire de forcer l'une des deux
 *   mais a tout hasard, la possibilite est laissee
 */
void pmSpring2(pm_edge *Root)
{
  pm_edge *Cur1, *Cur2, *Cur3;
  pm_vertex *Vtx;
  
  Cur1 = Root->next;
  while (Cur1 != Root){
    if (Cur1->oppo != NULL){
      if (Cur1->oppo->from->type != DONE &&
	  Cur1->from->type       != DONE){
	Cur2       = pmNewEdge(NULL,NULL,NULL,Cur1,INNER);
        Vtx        = pmNewVtx(Cur2);
	Vtx->type  = DONE;
	Cur3       = pmNewEdge(Vtx,NULL,NULL,Cur1->oppo,INNER);
	Cur1->oppo->type = INNER;
	Cur1->oppo->oppo = Cur3;
	Cur1->type       = INNER;
	Cur1->oppo       = Cur2;
	if (Cur1->type == FORCED){
	  Cur2->prev = Cur3;
	  Cur3->next = Cur2;
	  Cur2->next = pmNewEdge(Vtx,Cur2,Cur3,NULL,WHITE);
	  Cur3->prev = Cur2->next;
	}else if (Cur3->oppo->type == FORCED){
	  Cur2->next = Cur3;
	  Cur3->prev = Cur2;
	  Cur3->next = pmNewEdge(Vtx,Cur3,Cur2,NULL,WHITE);
	  Cur2->prev = Cur3->next;
	}else if (pmRandom(2)==1){
	  Cur2->prev = Cur3;
	  Cur3->next = Cur2;
	  Cur2->next = pmNewEdge(Vtx,Cur2,Cur3,NULL,WHITE);
	  Cur3->prev = Cur2->next;
	}else{
	  Cur2->next = Cur3;
	  Cur3->prev = Cur2;
	  Cur3->next = pmNewEdge(Vtx,Cur3,Cur2,NULL,WHITE);
	  Cur2->prev = Cur3->next;
	}	  
	Cur1=Cur3;
      }
      Cur1 = Cur1->oppo;
    }
    Cur1 = Cur1->next;
  }
}

/* arbres ternaires -> bourgeons pour non separables
 * 1 bourgeon avant chaque demi arete non terminale
 * Pas de choix.
 */ 

void pmSpring3(pm_edge *Root)
{
  pm_edge *Cur1, *Cur2;
  
  Cur1 = Root->next;
  while (Cur1 != Root){
    if (Cur1->oppo != NULL){
      Cur2       = pmNewEdge(Cur1->from,Cur1->prev,Cur1,NULL,WHITE);
      Cur1->prev->next = Cur2;
      Cur1->prev       = Cur2;
      Cur1->type       = VIRT;
      if (Cur1 == Cur1->from->root) Cur1->from->root=Cur2;
      Cur1->oppo->type = VIRT;
      Cur1             = Cur1->oppo;
    }
    Cur1 = Cur1->next;
  }
}

/* arbres ternaires -> bourgeons pour cubiques
 * deux bourgeons par sommets, diametralement opposes.
 * 2 choix -> possibiliter de forcer.
 */


void pmSpring4(pm_edge *Root)
{
  pm_edge *Cur1, *Cur2, *Cur3, *Cur4, *Cur5, *Cur6;
  pm_vertex *Vtx;
  pm_edge C;
  pm_vertex V;

  C.oppo = Root;
  C.from = &V;
  V.type = DONE;
  Cur1 = &C;

  while (Cur1 != Root){
    Vtx = Cur1->from;
    if (Vtx->type != DONE){
      if (Vtx->type == FORCED)
	Cur2 = Vtx->root;
      else {
	switch (pmRandom(2)){
	case 1: Cur2 = Cur1; break;
	case 2: Cur2 = Cur1->next; 
	}      
      }
      Vtx->type  = DONE;
      Vtx->root  = Cur2;
      Cur3       = pmNewEdge(Vtx,Cur2,NULL,NULL,WHITE);
      Cur4       = pmNewEdge(Vtx,Cur3,Cur2->prev,NULL,VIRT);
      Cur5       = pmNewEdge(NULL,Cur2->next->next,NULL,NULL,WHITE);
      Vtx        = pmNewVtx(Cur5);
      Vtx->type  = DONE;
      Cur6       = pmNewEdge(Vtx,Cur5,Cur2->next,Cur4,VIRT);
      Cur3->next       = Cur4;
      Cur5->next       = Cur6;
      Cur3->prev->next = Cur3;
      Cur4->next->prev = Cur4;
      Cur5->prev->next = Cur5;
      Cur6->next->prev = Cur6;
      Cur4->oppo       = Cur6;
      Cur5->prev->from = Vtx;
      Cur6->next->from = Vtx;
    }
    if (Cur1->oppo != NULL)
      Cur1 = Cur1->oppo;
    Cur1 = Cur1->next;
  }
}


/* arbre binaire -> bourgeons pour cartes bicubiques
 *   2 bourgeons par aretes internes
 *      3 possibilites au hasard
 */
void pmSpring5(pm_edge *Root)
{
  pm_edge *Cur1, *Cur2, *Cur3;
  pm_vertex *Vtx;
  
  Cur1 = Root->next;
  while (Cur1 != Root){
    if (Cur1->oppo != NULL){
      if (Cur1->oppo->from->type != DONE &&
	  Cur1->from->type       != DONE){
	Cur2       = pmNewEdge(NULL,NULL,NULL,Cur1,INNER);
        Vtx        = pmNewVtx(Cur2);
	Vtx->type  = DONE;
	Cur3       = pmNewEdge(Vtx,NULL,NULL,Cur1->oppo,INNER);
	Cur1->oppo->type = INNER;
	Cur1->oppo->oppo = Cur3;
	Cur1->type       = INNER;
	Cur1->oppo       = Cur2;
	switch (pmRandom(3)){
	case 1:
	  Cur2->prev = Cur3;
	  Cur3->next = Cur2;
	  Cur2->next = pmNewEdge(Vtx,Cur2,NULL,NULL,WHITE);
	  Cur3->prev = pmNewEdge(Vtx,Cur2->next,Cur3,NULL,WHITE);
	  Cur2->next->next = Cur3->prev;
	  break;
	case 2:
	  Cur2->next = pmNewEdge(Vtx,Cur2,Cur3,NULL,WHITE);
	  Cur3->prev = Cur2->next;
	  Cur3->next = pmNewEdge(Vtx,Cur3,Cur2,NULL,WHITE);
	  Cur2->prev = Cur3->next;
	  break;
	case 3:
	  Cur2->next = Cur3;
	  Cur3->prev = Cur2;
	  Cur3->next = pmNewEdge(Vtx,Cur3,NULL,NULL,WHITE);
	  Cur2->prev = pmNewEdge(Vtx,Cur3->next,Cur2,NULL,WHITE);
	  Cur3->next->next = Cur2->prev;
	  break;
	}  
	Cur1=Cur3;
      }
      Cur1 = Cur1->oppo;
    }
    Cur1 = Cur1->next;
  }
}


/* ici on fait la cloture des differents types d'arbres.
   les nouvelles aretes sont marquees OUTER.
   En meme temps, dans l'ordre postfixe, on elimine les aretes VIRT.
 */


pm_edge *pmBalance(pm_edge *Root)
{
  pm_edge *Cur1;
  pm_edge *Free = Root;
  long h = 0, min = 0;
  
  Cur1 = Root->next;
  while (Cur1 != Root){
    if (Cur1->oppo != NULL){
      Cur1 = Cur1->oppo;
    }else{
      switch (Cur1->type){
      case WHITE : h++; break;
      case BLACK : 
      case BLACK2: h--; break;
      default    :
	  break;
      }
      if (h<min){ min = h; Free = Cur1;}
    }
    Cur1 = Cur1->next;
  }
  return Free;
}

pm_edge *pmClosure(pm_edge *Free, pmStck *Stk)
{
  pm_edge *Cur1;
  pm_edge *Root;
  pm_vertex *Vtx;
  long deg=1;

  Free->oppo = pmNewEdge(NULL,NULL,NULL,Free,COFREE);
  if (Free->type == BLACK2)  
    Free->type = FREE2;
  else 
    Free->type = FREE;
  Root       = Free->oppo;
  Vtx        = pmNewVtx(Root);
  
  Cur1 = Free->next;
  while (Cur1 != Free){
    if (Cur1->oppo != NULL){
      Cur1 = Cur1->oppo;
      if (Cur1->type == VIRT){
	if (Cur1->oppo->type != VIRT){
	  Cur1->oppo->next->prev = Cur1->oppo->prev;
	  Cur1->oppo->prev->next = Cur1->oppo->next;
	  Cur1->next->prev = Cur1->prev;
	  Cur1->prev->next = Cur1->next;
	  Cur1->type = DLTD;
	}else
	  Cur1->type = DLTD;
      }
    }else{
      switch (Cur1->type){
      case WHITE : 
	pmStckIn(Cur1, Stk); break;
      case BLACK : 
      case BLACK2:
      case ROOT  :
	Cur1->oppo = pmStckOut(Stk); 
	if (Cur1->oppo == NULL){
	  if (Cur1->type == BLACK2)
	    Cur1->type = FREE2;
	  else 
	    Cur1->type = FREE;
	  Cur1->oppo = pmNewEdge(Vtx,NULL,Root,Cur1,COFREE);
	  Root->prev = Cur1->oppo;
	  Root       = Root->prev;
	  deg++;
	}
	else{
	  Cur1->type = OUTER;
	  Cur1->oppo->oppo = Cur1;
	  Cur1->oppo->type = OUTER;
	}
	break;
      default    : 
	  break;
      }
    }
    Cur1 = Cur1->next;
  }
  Free->oppo->next = Root;
  Root->prev       = Free->oppo;
  deg = pmRandom(deg);
  while (deg--) Root = Root->next;
  while (Root->oppo->type == FREE2) Root = Root->next;
  return Root;
}

/* dans certains cas on veut se debarrasser du sommet de degre 2 */

pm_edge *pmSuppress(pm_edge *Root)
{
  if (Root->next == Root->prev){
    Root->oppo->oppo = Root->next->oppo;
    Root->next->oppo->oppo = Root->oppo;
    Root->type = DLTD;
    Root->next->type = DLTD;
    Root->from->type = DLTD;
    Root = Root->next->oppo;
  }
  return Root;
}











