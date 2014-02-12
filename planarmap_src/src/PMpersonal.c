
/* fonctions d'interet personnel de l'auteur... */


/* affichage des aretes d'un arbre par parcours en profondeur.
   On marque toutes les aretes au cours du passage pour eviter
   les problemes de reutilisation de marques. 
   Les feuilles blanches sont marques o les noires *
*/
   
void pmPrintEdgeTree(pm_edge *Root){
  pm_edge *Cur1;
  long mark;
  
  mark = pmNewMark();
  
  Cur1       = Root->next;
  Cur1->mark = mark;
  while (Cur1 != Root){
    if (Cur1->oppo != NULL && Cur1->type != OUTER){
      Cur1 = Cur1->oppo;
      if (Cur1->mark != mark){  /* first time */
	if (Cur1->type == VIRT)
	  printf("[");
	else
	  printf("(");
	Cur1->mark = mark;
      }else 
	if (Cur1->type == VIRT)
	  printf("]");
	else
	  printf(")");
    }else{
      switch (Cur1->type){
      case OUTER : printf("-"); break; /* should not happen */
      case WHITE : printf("o");	break;
      case FREE  : printf("#");	break;
      case ROOT  : printf("@"); break;
      default    : printf("*");	
      }
    }
    Cur1       = Cur1->next;
    Cur1->mark = mark;
  }
  switch (Cur1->type){
  case WHITE : printf("o"); break;
  case FREE  : printf("#"); break;
  case ROOT  : printf("@"); break;
  default    : printf("*");
  }
  printf("\n");
}





void printCirclePacking(pm_edge *Root){
  pm_edge **List;
  pm_edge *Cur1;
  pm_vertex *Vtx1;
  
  long n=0, i=0, j, k=1;

  Vtx1 = Root->from;  
  while (Vtx1 != NULL){
    i++;
    for (Cur1 = Vtx1->root->next; Cur1 != Vtx1->root; Cur1 = Cur1->next)
      i++;
    Vtx1 = Vtx1->next;
    n++;
  }
  for(Cur1 = Root->Next; Cur1 != Root; Cur1=Cur1->Next) k++;
  
  printf("%ld %ld %ld 0\n", n,i,k);
  printf("%ld 0\n", Root->from->label);
  for(Cur1 = Root->Next; Cur1 != Root; Cur1=Cur1->Next)
    printf("%ld 2\n", Cur1->from->label);


  List = (pm_edge **)calloc(i,sizeof(pm_edge *));
  
  Vtx1 = Root->from;
  Cur1 = Vtx1->root; 
  while(Cur1 != Vtx1->root->prev){ 
    if (Cur1->label > 0) 
      List[Cur1->label-1] = Cur1;
    else 
      List[i+Cur1->label] = Cur1;
    Cur1 = Cur1->next; 
  }
  if (Cur1->label > 0) 
    List[Cur1->label-1] = Cur1;
  else 
    List[i+Cur1->label] = Cur1;
  while(Vtx1->next != NULL){ 
    Vtx1 = Vtx1->next; 
    Cur1 = Vtx1->root; 
    while(Cur1 != Vtx1->root->prev){ 
      if (Cur1->label > 0) 
	List[Cur1->label-1] = Cur1;
      else 
	List[i+Cur1->label] = Cur1;
      Cur1 = Cur1->next; 
    }
    if (Cur1->label > 0) 
      List[Cur1->label-1] = Cur1;
    else 
      List[i+Cur1->label] = Cur1;
  }
  printf("\n");
  printf("\n");
  printf("\n");
  for (j=0; j<i; j++){
    printf("%ld ", List[j]->face->label);
    printf("%ld ", List[j]->oppo->face->label);
    if (List[j]->oppo->label > 0) printf("%ld ", List[j]->oppo->label);
    else printf("%ld ", i+List[j]->oppo->label+1);
    if (List[j]->next->label > 0) printf("%ld ", List[j]->next->label);
    else printf("%ld ", i+List[j]->next->label+1);
    if (List[j]->prev->label > 0) printf("%ld ", List[j]->prev->label);
    else printf("%ld ", i+List[j]->prev->label+1);
    printf("\n");
  }
  free(List);
}

void printDualCirclePacking(pm_edge *Root){
  pm_edge **List;
  pm_edge *Cur1;
  pm_vertex *Vtx1;
  
  long n=0, i=0, j, k=1;

  Vtx1 = Root->from;  
  while (Vtx1 != NULL){
    i++;
    for (Cur1 = Vtx1->root->next; Cur1 != Vtx1->root; Cur1 = Cur1->next)
      i++;
    Vtx1 = Vtx1->next;
    n++;
  }
  for(Cur1 = Root->next; Cur1 != Root; Cur1=Cur1->next) k++;
  
  printf("%ld %ld %ld 0\n", i/2+2-n,i,k);
  printf("%ld 0\n", Root->face->label);
  for(Cur1 = Root->next; Cur1 != Root; Cur1=Cur1->next)
    printf("%ld 2\n", Cur1->face->label);


  List = (pm_edge **)calloc(i,sizeof(pm_edge *));
  
  Vtx1 = Root->from;
  Cur1 = Vtx1->root; 
  while(Cur1 != Vtx1->root->prev){ 
    if (Cur1->label > 0) 
      List[Cur1->label-1] = Cur1;
    else 
      List[i+Cur1->label] = Cur1;
    Cur1 = Cur1->next; 
  }
  if (Cur1->label > 0) 
    List[Cur1->label-1] = Cur1;
  else 
    List[i+Cur1->label] = Cur1;
  while(Vtx1->next != NULL){ 
    Vtx1 = Vtx1->next; 
    Cur1 = Vtx1->root; 
    while(Cur1 != Vtx1->root->prev){ 
      if (Cur1->label > 0) 
	List[Cur1->label-1] = Cur1;
      else 
	List[i+Cur1->label] = Cur1;
      Cur1 = Cur1->next; 
    }
    if (Cur1->label > 0) 
      List[Cur1->label-1] = Cur1;
    else 
      List[i+Cur1->label] = Cur1;
  }
  printf("\n");
  printf("\n");
  printf("\n");
  for (j=0; j<i; j++){
    printf("%ld ", List[j]->face->label);
    printf("%ld ", List[j]->oppo->face->label);
    if (List[j]->oppo->label > 0) printf("%ld ", List[j]->oppo->label );
    else printf("%ld ", i+List[j]->oppo->label+1);
    if (List[j]->Next->label > 0) printf("%ld ", List[j]->Next->label );
    else printf("%ld ", i+List[j]->Next->label+1);
    if (List[j]->Prev->label > 0) printf("%ld ", List[j]->Prev->label );
    else printf("%ld ", i+List[j]->Prev->label+1);
    printf("\n");
  }  
  free(List);
}
