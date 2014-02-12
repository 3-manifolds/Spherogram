#include <stdio.h>
#include <stdlib.h>


/******************************/
/* this is an example of program using the generator */
/* and supposingly providing a user friendly interface... */
/******************************/

#include "PMdef.h"
#include "stats.h"
#include "interface.h"
#include "PMplanmap.h"
#include "PMdisplay.h"

//#include "PMpersonal.c"

int main(int argc, char *argv[])
{
  pmSize Size;
  pmMethod Meth;
  pmOutput Outp;
  pmStats Stat;
  pmMemory Mem;

  pmCumul Cumul={NULL,NULL,NULL,NULL};
   

  pmMap Map;

  long i;

  if (! pmParseArgs(argc, argv,
		    &Size, &Meth, &Outp, &Stat))
    exit(2);

  if (! pmInitRND(&Meth))
    exit(17);
  
  if (! pmSetParameters(&Size, &Meth))
    exit(2);
 
  if (! pmMemoryInit(&Size, &Meth, &Mem))
    exit(19);
  if (! pmExtendMemory(&Size, &Meth, &Mem, Stat.core))
    exit(19);

  for (i=0; i < Stat.nb; i++){
    Map.i=i;
    pmPlanMap(&Size, &Meth, &Mem, &Map);
    pmPrint(&Map, &Outp);
    pmStatistic(&Map, &Stat, &Cumul);
    pmFreeMap(&Map);
  }

  exit(0);
}
