typedef struct cumul{
  long *allDist, *maxDist, *gauss, *maxgauss;
}pmCumul;

extern void pmStatistic(pmMap *Map, pmStats *Stat, pmCumul *Cumul);
extern long pmStatGauss(pmMap *Map); 
