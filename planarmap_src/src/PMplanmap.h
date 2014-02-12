extern int pmSetParameters(pmSize *Size, pmMethod *Meth);
extern int pmMemoryInit(pmSize *S, pmMethod *Meth, pmMemory *M);
extern int pmExtendMemory(pmSize *S, pmMethod *Meth, pmMemory *M, char OtherReason);
extern int pmPlanMap(pmSize *S, pmMethod *Meth, pmMemory *M, pmMap *Map);
extern int pmFreeMap(pmMap *Map);
