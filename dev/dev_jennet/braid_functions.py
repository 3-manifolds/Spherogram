from __future__ import print_function
"""
Initialize a link from a sage braid word.
"""
import spherogram
import sage.groups.braid as braid
def braidwordToCrossings(braidword):
    
    braidgens = list(braidword.parent().gens())
    strands = [spherogram.Strand('s'+repr(i)) for i in range(len(braidgens)+1)]
 
    s = [(x,0) for x in strands] # start                                                               
    l = [(x,1) for x in strands] # loose ends     
                                                     
    braidsylls = []
    for b in braidword.syllables():
        if b[1]>0:
            braidsylls = braidsylls+[b[0]]*abs(b[1])
        if b[1]<0:
            braidsylls = braidsylls+[b[0]**-1]*abs(b[1])
    xings = [0]*len(braidsylls)

    for i, b in enumerate(braidsylls): 
        # for each syllable, there is a single crossing               
        label = "x"+repr(i)
        xings[i] = spherogram.Crossing(label)
        for j, a in enumerate(braidgens): 
            # j tells us which two strands are crossing                

            if b == a: # if crossing is negative  
                print('-')
                xings[i][1] = l[j][0][l[j][1]]
                xings[i][0] = l[j+1][0][l[j+1][1]]
                l[j]   = (xings[i],2)
                l[j+1] = (xings[i],3)
                    
            if b**-1 == a: # if crossing is positive      
                print('+')
                xings[i][0] = l[j][0][l[j][1]]
                xings[i][3] = l[j+1][0][l[j+1][1]]
                l[j]   = (xings[i],1)
                l[j+1] = (xings[i],2)
                    
    for i in range(len(s)):
        s[i][0][s[i][1]] = l[i][0][l[i][1]]

    crossings = xings+strands
    return crossings
