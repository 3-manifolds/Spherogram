from __future__ import print_function
"""
Initialize a torus knot
"""
import spherogram
import sage.groups.braid as braid

def torus_knot(name, method='calc'):
    
    """Returns a (p,q)-torus knot, as an instance of the Link class.
    Example:
    >>> k=torus_knot(2,3)
    >>> k.determinant()
    3
    """

    comma=0
    for i, x in enumerate(name):
        if x == ',':
            comma = i

    p = int(name[2:comma])
    q = int(name[comma+1:-1])
        
    if method == 'braid':
        B = braid.BraidGroup(p-1)
        b = B([i for i in range(0,p-1)]*q)
        return spherogram.Link(b)

    else: 
        to_mirror=False
	
        if p==0 or q==0:
            raise Exception("Torus_knot(p,q) requires p, q, non-zero.")
		        
        #If q<0 or p<0 run torus_knot for |p|, |q| and then if p*q<0 take the mirror.
	if p<0 or q<0:
		p=abs(p)
		q=abs(q)
	if p<0 and q<0:
		p=abs(p)
		q=abs(q)
		to_mirror=True
	if p==2:
		our_crossings=list()
		for i in range (q):
			our_crossings.append(spherogram.Crossing(i))
		if q > 1:
			#set up conditions true for all two strand situations if q>1
			our_crossings[0][0]=our_crossings[q-1][1]
			our_crossings[0][3]=our_crossings[q-1][2]
			our_crossings[0][1]=our_crossings[1][0]
			our_crossings[0][2]=our_crossings[1][3]		
                        #set up in between crossings
			for i in range (1, q-1):
				our_crossings[i][1]=our_crossings[i+1][0]
				our_crossings[i][2]=our_crossings[i+1][3]
			if to_mirror:
				return Link(our_crossings).mirror()
			else:
				return Link(our_crossings)
        
		if q==1:
			our_crossings[0][0]=our_crossings[0][1]
			our_crossings[0][2]=our_crossings[0][3]
		
			for i in range (q):
				print(spherogram.Crossing(i).adjacent)
        	if to_mirror:
				return spherogram.Link(our_crossings).mirror()
        	else:
        		return Link(our_crossings)
    
	if p != 2:
		our_crossings=dict()
		for i in range (q):
			for j in range (p-1):
				our_crossings[(i,j)]=spherogram.Crossing((i,j))
        
		#set up connecting ends
		our_crossings[(0,0)][3]=our_crossings[(q-1,0)][2]
		our_crossings[(0,p-2)][0]=our_crossings[(q-1,p-2)][1]
		
		#middle strands of connecting ends
		for i in range (p-2):
			our_crossings[(0,i)][0]=our_crossings[(q-1,i+1)][2]
		
		if q>1:	
			#set up side connections
			for i in range (q-1):
				our_crossings[(i,0)][2]=our_crossings[(i+1,0)][3]
				our_crossings[(i,p-2)][1]=our_crossings[(i+1,p-2)][0]
			
			#set up connections between crossings
            
			for i in range (p-2): #was p-3
				for j in range (q):
					our_crossings[(j, i)][1]=our_crossings[(j, i+1)][3]
            
			for i in range (1, p-1): #was p-2
				for j in range (q-1):
					our_crossings[(j, i)][2]=our_crossings[(j+1, i-1)][0]	
        
		crossings_list=our_crossings.values()
		
        if to_mirror:
            return spherogram.Link(crossings_list).mirror()
        else:
            return spherogram.Link (crossings_list)
