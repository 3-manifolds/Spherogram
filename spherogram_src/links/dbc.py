def double_branched_cover(M):
    """
    The double branched cover of a (2,2)-tangle, given as the link complement of its braid closure (denominator closure) and axis.
    The axis must be the last component. 
    
    >>> L=Manifold('L10nXXX')  #find the right link  NOT THIS ONE!!!
    >>> M=double_branched_cover(L)
    >>> M.identify()

    """
    
    #later change input from M, link exterior (annular closure of T), to T itself and call annular closure function
    num_components=M.num_cusps()
    
    #assume last component is one we don't fill
    fill_list=[(2,0)]*(num_components-1)+[(0,0)]
    
    
    M.dehn_fill(fill_list)
    covers_list=M.covers(2)
    
    #Look for cover with 2 unfilled components and all others with filling (1,0)
    good_covers_list=[]
    
    for C in covers_list:
        
        cusp_pairs=C.cusp_info('filling')
        
        #number of (0,0) cusps
        zero_zero_count=0
        #number of (1,0) cusps
        one_zero_count=0
        
        #list of indices where (0,0) cusps occur
        zero_zero_indices=[]
        
        for i in range(len(cusp_pairs)):
            pair=cusp_pairs[i]
            if pair == (1,0):
                one_zero_count+=1
            elif pair == (0,0):
                zero_zero_count+=1
                zero_zero_indices.append(i)
                
            
        if zero_zero_count == 2 and one_zero_count == len(cusp_pairs)-2:
            good_covers_list.append([C,zero_zero_indices])
            
    if len(good_covers_list)>1:
        raise Exception("MULTIPLE GOOD COVERS")
    elif len(good_covers_list)==0:
        raise Exception("NO GOOD COVERS FOUND")
    else:
        
        #cover with two (0,0) cusps and other cusps (1,0)
        good_cover=good_covers_list[0][0]
        
        #location of (0,0) cusps of the good cover
        good_cover_unfilled_indices=good_covers_list[0][1]
        
        #fill one of the unfilled cusps to get double branched cover of the tangle 
        good_cover.dehn_fill((1,0),good_cover_unfilled_indices[0])
        
        #retriangulate as a one-cusped manifold
        good_cover_f=good_cover.filled_triangulation()
        
        return(good_cover_f)