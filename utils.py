def getIntersect(all_blocks: list):
    '''
    This function takes in a list of lists and groups sub-lists if their intersects is not NULL.
    '''
    # Base case
    if len(all_blocks) == 0:
        return None
    
    intersect = []
    tmp = []
    base = set()
    for cur in range(len(all_blocks)):
        # Does it intersect?
        if base.intersection(all_blocks[cur]) == set() and cur != 0:
            intersect.append(tuple(tmp))
            tmp.clear()
            tmp.append(all_blocks[cur])
            base = set(all_blocks[cur])
        else:
            tmp.append(all_blocks[cur])
            # Find new members
            for _idx in all_blocks[cur]:
                if _idx not in base:
                    base.add(_idx)
    intersect.append(tuple(tmp))
    
    return intersect

