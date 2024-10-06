import logging
logger = logging.getLogger('util')

def resolve_int_list(intlist,minVal=None,maxVal=None):
    ''' resolve integer list. for example, if the given string is
       "0-2,4-5,3-6,9-11", then this function will return 
       [0,1,2,3,4,5,6,9,10,11]. will return None if any errors are encountered. 
        arg intlist  : the 'integer list string'  to be resolved. '''
    spllist = intlist.split(",")
    ret = []
    for bb in spllist:
        if len(bb.strip())==0:
            continue
        b = bb.split('-')
        intblock = []
        for ba in b:
            if len(ba.strip())==0:
                continue
            try:
                intblock.append(int(ba))
            except ValueError:
                return None
        if len(intblock)==0:
            continue
        if len(intblock)==1:
            ret.append(intblock[0])
        else:
            intblock = [intblock[0],intblock[1]]
            intblock.sort()
            ndata = intblock[1]-intblock[0]+1
            for i in range(ndata):
                ret.append(i+intblock[0])
        if maxVal is not None:
            for i in ret:
                if i>maxVal:
                    return None
        if minVal is not None:
            for i in ret:
                if i<minVal:
                    return None
    ret = list(set(ret))
    ret.sort()
    return ret

