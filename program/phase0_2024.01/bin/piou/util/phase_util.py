import re
import logging
import os
logger = logging.getLogger("piou.util.phase_util")

def get_mass_list(input):
    ''' obtain the list of 'mass' in atomic units associated to each atom.
        will return None if any errors are encountered. '''
    element_mass = get_element_mass_dict(input)
    if element_mass is None:
        return None
    atom_table = input.get_root_entry().get_entry('structure.atom_list.atoms.table')
    if atom_table is None:
        logger.error("atom table is undefined.")
        return None
    natm = atom_table.get_num_data()
    mass = []
    for i in range(natm):
        elem = atom_table.get_data(i,'element')
        if elem is None:
            logger.error("'element' undefined in atom table.")
            return None
        if not elem in element_mass.keys(): 
            logger.error("element : "+elem+ " undefined in element table")
            return None
        mass.append(element_mass[elem])
    return mass

def get_element_mass_dict(input):
    ''' obtain a dictionary for the element and its mass. 
        Will return None if any errors are encountered. '''
    element_table = input.get_root_entry().get_entry('structure.element_list.table')
    if element_table is None:
        logger.error('element table undefined.')
        return None
    map={}
    nelem=element_table.get_num_data()
    for i in range(nelem):
        elem = element_table.get_data(i,'element')
        if elem is None:
            logger.error("'element' undefined in element table")
            return None
        mas = 51577.50
        mass = element_table.get_data(i,'mass')
        if mass is not None:
            mas = mass
        map[elem] = mas
    return map

def get_xctype_from_pp(input):
    ''' obtain the xctype defined in the pp file. Will return None 
        if any errors are encountered, or if the xctypes are inconsistent 
        among the defined pp files. '''
    ppfiles = get_pp_path_as_list(input)
    if ppfiles is None:
        return None
    xc=None
    xctype={}
    for ppfile in ppfiles:
        if len(ppfile)==0:
            continue
        if not os.path.exists(ppfile):
            logger.error("ppfile : "+ppfile+" does not exist")
            return None
        firstline=False
        f = open(ppfile,'r')
        for line in f:
            line = line.strip()
            if line.startswith("#") or line.startswith("*") or len(line)==0:
                continue
            if not firstline:
                firstline=True
                continue
            words=line.split(":")
            if len(words)<2:
                logger.error("second line of ppfile "+ppfile+" is not valid : "+line)
                return None
            xctype[ppfile]=words[0].strip()
            xc = words[0].strip()
            break
        f.close()

    for k0,v0 in xctype.items():
        for k1,v1 in xctype.items():
            if v0 != v1:
                logger.error("XC type of "+k0+" is "+v0+" while that of "+k1+" is "+v1)
                return None
    return xc

def get_natm2(input):
    ''' get the value of 'natm2', which represents the number of atoms 
        with the 'weight' attribute taken into account.  Will return None
        if any errors are encountered. '''
    atom_table = input.get_root_entry().get_entry('structure.atom_list.atoms.table')
    if atom_table is None:
        return None
    natm = atom_table.get_num_data()
    weight = get_attributes(atom_table,'weight')
    natm2 = 0
    for i in range(natm):
        wei = 1
        if weight is not None and weight[i]==2:
            wei=2
        natm2 += wei
    return natm2
    
def ppfiles_exist(input):
    ''' determine whether ppfiles exist for elements defined in the element_list block '''
    ppfiles = get_pp_path_as_list(input)
    elem_table = input.get_root_entry().get_entry('structure.element_list.table')
    nelem = elem_table.get_num_data()
    for i in range(nelem):
        element = elem_table.get_data(i,'element')
        ii = map_elemname_to_id(elem_table,element)
        ppfile = ppfiles[ii-1]
        if not os.path.exists(ppfile):
            logger.error("ppfile : "+ppfile+" does not exist")
            return False
    return True

def get_num_valence_electrons(input):
    ''' obtain the number of valence electrons.
        Will return None if any errors are encountered. '''
    ppfiles = get_pp_path_as_list(input)
    if ppfiles is None:
        logger.error("failed to obtain the number of valence electrons")
        return None
    atom_table = input.get_root_entry().get_entry('structure.atom_list.atoms.table')
    elem_table = input.get_root_entry().get_entry('structure.element_list.table')
    natm = atom_table.get_num_data()
    nelem = elem_table.get_num_data()
    valence = {}
    for i in range(nelem):
        element = elem_table.get_data(i,'element')
        ii = map_elemname_to_id(elem_table,element)
        ppfile = ppfiles[ii-1]
        if len(ppfile)==0:
            continue
        if not os.path.exists(ppfile):
            logger.error("ppfile : "+ppfile+" does not exist")
            return None
        f = open(ppfile,'r')
        for line in f:
            line = line.strip()
            if line.startswith("#") or line.startswith("*") or len(line)==0:
                continue
            words = line.split()
            if len(words)<10:
                logger.error("first line of "+ppfile+" is not valid : "+line)
                return None
            try:
                valence[element] = int(words[1])
                break
            except ValueError:
                logger.error("invalid first line : "+line)
                return None
        f.close()

    nv = 0
    for i in range(natm):
        wei = 1.0
        weight = atom_table.get_data(i,'weight')
        if weight is not None and isinstance(weight,int) and weight==2:
            wei = 2.0
        try:
            nv += wei * valence[atom_table.get_data(i,'element')] 
        except KeyError:
            return None

    return nv

def get_pp_path_as_list(input):
    ''' obtain a list of pp files. Will return None if atomic coordinates 
        and element info are  undefined in the input '''
    atom_table = input.get_root_entry().get_entry('structure.atom_list.atoms.table')
    elem_table = input.get_root_entry().get_entry('structure.element_list.table')

    if atom_table is None:
        logger.error('atom table is undefined.')
        return None

    nat = atom_table.get_num_data()
    elem_dict={}
    for i in range(nat):
        elem = atom_table.get_data(i,'element')
        if not elem in elem_dict.keys():
            elem_dict[elem] = False

    nelem = elem_table.get_num_data()
    ppfiles = []
    for i in range(nelem):
        elemname = elem_table.get_data(i,'element')
        if elemname is None:
            logger.error("element name undefined")
            return None

        if elemname in elem_dict.keys():
            elem_dict[elemname] = True
        ii = map_elemname_to_id(elem_table,elemname)
        if ii is None:
            return None
        ppfiles.append(input.get_filenames_data().get_file_name('F_POT('+str(ii)+')'))

    return ppfiles

def map_elemname_to_id(elem_table,elemname):
    ''' map the element name and its id. will return None  
        if 'elemname' cannot be found in elem_table.'''
    nelem = elem_table.get_num_data()
    for i in range(nelem):
        elem = elem_table.get_data(i,'element')
        if elem.lower() != elemname.lower():
            continue
        ii = elem_table.get_data(i,'id')
        if ii is None:
            ii = elem_table.get_data(i,'no')
        if ii is None:
            ii = i+1
        if ii<1 or ii>16:
            logger.error("invalid element id : "+str(ii))
            return None
        return ii
    return None

def get_attributes(thetable,attrname):
    ''' extract the specified attribute as a list from the specified table.'''
    if thetable is None:
        return None

    if not thetable.valid_identifier(attrname):
        #logger.error("attribute : "+attrname+" undefined in ["+thetable.get_long_name()+"]")
        return None

    n=thetable.get_num_data()
    ret=[]
    for i in range(n):
        ret.append(thetable.get_data(i,attrname))
    return ret

def get_atom_attributes(input,attrname):
    ''' obtain the specified attribute for the atoms as a list.
        will return None if an error is encountered, or if the 
        specified attribute is undefined. '''
    atom_table = input.get_root_entry().get_entry('structure.atom_list.atoms.table')
    if atom_table is None:
        return None
    return get_attributes(atom_table,attrname)

def is_md(method):
    ''' return True if the specified method corresponds to 
        a molecular dynamics simulation '''
    if method is None:
        return False
    return method=='velocity_verlet' or method=='temperature_control' or method=='velocity_scaling'

def get_default_ppfile_from_dir(dire,element_name,relative=False):
    ''' obtain the path to the default pp file corresponding to the specified element. 
        will return None if the pp file could not be resolved. if relative is set to True,
        will return the file name by relative path with respect to cwd. '''
    dirlist = os.listdir(dire)
    pps=[]
    if element_name is None:
        return None
    for file in dirlist:
        if file.startswith(element_name+"_"):
            ppfile_elem = file.split('.')[0].split('_')
            id = ppfile_elem[len(ppfile_elem)-1]
            gga_lda=file.split('_ggapbe_')
            if len(gga_lda)>=2 and file.endswith('.pp') and len(id)<=2:
                pps.append(file)
    if len(pps)>=1:
        return pps[len(pps)-1]
    return None

def contains_uspp(input):
    pps = get_pp_path_as_list(input)
    for pp in pps:
        if len(pp.split('_us_'))>1 or len(pp.split('_paw_'))>1:
            return True
    return False

class PPFile:
    xctype = None
    pptype = None
    path = None
    def __init__(self,path,xctype,pptype):
        self.path = path
        if xctype=='ggapbe' or xctype=='ldapw91':
            self.xctype = xctype
        if pptype == 'nc' or pptype=='us' or pptype=='paw':
            self.pptype = pptype

    def get_pp_path(self):
        return self.path

    def get_xc_type(self):
        return self.xctype

    def get_pp_type(self):
        return self.pptype

