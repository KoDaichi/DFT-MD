''' module for the PHASE input/output (F_DYNM) data. '''
# -*- coding: utf-8 -*-
import re

import unicodedata

import logging
logger = logging.getLogger("piou.config.phase")

import os
from os.path import join

import difflib

from piou.util import pyutil 
from piou.util import phase_util 

from piou.data.constants import *
from piou.data import elements
from piou.data import properties

from piou import config
from piou.config import functions

def resolve_unit(unit_name):
    ''' return the 'type' of the unit associated with the given unit name.
        for example, 'bohr' will return 'length'. will return None if
        no such unit exists.
        parameter : unit_name the name of the unit. '''
    if unit_name == None:
        return None
    for typ,units in unit_spec.items():
        if unit_name.lower() in units:
            return typ
    return None

def get_value(spec,val,index=0):
    ''' obtain converted value according to spec.
        if spec ==None, or if any errors are encountered,
        will return str(val). '''
    if spec==None or not val_type in spec:
        return str(val)

    vtyp = spec[val_type].split(',')[index].strip()
    if vtyp == val_type_float:
        try:
            return float(val)
        except ValueError:
            return str(val)
    if vtyp == val_type_int:
        try:
            return int(val)
        except ValueError:
            return str(val)
    if vtyp == val_type_bool:
        return pyutil.parse_bool(val)
    return str(val)

def get_value_in_default_units(inpelem):
    ''' when a value of a primitive entry is of type float, it can have units.
        this method will return the value converted to the default units, ie in hartree
        for energy, in bohr for length etc. inpelem must be an instance of PrimitiveEntry'''
    if not isinstance(inpelem,PrimitiveEntry):
        return None
    un  = inpelem.get_unit()
    unt = resolve_unit(un)
    val = inpelem.get_value()
    if unt is None:
        return val

    if isinstance(val,float):
        f = get_unit_factor(un,default_units[unt])
        if f is not None:
            return val * get_unit_factor(un,default_units[unt])
        else:
            return val

class Line:
    ''' the class which encapsulates a line of an input file'''
    def __init__(self,line,line_no):
        self._line=line
        self._line_no=line_no

    def get_line_no(self):
        ''' get the no. of the line '''
        return self._line_no

    def get_line(self):
        ''' get the line it self '''
        return self._line

    def set_line(self,line):
        ''' the setter for the line '''
        self._line=line

    def add_line(self,line):
        ''' the line may also be 'appended' '''
        self._line += line
    
    def __str__(self):
        return str(self._line_no)+":"+self._line

class Input(config.AtomConfigGenerator):
    ''' the main class of this module. 

        typical usage are as follows...

        inp = input.Input()

        inp.validate() -> input validation will be performed.
        root_block=inp.get_root_entry() -> obtain the 'root entry' of F_INP.
        control = root_block.get_entry('control') -> obtain the 'control' block
        control.is_block() -> True
        cpuma=control.get_entry('cpumax') -> get an instance of class PrimitiveEntry
        print cpuma.get_value()+cpuma.get_unit() -> get something like 1 day
        cpuma.set_value('2') -> the 'control.cpumax' parameter will become '2'
        ...
        ...
        ...
        input.save() -> saves to F_INP
        input.save_to(file) -> saves to file '''
    root_block=None
    curr_pblock=None
    ljust_len=70
    def __init__(self,specfile_name='phase.spec',specdir=None,coord_file='nfinp.data'):
        ''' parameter : specfile_name the filename of the specification file.
            if given, type conversions of primitive types will be performed.
            parameter : specdir the directory on which the spec file resides.
            defaults to psf.data.inputspec'''
        super().__init__()
        self.nwarn=0
        self.nerror=0
        self.inp_exists=False
        self.filename_data = FileNamesData(self,nfinp=coord_file)
        #if coord_file is None:
        #if not self.filename_data.file_names_data_exists():
        #    return    

        self.filename_data.debug_output()
        if coord_file is None:
            self.file_name = self.filename_data.get_file_name('F_INP')
        else:
            self.file_name = coord_file
        if self.file_name is None:
            return

        if not os.path.exists(self.file_name):
            self.inp_exists = False
        else:
            self.inp_exists = True

        self.root_block = BlockEntry("root",self,0)
        self.curr_pblock = self.root_block
        if specfile_name!=None:
            self.__init_spec(specfile_name,specdir)

        self.__parse_input()
        self.__dump()
        config.AtomConfigGenerator.__init__(self)

    def is_valid(self):
        #return self.inp_exists
        return True

    def get_name(self):
        return 'PHASE input format'

    def get_defaultname(self,args):
        return 'nfinp.data'

    def _gen_atomic_configuration(self):
        if self.get_root_entry().get_entry('structure.atom_list') is None or\
           self.get_root_entry().get_entry('structure.atom_list.atoms') is None:
            logger.error("structure.atom_list.atoms undefined. ")
            return

        coordsys = self.get_root_entry().get_entry('structure.atom_list.coordinate_system')
        cs = config.FRACTIONAL
        if coordsys is not None and coordsys.get_value()=='cartesian':
            cs = config.CARTESIAN

        uni_coords = config.BOHR
        if self.get_root_entry().get_entry('structure.atom_list.atoms').get_default_unit('length')\
           is not None:
            struni = self.get_root_entry().get_entry('structure.atom_list.atoms').get_default_unit('length')
            if struni.lower() == 'angstrom':
                uni_coords = config.ANGSTROM
            elif struni.lower()=='nm':
                uni_coords = config.NM
        
        cell = None
        a_vector = self.get_root_entry().get_entry('structure.unit_cell.a_vector')
        b_vector = self.get_root_entry().get_entry('structure.unit_cell.b_vector')
        c_vector = self.get_root_entry().get_entry('structure.unit_cell.c_vector')
        if a_vector is not None and b_vector is not None and c_vector is not None:
            cell = config.UnitCell(\
                          a_vector=a_vector.get_value(),\
                          b_vector=b_vector.get_value(),\
                          c_vector=c_vector.get_value())
        if cell is None:
            a = self.get_root_entry().get_entry('structure.unit_cell.a')
            b = self.get_root_entry().get_entry('structure.unit_cell.b')
            c = self.get_root_entry().get_entry('structure.unit_cell.c')
            alpha = self.get_root_entry().get_entry('structure.unit_cell.alpha')
            beta  = self.get_root_entry().get_entry('structure.unit_cell.beta')
            gamma = self.get_root_entry().get_entry('structure.unit_cell.gamma')
            if a is not None and b is not None and c is not None and \
               alpha is not None and beta is not None and gamma is not None:
                cell = config.UnitCell(\
                       a=a.get_value(),b=b.get_value(),c=c.get_value(),\
                       alpha=alpha.get_value(),beta=beta.get_value(),gamma=gamma.get_value())

        if cell is not None:
            celluni = config.BOHR 
            cellu = self.get_root_entry().get_entry('structure.unit_cell').get_default_unit('length')
            if cellu is not None and cellu.lower()!='bohr':
                if cellu.lower() == 'angstrom':
                    celluni = config.ANGSTROM
                elif cellu.lower()=='nm':
                    celluni = config.NM
            cell.set_length_unit(celluni) 

        atomconfig = config.AtomConfig(\
        coordinate_system=cs,unit_cell=cell,length_unit=uni_coords)
        
        atm_table = self.get_root_entry().get_entry('structure.atom_list.atoms.table')
        natm = atm_table.get_num_data()
        for i in range(natm):
            elem = atm_table.get_data(i,'element')
            rx = atm_table.get_data(i,'rx')
            ry = atm_table.get_data(i,'ry')
            rz = atm_table.get_data(i,'rz')
            mob = atm_table.get_data(i,'mobile')
            weight = atm_table.get_data(i,'weight')
            atomconfig.add_atom(elem,rx,ry,rz,mobile=mob,weight=weight)
        logger.debug(str(atomconfig))
        atom = atomconfig.get_atom_at(1)

        return [atomconfig]

    def inpfile_exists(self):
        return self.inp_exists

    def get_root_entry(self):
        ''' returns the 'root entry' of the F_INP file. '''
        return self.root_block

    structure_template=\
    "\
structure{\n\
    atom_list{\n\
        atoms{\n\
            #tag element rx ry rz weight mobile\n\
        }\n\
    }\n\
    unit_cell{\n\
    }\n\
    element_list{\n\
        #tag element atomicnumber mass zeta deviation\n\
    }\n\
}\n\
    "
    def replace_atomic_configuration(self,the_atomic_coordinates,frame_no=None,to_file=None):
        atomtable = self.get_root_entry().get_entry('structure.atom_list.atoms.table')
        blunitcell = self.get_root_entry().get_entry('structure.unit_cell')
        elem_table = self.get_root_entry().get_entry('structure.element_list.table')
        if atomtable is None or blunitcell is None or elem_table is None and to_file is not None:
            file=open(to_file,'w')
            file.write(self.structure_template)
            file.flush()
            file.close()
            self.__parse_input(to_file)
            atomtable = self.get_root_entry().get_entry('structure.atom_list.atoms.table')
            blunitcell = self.get_root_entry().get_entry('structure.unit_cell')
            elem_table = self.get_root_entry().get_entry('structure.element_list.table')
        atomtable.reset_data()
        atomic_coordinates = the_atomic_coordinates
        atmblock = self.get_root_entry().get_entry('structure.atom_list.atoms')
        if not isinstance(the_atomic_coordinates,config.AtomConfig):
            frame=0
            if frame_no is not None and frame_no >0:
                frame=frame_no
            elif frame_no is not None and frame_no<0:
                frame = len(the_atomic_coordinates)-1
            atomic_coordinates = the_atomic_coordinates[frame]

        mode = config.FRACTIONAL
        if not atomic_coordinates.is_cell_ready():
            mode=config.CARTESIAN
        elem_names=[]
        unit=config.BOHR
        for atom in atomic_coordinates:
            element = atom.get_element_name()
            if not element in elem_names:
                elem_names.append(element)
            pos = atom.get_pos(mode=mode,to_unit=config.BOHR)
            mob = atom.is_mobile()
            wei = atom.get_weight()
            row = {'element':element,'rx':pos[0],'ry':pos[1],'rz':pos[2],'mobile':mob,'weight':wei}
            atomtable.add_row_data(row)

        if not atomic_coordinates.is_cell_ready():
            mode=config.CARTESIAN
            atlist=self.get_root_entry().get_entry('structure.atom_list')
            atlist.add_entry(PrimitiveEntry('coordinate_system','cartesian'))
            #self.get_root_entry().get_entry('structure.atom_list.atoms').\
            #set_default_unit('angstrom')
        else:
            (avec, bvec, cvec) = atomic_coordinates.get_unit_cell().\
            get_lattice_vector(to_unit=config.BOHR)
            a = VectorEntry('a_vector')
            for av in avec:
                a.add_value(av)
            b = VectorEntry('b_vector')
            for bv in bvec:
                b.add_value(bv)
            c = VectorEntry('c_vector')
            for cv in cvec:
                c.add_value(cv)
            blunitcell.add_entry(a)
            blunitcell.add_entry(b)
            blunitcell.add_entry(c)

        ppdir = properties.Property().get_property('pp.default_pp_dir')
        if ppdir is None:
            logger.warn('failed to resolve the pp directory; \
            an inappropriate file_names.data file will be created')
        elem_table.reset_data()
        ii=1
        for elemname in elem_names:
            if elemname is None or len(elemname)==0:
                continue
            mas = elements.ElementInfo().get_element_attribute(elemname,'mass')
            an = elements.ElementInfo().get_element_attribute(elemname,'atomic_number')
            row = {'element':elemname,'atomicnumber':an,'mass':mas,'zeta':0,'deviation':1.83}
            elem_table.add_row_data(row)
            nam = elements.ElementInfo().get_element_attribute(elemname,'name')
            if nam is None or len(nam)==0:
                logger.error('failed to resolve element : '+elemname)
                continue
            if ppdir is None:
                continue
            pppath = phase_util.get_default_ppfile_from_dir(ppdir,nam)
            if pppath is None:
                logger.warn('failed to resolve pp file for element : '+elemname)
                continue
            logger.debug("relpath "+pyutil.relpath(os.getcwd(),ppdir))
            self.get_filenames_data().set_file_name("F_POT("+str(ii)+")",\
            pyutil.relpath(os.getcwd(),ppdir)+os.path.sep+pppath)
            ii += 1

    def export_atomic_configuration(self,the_atomic_coordinates,to_file=None,frame_no=None,all_frames=False):
        ''' generate & set the structure.atom_list.atoms table '''
        if self.root_block is None:
            self.root_block = BlockEntry("root",self,0)
            self.curr_pblock = self.root_block
        self.replace_atomic_configuration(the_atomic_coordinates,frame_no,to_file)
        if to_file is not None:
            f=open(to_file,mode='w')
            f.write(self.__dump())
            f.flush()
            f.close()

    def get_num_valence_electrons(self):
        return phase_util.get_num_valence_electrons(self)

    def get_filenames_data(self):
        ''' returns an instance of a file_names.data object '''
        return self.filename_data

    def get_file_name(self):
        ''' returns the file name of F_INP '''
        return self.file_name
    
    currdir='.'
    def set_curr_dir(self,currdir):
        ''' set the 'current directory' associated with the present input.
            parameter : currdir the 'current directory'. defaults to '.' '''
        self.currdir=currdir

    def get_curr_dir(self):
        ''' returns the 'current directory' associated with the present input.
            return : the 'current directory' associated with the present input. '''
        return self.currdir
        
    def save(self,logoffset=''):
        ''' save the input to filename.
            if direc is given, will save to the corresponding file under
            that directory'''
        fname = self.currdir + '/'+self.get_file_name()
        logger.info(logoffset+'dumping F_INP to : '+fname)
        try:
            fi = open(fname,'w')
            fi.write(str(self))
            fi.flush()
            fi.close()
        except IOError:
            logger.error(logoffset+'failed write to: '+fname)

    def increment_nerr(self):
        ''' increments the nerror parameter, which counts the number of errors 
            encountered during input validation. '''
        self.nerror += 1

    def increment_nwarn(self):
        ''' increments the nwarn parameter, which counts the number of warnings 
            encountered during input validation. '''
        self.nwarn += 1

    def reset_err_and_warn_count(self):
        ''' reset the 'nerror' and 'nwarn' attributes. '''
        self.nerror=0
        self.nwarn=0

    def __parse_input(self,inpfile=None):
        inpbuf = self.__get_preprocessed_input(inpfile)
        if inpbuf is None:
            return
        logger.debug('')
        logger.debug("parsing input file : " +self.file_name)
        self.__parse_core(inpbuf)
        logger.debug('parsed input')

    def __parse_core(self,inpbuf):
        parsing_table = False
        table = None
        for li in inpbuf:
            line = li.get_line()
            ent = None
            if line.endswith('{'):
                block = BlockEntry(line.split('{')[0],inp=None,line_no=li.get_line_no())
                self.curr_pblock.add_entry(block)
                self.curr_pblock = block
            elif line=='}':
                self.curr_pblock.close_block()
                if self.curr_pblock.get_parent() is None:
                    logger.error('parent block is None at line '+str(li.get_line_no())+'; probably due to excess \'}\'')
                    self.increment_nerr()
                else:
                    self.curr_pblock = self.curr_pblock.get_parent()
                parsing_table = False
                table = None
            elif line.startswith('#units'):
                line = line.replace('#units','')
                ll = line.split()
                for l in ll:
                    self.curr_pblock.set_default_unit(l)
            elif line.startswith('#tag') or line.startswith('#default') or parsing_table:
                parsing_table = True
                if line.startswith('#tag'):
                    if table==None:
                        table = TableEntry(line_no=li.get_line_no())
                    line = line.replace('#tag','').strip()
                    idents = line.split()
                    for ident in idents:
                        table.add_identifier(ident)
                    idents = table.get_identifiers()
                    defs = table.get_defaults()
                    keys = list(defs.keys())
                    for key in keys:
                        if not key in idents:
                            table.add_identifier(key)
                    self.curr_pblock.add_entry(table)
                elif line.startswith('#default'):
                    if table==None:
                        table = TableEntry(line_no=li.get_line_no())
                    line = line.replace('#default','')
                    defs = line.split(',')
                    for defa in defs:
                        defas = defa.strip().split('=')
                        if len(defas)==2:
                            table.add_default(defas[0].strip(),defas[1].strip())
                else:
                    vals = line.split()
                    iden = table.get_identifiers()
                    row = {}
                    count=0
                    defs = table.get_defaults()
                    for id in iden:
                        if len(vals)>count:
                            if vals[count] != '*':
                                row[id] = vals[count]
                            else:
                                if id in defs:
                                    row[id] = defs[id]
                                else:
                                    row[id] = '*'
                        else:
                            if id in defs:
                                row[id] = defs[id]
                            else:
                                row[id] = '*'
                        count += 1
                    table.add_row_data(row)
            else:
                vals = line.split('=')
                if len(vals)==2:
                    if vals[1].strip().startswith("\"") and vals[1].strip().endswith("\""):
                        ent = PrimitiveEntry(vals[0].strip(),line_no=li.get_line_no())
                        ent.set_value(vals[1].strip()[1:-1])
                        ent.set_double_quote(True)
                    else:
                        va = vals[1].strip().split()
                        if len(va)==1:
                            ent = PrimitiveEntry(vals[0].strip(),line_no=li.get_line_no())
                            ent.set_value(va[0].strip())
                        elif len(va) == 2:
                            ent = PrimitiveEntry(vals[0].strip(),line_no=li.get_line_no())
                            ent.set_value(va[0].strip())
                            ent.set_unit(va[1].strip())
                        elif len(va) >= 3:
                            #try:
                                #float(va[1].strip())
                            ent = VectorEntry(vals[0].strip(),line_no=li.get_line_no())
                            for v in va:
                                ent.add_value(v)
                            #except ValueError:
                            #    va[1].strip()
                            #    ent = PrimitiveEntry(vals[0].strip(),line_no=li.get_line_no())
                            #    ent.set_value(va[0].strip())
                            #    ent.set_unit(va[1].strip())
                    if ent != None:
                        self.curr_pblock.add_entry(ent)

    def get_input_as_list(self):
        ''' returns all input entries as a python list. '''
        li = []
        self.__dump(li)
        for l in li:
            if l.is_block() and not l.closed():
                logger.error('unclosed block : ['+l.get_long_name()+'] at line '+str(l.get_line_no()))
                self.increment_nerr()
            logger.debug(l.get_long_name())
        return li

    def validate(self,specfile_name='phase.spec',specdir=None):
        ''' validate the F_INP file according to the specification file.
            the specification files can be found under psf.data.inputspec directory.
            the contents of this file are hopefully self-explanatory...
            will raise an input.InputValidationError if any errors
            are encountered during the validation.

            parameter : specfile_name the name of the specification file.
            defaults to 'phase.spec' '''

#        self.__init_spec(specfile_name,specdir)

        logger.debug('')
        logger.info('validating input...')
        val = self.__validate()
        logger.info('')
        logger.info('input validation done.')
        logger.info('')
        nerr = self.get_nerror()
        #if nerr>0 or not self.filename_data.is_valid():
        if nerr>0:
            raise InputValidationError('input validation failed.')

    def __get_names(self,key,inf):
        ret=[key]
        if 'aliases' in inf:
            alts=inf['aliases']
            for alt in alts:
                ret.append(alt)
        return ret

    def __init_spec(self,specfile_name,sdir=None):
        inpspecdir=sdir
        if sdir is None:
            inpspecdir = join(pyutil.map_package_to_dir('piou'),'config')

        self.infos={}
        self.__init_spec_core(inpspecdir,specfile_name)

        logger.debug('--input specification--')
        keys = list(self.infos.keys())
        for key in keys:
            logger.debug('')
            logger.debug('name: '+key)
            info = self.infos[key]
            ik = list(info.keys())
            for i in ik:
                logger.debug(i+': '+str(info[i]))

    re_map={}
    block_long_name=[]
    nonblock_long_name=[]
    def __init_spec_core(self,dirname,fname):
        logger.debug('input specification dir: '+dirname)
        logger.debug('specification file name: '+fname)
        try:
            f=open(join(dirname,fname),'r',encoding='utf_8_sig')
        except TypeError:
            f=open(join(dirname,fname),'r')
        tmp={}
        alt_names=[]
        begin=False
        nam=''
        self.block_long_name=[]
        self.nonblock_long_name=[]
        for line in f:
            line=line.strip()
            if line.startswith("#"):
                continue
            if line==begin_tag:
                begin=True
                continue
            vals = line.split('=')
            if len(vals)>=2 and begin:
                vals[0]=vals[0].strip().lower()
                vals[1]=vals[1].strip().lower()
                if vals[0]!='name':
                    tmp[vals[0]] = vals[1]
                else:
                    nam=vals[1].replace("\\.",".").replace("$","")
                    self.infos[vals[1]] = tmp
                    self.re_map[vals[1]] = re.compile(vals[1])
                if vals[0]=='aliases':
                    tmp[vals[0]] = vals[1].split(',')
                    for r in tmp[vals[0]]:
                        self.re_map[r] = re.compile(r)
                if vals[0]=='file':
                    subfile=vals[1].strip()
                    self.__init_spec_core(dirname,subfile)

            if line==end_tag:
                if tmp['entry_type']=='block':
                    self.block_long_name.append(nam)
                else:
                    self.nonblock_long_name.append(nam)
                nam=''
                tmp={}
                continue

        f.close()

    def resolve_spec(self,lname):
        ''' returns the specification dictionary associated with
            the input entry whose 'longname' is lname.
            parameter : lname the 'longname' of the input entry. '''
        try:
            self.infos
        except AttributeError:
            return NOT_READY

        for key in self.infos:
            rkey = self.re_map[key]
            if rkey.match(lname):
                return self.infos[key]
            if 'aliases' in self.infos[key]:
                alts=self.infos[key]['aliases']
                for alt in alts:
                    ralt = self.re_map[alt]
                    if ralt.match(lname):
                        return self.infos[key]
        return None

    def __resolve_condition(self,stri):
        list = stri.split("__if__")
        namee = list[0].strip().split("__or__")
        if len(namee)==1:
            namee=list[0].strip().split(",")
        name=[]
        for na in namee:
            name.append(na.strip())
        bcondi=True
        the_value=None
        the_tag = None
        ltag=[]
        lval=[]
        band = False
        count = 0
        l0 = []
        if len(list)>=2:
            bcondi = False
            l0 = list[1].split("&&")
            if len(l0)>=2:
                band=True
            else:
                l0 = list[1].split("||")
            for l1 in l0:
                ll = l1.split("__eq__")
                if len(ll)>=2:
                    lll=ll[1].strip().split("__or__")
                    inp = self.get_root_entry().get_entry(ll[0].strip())
                    if inp is not None and isinstance(inp,PrimitiveEntry):
                        for llll in lll:
                            if str(inp.get_value()).strip().lower()==llll.strip():
                                the_value = llll.strip()
                                the_tag = ll[0].strip()
                                if band:
                                    count += 1
                                    ltag.append(the_tag)
                                    lval.append(the_value)
                                else:
                                    return (True,name,[the_tag],[the_value])
                            elif llll.strip() == '__defined__':
                                the_value = str(inp.get_value()).lower()
                                the_tag = ll[0].strip()
                                if band:
                                    count += 1
                                    ltag.append(the_tag)
                                    lval.append(the_value)
                                else:
                                    return (True,name,[the_tag],[the_value])
                    elif inp is None:
                        for llll in lll:
                            if llll.strip() == '__undefined__':
                                the_tag = ll[0].strip()
                                the_value = 'undefined'
                                if band:
                                    count += 1
                                    ltag.append(the_tag)
                                    lval.append(the_value)
                                else:
                                    return (True,name,[the_tag],[the_value])

  
        if band and count==len(l0):
            bcondi=True

        return (bcondi,name,ltag,lval)
                            
    def _on_off_if_bool(self,boolval):
        if boolval is None:
            return 'off'
        if str(boolval).strip().lower()=='true':
            return 'on'
        if str(boolval).strip().lower()=='false':
            return 'off'
        return str(boolval).strip().lower()

    def __validate(self):
        li = self.get_input_as_list()

        #check for required/recommended entries
        logger.info('')
        logger.info('checking if required/recommended entries exist...')
        logger.info('')
        speckeys = list(self.infos.keys())
        boo = True
        for key in speckeys:
            inf=self.infos[key]
            if importance in inf:
                boo = self.__validate_importance(li,key,inf)
            count=1
            while True:
                if importance+str(count) in inf:
                    boo = self.__validate_importance(li,key,inf,str(count))
                    count+=1
                else:
                    break

        #check validity of each entry.
        logger.info('')
        logger.info('checking whether each entires are valid...')
        logger.info('')
        for inputelem in li:
            ln = inputelem.get_long_name()
            keys = list(self.infos.keys())
            tmpstr = "line "+str(inputelem.get_line_no()).rjust(6)+': '
            found=False
            done_iter=False
            for key in keys:
                if done_iter:
                    break
                inf = self.infos[key]
                nams = self.__get_names(key,inf)
                for nam in nams:
                    rnam = self.re_map[nam]
                    if rnam.match(ln.lower().strip()):
                        found=True
                        logger.debug('')
                        logger.debug(tmpstr+'found ['+ln+'] in specification')
                        inputelem.set_specification(inf)
                        logger.debug(' entry type : '.ljust(20)+inf[entry_type])

                        if description in inf and len(inf[description])!=0:
                            logger.debug('  description:'.ljust(20)+inf[description]) #spec[description])

                        ok = self.__validate_entry(inputelem,inf)

                        pref=tmpstr+' entry ['+ln+']'
                        if ok:
                            if not inputelem.is_block():
                                logger.info(pref.ljust(self.ljust_len)+' matches the specification.')
                        else:
                            logger.error(pref.ljust(self.ljust_len)+' does not match the specification.')
                            logger.error("")
                            boo=False
                        done_iter=True
                        break
            if not found:
                pref = tmpstr+' could not find   ['+ln+']'
                logger.warn('')
                logger.warn(pref.ljust(self.ljust_len)+' in specification.')
                self.nwarn += 1
                boo=False
                close=[]
                if inputelem.get_type()==BLOCK:
                    close = difflib.get_close_matches(ln,self.block_long_name,10,0.85)
                else:
                    close = difflib.get_close_matches(ln,self.nonblock_long_name,10,0.85)

                if len(close)==0:
                    logger.warn("              no candidate entries")
                for cl in close:
                    logger.warn("              candidate entry: ["+cl+"]")
                logger.warn('')
        return boo

    def __validate_entry(self,inputelem,inf):
        ok = False
        etype = inf[entry_type]
        if etype == entry_type_block:
            ok=self.__validate_block(inputelem,inf)
        elif etype == entry_type_primitive:
            ok=self.__validate_primitive(inputelem,inf)
        elif etype == entry_type_table:
            ok=self.__validate_table(inputelem,inf)
        elif etype == entry_type_vector:
            ok=self.__validate_vector(inputelem,inf)

        o = self.__validate_aux('required_attribute',\
                                inputelem,inf,self.__validate_attribute)
        if not o:
            ok = False

        o = self.__validate_aux('recommended_attribute',\
                                inputelem,inf,self.__validate_attribute)
        if not o:
            ok = False

        o = self.__validate_aux('check_existence',\
                                inputelem,inf,self.__validate_file_existence)
        if not o:
            ok = False

        o = self.__validate_aux('invalid_choice',\
                                inputelem,inf,self.__validate_invalid_choice)
        if not o:
            ok = False

        o = self.__validate_aux('restriction',\
                                inputelem,inf,self.__validate_restriction)
        if not o:
            ok = False

        o = self.__validate_aux('function',\
                                inputelem,inf,self.__validate_function)
        if not o:
            ok = False

        return ok

    def __validate_aux(self,name,inpelem,inf,func):
        ret = True
        if name in inf:
            o=func.__call__(name,inpelem,inf)
            if not o:
                ret = o
        count = 1
        while True:
            if name+str(count) in inf:
                o = func.__call__(name+str(count),inpelem,inf)
                if not o:
                    ret = False
                count += 1 
            else:
                break
        return ret

    def __validate_attribute(self,nam,inpelem,inf):
        (bcondi,iden,tag,val) = self.__resolve_condition(inf[nam])
        ok = True
        breq=True
        t = "required"
        if nam.startswith('recommended'):
            breq=False
            t = "recommended"
        wh=''
        if len(tag)!=0:
            wh = ' when '
        if bcondi and not inpelem.valid_identifier(iden[0]):
            msgstr = 'attribute ['+iden[0]+'] for table entry ['+inpelem.get_long_name()+'] is '+t+wh
            for ii in range(len(tag)):
                if ii<len(tag)-1:
                    msgstr += "["+tag[ii]+"]"+" is '"+self._on_off_if_bool(val[ii])+"' and "
                else:
                    msgstr += "["+tag[ii]+"]"+" is '"+self._on_off_if_bool(val[ii])+"'"
            if breq:
                logger.error(msgstr)
                self.increment_nerr()
                if ok:
                    ok=False
            else:
                logger.warn(msgstr)
                self.increment_nwarn()
        return ok

    def __validate_importance(self,li,key,inf,index=''):
        nams = self.__get_names(key,inf)
        condistr = "."
        (bcondi,impor,tag,val) = self.__resolve_condition(inf[importance+index])
        imp = impor[0]
        if len(val)>0:
            condistr = ' when '
            for ii in range(len(tag)):
                if ii<len(tag)-1:
                    condistr += "["+tag[ii]+"]"+" is '"+self._on_off_if_bool(val[ii])+"' and "
                else:
                    condistr += "["+tag[ii]+"]"+" is '"+self._on_off_if_bool(val[ii])+"'"

        found = False
        if imp == importance_required and bcondi:
            for inputelem in li:
                for nam in nams:
                    rnam = self.re_map[nam]
                    if rnam.match(inputelem.get_long_name()):
                        found = True
                        stri ='found ['+inputelem.get_long_name()+'] at line '+str(inputelem.get_line_no())+','
                        logger.info(stri.ljust(self.ljust_len)+' a required entry'+condistr)
                        break
            if not found:
                stri ='could not find ['+key.replace('\\.','.').replace('$','')+'],'
                logger.error(stri.ljust(self.ljust_len)+' a required entry'+condistr)
                self.nerror += 1
                return False
        elif imp == importance_recommended and bcondi:
            for inputelem in li:
                for nam in nams:
                    rnam = self.re_map[nam]
                    if rnam.match(inputelem.get_long_name()):
                        found = True
                        stri = 'found ['+inputelem.get_long_name()+'] at line '+str(inputelem.get_line_no())+','
                        logger.info(stri.ljust(self.ljust_len)+' a recommended entry'+condistr)
                        break
            if not found:
                stri ='could not find ['+key.replace('\\.','.').replace('$','')+'],'
                logger.warn(stri.ljust(self.ljust_len)+' a recommended entry'+condistr)
                self.nwarn+=1
        return True

    def __validate_file_existence(self,name,inputelem,inf):
        (bcondi,files,thetag,theval)=self.__resolve_condition(inf[name])
        ok=True
        if bcondi and inputelem.get_type()==PRIMITIVE:
            for fi in files:
                fname=self.filename_data.get_file_name(fi)
                if fname is None or not os.path.exists(fname):
                    st = "file : "+fi.upper()+" = '"+fname+"' must exist when "
                    for ii in range(len(thetag)):
                        if ii==len(thetag)-1:
                            st += '['+thetag[ii]+']'+" is '"+self._on_off_if_bool(theval[ii])+"'"
                        else:
                            st += '['+thetag[ii]+']'+" is '"+self._on_off_if_bool(theval[ii])+"' and "
                    logger.error(st)
                    self.nerror+=1
                    ok=False
        return ok

    def __validate_invalid_choice(self,name,inputelem,inf):
        (bcondi,vnames,thetag,theval)=self.__resolve_condition(inf[name])
        if bcondi and inputelem.get_type()==PRIMITIVE:
            v = str(inputelem.get_value()).strip().lower()
            for vname in vnames:
                if v==vname.strip().lower():
                    st = "["+inputelem.get_long_name()+"] cannot be '"+self._on_off_if_bool(v)+"' when "
                    for ii in range(len(thetag)):
                        if ii==len(thetag)-1:
                            st += '['+thetag[ii]+']'+" is '"+self._on_off_if_bool(theval[ii])+"'"
                        else:
                            st += '['+thetag[ii]+']'+" is '"+self._on_off_if_bool(theval[ii])+"' and "
                    logger.error(st)
                    self.nerror+=1
                    return False
        return True

    def __validate_restriction(self,name,inputelem,inf):
        (bcondi,vnamee,thetag,theval)=self.__resolve_condition(inf[name])
        if bcondi and inputelem.get_type()==PRIMITIVE:
            rfound=False
            for vname in vnamee:
                if vname.strip()=='__undefined__' :
                    rfound=False
                    break
            for vname in vnamee:
                if vname.strip()==str(inputelem.get_value()).strip().lower():
                    rfound=True
                    break
            if not rfound:
                vnamestr = "["+inputelem.get_long_name()+"]"
                for ii in range(len(vnamee)):
                    if vnamee[ii]=='__undefined__':
                        vnamestr += " should not be defined "
                    else:
                        vnamestr += " should be restricted to '"+self._on_off_if_bool(vnamee[ii])+"' "
                    if ii==len(vnamee)-1:
                        vnamestr += "when "
                    else:
                        vnamestr += "or"
                for ii in range(len(thetag)):
                    if ii==len(thetag)-1:
                        vnamestr += '['+thetag[ii]+']'+" is '"+self._on_off_if_bool(theval[ii])+"'"
                    else:
                        vnamestr += '['+thetag[ii]+']'+" is '"+self._on_off_if_bool(theval[ii])+"' and "
                logger.error("")
                logger.error(vnamestr)
                self.nerror+=1
                return False
        return True

    def __validate_function(self,name,inputelem,inf):
        func = inf[name]
        (bcondi,funcname,thetag,theval)=self.__resolve_condition(func)
        if bcondi:
            fobj = getattr(functions,funcname[0])
            if fobj != None and callable(fobj):
                ok = fobj.__call__({'input':self,'currelem':inputelem})
                return ok
        return True

    def __validate_block(self,inputelem,spec):
        if inputelem.get_type() != BLOCK:
            logger.error("  the type of "+inputelem.get_long_name()+" is not BLOCK")
            self.nerror+=1
            return False
        nam = inputelem.get_long_name()
        logger.debug('  has table  : '.ljust(20)+spec[has_table])

        ents = inputelem.get_entries()
        spechastable = pyutil.parse_bool(spec['has_table'])
        hastable = False
        for ent in ents:
            if not ent.is_block() and ent.get_type()==TABLE:
                hastable = True
                break
        if spechastable==hastable:
            if spechastable:
                logger.debug('  tabular data found under ['+nam+'].')
            else:
                logger.debug('  tabular data not found under ['+nam+'].')
        else:
            if spechastable:
                logger.error('  tabular data should exist under ['+nam+'].')
            else:
                logger.error('  tabular data should not exist under ['+nam+'].')
            self.nerror+=1

        return spechastable==hastable

    def __validate_primitive(self,inputelem,spec):
        if inputelem.get_type()!=PRIMITIVE:
            logger.error("  the type of "+inputelem.get_long_name()+" is not PRIMITIVE")
            logger.error("  the type of "+inputelem.get_long_name()+" : "+str(inputelem.get_type()))
            self.nerror+=1
            return False
        logger.debug('  val_type   : '.ljust(20)+spec[val_type])
        if unit_type in spec:
            logger.debug('  unit_type  :'.ljust(20)+spec[unit_type])

        choices=[]
        if choice in spec:
            logger.debug('  choices    :'.ljust(20)+spec[choice])
            choices = spec[choice].split(',')
        vtyps = spec['val_type'].split('__or__')
        val = inputelem.get_str_value().lower()
        un  = inputelem.get_unit().lower()

        vtyp_passed = self.__check_val_type(vtyps,val)
        vra_passed=True
        if val_range in spec:
            vra_passed  = self.__check_val_range(spec[val_range],val)

        unitok = True
        if unit_type in spec and len(un)!=0:
            unitok=False
            ul = unit_spec[spec[unit_type]]
            if not un in ul:
                logger.error('  the unit ['+un+'] is not type ['+spec[unit_type]+']')
                self.nerror+=1
            else:
                unitok=True
                logger.debug('  unit_type is OK.')

        choiok=True
        if len(choices)!=0:
            choiok=False
            for choi in choices:
                if choi.strip() == val.lower():
                    choiok=True
                    break
        if len(choices)!=0 and choiok:
            logger.debug('  found '+val+' in choice.')
        elif len(choices)!=0 and not choiok:
            logger.error('  could not find     '+val+' in choice.')
            self.nerror+=1
            close = difflib.get_close_matches(val,choices,1,0.8)
            if close!=None and len(close)==1:
                logger.error('  candidate choice : '+close[0])
            else:
                logger.error('  no candidate choice')
        return vtyp_passed and unitok and choiok and vra_passed

    def __validate_table(self,inputelem,spec):
        if(inputelem.get_type()!=TABLE):
            logger.error("  the type of "+inputelem.get_long_name()+" is not TABLE")
            self.nerror+=1
            return False
        logger.debug('  num. data  :'.ljust(20)+str(inputelem.get_num_data()))
        colnameok = self.__validate_column_name(inputelem,spec)
        dataok=True
        if colnameok:
            if 'max_nrow' in spec:
                try:
                    dataok = int(spec['max_nrow'])>=inputelem.get_num_data()
                except ValueError:
                    pass
                if not dataok:
                    logger.error("too much data defined ("+str(inputelem.get_num_data())\
                                 +") in ["+inputelem.get_long_name()+"]. the maximum is "+spec['max_nrow']+".")
                    self.nerror+=1
            if dataok:
                dat = inputelem.get_str_table_data()
                for row in dat:
                    logger.debug('  validating row data no. : '+str(dat.index(row)+1))
                    if not self.__validate_row_data(row,inputelem.get_identifiers(),spec):
                        dataok=False
        return colnameok and dataok

    def __validate_column_name(self,inputelem,spec):
        all_columns_required=True
        noerr = True
        if 'all_columns_required' in spec and spec['all_columns_required']=='false':
            all_columns_required=False
        idents = inputelem.get_identifiers()
        colnames = spec[columns].split(',')
        if len(idents) != len(set(idents)):
            logger.error('duplicate entry in column names')
            for ident in idents:
                logger.error(ident)
            self.increment_nerr()
            noerr = False
        vtyps  = spec[val_type].split(',')
        for iden in idents:
            if not iden in colnames:
                logger.error('  could not find '+iden+' in column names.')
                self.nerror+=1
                close = difflib.get_close_matches(iden,colnames,1,0.8)
                if close!=None and len(close)==1:
                    logger.error('  candidate choice : '+close[0])
                else:
                    logger.error('  no candidate choice')
                noerr = False
            else:
                choi = None
                choid=choice+str(colnames.index(iden)+1)
                if choid in spec:
                    choi = spec[choid]
                desid=description+str(colnames.index(iden)+1)
                des = None
                if desid in spec:
                    des = spec[desid]
                infostr ='  found '+iden.ljust(11)+' in column names; val_type : ['+vtyps[colnames.index(iden)]+']'
                if choi != None:
                    infostr += '; choices: '+choi
                if des != None and len(des)!=0 :
                    infostr += '; description: '+des
                logger.debug(infostr)
        if all_columns_required:
            for col in colnames:
                if not col in idents:
                    logger.error('  undefined attribute : '+col)
                    noerr = False
        return noerr

    def __validate_row_data(self,rowdata,colnames,spec):
        scolname = spec[columns].split(',')
        vtyps = spec[val_type].split(',')
        ret = True
        for name in colnames:
            index = scolname.index(name)
            if not name in rowdata:
                ret = False
                continue
            if not self.__check_val_type([vtyps[index]],str(rowdata[name]),name,True):
                ret = False
            choi = None
            choid=choice+str(scolname.index(name)+1)
            if choid in spec:
                choi = spec[choid].split(',')
                if str(rowdata[name].strip())=='*':
                    continue
                if str(rowdata[name]).lower().strip() in choi:
                    logger.debug( '    found ['+str(rowdata[name])+'] in choice.')
                else:
                    logger.error('   could not find ['+str(rowdata[name])+ '] in ['+name+']')
                    self.nerror+=1
                    close = difflib.get_close_matches(str(rowdata[name]).lower().strip(),choi,1,0.8)
                    if close!=None and len(close)==1:
                        logger.error('   candidate choice : '+close[0])
                    else:
                        logger.error('   no candidate choice')
                    ret = False
            raid = val_range+str(scolname.index(name)+1)
            if raid in spec:
                if str(rowdata[name].strip())=='*':
                    continue
                if not self.__check_val_range(spec[raid],rowdata[name]):
                    ret = False
        return ret

    def __validate_vector(self,inputelem,spec):
        nam = inputelem.get_long_name()
        if inputelem.get_type()!=VECTOR:
            logger.error("  the type of "+nam+" is not VECTOR")
            self.nerror+=1
            return False
        logger.debug('  num_elements:'.ljust(20)+spec[num_elements])
        logger.debug('  val_type    :'.ljust(20)+spec[val_type])
        numelem = int(spec[num_elements])
        val = inputelem.get_value()
        cknumelem=False
        try:
            if val ==None or len(val) != numelem:
                logger.error('  the number of elements is not '+spec[num_elements])
                self.nerror+=1
            else:
                logger.debug( '  the number of elements is OK.')
                cknumelem = True
        except TypeError:
            logger.error(str(val)+" is not a sized object")
            return False

        valtyp = spec[val_type].split(',')
        vtyp = True
        for i in range(numelem):
            if i <= len(val)-1:
                b = self.__check_val_type([valtyp[i]],val[i])
                if not b:
                    vtyp=False

        return cknumelem and vtyp

    def __check_val_range(self,ra,val):
        rang = ra.strip().split(',')
        if len(rang)>=2 and len(rang[1])!=0:
            try:
                mi = float(rang[0].strip())
                ma = float(rang[1].strip())
                if float(val) >= mi and float(val) <= ma:
                    logger.debug( '    value is within range : ['+str(mi)+','+str(ma)+']')
                    return True
                else:
                    logger.error('    value is not within range : ['+str(mi)+','+str(ma)+']')
                    self.nerror+=1
            except ValueError:
                logger.error('invalid spec')
                self.nerror+=1
                return False
        elif ra.strip().startswith(','):
            try:
                ma = float(rang[0].strip())
                if float(val)<=ma:
                    logger.debug( '    value is smaller than : ['+str(ma)+']')
                    return True
                else:
                    logger.error('    value is larger than : ['+str(ma)+']')
                    self.nerror+=1
                    return False
            except ValueError:
                logger.error('invalid spec')
                self.nerror+=1
                return False
        elif ra.strip().endswith(','):
            try:
                mi = float(rang[0].strip())
                if float(val)>=mi:
                    logger.debug( '    value is larger than : ['+str(mi)+']')
                    return True
                else:
                    logger.error('    value is smaller than : ['+str(mi)+']')
                    self.nerror+=1
                    return False
            except ValueError:
                logger.error('invalid spec')
                self.nerror+=1
                return False
        return False

    fractre = re.compile('\d+$|\d+/\d+$')
    def __check_val_type(self,vtyps,val,tag='',tabular=False):
        vtyp_passed=False
        if tabular and val=='*':
            if len(tag) != 0:
                stri = '['+tag+' = '+val+']'
                logger.debug('  val_type of '+stri.ljust(10)+' is OK (default value).')
            else:
                logger.debug('  val_type is OK (default value).')
            return True
        for vtyp in vtyps:
            vtyp = vtyp.strip()
            if vtyp==val_type_fract:
                msg = 'a [fraction].'
                if self.fractre.match(val):
                    vtyp_passed=True
                
            if vtyp==val_type_int:
                msg = 'an [int].'
                try:
                    int(val)
                    vtyp_passed = True
                except ValueError:
                    pass

            if vtyp==val_type_float:
                msg = 'a [float].'
                try:
                    float(str(val).lower().replace('d','e'))
                    vtyp_passed = True
                except ValueError:
                    pass
            if vtyp == val_type_bool:
                msg = 'a [bool].'
                if not pyutil.is_bool(val):
                    pass
                else:
                    vtyp_passed = True
            if vtyp==val_type_string:
                vtyp_passed=True
        if not vtyp_passed:
            logger.error('  the value ['+val+'] is not '+msg)
            self.nerror+=1
        else:
            if len(tag) != 0:
                stri = '['+tag+' = '+val+']'
                logger.debug('  val_type of '+stri.ljust(10)+' is OK.')
            else:
                logger.debug('  val_type is OK.')
        return vtyp_passed

    def __dump(self,li=None):
        if self.root_block==None:
            return None
        entries = self.root_block.get_entries()
        inpstr=''
        for ent in entries:
            if ent.is_block():
                if li!=None:
                    li.append(ent)
                tmpstr=self.__dump_block(ent,li)
                if tmpstr!=None:
                    inpstr += tmpstr
        if li==None:
            logger.debug('input :')
            logger.debug('\n'+inpstr)
        return inpstr

    def get_nerror(self):
        return self.nerror

    def get_nwarn(self):
        return self.nwarn

    def __str__(self):
        return self.__dump()

    def __dump_block(self,bloc,li=None):
        inpstr=''
        if li==None:
            logger.debug('long name: '+bloc.get_long_name())
        ents = bloc.get_entries()
        if len(ents)==0:
            return
        depth = bloc.get_block_depth()
        ws = pyutil.get_ws(depth)
        inpstr += ws+bloc.get_name()+'{\n'
        units = bloc.get_default_units()
        if len(units)!=0:
            inpstr += ws+'  #units'
            for k,v in units.items():
                inpstr += ' '+v
            inpstr += '\n'
        for entry in ents:
            if li!=None:
                li.append(entry)
            if entry.is_block():
                tmpstr=self.__dump_block(entry,li)
                if tmpstr!=None:
                    inpstr += tmpstr
            else:
                if li==None:
                    logger.debug('long name: '+entry.get_long_name())
                entry.set_ws(ws)
                inpstr += str(entry)+'\n'
        inpstr+= ws+'}\n'
        return inpstr

    def __get_preprocessed_input(self,inp=None):
        inpfile = None
        if inp is None and self.file_name is None and not \
            os.path.exists(self.filename_data.get_file_name('F_INP')):
            logger.error('F_INP does not exist.')
        try:
            if self.file_name is None and self.filename_data is not None:
                try:
                    inpfile = open(self.filename_data.get_file_name('F_INP'),encoding='utf_8_sig')
                except TypeError:
                    inpfile = open(self.filename_data.get_file_name('F_INP'))
            elif inp is None:
                try:
                    inpfile = open(self.file_name,encoding='utf_8_sig')
                except TypeError:
                    inpfile = open(self.file_name)
            elif inp is not None:
                try:
                    inpfile = open(inp,encoding='utf_8_sig')
                except TypeError:
                    inpfile = open(inp)
        except IOError:
            #logger.error(' input open error.')
            return None
        inpbuf=[]
#        for line in inpfile:
        icount=0
        for ll in inpfile:
            icount = icount+1
            try:
                l = ll.decode('utf8')
            except AttributeError:
                l = ll
            if self.contains_zenkaku(l):
                logger.error('detected ZENKAKU character at line '+str(icount))
                self.increment_nerr()
            line = Line(ll,icount)
            line.set_line(line.get_line().strip())
            line.set_line(line.get_line().replace('!#','#'))
            if line.get_line().startswith('!') or line.get_line().startswith('//'):
                continue
            lines = line.get_line().split('!')
            lines = lines[0].split('//')
            if(not lines[0].startswith("#")):
                lines = lines[0].split(',')
            for l in lines:
                l=l.strip()
                if len(l)==0:
                    continue
                inpbuf.append(Line(l,line.get_line_no()))
        inpbuf2=[]
        for line in inpbuf:
            line.set_line(line.get_line().strip())
            lines = line.get_line().split('{')
            if len(lines)==1:
                inpbuf2.append(Line(lines[0],line.get_line_no()))
            else:
                endswithlb = line.get_line().endswith('{')
                for ll in lines:
                    ll=ll.strip()
                    if len(ll)==0:
                        continue
                    inpbuf2.append(Line(ll+'{',line.get_line_no()))
                if not endswithlb:
                    inpbuf2[len(inpbuf2)-1].set_line(inpbuf2[len(inpbuf2)-1].get_line()[0:-1])

        inpbuf = []
        for line in inpbuf2:
            line.set_line(line.get_line().strip())
            lines = line.get_line().split('}')
            if len(lines)==1:
                inpbuf.append(Line(lines[0],line.get_line_no()))
            else:
                for i in range(len(lines)-1):
                    ll = lines[i]
                    ll=ll.strip()
                    if len(ll)==0:
                        inpbuf.append(Line('}',line.get_line_no()))
                    else:
                        inpbuf.append(Line(ll,line.get_line_no()))
                        inpbuf.append(Line('}',line.get_line_no()))

        logger.debug('preprocessed input : ')
        for line in inpbuf:
            logger.debug(line)

        inpfile.close()
        return inpbuf

    def contains_zenkaku(self,l):
        for c in l:
            if unicodedata.east_asian_width(c) in 'WFA':
                return True
        return False    

class InputEntry:
    ''' the base class for input entries '''
    def __init__(self,name,line_no=-1):
        ''' parameters : name the name for this entry '''
        self.name = name
        self.type = -1
        self.ws = ''
        self.parent=None
        self.spec=None
        self.line_no=line_no

    def set_line_no(self,line_no):
        ''' set the line no associated to this input entry.'''
        self.line_no=line_no

    def get_line_no(self):
        ''' get the line no associated to this input entry.'''
        return self.line_no

    def get_name(self):
        ''' returns the name of this entry '''
        return self.name

    def set_specification(self,spe):
        ''' set the dictionary object which specifies
            the present input entry. '''
        self.spec = spe

    def get_specification(self):
        ''' returns the specification dictionary.
            will return None if no specification is set. '''
        return self.spec

    def get_type(self):
        ''' returns the 'type' of this entry;
            it will be one of BLOCK, TABLE, PRIMITIVE or VECTOR
            subclass MUST implement this method. '''
        raise NotImplementedError()

    def get_parent(self):
        ''' get the reference to the parent block '''
        return self.parent

    def _set_parent(self,parent):
        self.parent = parent

    def get_long_name(self):
        ''' get the 'long name' of this input entry.
            examples of a 'long name' :
            control.cpumax
            accuracy.ksampling
            accuracy.ksampling.method
            structure.atom_list.atoms.table
            ...
            ... '''
        ret=[]
        ent=self
        while 1:
            if ent==None:
                break
            ret.append(ent.get_name())
            ent = ent.get_parent()
        ret.reverse()
        retret=''
        for re in ret:
            retret += re+'.'
        retret = retret[0:-1].replace('root.','').lower()
        return retret

    def is_block(self):
        ''' returns True if this input entry is of type BLOCK '''
        return self.get_type()==BLOCK

    def set_ws(self,ws):
        ''' sets the amount of white space for printing
            parameter : ws the white-space string '''
        self.ws = ws

    def delete_me(self):
        ''' remove this entry from its parent block
            (this method will not work for the root block,
            whose parent is None). '''
        if self.parent != None:
            self.parent.remove_entry(self)

class BlockEntry(InputEntry):
    ''' a class which represents a BLOCK entry '''
    def __init__(self,name,inp=None,line_no=-1):
        InputEntry.__init__(self,name,line_no)
        self.default_unit = {}
        self.entries = []
        self.input=None
        self.__closed = False
        if inp != None:
            self.input = inp

    def close_block(self):
        self.__closed = True 

    def closed(self):
        return self.__closed

    def get_input_obj(self):
        if self.input!=None:
            return self.input
        ent=self
        while 1:
            if ent==None:
                return None
            if ent.input!=None:
                return ent.input
            ent = ent.get_parent()

    def get_type(self):
        return BLOCK

    def set_default_units(self,du):
        ''' set the 'default unit' associated with this block.
            parameter: du a Python dictionary;
            example: du[unit_type_force] >> 'hartree/bohr' '''
        self.default_unit = du

    def set_default_unit(self,uval,typ=None):
        ''' set the 'default_unit' the unit type will be
            automatically resolved unless you give it
            explicitly.
            parameter : uval the value of the unit
            typ : the unit type (optional) '''
        if typ==None:
            typ = resolve_unit(uval)
            if typ != None:
                self.default_unit[typ]=str(uval)
        else:
            self.default_unit[typ] = str(uval)

    def get_default_unit(self,utyp):
        ''' get the 'default unit'
            will return None if the specified unit-type
            does not exist.
            parameter : utyp the type of the unit to be obtained. '''
        if not utyp in self.default_unit:
            return None
        return self.default_unit[utyp]

    def get_default_units(self):
        ''' get the 'default unit' dictionary itself. '''
        return self.default_unit

    def add_entry(self,entry):
        ''' add a new entry to this block. if an entry of the same name
            already exist, it will be replaced.
            parameter : entry the entry to be added.
            it must be an instance of a derived class of InputEntry. '''
        if not isinstance(entry,InputEntry):
            logger.error(type(entry))
            raise TypeError('only an instance of InputEntry can be added.')

        ents = self.get_entries()
        nam=entry.get_name().lower()
        exi=False
        for ent in ents:
            if ent.get_name().lower()==nam:
                self.entries[self.entries.index(ent)]=entry
                exi=True
                break
        if not exi:
            self.entries.append(entry)
        entry._set_parent(self)
        lname = entry.get_long_name()
        inp=self.get_input_obj()
        sp=None
        if inp != None:
            sp = inp.resolve_spec(lname)
        else:
            sp=NOT_READY

        if sp==None and sp!=NOT_READY:
#            logger.error('unknown entry : '+lname)
            pass
        elif sp==NOT_READY:
            return
        else:
            entry.set_specification(sp)

    def remove_entry(self,entry):
        ''' remove the specified entry (if it exists)
            parameter : entry the entry to be removed. '''
        if entry in self.entries:
            self.entries.remove(entry)

    def get_entries(self):
        ''' get all the entries contained in this block. '''
        return self.entries

    def get_entry(self,entry_name):
        ''' get the entry
            parameter : entry_name the name of the entry to be obtained;
            (case insensitive). if entry_name is delimited by a '.', then
            get_entry will recursively call itself assuming that the list obtained by
            entry_name.split('.')[:1] represent a block of the corresponding depth.
            will return None if no such entry exists. '''
        entnam = entry_name.split(".")
        if len(entnam)<=1:
            for ent in self.entries:
                if ent.get_name().lower()==entry_name.lower():
                    return ent
        else:
            ent = self.get_entry(entnam[0])
            if ent == None:
                return None
            str = ".".join(entnam[1:])
            return ent.get_entry(str)
        return None

    def entry_count(self):
        ''' get the number of entries defined in this block '''
        return len(self.entries)

    def get_block_depth(self):
        ''' get the 'depth' of the block '''
        parent = self
        depth = 0
        while 1:
            parent = parent.get_parent()
            if parent == None:
                return depth
            depth += 1
    def __str__(self):
        return self.get_name()

class PrimitiveEntry(InputEntry):
    ''' this class represents a 'primitive entry'.
        a primitive entry is specified by
        name = value unit
        a value can be an int, float, bool or string;
        only a float can have units. '''

    def __init__(self,name,val='',uni='',dq=False,line_no=-1):
        ''' parameter : name the name of this entry '''
        InputEntry.__init__(self,name,line_no)
        self.unit = str(uni)
        self.value = str(val)
        self.dq = dq

    def get_type(self):
        return PRIMITIVE

    def set_unit(self,unit):
        ''' set the unit for this entry
            paremeter : unit the unit; it must be a string '''
        self.unit = str(unit)

    def set_value(self,value):
        ''' set the value for this entry
            parameter : value the value; it must be a string. '''
        tmpstr=str(value)
        if str(value)=='False':
            tmpstr='off'
        elif str(value)=='True':
            tmpstr='on'
        self.value = tmpstr

    def get_unit(self):
        ''' get the unit '''
        return self.unit

    def get_value(self):
        ''' get the value. if the specification is given,
            will return a suitably converted value. '''
        return get_value(self.spec,self.value)

    def get_str_value(self):
        ''' get the value in string. '''
        return str(self.value)

    def set_double_quote(self,dq):
        ''' if set to True, will double quote the value. '''
        self.dq = dq

    def is_double_quoted(self):
        ''' returns True if the value should be double-quoted on save. '''
        return self.dq

    def __str__(self):
        if not self.dq:
            return self.ws+'  '+self.name+' = '+str(self.value)+' '+self.unit
        else:
            return self.ws+'  '+self.name+' = \"'+str(self.value)+'\"'

class VectorEntry(InputEntry):
    ''' this class represents a 'vector entry'. note that
        the only usage of this data structure is for the unit cell. '''
    def __init__(self,name,line_no=-1):
        InputEntry.__init__(self,name,line_no)
        self.value=[]

    def set_value(self,value):
        ''' set the value.
            parameter : value the value of
            the vector entry. it must be a Python list. '''
        self.value = value

    def add_value(self,va):
        ''' add values to the vector
            parameter : va the value to be added '''
        self.value.append(str(va))

    def get_value(self):
        ''' get the value.
            returns a python list. '''
        if self.spec==None:
            return self.value
        ret=[]
        inde=0
        for val in self.value:
            ret.append(get_value(self.spec,val,inde))
            inde=inde+1
        return ret

    def get_type(self):
        return VECTOR

    def __str__(self):
        ret = self.name+' ='
        for val in self.value:
            ret += ' '+str(val)
        return self.ws+'  '+ret

class TableEntry(InputEntry):
    ''' this class represents a tabular data.
        refer to the PHASE manual for details of
        this data atructure'''
    def __init__(self,line_no=-1):
        InputEntry.__init__(self,'table',line_no)
        self.identifiers = []
        self.data = []
        self.default = {}

    def reset_data(self):
        self.data = []

    def get_num_data(self):
        ''' returns the number of row data included in this table '''
        return len(self.data)

    def get_type(self):
        return TABLE

    def add_default(self,iden,val):
        ''' a table can have 'default values'.
            the default value will be stored in a Python dictionary
            parameter : iden the identifier for the default value.
            parameter : val  the actual value of the default value. '''
        logger.debug(iden)
        logger.debug(val)
        self.default[iden] = str(val)

    def get_defaults(self):
        return self.default

    def add_identifier(self,ident):
        ''' add table identifiers
            (note that the table identifiers are stored in lists). '''
        self.identifiers.append(ident)

    def set_identifiers(self,ident):
        self.identifiers = ident

    def get_identifiers(self):
        ''' returns a list of column identifiers. '''
        return self.identifiers

    def valid_identifier(self,iden):
        ''' return True if iden is included in the list of identifiers, False otherwise '''
        return iden in self.identifiers

    def get_table_data(self):
        ''' returns the table data it self.
            table data are stored in a list of Python dictionaries.
            for instance, get_table_data()[0]['element'] should return
            the element of the 0-th atom. '''
        if self.spec==None:
            return self.data
        ret=[]
        for k in range(len(self.data)):
            self.data[k] = self.get_row_data(k)
        return self.data

    def get_str_table_data(self):
        return self.data

    def get_corresponding_data(self,key,val,ident):
        ''' similar to get_data, except that the row index is resolved differently.
            the row whose 'key' column equals 'val' will be adopted.
            will return None if 'key' is unfound in the identifiers,
            or if no  equals 'val' '''
        if not key in self.identifiers:
            return None
        for i in range(len(self.data)):
            rowdata = self.data[i]
            if not key in list(rowdata.keys()):
                return None
            coldat = rowdata[key]
            if coldat == val:
                return get_data(i,ident)
        return None
        
    def get_data(self,row,ident):
        ''' returns the specified data
            parameter : row the row index of the data.
            ident : the identifier for the data.
            if the specification file is given, a suitably converted
            value will be returned. '''

        if row>=len(self.data):
            return None
        if not ident in self.data[row]:
            return None

        if self.spec==None:
            return self.data[row][ident]
        cols=self.spec[columns].split(',')
        if not ident.lower() in cols:
            return self.data[row][ident]
        inde = cols.index(ident)
        return get_value(self.spec,self.data[row][ident],inde)

    def get_row_data(self,row):
        ''' get the specified row-data
            parameter : row the row-index of the data
            if the specification file is given, a suitably converted
            value will be returned.'''
        if row>=len(self.data):
            return None

        if self.spec==None:
            return self.data[row]

        keys=list(self.data[row].keys())
        for key in keys:
            self.data[row][key]=self.get_data(row,key)
        return self.data[row]

    def get_str_row_data(self,row):
        ''' get the row data as it is.'''
        if row>=len(self.data):
            return None
        return self.data[row]

    def add_row_data(self,row):
        ''' add row data.
            remark: dat must be a python dictionary '''
        rkeys = list(row.keys())
        for rkey in rkeys:
            row[rkey] = str(row[rkey])

        for rkey in rkeys:
            if not rkey in self.identifiers:
                self.identifiers.append(rkey)
        self.data.append(row)

    def remove_row_data(self,rowdata):
        ''' remove the specified row. will do nothing if the
            specified row-data does not exist.
            parameter : the row data to be removed. '''
        if rowdata in self.data:
            self.data.remove(rowdata)

    def remove_row_data_by_index(self,index):
        if index>=len(self.data):
            return
        self.data.pop(index)

    def set_data(self,row,ident,value):
        if row>=len(self.data):
            logger.warn('at set_data : row index is larger than data.')
            logger.warn('nothing to do')
            return
        self.data[row][ident] = str(value)

    def set_row_data(self,row,dat):
        ''' remark: dat must be a python dictionary '''
        if row>=len(self.data):
            logger.warn('at set_data : row index is larger than data.')
            logger.warn('nothing to do')
            return
        for da in dat:
            dat[da] = str(dat[da])
        self.data[row] = dat

    def __str__(self):
        ret = self.ws+'  #tag'
        for iden in self.identifiers:
            ret += '    '+iden
        for row in self.data:
            ret += '\n'+self.ws+'    '
            for iden in self.identifiers:
                if iden in row:
                    tmpstr=str(row[iden])
                    if str(row[iden])=='False':
                        tmpstr='off'
                    elif str(row[iden])=='True':
                        tmpstr='on'
                    ret += '    '+tmpstr
                else:
                    ret += '    *'
        return ret

class FileNamesData:
    fnames_keys = [
     'F_INP'
    ,'F_POT(1)'
    ,'F_POT(2)'
    ,'F_POT(3)'
    ,'F_POT(4)'
    ,'F_POT(5)'
    ,'F_POT(6)'
    ,'F_POT(7)'
    ,'F_POT(8)'
    ,'F_POT(9)'
    ,'F_POT(10)'
    ,'F_POT(11)'
    ,'F_POT(12)'
    ,'F_POT(13)'
    ,'F_POT(14)'
    ,'F_POT(15)'
    ,'F_POT(16)'
    ,'F_CHGT'
    ,'F_VLC'
    ,'F_CHR'
    ,'F_WANNIER'
    ]

    neb_file_idents={
      'F_NEB_OUT':'output_neb'
     ,'F_IMAGE(0)':'endpoint0.data'
     ,'F_IMAGE(1)':''
     ,'F_IMAGE(2)':''
     ,'F_IMAGE(3)':''
     ,'F_IMAGE(4)':''
     ,'F_IMAGE(5)':''
     ,'F_IMAGE(6)':''
     ,'F_IMAGE(7)':''
     ,'F_IMAGE(8)':''
     ,'F_IMAGE(9)':''
     ,'F_IMAGE(10)':''
     ,'F_IMAGE(11)':''
     ,'F_IMAGE(12)':''
     ,'F_IMAGE(13)':''
     ,'F_IMAGE(14)':''
     ,'F_IMAGE(15)':''
     ,'F_IMAGE(16)':''
     ,'F_IMAGE(17)':''
     ,'F_IMAGE(18)':''
     ,'F_IMAGE(19)':''
     ,'F_IMAGE(20)':''
     ,'F_IMAGE(21)':''
     ,'F_IMAGE(22)':''
     ,'F_IMAGE(23)':''
     ,'F_IMAGE(24)':''
     ,'F_IMAGE(25)':''
     ,'F_IMAGE(-1)':'endpoint1.data'
     ,'F_NEB_STOP':'./nfnebstop.data'
     ,'F_NEB_CNTN':'./neb_continue.data'
     ,'F_NEB_ENF':'./nfnebenf.data'
     ,'F_NEB_DYNM':'./nfnebdynm.data'
    }

    fname_idents={
     'F_INP':'./nfinp.data'
    ,'F_POT(1)':'pot.01'
    ,'F_POT(2)':'pot.02'
    ,'F_POT(3)':'pot.03'
    ,'F_POT(4)':'pot.04'
    ,'F_POT(5)':'pot.05'
    ,'F_POT(6)':'pot.06'
    ,'F_POT(7)':'pot.07'
    ,'F_POT(8)':'pot.08'
    ,'F_POT(9)':'pot.09'
    ,'F_POT(10)':'pot.10'
    ,'F_POT(11)':'pot.11'
    ,'F_POT(12)':'pot.12'
    ,'F_POT(13)':'pot.13'
    ,'F_POT(14)':'pot.14'
    ,'F_POT(15)':'pot.15'
    ,'F_POT(16)':'pot.16'
    ,'F_PKB':'./vkb.data'
    ,'F_PD':'./vd.data'
    ,'F_PPC':'./vpc.data'
    ,'F_STOP':'./nfstop.data'
    ,'F_CNST':'./nfcnst.data'
    ,'F_OPGR':'./opgr.data'
    ,'F_MATBP':'./matrix.BP'
    ,'F_KPOINT':'./kpoint.data'
    ,'F_KINDEX':'./f.kp0'
    ,'F_SPG':'./bnprp4.i5'
    ,'F_KPGN':'./bnkpgn.i5'
    ,'F_OTP':'./nfotp.data'
    ,'F_DYNM':'./nfdynm.data'
    ,'F_EGRID':'./nfegrid.data'
    ,'F_LDOS':'./nfldos.data'
    ,'F_ENF':'./nfefn.data'
    ,'F_GPT':'./nfgpt.data'
    ,'F_CHGT':'./nfchgt.data'
    ,'F_CHGO':'./nfchgo.data'
    ,'F_CHGU':'./nfchgu.data'
    ,'F_ENERG':'./nfenergy.data'
    ,'F_CNTN':'./continue.data'
    ,'F_CNTN_BIN':'./continue_bin.data'
    ,'F_ZAJ':'./zaj.data'
    ,'F_VLC':'./nfvlc.cube'
    ,'F_CNTN_BIN_STM':'./continue_bin_stm.data'
    ,'F_DOS':'./dos.data'
    ,'F_CHR':'./nfchr.cube'
    ,'F_WFk':'./nfwfk.data'
    ,'F_WANNIER':'./nfwannier.cube'
    ,'F_CNTN_WAN':'./nfcontinue_wannier.data'
    ,'F_POT_WAN':'./nfpotential_wannier.data'
    ,'F_ELF':'./nfelf.data'
    ,'F_BERRY':'./berry.data'
    ,'F_EFFCHG':'./effchg.data'
    ,'F_FORCE':'./force.data'
    ,'F_MODE':'./mode.data'
    ,'F_EPSILON':'./epsilon.data'
    ,'F_EPSOUT':'./eps.data'
    ,'F_NLO':'./nlo.data'
    ,'F_MAGOPT':'./magopt.data'
    ,'F_PHDOS':'./phdos.data'
    ,'F_PHDOS':'./phdos.data'
    ,'F_OCCMAT':'./occmat.data'
    ,'F_STRFRC':'./strfrc.data'
    ,'F_PSTRN':'./positron.cube'
    ,'F_VELEC':'./electron.cube'
    ,'F_EPPAIR':'./ep_pair.cube'
    ,'F_VELEC_GRAD':'./electron_grad.cube'
    ,'F_EFERMI':'./nfefermi.data'

    ,'F_ZAJ_filetype':'unified'
    ,'F_CHGT_filetype':'unified'
    ,'F_CNTN_BIN_filetype':'unified'
    ,'F_CNTN_filetype':'unified'

    ,'F_ZAJ_in_filetype':'nogiven'
    ,'F_CHGT_in_filetype':'nogiven'
    ,'F_CNTN_BIN_in_filetype':'nogiven'
    ,'F_CNTN_in_filetype':'nogiven'

    ,'F_CNTN_BIN_PAW':'./continue_bin_paw.data'
    ,'F_CNTN_BIN_RSNLPP':'./continue_bin_rsnlpp.data'

    ,'F_HYBRIDINFO':'./nfhybridinfo.data'
    ,'F_SCF_ZAJ':'./zaj.data'
    
    ,'F_POS':'./nfinp.data'
    ,'F_COORD_ATTR':'./nfinp.data'

    ,'F_CHKPNT':'chkpnt'
    ,'F_LATCONST':'./nflatconst.data'
    ,'F_DFTD3PAR': './dftd3par.data'
    }

    mpifiletypes_idents = {
    'F_ZAJ_filetype' : 'unified'
    ,'F_CHGT_filetype' : 'unified'
    ,'F_CNTN_BIN_filetype' : 'unified'
    ,'F_CNTN_filetype' : 'unified'

    ,'F_ZAJ_in'      : './zaj.data'
    ,'F_CNTN_in'     : './continue.data'
    ,'F_CNTN_BIN_in' : './continue_bin.data'
    ,'F_CHGT_in'     : './nfchgt.data'

    ,'F_ZAJ_bak'      : './zaj.data'
    ,'F_CNTN_bak'     : './continue.data'
    ,'F_CNTN_BIN_bak' : './continue_bin.data'
    ,'F_CHGT_bak'     : './nfchgt.data'

    ,'F_ZAJ_in_filetype' : 'nogiven'
    ,'F_CHGT_in_filetype' : 'nogiven'
    ,'F_CNTN_BIN_in_filetype' : 'nogiven'
    ,'F_CNTN_in_filetype' : 'nogiven'
    }

    f_param_name_idents = {
    'F_PRM':'./m_ArraySize_Parameters.f90'
    }

    default_identifiers ={
      '&fnames':fname_idents
     ,'&nebfiles':neb_file_idents
     ,'&mpifiletypes':mpifiletypes_idents
     ,'&f_param_name':f_param_name_idents
    }

    def __init__(self,inp,nfinp=None):
        self.fnames_data_exists=False
        self.inp = inp
        self.file_name = join('.','file_names.data')
        if not os.path.exists(self.file_name):
            self.fnames_data_exists=False
        else:
            self.fnames_data_exists=True
            #return
            #if nfinp is None:
            #    logger.debug('file_names.data does not exist.')
            #    self.valid_fnames=False
            #    return
            #else:
            #    f=open(self.file_name,'w')
            #    f.write('&fnames\n')
            #    f.write("F_INP='"+nfinp+"'\n")
            #    f.write("/")
            #    f.flush()
            #    f.close()
        self.entry = self.default_identifiers
        #self.entry["&fnames"]=self.default_identifiers
        self.valid_fnames = False
        if self.fnames_data_exists:
            self.valid_fnames = self._ini()

    def is_valid(self):
        return self.valid_fnames

    def get_file_name(self,ident,file_set="&fnames",to_upper=True):
        ''' get the file name corresponding to the specified identifier.
            parameter : ident the identifier '''
        ide = str(ident)
        try:
            if to_upper:
                ide = str(ident).upper()
            return self.entry[file_set][ide]
        except:
            return None

    def set_file_name(self,ident,file_name,file_set="&fnames"):
        ''' set the file name corresponding to the specified identifier.
            parameter : ident the identifier.
            parameter : file_name the filename for the file corresponding to ident. '''
        f = self.entry[file_set]
        f[ident] = file_name
        try:
            logger.debug(self.fnames_keys.index(ident))
        except ValueError:
            logger.debug("added "+ident+" to the list of keys")
            self.fnames_keys.append(ident)

    def save(self,logoffset=''):
        ''' save the file_names.data file to disk.
            if direc is given, the method will save to the file_names.data
            under dir. '''
        fname = join(self.inp.get_curr_dir(),'file_names.data')

        logger.info(logoffset+"dumping file_names.data to : "+fname)

        try:
            f=open(fname,"w")
            string = self.__str__()
            f.write(string)
            f.flush()
            f.close()
            ekeys = list(self.entry.keys())
        except IOError:
            logger.error(logoffset+"failed write to "+self.file_name)

    def __str__(self):
        ret=''
        ekeys = list(self.entry.keys())
        for ekey in ekeys:
            ret += ekey+'\n'
            iden = self.entry[ekey]
            if ekey=='&fnames':
                for key in self.fnames_keys:
                    if key in iden and len(iden[key])!=0:
                        ret += key+' = \''+iden[key]+'\'\n'
            else:
                ikeys = list(iden.keys())
                for ikey in ikeys:
                    if ikey in iden:
                        ret += ikey+' = \''+iden[ikey]+'\'\n'
            ret += '/\n'
        return ret

    def file_names_data_exists(self):
        return self.fnames_data_exists

    def _ini(self):
        logger.debug("")
        logger.debug("checking the file_names.data file...")
        boo = True
        try:
            fi = open(self.file_name,"r",encoding='utf_8_sig')
        except TypeError:
            fi = open(self.file_name,"r")
        currident = ''
        for line in fi:
            line = line.strip()
            if line.startswith("&"):
                if line in list(self.default_identifiers.keys()):
                    self.entry[line] = self.default_identifiers[line]
                else:
                    self.entry[line] = {}
                currident = line
            elif line.startswith("/"):
                currident = ''
            elif line.startswith("!"):
                continue
            else:
                if not len(currident)==0:
                    bb = line.split("=")
                    ident = ''
                    val = ''
                    if len(bb)>1:
                        ident = bb[0].strip()
                        if len(bb[1])<1:
                            continue
                        if ident in list(self.entry[currident].keys()):
                            val = bb[1].strip()[1:len(val)-1]
                            self.entry[currident][ident] = val
                        else:
                            logger.error("  unknown file pointer : "+ident.rjust(15)+" in "+currident)
                            boo=False
                            self.inp.increment_nerr()
        if boo:
            logger.debug("no errors were found in the file_names.data file")
        else:
            logger.error("encountered errors in the file_names.data file")

        return boo

    def debug_output(self):
        logger.debug(self.file_name)
        ekeys = list(self.entry.keys())
        for ekey in ekeys:
            logger.debug('entry : '+ekey)
            iden = self.entry[ekey]
            if ekey=='&fnames':
                for key in self.fnames_keys:
                    if key in iden:
                        logger.debug('  '+key.ljust(10)+' : '+iden[key])
            else:
                ikeys = list(iden.keys())
                for ikey in ikeys:
                    if ikey in iden:
                        logger.debug('  '+ikey.ljust(14)+' : '+iden[ikey])

    def convert_pp_to_abspath(self):
        ''' convert the path to the pseudo-potential file to absolute'''
        if not '&fnames' in self.entry:
            logger.warn(' no &fname entry in file_names.data')
            return
        fnams = self.entry['&fnames']
        idents = list(fnams.keys())
        for ident in idents:
            if ident.startswith('F_POT('):
                fnams[ident] = os.path.abspath(fnams[ident])

    def convert_to_abspath(self,ident):
        ''' convert the path to file specified by ident. '''
        if not '&fnames' in self.entry:
            logger.warn(' no &fname entry in file_names.data')
            return
        fnams = self.entry['&fnames']
        if not ident in fnams:
            logger.warn('could not find '+ident+' in file_names.data')
            return
        fnams[ident] = os.path.abspath(fnams[ident])
        
class InputValidationError(Exception):
    ''' exception raised on input validation error. '''
    def __init__(self, value):
         self.value = value

    def __str__(self):
         return repr(self.value)

class F_DYNM(config.AtomConfigGenerator):
    ''' representation of the F_DYNM file, which is the 
        standard corodinate file for PHASE. '''
    #valid = False
    IN_HEADER=0
    IN_DATA=2

    def __init__(self,coord_file='nfdynm.data'):
        super().__init__()
        self.dynm = coord_file
        #if os.path.exists(self.dynm):
        #    self.valid = True
        self.valid=True
        config.AtomConfigGenerator.__init__(self) 
    
    def get_name(self):
        return 'PHASE output format (F_DYNM format)'

    def get_defaultname(self,args):
        return 'nfdynm.data'

    def _gen_atomic_configuration(self):
        dynmfile = open(self.dynm)
        ret = []
        unitcell = None
        ntyp=0
        natm=0
        avec=[]
        bvec=[]
        cvec=[]
        lineno=0
        atm_typ=[]
        species_map = {}
        symbol = []
        mode = self.IN_HEADER
        atom_config=None
        atm_count=0
        logger.info("parsing F_DYNM file : "+self.dynm)
        for line in dynmfile:
            line = line.strip()
            lineno+=1
            if mode==self.IN_HEADER and line.startswith('cps'):
                for at in atm_typ:
                    symbol.append(species_map[at])
                for sym in symbol:
                    logger.debug("species : "+sym)
                unitcell = config.UnitCell(a_vector=avec,b_vector=bvec,c_vector=cvec)    
                mode = self.IN_DATA
                atm_count=0

            if line.startswith('cps'):
                atom_config = config.AtomConfig(coordinate_system=config.CARTESIAN, \
                             length_unit=config.BOHR, unit_cell=unitcell)
                ret.append(atom_config)
                atm_count=0
                continue

            if line.startswith('#'):
                mode=self.IN_HEADER
                words = line.split('#')
                if len(words)<2:
                    continue
                words = line.split('#')[1].split()
                if len(words)==0:
                    continue
                word0 = words[0]
                try:
                    if word0=='a_vector':
                        avec=[]
                        bvec=[]
                        cvec=[]
                        atm_typ=[]
                        species_map={}
                        symbol=[]
                        avec.append(float(words[2]))
                        avec.append(float(words[3]))
                        avec.append(float(words[4]))
                    elif word0=='b_vector':
                        bvec.append(float(words[2]))
                        bvec.append(float(words[3]))
                        bvec.append(float(words[4]))
                    elif word0=='c_vector':
                        cvec.append(float(words[2]))
                        cvec.append(float(words[3]))
                        cvec.append(float(words[4]))
                    elif word0=='ntyp':
                        ntyp = int(words[2])
                        natm = int(words[5])
                    elif word0=='(natm->type)':
                        for word in words:
                            if word==word0:
                                continue
                            atm_typ.append(int(word))
                    elif word0.startswith('(speciesname)'):
                        species_map[int(words[1])] = words[3]
                except ValueError:
                    logger.error("encountered error at line "+str(lineno)+" in file "+self.dynm)
                    self.valid=False

            if mode==self.IN_DATA:
                words = line.split()
                try:
                    geta = 0
                    if len(words)==10: # molecular-dynamics case
                        geta += 3
                    fx = float(words[4+geta])
                    fy = float(words[5+geta])
                    fz = float(words[6+geta])
                    #atm = config.Atom(atom_config,element=species_map[atm_typ[atm_count]],rx=float(words[1]),ry=float(words[2]),rz=float(words[3]),mobile=None,weight=None,mode=config.CARTESIAN)
                    atom_config.add_atom(element=species_map[atm_typ[atm_count]],\
                    rx=float(words[1]),ry=float(words[2]),rz=float(words[3]),fx=fx,fy=fy,fz=fz)
                    #atom_config.add_atom(atom = atm)
                    #atm.set_attribute('forcex',forx)
                    #atm.set_attribute('forcey',fory)
                    #atm.set_attribute('forcez',forz)
              
                    atm_count+=1
                except ValueError:
                    logger.error("encountered error at line "+str(lineno)+" in file "+self.dynm)
                    self.valid=False
        logger.info("done.")
        logger.info("number of frames in this F_DYNM file : "+str(len(ret)))
        return ret
    
    def _write_dynm_header(self,file,coords):
        if not coords.is_cell_ready():
            return
        unitcell = coords.get_unit_cell()
        (avec,bvec,cvec) = unitcell.get_lattice_vector(to_unit=config.BOHR)
        atomlist = coords.get_atom_list()
        
        file.write('#\n')
        file.write('#   a_vector = '+str(avec[0])+"     "+str(avec[1])+"     "+str(avec[2])+'\n')
        file.write('#   b_vector = '+str(bvec[0])+"     "+str(bvec[1])+"     "+str(bvec[2])+'\n')
        file.write('#   c_vector = '+str(cvec[0])+"     "+str(cvec[1])+"     "+str(cvec[2])+'\n')
        elements = coords.get_elem_list()
        file.write('#   ntyp = '+str(len(elements))+' natm = '+str(coords.get_num_atom())+'\n')
        prefix =   '# (natm->type) '
        row=prefix
        for i in range(len(atomlist)):
            atom=atomlist[i]
            ind = elements.index(atom.get_element_name())+1
            row += '    '+str(ind)
            if (i+1)%10==0:
                file.write(row+'\n')
                row=prefix
        if len(atomlist)%10!=0:
            file.write(row+'\n')

        for i in range(len(elements)):
            file.write('# (speciesname)     '+str(i+1)+' :    '+elements[i]+'\n')
        file.write('#'+'\n')

    def export_atomic_configuration(self,atomic_coordinates,to_file='nfdynm.data',\
        frame_no=None,all_frames=True):
        if to_file is None: 
            return # generating the f_dynm instance does not make sense if you're not outputting it to a file

        frameno = len(atomic_coordinates)-1
        if not all_frames and frame_no>=0 and frame_no<frameno:
            frameno = frame_no

        file=open(to_file,'w')
        atomic_coordinate = atomic_coordinates[0]
        self._write_dynm_header(file,atomic_coordinate)
        unitcell0 = atomic_coordinate.get_unit_cell()
        #com = atomic_coordinate.get_center_of_mass()
        for icoord in range(len(atomic_coordinates)):
            if not all_frames and icoord != frameno:
                continue
            coord = atomic_coordinates[icoord]
            #com_now = coord.get_center_of_mass()
            #coord.shift((com[0]-com_now[0],com[1]-com_now[1],com[2]-com_now[2]))
            unitcell=coord.get_unit_cell()
            if unitcell0!=unitcell:
                self._write_dynm_header(file,coord)
            unitcell0 = unitcell
            file.write('  cps cpd and forc at (iter_ion, iter_total =     *      * )\n')
            atomlist=coord.get_atom_list()
            for index in range(len(atomlist)):
                atom = atomlist[index]
                pos = atom.get_pos(mode=config.CARTESIAN,to_unit=config.BOHR)
                forcx = atom.get_attribute('forcex')
                forcy = atom.get_attribute('forcey')
                forcz = atom.get_attribute('forcez')
                if forcx is None or forcy is None or forcz is None:
                    forcx = atom.fx
                    forcy = atom.fy
                    forcz = atom.fz
                forc = '0.0 0.0 0.0'
                if forcx is not None and forcy is not None and forcz is not None:
                    forc = str(forcx[0])+' '+str(forcy[0])+' '+str(forcz[0])
                veloc = ''
                velocx = atom.get_attribute('velocx')
                velocy = atom.get_attribute('velocy')
                velocz = atom.get_attribute('velocz')
                if velocx is not None and velocy is not None and velocz is not None:
                    veloc = str(velocx[0])+' '+str(velocy[0])+' '+str(velocz[0])+' '

                file.write(str(index+1)+' '+str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+' '+veloc+forc+'\n')

        file.flush()
        file.close()

    def is_valid(self):
        return self.valid

    def supports_frame(self):
        return True

