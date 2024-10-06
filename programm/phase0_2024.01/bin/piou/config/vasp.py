''' maipulate I/O files for the VASP code '''
import logging
logger = logging.getLogger('piou.config.vasp')
import os
from piou import config
from piou.data import properties

class VASP_input(config.AtomConfigGenerator):
    inpfile = None
    elems = None

    def get_defaultname(self,args):
        return 'POSCAR'

    def __init__(self,coord_file):
        super().__init__()
        if coord_file is None:
            self.set_is_valid(False)
            return
        self.inpfile = coord_file

    def get_elem_map(self):
        return self.elems

    def _gen_atomic_configuration(self):
        f=open(self.inpfile,'r')
        line=[]
        for lin in f:
            line.append(lin.strip())
        f.close()

        if len(line)<8:
            logger.error('insufficient number of lines in POSCAR')
            return None

        pot = 'POTCAR'
        dir = os.path.dirname(self.inpfile)
        if len(dir)>0:
            pot = dir+os.sep+'POTCAR'
        if os.path.exists(pot):
            f=open(pot,'r')
            key = 'VRHFIN'
            elems=[]
            for li in f:
                li = li.strip()
                if li.startswith(key):
                    w0 = li.split('=')
                    if len(w0)>=2:
                        elems.append(w0[1].strip().split(':')[0])
            msg = 'elements resolved from the POTCAR file :'
            for ele in elems:
                msg += ' '+ele
            print(msg)

        factor=1.0
        try:
            factor = float(line[1].split()[0])
        except ValueError:
            logger.error("the 1st line is invalid : "+line[1])
            return None

        try:
            av=line[2].split()
            av[0] = float(av[0])*factor
            av[1] = float(av[1])*factor
            av[2] = float(av[2])*factor
            bv=line[3].split()
            bv[0] = float(bv[0])*factor
            bv[1] = float(bv[1])*factor
            bv[2] = float(bv[2])*factor
            cv=line[4].split()
            cv[0] = float(cv[0])*factor
            cv[1] = float(cv[1])*factor
            cv[2] = float(cv[2])*factor
        except ValueError:
            logger.error('failed to read the unit cell')
            return None
        except KeyError:
            logger.error('failed to read the unit cell')
            return None

        unitcell = config.UnitCell(\
        a_vector=[av[0],av[1],av[2]],\
        b_vector=[bv[0],bv[1],bv[2]],\
        c_vector=[cv[0],cv[1],cv[2]])

        nat_per_elem_line = line[5].split()
        ver5=False
        ioff = 0
        try:
            iii = int(nat_per_elem_line[0])
        except ValueError:
            ver5=True
        if ver5:
            logger.info("detected POSCAR for VASP ver.5 (and over)")
            elems = line[5].split()
            nat_per_elem_line  = line[6].split()
            ioff = 1
        else:
            logger.info("detected POSCAR for VASP ver.4")

        nat_per_elem = []
        for i in range(len(nat_per_elem_line)):
            try:
                ii = int(nat_per_elem_line[i])
                nat_per_elem.append(ii)
            except ValueError:
                if not ver5:
                    logger.error('invalid specification for the number of atoms : '+str(line[5+ioff]))
                    return None
                else:
                    continue
        if len(nat_per_elem)!=len(elems):
            logger.error('inconsistent number of elements defined in the POSCAR and POTCAR files')
            logger.error(str(len(nat_per_elem))+' '+str(len(elems)))
            return None

        nextind=6+ioff
        if line[6+ioff].lower().startswith('s'):
            nextind=7+ioff

        atomcount=0
        for i in range(len(nat_per_elem)):
            atomcount += nat_per_elem[i]
        if atomcount+nextind>len(line):
            logger.error('insufficient number of lines in POSCAR')
            return None

        mode=config.FRACTIONAL
        factatom=1.0
        if line[nextind].lower().startswith('c') or \
            line[nextind].lower().startswith('k'):
            mode=config.CARTESIAN
            factatom=factor
        atomconf=config.AtomConfig(coordinate_system=mode,\
        length_unit=config.ANGSTROM,unit_cell=unitcell)
        unitcell.set_length_unit(config.ANGSTROM)

        self.elems=[]
        atomcount = 0
        for ielem in range(len(nat_per_elem)):
            elem = elems[ielem]
            for i in range(nat_per_elem[ielem]):
                self.elems.append(elem)
                atomcount += 1
                index = atomcount+nextind
                pos = line[index].split()
                for i in range(3):
                    try:
                        pos[i] = float(pos[i])*factatom
                    except ValueError:
                        logger.error('failed to read the atomic coordinates from the POSCAR file')
                        return None
                    except IndexError:
                        logger.error('failed to read the atomic coordinates from the POSCAR file')
                        return None
                atomconf.add_atom(element=elem,rx=pos[0],ry=pos[1],rz=pos[2])

        return [atomconf]

    mobile_atms = None
    def set_mobile_atoms(self,matoms):
        self.mobile_atms = matoms

    def export_atomic_configuration(self,atomic_coordinates,to_file=None,frame_no=None,all_frames=False):
        frameno=config.resolve_frame_no(frame_no,len(atomic_coordinates))
        atomic_coordinate = atomic_coordinates[frameno]
        if not atomic_coordinate.is_cell_ready():
            logger.error('the unit cell is not set')
            return

        file=open(to_file,'w')
        file.write('first line is a comment \n')
        file.write('1\n')
        (avec,bvec,cvec) = atomic_coordinate.get_unit_cell().get_lattice_vector(to_unit=config.ANGSTROM)
        file.write(str(avec[0])+" "+str(avec[1])+" "+str(avec[2])+'\n')
        file.write(str(bvec[0])+" "+str(bvec[1])+" "+str(bvec[2])+'\n')
        file.write(str(cvec[0])+" "+str(cvec[1])+" "+str(cvec[2])+'\n')

        elements = atomic_coordinate.get_elem_list()
        for element in elements:
            file.write(element+' ')
        file.write('\n')
        for element in elements:
            file.write(str(self._get_num_elem(element,atomic_coordinate))+' ')
        file.write('\n')
        if self.mobile_atms is not None:
            file.write('Selective dynamics\n')
        file.write('cartesian\n')

        for element in elements:
            for atom in atomic_coordinate:
                if atom.get_element_name() == element:
                    pos = atom.get_pos(mode=config.CARTESIAN, to_unit=config.ANGSTROM)
                    if self.mobile_atms is not None:
                        if atom.idi in self.mobile_atms:
                            file.write(str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+' T T T\n')
                        else:
                            file.write(str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+' F F F\n')
                    else:
                        file.write(str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+'\n')

        file.flush()
        file.close()

        # ppfile
        ppdir = properties.Property().get_property('pp.vasp_pp_dir')
        if ppdir is not None and os.path.exists(ppdir) and os.path.isdir(ppdir):
            potcar = open('POTCAR','w')
            potfound=0
            for element in elements:
                dirs = os.listdir(ppdir)
                for dir in dirs:
                    if str(dir).startswith(element):
                        pot = open(ppdir+os.sep+dir+os.sep+'POTCAR',mode='r')
                        for l in pot:
                            potcar.write(l)
                        pot.close()
                        potfound += 1
                        break

            if potfound<len(elements):
                logger.warn('failed to resolve pp files for some of the elements')

    def _get_num_elem(self,elemname,coord):
        ret = 0
        for atom in coord:
            if elemname==atom.get_element_name():
                ret+=1
        return ret

    def get_name(self):
        return 'VASP input format'

class VASP_output(config.AtomConfigGenerator):
    def __init__(self,coord_file):
        super().__init__()
        if coord_file is None:
            self.set_is_valid(False)
            return
        self.inpfile = coord_file

    def supports_frame(self):
        return True

    def supports_export(self):
        return False

    def _gen_atomic_configuration(self):
        f = open(self.inpfile,'r')
        title_read = False
        fac_read = False
        avec = None
        bvec = None
        cvec = None
        cell = None
        nat_per_elem_read = False
        elemarray = []
        elems = None
        fac = 0.0
        frame = []
        iat = 0
        for line in f:
            line = line.strip()
            if not title_read:
                title_read = True
                continue
            if not fac_read:
                fac = float(line)  
                fac_read = True
                continue
            words = line.split()
            if avec is None:
                avec = (float(words[0])*fac,float(words[1])*fac,float(words[2])*fac) 
                continue
            if bvec is None:
                bvec = (float(words[0])*fac,float(words[1])*fac,float(words[2])*fac) 
                continue
            if cvec is None:
                cvec = (float(words[0])*fac,float(words[1])*fac,float(words[2])*fac) 
                cell = config.UnitCell(a_vector=[avec[0],avec[1],avec[2]],\
                                       b_vector=[bvec[0],bvec[1],bvec[2]],\
                                       c_vector=[cvec[0],cvec[1],cvec[2]])
                cell.set_length_unit(config.ANGSTROM)
                continue
            if line.lower().startswith('direct'):
                curr_config = config.AtomConfig(\
                coordinate_system=config.FRACTIONAL,length_unit=config.ANGSTROM,\
                unit_cell=cell)
                frame.append(curr_config) 
                iat = 0
                continue
            if elems is None:
                elems = words
                continue
            if not nat_per_elem_read:
                for i in range(len(words)):
                    for j in range(int(words[i])):
                        elemarray.append(elems[i])
                nat_per_elem_read = True
                continue
            xx = float(words[0])
            yy = float(words[1])
            zz = float(words[2])
            if xx>0.95:
                xx -= 1
            if yy>0.95:
                yy -= 1
            if zz>0.95:
                zz -= 1
            curr_config.add_atom(element=elemarray[iat],rx = xx,ry = yy,rz = zz)
            iat += 1
        return frame

    def export_atomic_configuration(self,atomic_coordinates,to_file=None,frame_no=None,all_frames=False):
        pass

    def get_name(self):
        return 'VASP output format'

    def get_defaultname(self,args):
        return 'XDATCAR'

