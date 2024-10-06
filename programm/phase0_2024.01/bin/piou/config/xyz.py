import logging
logger = logging.getLogger('piou.config.xyz')
import os

from piou import config
from piou.config import cif

class XYZ(config.AtomConfigGenerator):
    xyzfile=None
    def __init__(self,coord_file):
        super().__init__()
        self.xyzfile = coord_file

    def get_defaultname(self,args):
        return 'coord.xyz'

    def get_name(self):
        return 'XYZ format'

    def export_atomic_configuration(self,atomic_coordinates,to_file=None,frame_no=None,all_frames=True):
        frameno = len(atomic_coordinates)-1
        if not all_frames and frame_no>=0 and frame_no<frameno:
            frameno = frame_no
        file = open(to_file,mode='w')
        for icoord in range(len(atomic_coordinates)):
            if not all_frames and icoord != frameno:
                continue
            coord = atomic_coordinates[icoord]
            unitcell = coord.get_unit_cell()
            file.write(str(coord.get_num_atom())+'\n')
            if coord.is_cell_ready():
                (a,b,c,alpha,beta,gamma) = unitcell.get_lattice_constant(to_unit=config.ANGSTROM)
                file.write(str(a)+" "+str(b)+" "+str(c)+" "+str(alpha)+" "+str(beta)+" "+str(gamma)+"\n")
            else:
                file.write('\n')
            for atom in coord:
                elem = atom.get_element_name()
                pos = atom.get_pos(mode=config.CARTESIAN, to_unit=config.ANGSTROM)
                file.write(elem+" "+str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+'\n')
        file.flush()
        file.close()

    def _gen_atomic_configuration(self):
        mdprefix = os.path.basename(self.xyzfile).split('.')[0]
        dir = os.path.dirname(self.xyzfile)
        mol2='grid.mol2'
        ci = mdprefix+'.cif'
        if len(dir)>0:
            mol2 = dir+os.sep+'grid.mol2'
        unit_cell = None
        if os.path.exists(mol2):
            f = open(mol2)
            in_grid=False
            cvec=[]
            for line in f:
                line = line.strip()
                if line=='grid file':
                    in_grid=True
                    continue
                if in_grid:
                    if line=='@<TRIPOS>BOND':
                        break
                    if line=='@<TRIPOS>ATOM':
                        continue
                    words=line.split()
                    if len(words)>=5:
                        try:
                            cvec.append(\
                            (float(words[2]),float(words[3]),float(words[4])))
                        except ValueError:
                            logger.warn('invalid grid.mol2 file')
                            continue
                    if line.startswith('8'):
                        break
            f.close()
            if len(cvec)>=5:
                avec=[cvec[4][0]-cvec[0][0],cvec[4][1]-cvec[0][1],cvec[4][2]-cvec[0][2]]
                bvec=[cvec[2][0]-cvec[0][0],cvec[2][1]-cvec[0][1],cvec[2][2]-cvec[0][2]]
                cvec=[cvec[1][0]-cvec[0][0],cvec[1][1]-cvec[0][1],cvec[1][2]-cvec[0][2]]
                unit_cell=config.UnitCell(a_vector=avec,b_vector=bvec,c_vector=cvec)
                unit_cell.set_length_unit(config.ANGSTROM)
        elif os.path.exists(ci):
            cifcoord=cif.CIF(coord_file=ci)
            ac = cifcoord.get_atomic_configuration()
            unit_cell = ac.get_unit_cell()

        f=open(self.xyzfile,'r')
        nat = None
        a=None
        b=None
        c=None
        alpha=None
        beta=None
        gamma=None
        in_header=True
        currcoords=None
        fra = []
        header_count=0

        for line in f:
            line=line.strip()

            if not in_header:
                if nat is None:
                    logger.error("invalid specification for the number of atoms")
                    return NOne
                words=line.split()
                if len(words)>=4:
                    try:
                        elem=words[0]
                        rx=float(words[1])
                        ry=float(words[2])
                        rz=float(words[3])
                        currcoords.add_atom(element=elem,rx=rx,ry=ry,rz=rz)
                    except ValueError:
                        continue
                if currcoords.get_num_atom()==nat:
                    in_header=True
                    header_count=0

            else:
                header_count += 1
                if header_count ==1:
                    try:
                        nat=int(line)
                    except ValueError:
                        logger.error(\
                        'the first line of an xyz-file must be the number of atoms')
                if header_count==2:
                    in_header=False
                    try:
                        words = line.split()
                        if len(words)>=6:
                            a = float(words[0])
                            b = float(words[1])
                            c = float(words[2])
                            alpha = float(words[3])
                            beta = float(words[4])
                            gamma = float(words[5])
                            unit_cell = config.UnitCell(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)
                            unit_cell.set_length_unit(config.ANGSTROM)
                            currcoords.set_unit_cell(cell=unit_cell)
                    except:
                        pass
                    currcoords=config.AtomConfig(\
                    coordinate_system=config.CARTESIAN,length_unit=config.ANGSTROM,unit_cell=unit_cell)
                    fra.append(currcoords)

        f.close()

        logger.info("the number of frames defined in this XYZ file : "+str(len(fra)))

        return fra

    def supports_frame(self):
        return True

