''' manipulate OpenMX I/O files '''
import logging
logger = logging.getLogger('piou.config.openmx')

from piou import config
from piou.config import xyz

class OMX_input(config.AtomConfigGenerator):
    inpfile=None
    def __init__(self,coord_file):
        super().__init__()
        if coord_file is None:
            self.set_is_valid(False)
            return
        self.inpfile=coord_file

    def _gen_atomic_configuration(self):
        coord_key = 'Atoms.SpeciesAndCoordinates'
        cell_key  = 'Atoms.UnitVectors'
        file = open(self.inpfile,'r')
        inp={}
        in_table=False
        table_name=None
        for line in file:
            line=line.strip()
            if line.startswith('#'):
                continue
            line=line.split('#')[0].strip()
            if not (line.startswith('<') or line.startswith('>')):
                words=line.split()
                if len(words)>1:
                    inp[words[0]] = words[1]
            elif line.startswith('<'):
                table_name=line[1:]
                inp[table_name] = []
                in_table=True
            elif line.endswith('>'):
                in_table=False
            if in_table and not (line.startswith('<') or line.startswith('>')):
                inp[table_name].append(line)
        file.close()

        logger.debug("the OpenMX input file dict")
        logger.debug(str(inp))

        coord_unit = 'Ang'
        if coord_key+'.Unit' in inp:
            coord_unit = inp[coord_key+'.Unit']

        cell_unit='Ang'
        if cell_key+'.Unit' in inp:
            cell_unit = inp[cell_key+'.Unit']

        coord_unit_p = config.ANGSTROM
        coord_sys_p = config.CARTESIAN
        if coord_unit=='AU':
            coord_unit_p = config.BOHR
        if coord_unit=='FRAC':
            coord_sys_p = config.FRACTIONAL
        conf = config.AtomConfig(coordinate_system=coord_sys_p, length_unit=coord_unit_p)

        cell_unit_p = config.ANGSTROM
        if cell_unit=='AU':
            cell_unit_p = config.BOHR
        ucell = config.UnitCell()
        ucell.set_length_unit(cell_unit_p)
        conf.set_unit_cell(ucell)
        if cell_key in inp:
            cellvecs = inp[cell_key]
            if len(cellvecs)<3:
                logger.error('invalid cell : '+str(cellvecs))
            else:
                avec=cellvecs[0].split()
                bvec=cellvecs[1].split()
                cvec=cellvecs[2].split()
                if len(avec)>=3 and len(bvec)>=3 and len(cvec)>=3:
                    try:
                        ucell.set_cell(\
                        a_vector=[float(avec[0]),float(avec[1]),float(avec[2])],\
                        b_vector=[float(bvec[0]),float(bvec[1]),float(bvec[2])],\
                        c_vector=[float(cvec[0]),float(cvec[1]),float(cvec[2])])
                    except ValueError:
                        logger.error('invalid cell vector : '+\
                        str(cellvecs[0])+" "+str(cellvecs[1])+" "+str(cellvecs[2]))

        if coord_key in inp:
            atomlist = inp[coord_key]
            for atom in atomlist:
                ats = atom.split()
                if len(ats)<5:
                    continue
                elem = ats[1]
                try:
                    rx = float(ats[2])
                    ry = float(ats[3])
                    rz = float(ats[4])
                    conf.add_atom(element=elem,rx=rx,ry=ry,rz=rz)
                except ValueError:
                    logger.error('invalid atom specification : '+str(atom))
                
        return [conf]

    def export_atomic_configuration(self,atomic_coordinates,to_file=None,frame_no=None,all_frames=False):
        frameno=config.resolve_frame_no(frame_no,len(atomic_coordinates))
        atomic_coordinate = atomic_coordinates[frameno]
        file=open(to_file,mode='w')
        nat=atomic_coordinate.get_num_atom()
        file.write('Atoms.Number '+str(nat)+'\n')
        file.write('Atoms.SpeciesAndCoordinates.Unit Ang\n')
        file.write('<Atoms.SpeciesAndCoordinates\n')
        atom_list=atomic_coordinate.get_atom_list()
        for i in range(nat):
            atom = atom_list[i]
            elemname = atom.get_element_name()
            [rx,ry,rz] = atom.get_pos(mode=config.CARTESIAN, to_unit=config.ANGSTROM)
            file.write(str(i+1)+' '+elemname+' '+str(rx)+' '+str(ry)+' '+str(rz)+'\n')
        file.write('Atoms.SpeciesAndCoordinates>\n')
        if atomic_coordinate.is_cell_ready():
            file.write('Atoms.UnitVectors.Unit Ang\n')
            (avec,bvec,cvec) = atomic_coordinate.get_unit_cell().get_lattice_vector(to_unit=config.ANGSTROM)
            file.write('<Atoms.UnitVectors\n')
            file.write(str(avec[0])+' '+str(avec[1])+' '+str(avec[2])+'\n')
            file.write(str(bvec[0])+' '+str(bvec[1])+' '+str(bvec[2])+'\n')
            file.write(str(cvec[0])+' '+str(cvec[1])+' '+str(cvec[2])+'\n')
            file.write('Atoms.UnitVectors>\n')
        file.flush()
        file.close()

    def get_name(self):
        return 'OpenMX input format'

    def get_defaultname(self,args):
        return 'omx.dat'

import os 
from piou.config import cif
class OMX_output(config.AtomConfigGenerator):
    mdfile=None
    mdprefix=None

    def supports_export(self):
        return False

    def __init__(self,coord_file):
        if coord_file is None:
            self.set_is_valid(False)
            return
        self.mdfile = coord_file
        self.mdprefix = os.path.basename(coord_file).split('.')[0]

    def supports_frame(self):
        return True

    def _gen_atomic_configuration(self):
        ci = self.mdprefix+'.cif'
        unitcell=None
        if os.path.exists(ci):
            cifcoord=cif.CIF(coord_file=self.mdprefix+".cif")
            ac = cifcoord.get_atomic_configuration()
            unitcell = ac.get_unit_cell()
        fra = []
        f = open(self.mdfile,'r')
        natm=0
        for line in f:
            line = line.strip()
            try:
                natm = int(line)
                conf = config.AtomConfig(coordinate_system=config.CARTESIAN,\
                length_unit=config.ANGSTROM,unit_cell=unitcell)
                fra.append(conf)
                continue
            except ValueError:
                pass
            if line.startswith('time'):
                continue
            words = line.split()
            if len(words)>=4:
                elem = words[0]
                try:
                    rx=float(words[1])
                    ry=float(words[2])
                    rz=float(words[3])
                    conf.add_atom(element=elem,rx=rx,ry=ry,rz=rz)
                except ValueError:
                    logger.error('invalid line  : '+line)
                    continue
        f.close()

        return fra

    def export_atomic_configuration(self,atomic_coordinates,to_file=None,frame_no=None,all_frames=True):
        ''' the xxx.md file seems to be equivalent to the xyz format... '''
        xyzfile=xyz.XYZ(to_file)
        xyzfile.export_atomic_configuration(atomic_coordinates,to_file,frame_no,all_frames)

    def get_name(self):
        return 'OpenMX output format'

    def get_defaultname(self,args):
        return 'omx.md'

