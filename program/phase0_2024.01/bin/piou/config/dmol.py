import logging
logger = logging.getLogger('piou.config.dmol')
from piou import config

class DMol(config.AtomConfigGenerator):
    inpfile=None
    def __init__(self,coord_file):
        super().__init__()
        self.inpfile=coord_file

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
            if coord.is_cell_ready():
                (avec,bvec,cvec) = unitcell.get_lattice_vector(to_unit=config.BOHR)
                file.write('$cell vectors\n')
                file.write('    '+str(avec[0])+' '+str(avec[1])+' '+str(avec[2])+'\n')
                file.write('    '+str(bvec[0])+' '+str(bvec[1])+' '+str(bvec[2])+'\n')
                file.write('    '+str(cvec[0])+' '+str(cvec[1])+' '+str(cvec[2])+'\n')
            file.write('$coordinates\n')
            for atom in coord:
                pos = atom.get_pos(mode=config.CARTESIAN, to_unit=config.BOHR)
                elem = atom.get_element_name()
                file.write(elem+" "+str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+'\n')
            file.flush()
            file.close()

    def _gen_atomic_configuration(self):
        file=open(self.inpfile)
        incell=False
        incoords=False
        cellvec=[]
        atms = []
        for line in file:
            line=line.strip()
            words=line.split()
            if line.startswith('$cell vectors'):
                incell=True
                continue
            if line.startswith('$coordinates'):
                incoords=True
            if incell and not incoords and len(words)>=3:
                try:
                    vec = (float(words[0]),float(words[1]),float(words[2]))
                    cellvec.append(vec)
                    continue
                except ValueError:
                    logger.error('invalid cell : '+line)
                    continue
            if incoords and len(words)>=4:
                try:
                    elem=words[0]
                    rx = float(words[1])
                    ry = float(words[2])
                    rz = float(words[3])
                    atms.append((elem,rx,ry,rz))
                except ValueError:
                    logger.error('invalid atomic coordinate : '+line)
                    continue
        conf = config.AtomConfig(coordinate_system=config.CARTESIAN, length_unit=config.BOHR, \
        a_vector=cellvec[0],b_vector=cellvec[1],c_vector=cellvec[2])
        for atm in atms:
            elem=atm[0]
            rx=atm[1]
            ry=atm[2]
            rz=atm[3]
            conf.add_atom(element=elem,rx=rx,ry=ry,rz=rz)

        return [conf]

    def get_name(self):
        return 'DMol3 format'

    def get_defaultname(self,args):
        return 'coord.dmol'

