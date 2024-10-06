import logging
logger = logging.getLogger('piou.config.cif')
from piou import config
from piou.config import AtomConfigGenerator

class CIF(AtomConfigGenerator):
    coord_file=None
    valid = True
    def __init__(self,coord_file=None):
        super().__init__()
        if coord_file is None:
            #logger.error("invalid file")
            return
        self.coord_file = coord_file

    def supports_frame(self):
        return True

    def _gen_atomic_configuration(self):
        coords = []
        ciffile = open(self.coord_file,'r')
        frameno=-1
        inComment=False
        coordsys=config.FRACTIONAL
        unit=config.ANGSTROM
        atoms=[]
        a=None
        b=None
        c=None
        alpha=None
        beta=None
        gamma=None
        inAtomSiteHeader=False
        inAtomSite=False
        for line in ciffile:
            line=line.strip()
            if line.startswith('#'):
                continue
            if line.startswith(';') and not inComment:
                inComment=True
            elif line.startswith(';') and inComment:
                inComment=False
            if inComment:
                continue
            if len(line)==0:
                continue
            words = line.split()
            if line.startswith('data_'):
                ll = line.split('data_')
                if len(ll)>1:
                    logger.info('procesing data : '+ll[1])
                frameno += 1
                currconf = config.AtomConfig(\
                coordinate_system=config.FRACTIONAL, length_unit=config.ANGSTROM)
                coords.append(currconf)
                a=None
                b=None
                c=None
                alpha=None
                beta=None
                gamma=None
                cell_ready=False
                atom_site_key=[]
                inAtomSiteHeader=False
                inAtomSite=False

            a=self._try_cell('_cell_length_a',words,a)
            b=self._try_cell('_cell_length_b',words,b)
            c=self._try_cell('_cell_length_c',words,c)
            alpha=self._try_cell('_cell_angle_alpha',words,alpha)
            beta=self._try_cell('_cell_angle_beta',words,beta)
            gamma=self._try_cell('_cell_angle_gamma',words,gamma)
            if currconf is not None and a is not None and b is not None and c is not None and\
                alpha is not None and beta is not None and gamma is not None:
                currconf.set_unit_cell(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)
                currconf.get_unit_cell().set_length_unit(config.ANGSTROM)
                cell_ready = True
            if line.startswith('_atom_site') and cell_ready:
                atom_site_key.append(line)
                inAtomSiteHeader = True
                inAtomSite = False
            elif line.startswith('_atom_site') and not cell_ready:
                atom_site_key.append(line)
                inAtomSiteHeader = True
                inAtomSite = False
            if inAtomSiteHeader and not line.startswith('_atom_site'):
                inAtomSiteHeader=False
                if '_atom_site_type_symbol' in atom_site_key and \
                    ('_atom_site_fract_x' in atom_site_key and \
                     '_atom_site_fract_y' in atom_site_key and \
                     '_atom_site_fract_z' in atom_site_key) or\
                    ('_atom_site_Cartn_x' in atom_site_key and \
                     '_atom_site_Cartn_y' in atom_site_key and \
                     '_atom_site_Cartn_z' in atom_site_key):
                    inAtomSite=True
                else:
                    inAtomSite=False
            if inAtomSite:
                if line.startswith('_') or line.startswith('loop_'):
                    inAtomSite=False
                    atom_site_key=[]
                    a=None
                    b=None
                    c=None
                    alpha=None
                    beta=None
                    gamma=None
                    cell_ready=False
                    continue
                label_x=None
                label_y=None
                label_z=None

                if '_atom_site_fract_x' in atom_site_key:
                    label_x=atom_site_key.index('_atom_site_fract_x')
                elif '_atom_site_Cartn_x' in atom_site_key:
                    label_x=atom_site_key.index('_atom_site_Cartn_x')
                    currconf.set_coordinate_system(config.CARTESIAN)

                if '_atom_site_fract_y' in atom_site_key:
                    label_y=atom_site_key.index('_atom_site_fract_y')
                elif '_atom_site_Cartn_y' in atom_site_key:
                    label_y=atom_site_key.index('_atom_site_Cartn_y')

                if '_atom_site_fract_z' in atom_site_key:
                    label_z=atom_site_key.index('_atom_site_fract_z')
                elif '_atom_site_Cartn_z' in atom_site_key:
                    label_z=atom_site_key.index('_atom_site_Cartn_z')

                if label_x is None or label_y is None or label_z is None:
                    logger.error('invalid coordinate : '+line)
                    continue
                label_sym=atom_site_key.index('_atom_site_type_symbol')
                if label_x>=0 and label_y>=0 and label_z>=0 and label_sym>=0 and len(words)>=4:
                    try:
                        rx=float(words[label_x].split(')')[0].split('(')[0])
                        ry=float(words[label_y].split(')')[0].split('(')[0])
                        rz=float(words[label_z].split(')')[0].split('(')[0])
                    except ValueError:
                        logger.error("invalid line : "+line)
                        self.valid=False
                        continue
                    currconf.add_atom(element=words[label_sym],rx=rx,ry=ry,rz=rz)

        ciffile.close()

        for conf in coords:
            if conf.get_num_atom()==0:
                coords.remove(conf)

        return coords

    def _try_cell(self,tag,words,myself):
        if myself is not None:
            return myself
        if words is None or len(words)==0:
            return None
        if words[0].startswith(tag) and len(words)>1:
            try:
                aa = float(words[1].split(')')[0].split('(')[0])
                return aa
            except ValueError:
                logger.error('invalid line : '+line)
        return None

    def export_atomic_configuration(self,atomic_coordinates,to_file=None,frame_no=None,all_frames=True):
        frameno = len(atomic_coordinates)-1
        if not all_frames and frame_no>=0 and frame_no<frameno:
            frameno = frame_no
        file = open(to_file,mode='w')
        for icoord in range(len(atomic_coordinates)):
            if not all_frames and icoord != frameno:
                continue
            coord = atomic_coordinates[icoord]
            cellok=True
            if not coord.is_cell_ready():
                logger.error('invalid unit cell')
                logger.error('only the atomic coordinates will be exported')
                cellok=False
                #continue
            file.write('data_frame'+str(icoord)+'\n')
            if cellok:
                unitcell = coord.get_unit_cell()
                vol = unitcell.get_volume(to_unit=config.ANGSTROM)
                file.write('_symmetry_cell_setting triclinic\n')
                file.write("_symmetry_space_group_name_H-M  'P1'\n")
                file.write("_symmetry_Int_Tables_number  1\n")
                (a,b,c,alpha,beta,gamma) = unitcell.get_lattice_constant(to_unit=config.ANGSTROM)
                file.write("_cell_length_a "+str(a)+'\n')
                file.write("_cell_length_b "+str(b)+'\n')
                file.write("_cell_length_c "+str(c)+'\n')
                file.write("_cell_angle_alpha "+str(alpha)+'\n')
                file.write("_cell_angle_beta "+str(beta)+'\n')
                file.write("_cell_angle_gamma "+str(gamma)+'\n')
                file.write("_cell_volume  "+str(vol)+'\n')
            atmtag='    _atom_site_fract_'
            if not cellok:
                atmtag = '    _atom_site_Cartn_'
            file.write("_cell_formula_units_Z "+str(coord.get_num_atom())+'\n')
            file.write("loop_\n")
            file.write("    _atom_site_label\n")
            file.write(atmtag+"x\n")
            file.write(atmtag+"y\n")
            file.write(atmtag+"z\n")
            file.write("    _atom_site_type_symbol\n")
            file.write("    _atom_site_occupancy\n")
            for atom in coord:
                elem = atom.get_element_name()
                pos = atom.get_pos(mode=config.FRACTIONAL)
                if not cellok:
                    pos=atom.get_pos(mode=config.CARTESIAN,to_unit=config.ANGSTROM)
                file.write(elem+str(atom.index()+1)+" "+str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+" "+elem+" 1\n")

        file.flush()
        file.close()

    def get_name(self):
        return 'CIF (crystallographic information file) format'

    def get_defaultname(self,args):
        return 'coord.cif'

