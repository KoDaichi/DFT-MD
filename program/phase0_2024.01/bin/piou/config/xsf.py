''' manipulate the XSF (XCrysDen Structure File) format file.
    supports only the 'CRYSTAL' section for the time being... '''

import logging
logger = logging.getLogger('piou.config.xsf')

from piou import config
from piou.data import elements

class XSF(config.AtomConfigGenerator):
    xsffile=None
    def get_defaultname(self,args):
        return 'coord.xsf'

    def get_name(self):
        return 'XSF (XCrysDen Structure File) format'

    def __init__(self,coord_file):
        super().__init__()
        #if coord_file is None:
        #    self.set_is_valid(False)
        #    logger.error('invalid file')
        #    return
        self.xsffile = coord_file

    def supports_volumetric_data(self):
        return True

    def supports_frame(self):
        return True

    def _gen_atomic_configuration(self):
        frame = []
        f=open(self.xsffile,mode='r')
        tag_crystal='CRYSTAL'
        tag_primvec='PRIMVEC'
        tag_convvec='CONVVEC'
        tag_primcoord='PRIMCOORD'
        tag_atoms='ATOMS'
        tag_animsteps='ANIMSTEPS'
        tag_begin_grid='BEGIN_BLOCK_DATAGRID'
        tag_datagrid='DATAGRID_3D'
        tag_begin_datagrid='BEGIN_DATAGRID_3D'
        tag_end_datagrid='END_DATAGRID_3D'
        tag_end_grid='END_BLOCK_DATAGRID'

        in_crystal=False
        in_primvec=False
        in_primcoord=False
        in_atom=False
        in_anim=False
        in_datagrid=False
        in_datagrid_data=False
        nanim=-1
        avec=None
        bvec=None
        cvec=None
        currconf=None
        nat=-1

        voldata=None
        vdata=None
        ndiv=None
        div=[]
        origin=None
        for line in f:
            line=line.strip()
            words=line.split()

            if line.startswith('#'):
                continue

            if len(line)==0:
                continue

            if words[0].startswith(tag_begin_grid):
                in_datagrid=True
                in_datagrid_data=False
                continue
            
            if in_datagrid and \
                (words[0].startswith(tag_begin_datagrid) or words[0].startswith(tag_datagrid)):
                in_datagrid_data=True
                continue

            if in_datagrid_data and words[0].startswith(tag_end_datagrid):
                in_datagrid_data=False
                if voldata is None:
                    voldata=[]
                avec=(div[2][0],div[2][1],div[2][2])
                bvec=(div[1][0],div[1][1],div[1][2])
                cvec=(div[0][0],div[0][1],div[0][2])
                unitcell = config.UnitCell(a_vector=avec,b_vector=bvec,c_vector=cvec)
                unitcell.set_length_unit(config.ANGSTROM)
                vvdata = config.VolumetricData(unitcell=unitcell,n1=ndiv[2],n2=ndiv[1],n3=ndiv[0],vdata=vdata)
                if origin is not None:
                    vvdata.set_origin(origin)
                currconf.add_volumetric_data(vvdata)
                voldata.append(vvdata)
                vdata=[]
                continue

            if in_datagrid and words[0].startswith(tag_end_grid):
                in_datagrid=False
                div=[]
                ndiv = None
                origin = None
                continue

            if in_datagrid and not in_datagrid_data:
                continue

            if in_datagrid_data and ndiv is None and len(words)>=3:
                try:
                    ndiv = []
                    ndiv.append(int(words[0]))
                    ndiv.append(int(words[1]))
                    ndiv.append(int(words[2]))
                    continue
                except ValueError:
                    logger.error('invalid data grid')
                    continue
             
            if in_datagrid_data and origin is None and len(words)>=3:
                origin=[]
                origin.append(float(words[0]))
                origin.append(float(words[1]))
                origin.append(float(words[2]))
                continue

            if in_datagrid_data and len(div)<3 and ndiv is not None and len(ndiv)==3 and len(words)>=3:
                try:
                    tmp=[]
                    tmp.append(float(words[0]))
                    tmp.append(float(words[1]))
                    tmp.append(float(words[2]))
                    div.append(tmp)
                    continue
                except ValueError:
                    logger.error('invalid spanning vector')
                    continue

            if in_datagrid_data and ndiv is not None and len(ndiv)==3 and len(div)==3:
                if vdata is None:
                    vdata=[]
                try:
                    for word in words:
                        vdata.append(float(word))
                except ValueError:
                    logger.error('non-float value found in grid data')
                    continue

            if in_datagrid:
                continue

            if words[0]==tag_primvec:
                in_primvec=True
                in_primcoord=False
                in_atom=False
                continue

            if words[0]==tag_animsteps:
                in_anim=True
                in_atom=False
                if len(words)>1:
                    try:
                        nanim=int(words[1])
                    except ValueError:
                        pass

            if words[0]==tag_primcoord:
                in_primcoord=True
                in_atom=False
                nat=-1
                continue

            if words[0].startswith('BEGIN'):
                in_atom=False
                in_primcoord=False
                continue

            if words[0]==tag_atoms:
                in_atom=True
                currconf=config.AtomConfig(coordinate_system=config.CARTESIAN,\
                length_unit=config.ANGSTROM)
                frame.append(currconf)
                continue

            if in_atom:
                if len(words)<4:
                    continue
                try:
                    elem = self._resolve_elem(words[0])
                    rx = float(words[1])
                    ry = float(words[2])
                    rz = float(words[3])
                    currconf.add_atom(element=elem,rx=rx,ry=ry,rz=rz)
                    if currconf.get_num_atom()==nat:
                        nat=-1
                    continue
                except ValueError:
                    logger.error('the atom line : '+line+' is not valid')
                    return None

            if in_primcoord:
                if words[0].startswith(tag_crystal):
                    continue
                if nat<0:
                    try:
                        nat = int(words[0])
                    except ValueError:
                        logger.error('the number of atoms is not given properly')
                        return None
                    currconf=config.AtomConfig(coordinate_system=config.CARTESIAN,\
                    length_unit=config.ANGSTROM,a_vector=avec,b_vector=bvec,c_vector=cvec)
                    frame.append(currconf)
                    continue
                if len(words)<4:
                    logger.error(line)
                    logger.error('the primcoord specification must have atleast 4 words')
                    return None
                try:
                    elem = self._resolve_elem(words[0])
                    rx = float(words[1])
                    ry = float(words[2])
                    rz = float(words[3])
                    currconf.add_atom(element=elem,rx=rx,ry=ry,rz=rz)
                    if currconf.get_num_atom()==nat:
                        nat=-1
                    continue
                except ValueError:
                    logger.error('the primcoord line : '+line+' is not valid')
                    return None

            if in_primvec:
                if avec is not None and bvec is not None and cvec is not None:
                    avec=None
                    bvec=None
                    cvec=None
                if avec is None:
                    try:
                        avec=[float(words[0]),float(words[1]),float(words[2])]
                    except IndexError:
                        logger.error('invalid primvec specification : '+line)
                        return None
                    except ValueError:
                        logger.error('invalid primvec specification : '+line)
                        return None
                elif bvec is None:
                    try:
                        bvec=[float(words[0]),float(words[1]),float(words[2])]
                    except IndexError:
                        logger.error('invalid primvec specification : '+line)
                        return None
                    except ValueError:
                        logger.error('invalid primvec specification : '+line)
                        return None
                elif cvec is None:
                    try:
                        cvec=[float(words[0]),float(words[1]),float(words[2])]
                    except IndexError:
                        logger.error('invalid primvec specification : '+line)
                        return None
                    except ValueError:
                        logger.error('invalid primvec specification : '+line)
                        return None
                if avec is not None and bvec is not None and cvec is not None:
                    in_primvec=False
                    continue

        logger.info('the number of frames in this XSF file : '+str(len(frame)))

        if voldata is not None:
            for fra in frame:
                for vol in voldata:
                    fra.add_volumetric_data(vol)
        return frame

    def _resolve_elem(self,elem):
        try:
            ielem = int(elem)
            return elements.ElementInfo().get_element_name_from_atomic_number(ielem)
        except ValueError:
            return str(elem)

    def export_atomic_configuration(self,atomic_coordinates,to_file=None,frame_no=None,all_frames=True):
        file=open(to_file,mode='w')
        nfra = len(atomic_coordinates)
        if nfra>1:
            file.write('ANIMSTEPS '+str(nfra)+'\n')
        for ifra in range(len(atomic_coordinates)):
            frameid=''
            if nfra>1:
                frameid+= ' '+str(ifra+1)
            atomcoord = atomic_coordinates[ifra]
            vdata = atomcoord.get_volumetric_data()
            if vdata is None:
                if atomcoord.is_cell_ready():
                    file.write('CRYSTAL\n')
                    primvec = 'PRIMVEC'+frameid
                    file.write(primvec+'\n')
                    (avec,bvec,cvec) = atomcoord.get_unit_cell().get_lattice_vector(to_unit=config.ANGSTROM)
                    file.write(str(avec[0])+' '+str(avec[1])+' '+str(avec[2])+'\n')
                    file.write(str(bvec[0])+' '+str(bvec[1])+' '+str(bvec[2])+'\n')
                    file.write(str(cvec[0])+' '+str(cvec[1])+' '+str(cvec[2])+'\n')
                    primcoord='PRIMCOORD'+frameid
                    file.write(primcoord+'\n')
                    file.write(str(atomcoord.get_num_atom())+' 1\n')
                else:
                    atom = 'ATOMS'+frameid
                    file.write(atom+'\n')
                    #file.write(str(atomcoord.get_num_atom())+' \n')
                for atom in atomcoord:
                    elem = atom.get_element_name()
                    pos = atom.get_pos(mode=config.CARTESIAN,to_unit=config.ANGSTROM)
                    file.write(elem+' '+str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+'\n')
            else:
                atom = 'ATOMS'+frameid
                file.write(atom+'\n')
                for atom in atomcoord:
                    elem = atom.get_element_name()
                    pos = atom.get_pos(mode=config.CARTESIAN,to_unit=config.ANGSTROM)
                    file.write(elem+' '+str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+'\n')
                for vdat in vdata:
                    if not vdat.has_enough_data():
                        logger.error('not a valid grid data')
                        continue
                    file.write('BEGIN_BLOCK_DATAGRID_3D\n')
                    file.write('3D_grid_data\n')
                    file.write('DATAGRID_3D_G98CUBE\n')

                    ndiv=vdat.get_ndiv()
                    (avec,bvec,cvec) = atomcoord.get_unit_cell().get_lattice_vector(to_unit=config.ANGSTROM)
                    strndiv = str(int(ndiv[2]))+' '+str(int(ndiv[1]))+' '+str(int(ndiv[0]))
                    file.write(strndiv+'\n')
                    orig = vdat.get_origin()
                    file.write(str(orig[0])+' '+str(orig[1])+' '+str(orig[2])+'\n')
                    file.write(str(cvec[0])+' '+str(cvec[1])+' '+str(cvec[2])+'\n')
                    file.write(str(bvec[0])+' '+str(bvec[1])+' '+str(bvec[2])+'\n')
                    file.write(str(avec[0])+' '+str(avec[1])+' '+str(avec[2])+'\n')
                    vvdat = vdat.get_vdata()
                    count = 0
                    nrow=6
                    for v in vvdat:
                        file.write(str(v)+' ')
                        count += 1
                        if count%nrow==0:
                            file.write('\n')
                            count=0
                    file.write('\n')

                    file.write('END_DATAGRID_3D\n')
                    file.write('END_BLOCK_DATAGRID_3D\n')
        file.flush()
        file.close()

