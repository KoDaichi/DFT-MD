''' module for selecting and instantiating atomic configurations '''

import os
import sys
import logging
logger = logging.getLogger('piou.config.genconfig')

from piou.config import phase,openmx,cif,cube,xyz,vasp,xsf,dmol,obabel,lammps

from piou.util import pyutil

coordinate_file_type=(\
                  'phase_input','phase_output',\
                  'VASP_input','VASP_output',\
                  'OpenMX_input','OpenMX_output',\
                  'XSF','xyz','cube','cif','dmol',\
                  'LAMMPS_input', 'LAMMPS_output', \
                  'OpenBabel_interface'\
                  )

coord_generators={\
                  'phase_input':phase.Input,'phase_output':phase.F_DYNM,\
                  'VASP_input':vasp.VASP_input,'VASP_output':vasp.VASP_output,\
                  'OpenMX_input':openmx.OMX_input,'OpenMX_output':openmx.OMX_output,\
                  'LAMMPS_input':lammps.LAMMPS_input, 'LAMMPS_output':lammps.LAMMPS_output, \
                  'XSF':xsf.XSF,'xyz':xyz.XYZ,\
                  'cube':cube.Cube,'cif':cif.CIF,\
                  'dmol':dmol.DMol,'OpenBabel_interface':obabel.OpenBabelInterface\
                  }

def get_all_supported_generators(mode='w'):
    ''' return a list of all supported AtomicConfiguration generators '''
    tmpgenerators={}
    typetype=[]
    for ty in coordinate_file_type:
        tmpgen = coord_generators[ty](None)
        tmpgenerators[ty]=tmpgen
        if tmpgen.check_support() and not (mode=='w' and not tmpgen.supports_export()) and \
            not (mode=='r' and not tmpgen.supports_import()):
            typetype.append(ty)
    ret = []
    for ft in typetype:
        tmp=tmpgenerators[ft]
        stypes=tmp.get_all_subtypes(mode=mode)
        if stypes is not None:
            for stype in stypes:
                tmpgen = coord_generators[ft](None)
                tmpgen.set_subtype(stype)
                ret.append(tmpgen)
        else:
            ret.append(tmp)
    return ret

def instantiate_coords_generator(coordtyp,coordinate_file):
    coords_generator=None
    if coordtyp in coordinate_file_type:
        coordgen = coord_generators[coordtyp]
        if coordgen is None:
            logger.error('invalid coord-generator : '+coordtyp)
        else:
            coords_generator = coordgen(coord_file=coordinate_file)

    if coords_generator is None or not coords_generator.is_valid():
        logger.error("failed to instantiate the coords generator. ")
        return None

    return coords_generator

def select_coordinate_file_interactively(msg,default_typ,typ_lambda=None,msg_frame=None,must_exist=True,\
    mustnot_exist=False,ask_frame_no=True,ask_mode=False,extended_frame_spec=False,check_export_support=False,
    mode=None):
    ''' select the coordinate file interactively. 
        will return a tuple with the following elements : 
        (coord_generator,coordinate_type,coordinate_file,frameno) '''
    lam = lambda f:os.path.exists(f)
    invalid='file does not exist: '
    if not must_exist:
        lam=lambda f:True
    if mustnot_exist:
        lam=lambda f:not os.path.exists(f)
        invalid='file exists: '
    typetype=[]
    tmpgenerators={}
    for ty in coordinate_file_type:
        tmpgen = coord_generators[ty](None)
        tmpgenerators[ty]=tmpgen
        if not tmpgen.check_support() or\
            (mode == 'w' and not tmpgen.supports_export()) or\
            (mode == 'r' and not tmpgen.supports_import()):
            continue
        typetype.append(ty)
        
    while True:
        coordinate_type = pyutil.interactive(msg='select the type of '+msg,\
        condition=typ_lambda,default=default_typ,choice=typetype,return_index=False,\
        ask_mode=ask_mode)
        if check_export_support and not tmpgenerators[coordinate_type].supports_export():
            logger.error("output support unavailable for coordinate type '"+coordinate_type+"'")
        else:
            break

    subtype = tmpgenerators[coordinate_type].select_subtype(mode=mode) 

    defname = tmpgenerators[coordinate_type].get_defaultname(subtype)
    coordinate_file = pyutil.interactive(msg="the name "+msg, condition=lam,typ=str,\
                      default=defname,msg_invalid_input=invalid, ask_mode=ask_mode)

    coords_generator = instantiate_coords_generator(coordinate_type,coordinate_file)
    coords_generator.set_subtype(subtype)

    frame_no = -1
    iframe   = -1
    start_end_skip=None
    if coords_generator.supports_frame() and ask_frame_no:
        if msg_frame is None:
            msg_frame = "the frame no. (the last frame will be adopted if a negative value is specified)"
        while True:
            frame_no = pyutil.interactive(msg=msg_frame, condition=lambda ii : True,typ=str,default="-1",\
            ask_mode=ask_mode)
            try:
                iframe = int(frame_no)
                break
            except ValueError:
                pass
            start_end_skip=[]
            ff = frame_no.strip().split(',')
            if len(ff)<3:
                print('enter an int')
                continue
            else:
                try:
                    start_end_skip.append(int(ff[0]))
                    start_end_skip.append(int(ff[1]))
                    start_end_skip.append(int(ff[2]))
                    break
                except ValueError:
                    print('enter the start frame, end frame and the number of frames to skip '+\
                    'by integers delimited by commas')
                    continue

    if not extended_frame_spec:
        return (coords_generator,coordinate_type,coordinate_file,iframe)
    else:
        return (coords_generator,coordinate_type,coordinate_file,iframe,start_end_skip)

