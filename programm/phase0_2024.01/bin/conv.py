#!/usr/bin/env python3

import logging
logger = logging.getLogger('conv')

import piou

from piou.util import pyutil
from piou import config
from piou.config import genconfig

class Converter(piou.PHASE_IO_utility):
    from_coords_generator=None
    from_coordinate_type=None
    from_coordinate_file=None
    from_frame_no=None

    to_coords_generator=None
    to_coordinate_type=None
    to_coordinate_file=None
    to_frame_no=None

    start_end_skip = None

    ibtype=None
    obtype=None

    doall=False
    pack=False

    na = 1
    nb = 1
    nc = 1

    def __init__(self):
        pass

    def fill_opts(self,parser):
        parser.add_option("","--input_type",type="choice",choices=genconfig.coordinate_file_type,\
        dest="input_typ", help="specify the type of the input atomic coordinates")
        parser.add_option("","--output_type",type="choice",choices=genconfig.coordinate_file_type,\
        dest="output_typ", help="specify the type of the output atomic coordinates")
        parser.add_option("","--input_file",dest="input_file",\
        help="specify the type of the input atomic coordinate file",default=None)
        parser.add_option("","--output_file",dest="output_file",\
        help="specify the type of the output atomic coordinate file",default=None)
        parser.add_option("-f","--frameno",dest="frame_no",type=int,\
        help="specify the target frame no. (enter a negative value in order to output all frames when possible)",\
        default=-1)
        parser.add_option("","--ibtype",dest="ibtype",type=str,\
        help="specify the type for the babel program (for input)",default="")
        parser.add_option("","--obtype",dest="obtype",type=str,\
        help="specify the type for the babel program (for output)",default="")
        parser.add_option("","--babel_opts",dest="babel_opts",type=str,\
        help="options to be passed to the babel program",default="")
        parser.add_option("-a","--all",dest="doall",action='store_true',\
        help="specify this option if you want to do conversion against all supported file formats.",default=False)
        parser.add_option("-p","--pack",dest="pack",action='store_true',\
        help="specify this option if you want to 'pack' the atomic coordinates into the unit cell.",default=False)
        parser.add_option("","--na",dest="na",type=int,\
        help="specify supercell size for the a-axis",default=1)
        parser.add_option("","--nb",dest="nb",type=int,\
        help="specify supercell size for the b-axis",default=1)
        parser.add_option("","--nc",dest="nc",type=int,\
        help="specify supercell size for the c-axis",default=1)

    def initialize(self,options):
        if options.batch:
            self.from_coordinate_type = options.input_typ
            self.to_coordinate_type = options.output_typ
            self.from_frame_no = options.frame_no
            self.from_coordinate_file = options.input_file
            self.to_coordinate_file = options.output_file
            self.ibtype=options.ibtype
            self.obtype=options.obtype
            self.pack=options.pack
            self.na=options.na
            self.nb=options.nb
            self.nc=options.nc
            
            self.from_coords_generator = genconfig.instantiate_coords_generator(\
            self.from_coordinate_type,self.from_coordinate_file)
            self.from_coords_generator.set_subtype(self.ibtype)
            self.to_coords_generator = genconfig.instantiate_coords_generator(\
            self.to_coordinate_type,self.to_coordinate_file)
            self.to_coords_generator.set_subtype(self.obtype)

        else:
            self.pack = options.pack
            self.na=options.na
            self.nb=options.nb
            self.nc=options.nc
            (self.from_coords_generator,self.from_coordinate_type,\
            self.from_coordinate_file,self.from_frame_no,self.start_end_skip)=\
            genconfig.select_coordinate_file_interactively(\
                   msg="the input atomic coordinate file",default_typ=0,\
                   msg_frame='the frame no. (enter a negative value in order to output all frames when possible)',\
                   extended_frame_spec=True,mode='r')
            if not options.doall:
                (self.to_coords_generator,self.to_coordinate_type,self.to_coordinate_file,self.to_frame_no)=\
                genconfig.select_coordinate_file_interactively(\
                   msg="the output atomic coordinate file",default_typ=1,\
                   msg_frame="the frame no. (enter a negative value in order to output all frames when possible)",\
                   must_exist=False,mustnot_exist=True,ask_frame_no=False,ask_mode=True,check_export_support=True,mode='w')
            else:
                self.doall=True

    def get_name(self):
        return 'converter'

    def get_description(self):
        return \
        'atomic configuration converter utility.'

    def run(self):

        #self.from_coords_generator = genconfig.instantiate_coords_generator(\
        #self.from_coordinate_type,self.from_coordinate_file)

        #self.to_coords_generator = genconfig.instantiate_coords_generator(\
        #self.to_coordinate_type,self.to_coordinate_file)

        atom_coords = self.from_coords_generator.get_frame()
        if atom_coords is None or len(atom_coords)==0:
            logger.error('failed to get an atomic coordinates instance from '+self.from_coordinate_file)
            return

        if not self.doall:
            logger.info("converting input file '"+self.from_coordinate_file+"' ("+self.from_coordinate_type+' file)'+\
            " to output file '"+self.to_coordinate_file+"' ("+self.to_coordinate_type+' file)')
            thecoords=self.get_target_coords(atom_coords,self.to_coords_generator.supports_frame())
            if self.na>1 or self.nb>1 or self.nc>1:
                for coor in thecoords:
                    coor.to_super(self.na,self.nb,self.nc)
            if self.pack:
                for coor in thecoords:
                    coor.pack()
            self.to_coords_generator.export_atomic_configuration(\
            thecoords,to_file=self.to_coordinate_file,all_frames=True)
        else:
            generators = genconfig.get_all_supported_generators()
            for generator in generators:
                logger.info('exporting to : '+generator.get_name())
                thecoords=self.get_target_coords(atom_coords,generator.supports_frame())
                generator.export_atomic_configuration(thecoords,\
                to_file=generator.get_defaultname(None),all_frames=True)

        logger.info('... done')
        logger.info('')

    def get_target_coords(self,origcoords,supports_frame):
        retcoords=origcoords
        if self.start_end_skip is not None:
            istart = self.start_end_skip[0]
            iend = self.start_end_skip[1]
            iskip = self.start_end_skip[2]
            if istart>=len(origcoords):
                istart=0
            if iend<0:
                iend = len(origcoords)-1
            if iend<istart:
                iend=istart
            if iskip>=len(origcoords):
                iskip=1
            tmp=[]
            for i in range(len(origcoords)):
                if i<istart:
                    continue
                if i>iend:
                    break
                if i%iskip!=0:
                    continue
                tmp.append(origcoords[i])
            retcoords = tmp
        else:
            if self.from_frame_no>=0 and self.from_frame_no < len(origcoords):
                retcoords = [origcoords[self.from_frame_no]]

            if (self.from_frame_no<0 or self.from_frame_no>=len(origcoords)) and not supports_frame:
                retcoords = [origcoords[len(origcoords)-1]]

        logger.info('the number of frames to output : '+str(len(retcoords)))
        return retcoords

if __name__ == '__main__':
    piou.run_main(Converter)

