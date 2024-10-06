#!/usr/bin/env python3

import logging
import os
import errno
import math
import sys

import piou

from piou.util import pyutil
from piou.util import phase_util
from piou.util import vector_op

from piou.data import properties

from piou import config

from piou.config import genconfig
from piou.config import phase

logger = logging.getLogger("geninp")

phase_input_template=\
" control{\n\
    condition = initial\n\
    cpumax = 1 day\n\
    max_iteration = 1000000000\n\
}\n\
accuracy{\n\
}\n\
structure{\n\
    atom_list{\n\
        atoms{\n\
            #tag element rx ry rz mobile\n\
        }\n\
    }\n\
    unit_cell{\n\
    }\n\
    element_list{\n\
        #tag element atomicnumber mass zeta deviation\n\
    }\n\
}\n\
structure_evolution{\n\
}\n\
wavefunction_solver{\n\
    solvers{\n\
        #tag sol till_n dts dte itr var prec cmix submat\n\
    }\n\
}\n\
charge_mixing{\n\
    mixing_methods{\n\
        #tag no method rmxs rmxe itr var prec istr nbmix update\n\
    }\n\
}\n\
postprocessing{\n\
}\n\
printoutlevel{\n\
    base=1\n\
}\n\
"

file_names_data_template= "\
&fnames\n\
F_INP='./nfinp.data'\n\
/\n\
"

yes_no=('yes','no')

calc_type=('static','stropt','phonon','dos','md_nvt','md_nve','neb')
convergence_choice=('very_conservative','conservative','normal','aggressive','very_aggressive')
accuracy_choice=('low','normal','high')


class InputGenerator(piou.PHASE_IO_utility):
    frame_no = -1
    coordinate_file = 'nfinp.data'
    coordinate_type = 'phase_input'
    atomic_coordinates = None
    outdir = 'stropt'
    input_type = 'stropt'
    symmetry=True
    spin=False
    valid=None
    conver=3
    accuracy=1
    temperature=100

# for neb
    nreplica=5
    coordinate_file_end0 = 'nfdynm.data'
    coordinate_file_end1 = 'nfdynm.data'
    coordinate_type_end0 = 'phase_output'
    coordinate_type_end1 = 'phase_output'
    atomic_coordinate_end0 = None
    atomic_coordinate_end1 = None
    frame_no_end0 = -1
    frame_no_end1 = -1
    options = None

# OpenBabel file type
    btype = None

    def fill_opts(self,parser):
        parser.add_option("-i","--input",dest="input",help="specify the atomic coordinate file ",\
        default="nfinp.data")

        parser.add_option("-t","--type",type="choice",choices=genconfig.coordinate_file_type,\
        dest="type",help="specify the type of the atomic coordinate file, "+\
        " ie 'phase_input' for PHASE input file",default="phase_input")

        parser.add_option("-f","--frameno",type=int,dest="frameno",\
        help="specify the frame number (enter a negative value in order to use the last frame)",default=-1)

        parser.add_option("-d","--dest",dest="result",\
        help="specify the destination directory to where the results are stored. defaults to 'stropt'"\
       ,default="stropt")

        parser.add_option("-c","--calc_type",type="choice",choices=calc_type,\
                help="specify the calculation type of PHASE input.",dest="pinp_type",default="stropt")

        parser.add_option("","--temperature",type=float,dest="tempera",\
        help="specify the target temperature in Kelvin units. makes sence only when calc_type is md_nvt or md_nve",\
        default=100)

        parser.add_option("","--spin",action="store_true",dest="spin",default=False,\
        help="set this option in order to perform spin-polarized calculations.")

        parser.add_option("","--nosymm",action="store_true",dest="nosymm",default=False,\
        help="set this option in order to ignore crystal symmetry.")

        parser.add_option("-v","--convergence",type="choice", choices=convergence_choice,\
        help="specify how aggressive the solvers and mixers should be configured.",\
        dest="conver",default='aggressive')

        parser.add_option("-a","--accuracy",type="choice",\
        choices=accuracy_choice, help="specify the expected accuracy of the calculation.",\
        dest="accura",default='normal')

        parser.add_option("","--end0_coord_type",dest="end0coord_type",\
        help="specify the corodinate type for the initial state (for NEB)",default="phase_output")

        parser.add_option("","--end0_coord_file",dest="end0coord_file",\
        help="specify the corodinate file for the initial state (for NEB)",default="nfdynm.data")

        parser.add_option("","--end0_frameno",type=int,dest="end0frameno",\
        help="specify the frame number for the initial state (for NEB)",default=-1)

        parser.add_option("","--end1_coord_type",dest="end1coord_type",\
        help="specify the corodinate type for the final state (for NEB)",default="phase_output")

        parser.add_option("","--end1_coord_file",dest="end1coord_file",\
        help="specify the corodinate file for the final state (for NEB)",default="nfdynm.data")

        parser.add_option("","--end1_frameno",type=int,dest="end1frameno",\
        help="specify the frame number for the final state (for NEB)",default=-1)

        parser.add_option("","--btype",type=str,dest="btype",\
        help="specify the file type for the Open Babel program (makes sense only when --type=OpenBabel_Interface)",\
        default="")

    def index_of_tuple(self,targ,tup):
        for i in range(len(tup)):
            if tup[i]==targ:
                return i
        return -1
    
    def set_options(opts):
        ''' set the option instance '''
        self.options = opts

    def __init__(self):
        pass

    def initialize(self,options):
        coords_generator=None
        if options.batch:
            self.coordinate_file = options.input
            self.coordinate_type = options.type
            coords_generator = genconfig.instantiate_coords_generator(\
                               self.coordinate_type,self.coordinate_file)
            if coords_generator is None:
                return
            self.outdir = options.result
            self.input_type = options.pinp_type
            self.spin = options.spin
            self.conver = self.index_of_tuple(options.conver,convergence_choice)
            self.accuracy = self.index_of_tuple(options.accura,accuracy_choice)
            self.temperature = int(options.tempera)
            if (options.nosymm and self.input_type!='phonon') or self.input_type.startswith('md'):
                self.symmetry = False
            self.frame_no = options.frameno
            coords_generator.set_subtype(options.btype)
            if self.input_type=='neb':
                self.coordinate_type_end0 = options.end0coord_type
                self.coordinate_type_end1 = options.end1coord_type
                self.coordinate_file_end0 = options.end0coord_file
                self.coordinate_file_end1 = options.end1coord_file
                self.frame_no_end0 = options.end0frameno
                self.frame_no_end1 = options.end1frameno
                self.atomic_coordinate_end0 =\
                coords_generator_end0.get_atomic_configuration(nframe=self.frame_no_end0)
                self.atomic_coordinate_end1 =\
                coords_generator_end1.get_atomic_configuration(nframe=self.frame_no_end1)
        else:
            (coords_generator,self.coordinate_type,self.coordinate_file,self.frame_no)=\
            genconfig.select_coordinate_file_interactively(\
                   msg="the atomic coordinate file",default_typ=0,mode='r')

            tmp = pyutil.interactive(msg="select calculation type", \
                              default=1,choice=calc_type,return_index=False,skippable=True)
            if tmp is not None:
                self.input_type = tmp
                self.outdir = tmp
            else:
                self.__init__sub(coords_generator)
                return

            tmp = pyutil.interactive("the name of the output directory",\
                                      lambda f: len(f)!=0,str,self.input_type,skippable=True)

            if tmp is not None:
                self.outdir = tmp
            else:
                self.__init__sub(coords_generator)
                return

            if self.input_type=='neb':
                (coords_generator_end0,self.coordinate_type_end0,\
                 self.coordinate_file_end0,self.frame_no_end0)=\
                genconfig.select_coordinate_file_interactively(\
                msg="the atomic coordinate file for the initial state",default_typ=1,mode='r')

                self.atomic_coordinate_end0 =\
                coords_generator_end0.get_atomic_configuration(nframe=self.frame_no_end0)

                (coords_generator_end1,self.coordinate_type_end1,\
                 self.coordinate_file_end1,self.frame_no_end1)=\
                genconfig.select_coordinate_file_interactively(\
                msg="the atomic coordinate file for the final state",default_typ=1,mode='r')

                self.atomic_coordinate_end1 =\
                coords_generator_end1.get_atomic_configuration(nframe=self.frame_no_end1)

                tmp = pyutil.interactive('the number of replica',\
                                                    lambda nrep: nrep>=1,int,"5",skippable=True)
                if tmp is not None:
                    self.nreplica = tmp
                else:
                    self.__init__sub(coords_generator)
                    return
            elif self.input_type=='md_nvt':
                tmp = pyutil.interactive('the target temperature in Kelvin units',\
                                                    lambda temp: temp>0,int,"100",skippable=True)
                if tmp is not None:
                    self.temperature = tmp
                else:
                    self.__init__sub(coords_generator)
                    return
            elif self.input_type=='md_nve':
                tmp = pyutil.interactive('the initial temperature in Kelvin units',\
                                                    lambda temp: temp>=0,int,"100",skippable=True)
                if tmp is not None:
                    self.temperature = tmp
                else:
                    self.__init__sub(coords_generator)
                    return
            elif not self.input_type=='phonon':
                tmp = pyutil.interactive('take symmetry into account?',\
                      default=0,choice=yes_no,return_index=False,skippable=True)
                if tmp is not None:
                    self.symmetry = pyutil.parse_bool(tmp)
                else:
                    self.__init__sub(coords_generator)
                    return

            tmp = pyutil.interactive('take spin into account?',\
                  default=1,choice=yes_no,return_index=False,skippable=True)
            if tmp is not None:
                self.spin = pyutil.parse_bool(tmp)
            else:
                self.__init__sub(coords_generator)
                return
            tmp=pyutil.interactive('select the accuracy of the calculation.',\
                                           default=1,choice=accuracy_choice,skippable=True)
            if tmp is not None:
                self.accuracy = tmp
            else:
                self.__init__sub(coords_generator)
                return

            self.conver=pyutil.interactive(\
                        'select how aggressive the solvers and charge-mixing should be configured.',\
                         default=3,choice=convergence_choice)
     
        self.valid=True
        self.atomic_coordinates = coords_generator.get_atomic_configuration()
        if self.atomic_coordinates is None:
            logger.error('failed to get a valid atomic configuration from file : '+self.coordinate_file)
            sys.exit()
        logger.info("the coordinates defined in this file : "+str(self.atomic_coordinates))
            
    def get_description(self):
        return \
        "input data generation utility for PHASE"

    def get_name(self):
        return 'The input generator utility'

    def __init__sub(self,coords_generator):
        self.valid=True
        self.atomic_coordinates = coords_generator.get_atomic_configuration(nframe=self.frame_no)
        if self.atomic_coordinates is None:
            logger.error('failed to get a valid atomic configuration from file : '+self.coordinate_file)
            sys.exit()
        logger.info("the coordinates defined in this file : "+str(self.atomic_coordinates))

    def generate_input(self):
        if not self.valid:
            return
        logger.info("generating the input file for PHASE under directory : "+self.outdir)
        try:
            os.makedirs(self.outdir)
        except OSError as exc:
            if exc.errno == errno.EEXIST:
                pass
            else: 
                logger.error("failed to make directory : "+self.outdir)
                return

        os.chdir(self.outdir)

        fnames = open('file_names.data','w')
        fnames.write(file_names_data_template)
        fnames.flush()
        fnames.close()

        nfinp = open('nfinp.data','w')
        nfinp.write(phase_input_template)
        nfinp.flush()
        nfinp.close()

        target_input = phase.Input(coord_file='nfinp.data')

        self.do_structure(target_input)
        self.do_nbands(target_input)
        self.do_ksamp(target_input)
        self.do_cutoff(target_input) 
        self.do_initial_wf_and_chg(target_input)
        self.do_scf_convergence(target_input)
        self.do_str_evol(target_input)
        self.do_wf_solver(target_input)
        self.do_cd_mixing(target_input)
        self.do_postprocessing(target_input)
        self.do_neb(target_input)

        target_input.get_filenames_data().save()
        target_input.save()

        if not self.atomic_coordinates.is_cell_ready():
            logger.warn('the unit cell was not properly specified in : '+\
            self.coordinate_file+", ")
            logger.warn('so you will have to manually specify '+\
            'the unit cell and k-point sampling to the generated input.')

    def do_structure(self,target_input):
        for atom in self.atomic_coordinates:
            if self.input_type=='static' or self.input_type=='phonon' or self.input_type=='dos':
                atom.set_mobile(False)
            else:
                atom.set_mobile(True)

        target_input.replace_atomic_configuration(\
        self.atomic_coordinates)

        if self.spin:
            target_input.get_root_entry().get_entry('structure').\
            add_entry(phase.PrimitiveEntry('magnetic_state',"ferro"))

        if self.symmetry and not (self.input_type=='md_nve' or self.input_type=='md_nvt'):
            target_input.get_root_entry().get_entry('structure').\
            add_entry(phase.BlockEntry('symmetry'))
            target_input.get_root_entry().get_entry('structure.symmetry').\
            add_entry(phase.PrimitiveEntry('method','automatic'))
            atmtable=target_input.get_root_entry().get_entry('structure.atom_list.atoms.table')
            ats = phase_util.get_attributes(atmtable,'weight')
            if ats is not None:
                found_weight_2=False
                for at in ats:
                    if at==2:
                        found_weight_2=True
                        break
                if found_weight_2:
                    inv = phase.PrimitiveEntry('sw_inversion','on')
                    target_input.get_root_entry().get_entry('structure.symmetry').add_entry(inv)
            inv=self.atomic_coordinates.has_inversion_symmetry()
            if inv:
                logger.info('the system has inversion symmetry')
                inv = phase.PrimitiveEntry('sw_inversion','on')
                target_input.get_root_entry().get_entry('structure.symmetry').add_entry(inv)

    def do_nbands(self,target_input):
        nv = target_input.get_num_valence_electrons()
        natm = target_input.get_atomic_configuration().get_num_atom()
        sfactor=1.0
        if self.spin:
            sfactor=1.3
        if nv is not None:
            logger.info("number of valence electrons : "+str(nv))
            fact = properties.Property().get_property("geninp.band_tolerance_factor",float)
            if fact is not None:
                nbmin = int(math.ceil(nv * 0.5 * fact * sfactor))
                nbm = int(math.ceil(nv*0.5+natm*0.5))
                if nbm>nbmin:
                    nbmin = nbm
                if nbmin%2!=0:
                    nbmin += 1
                target_input.get_root_entry().get_entry('accuracy').\
                             add_entry(phase.PrimitiveEntry('num_bands',nbmin))
        else:
            logger.error("failed to get the number of valence electrons.")

    def do_ksamp(self,target_input):
        if self.atomic_coordinates.is_cell_ready():
            kmesh_dk = properties.Property().get_property("inpcheck.kmesh_dk",float)
            fact=1.0
            if self.accuracy==2:
                fact = 2
            elif self.accuracy==0:
                fact = 0.5
            if kmesh_dk is not None:
                (rb1,rb2,rb3) = \
                config.get_reference_kmesh(self.atomic_coordinates,kmesh_dk)
                if rb1 is None:
                    return
                rb1 = int(round(rb1*fact))
                rb2 = int(round(rb2*fact))
                rb3 = int(round(rb3*fact))
                if rb1 == 0:
                    rb1 = 1
                if rb2 == 0:
                    rb2 = 1
                if rb3 == 0:
                    rb3 = 1

                if self.input_type=='dos':
                    if rb1==1:
                        rb1=2
                    if rb2==1:
                        rb2=2
                    if rb3==1:
                        rb3=2

                gamma=False
                if rb1==1 and rb2==1 and rb3==1:
                    gamma=True
                blksamp = phase.BlockEntry('ksampling')
                if not gamma:
                    blmesh = phase.BlockEntry('mesh')
                    blksamp.add_entry(blmesh)
                    blmesh.add_entry(phase.PrimitiveEntry('nx',rb1))
                    blmesh.add_entry(phase.PrimitiveEntry('ny',rb2))
                    blmesh.add_entry(phase.PrimitiveEntry('nz',rb3))
                    if self.input_type=='dos':
                        sampm = phase.PrimitiveEntry('method','mesh')
                        blksamp.add_entry(sampm)
                        tetra = phase.PrimitiveEntry('method','tetrahedral')
                        blsmear = phase.BlockEntry('smearing')
                        blsmear.add_entry(tetra)
                        target_input.get_root_entry().get_entry('accuracy').\
                                     add_entry(blsmear)
                else:
                    sampm = phase.PrimitiveEntry('method','gamma')
                    red = phase.PrimitiveEntry('base_reduction_for_GAMMA','on')
                    sym = phase.PrimitiveEntry('base_symmetrization_for_GAMMA','on')
                    blksamp.add_entry(sampm)
                    blksamp.add_entry(red)
                    blksamp.add_entry(sym)
                target_input.get_root_entry().get_entry('accuracy').\
                             add_entry(blksamp)

    def do_scf_convergence(self,target_input):
        scfconv=None
        succession = None
        if self.input_type=='static' or self.input_type=='dos':
            scfconv = properties.Property().get_property('geninp.scf_criteria.static',float)
        elif self.input_type=='stropt' or self.input_type=='neb':
            scfconv = properties.Property().get_property('geninp.scf_criteria.stropt',float)
            succ = properties.Property().get_property('geninp.scf_succession.stropt',int)
            succession = phase.PrimitiveEntry('succession',succ)
        elif self.input_type=='md_nvt' or self.input_type=='md_nve' or self.input_type=='md':
            scfconv = properties.Property().get_property('geninp.scf_criteria.md',float)
            succ = properties.Property().get_property('geninp.scf_succession.md',int)
            succession = phase.PrimitiveEntry('succession',succ)
        elif self.input_type=='phonon':
            scfconv = properties.Property().get_property('geninp.scf_criteria.phonon',float)
        if scfconv is not None:
            criteriablock = phase.BlockEntry('scf_convergence')
            scfconventry = phase.PrimitiveEntry('delta_total_energy',scfconv,'hartree')
            criteriablock.add_entry(scfconventry)
            if succession is not None:
                criteriablock.add_entry(succession)
            target_input.get_root_entry().get_entry('accuracy').\
                         add_entry(criteriablock)

    def do_initial_wf_and_chg(self,target_input):
        target_input.get_root_entry().get_entry('accuracy').add_entry(\
        phase.PrimitiveEntry('initial_wavefunctions','atomic_orbitals'))

        target_input.get_root_entry().get_entry('accuracy').add_entry(\
        phase.PrimitiveEntry('initial_charge_density','atomic_charge_density'))


    def do_cutoff(self,target_input):
        cutoff_wf = properties.Property().get_property("geninp.cutoff",float)
        factacc = 1.0
        if self.accuracy==2:
            factacc = 1.44
        elif self.accuracy==0:
            factacc = 0.6
        cutoff_wf *= factacc
        if cutoff_wf is not None:
            wfentry = phase.PrimitiveEntry('cutoff_wf',cutoff_wf,uni='hartree')
            target_input.get_root_entry().get_entry('accuracy').\
                         add_entry(wfentry)
            factor = 4.0
            if phase_util.contains_uspp(target_input):
                fac = properties.Property().get_property('geninp.aug_factor',float)
                if fac is not None:
                    factor = fac
            cdentry = phase.PrimitiveEntry('cutoff_cd',cutoff_wf*factor,uni='hartree')
            target_input.get_root_entry().get_entry('accuracy').\
                         add_entry(cdentry)

    def do_wf_solver(self,target_input):
        wftable = target_input.get_root_entry().get_entry('wavefunction_solver.solvers.table')
        wfs=[]
        if self.conver>=2:
            nb=properties.Property().get_property('geninp.david_lmm',int)
            nid=2
            if nb is None:
                nb=256
            if target_input.get_num_valence_electrons()<nb:
                wfs.append({'no':1,'sol':'msd','till_n':1,'dts':'*','dte':'*','itr':'*','var':'linear',\
                   'prec':False,'cmix':'1','submat':False})
                wfs.append({'no':2,'sol':'davidson','till_n':2,'dts':'*','dte':'*','itr':'*','var':'linear',\
                   'prec':False,'cmix':'1','submat':False})
                nid=3
            else:
                wfs.append({'no':1,'sol':'lm+msd','till_n':1,'dts':'*','dte':'*','itr':'*','var':'linear',\
                   'prec':True,'cmix':'1','submat':True})
            wfs.append({'no':nid,'sol':'rmm3','till_n':-1,'dts':'*','dte':'*','itr':'*','var':'linear',\
                   'prec':True,'cmix':'1','submat':True})
        elif self.conver==1:
            wfs.append({'no':1,'sol':'lm+msd','till_n':-1,'dts':'*','dte':'*','itr':'*','var':'linear',\
                   'prec':True,'cmix':'1','submat':True})
        elif self.conver==0:
            wfs.append({'no':1,'sol':'davidson','till_n':-1,'dts':'*','dte':'*','itr':'*','var':'linear',\
                   'prec':False,'cmix':'1','submat':False})
     
        for wf in wfs:
            wftable.add_row_data(wf)

        if self.conver>=2:
            edelta = 1e-2
            if self.conver==3:
                edelta=5e-3
            if self.conver==2:
                edelta=1e-4
            rmmblock = phase.BlockEntry('rmm')
            rmmblock.add_entry(phase.PrimitiveEntry('edelta_change_to_rmm',float(edelta),'hartree'))
            target_input.get_root_entry().get_entry('wavefunction_solver').add_entry(rmmblock)

            submatblock = phase.BlockEntry('submat')
            submatblock.add_entry(phase.PrimitiveEntry('before_renewal','on'))
            target_input.get_root_entry().get_entry('wavefunction_solver').add_entry(submatblock)

    def do_cd_mixing(self,target_input):
        mix = 0.4
        nbmix=20 
        istr=3 
        if self.conver==0:
            mix = 0.2
            nbmix=5
            istr=5
        elif self.conver==1 or self.conver==2:
            mix = 0.3
            nbmix=5
            istr=3
        elif self.conver==3:
            mix = 0.4
            nbmix=15
            istr=3
        method = 'broyden2'
        if self.spin and self.conver<2:
            method = 'simple'
            mix=0.2
        elif self.spin and self.conver>=2:
            mix=0.1
        cmixtable = target_input.get_root_entry().get_entry('charge_mixing.mixing_methods.table')
        cmix0={'no':1,'method':method,'rmxs':mix,'rmxe':mix,'itr':40,'var':'linear','prec':True,\
                 'nbmix':nbmix,'istr':istr,'update':'renew'}
        cmixtable.add_row_data(cmix0)
        if self.spin and self.conver>=2:
            target_input.get_root_entry().get_entry('charge_mixing').\
            add_entry(phase.PrimitiveEntry('sw_recomposing','on'))
            target_input.get_root_entry().get_entry('charge_mixing').\
            add_entry(phase.PrimitiveEntry('spin_density_mixfactor',4))

    def do_str_evol(self,target_input):
        taufac = properties.Property().get_property("geninp.tau_mass_factor",float)
        if taufac is None:
            tafac=8.0
        ndiv = properties.Property().get_property("geninp.tau_ndiv",float)
        if ndiv is None:
            ndiv=20.0
        if self.input_type=='stropt':
            forc = properties.Property().get_property('geninp.forcmx',float)
            if forc is not None:
                forceconvblock=phase.BlockEntry('force_convergence')
                forconv = phase.PrimitiveEntry('max_force',forc,'hartree/bohr')
                forceconvblock.add_entry(forconv)
                target_input.get_root_entry().get_entry('accuracy').\
                             add_entry(forceconvblock)
            strmethod = properties.Property().get_property('geninp.stropt.method')
            if strmethod is not None:
                strevlblock = phase.BlockEntry('structure_evolution')
                methodent = phase.PrimitiveEntry('method',strmethod)
                strevlblock.add_entry(methodent)

                dt = config.get_dt_estimate(self.atomic_coordinates,taufac,ndiv)
                dtentry = phase.PrimitiveEntry('dt',dt,'au_time')
                strevlblock.add_entry(dtentry)
                if strmethod=='gdiis' or strmethod=='bfgs':
                    gdiisblock = phase.BlockEntry('gdiis')
                    inimethod = properties.Property().get_property('geninp.stropt.gdiis.initial')
                    if inimethod is not None:
                        gdiisblock.add_entry(phase.PrimitiveEntry('initial_method',inimethod))
                    forc2 = properties.Property().get_property('geninp.stropt.gdiis.forc2gdiis',float)
                    if forc2 is not None:
                        gdiisblock.add_entry(phase.PrimitiveEntry('c_forc2gdiis',forc2,'hartree/bohr'))
                    if gdiisblock.entry_count()!=0:
                        strevlblock.add_entry(gdiisblock)
                target_input.get_root_entry().add_entry(strevlblock)

        elif self.input_type=='phonon':
            phblock = phase.BlockEntry('phonon')
            phblock.add_entry(phase.PrimitiveEntry('sw_phonon','on'))
            phblock.add_entry(phase.PrimitiveEntry('sw_vibrational_modes','on'))
            phblock.add_entry(phase.PrimitiveEntry('sw_calc_force','on'))
            phblock.add_entry(phase.PrimitiveEntry('displacement',0.1,'bohr'))
            target_input.get_root_entry().add_entry(phblock)
            if not self.symmetry:
                target_input.get_root_entry().get_entry('structure').\
                add_entry(phase.BlockEntry('symmetry'))
                target_input.get_root_entry().get_entry('structure.symmetry').\
                add_entry(phase.PrimitiveEntry('method','automatic'))
            
        elif self.input_type=='md_nvt' or self.input_type=='md_nve' or self.input_type=='md':
            strevlblock = target_input.get_root_entry().get_entry('structure_evolution')
            dt = config.get_dt_estimate(self.atomic_coordinates,taufac,ndiv)
            if dt is None:
                return
            dtentry = phase.PrimitiveEntry('dt',dt,'au_time')
            strevlblock.add_entry(dtentry)

            method = 'temperature_control'
            if self.input_type=='md_nve':
                method = 'velocity_verlet'
            strevlblock.add_entry(phase.PrimitiveEntry('method',method))
    
            taufac = properties.Property().get_property("geninp.tau_mass_factor",float)
            if taufac is None:
                return
            nat = self.atomic_coordinates.get_natm2()
            qmass = config.get_qmass_estimate(\
            self.atomic_coordinates,taufac,nat,self.temperature)
            tempcontblock = phase.BlockEntry('temperature_control')
            tempcontblock.add_entry(phase.PrimitiveEntry('set_initial_velocity','on'))
            strevlblock.add_entry(tempcontblock)
            thermoblock = phase.BlockEntry('thermostat')
            tempcontblock.add_entry(thermoblock)
            thermotable = phase.TableEntry()
            thermoblock.add_entry(thermotable)
            row_data = {'no':'1','qmass':qmass,'temp':self.temperature}
            thermotable.add_row_data(row_data)

            target_input.get_root_entry().get_entry('printoutlevel').add_entry(\
                                  phase.PrimitiveEntry('iprivelocity',2))

            if self.input_type=='md_nvt':
                atmtable = target_input.get_root_entry().get_entry(\
                           'structure.atom_list.atoms.table')
                atmtable.add_identifier('thermo_group')
                nat = atmtable.get_num_data()
                for i in range(nat):
                    atmtable.set_data(i,'thermo_group','1')

    def do_postprocessing(self,target_input):
        if self.input_type=='dos':
            bldos = phase.BlockEntry('dos')
            target_input.get_root_entry().get_entry('postprocessing').add_entry(bldos)
            swdos = phase.PrimitiveEntry('sw_dos','on')
            method = phase.PrimitiveEntry('method','tetrahedral')
            bldos.add_entry(swdos)
            bldos.add_entry(method)

    def do_neb(self,target_input):
        if self.input_type != 'neb':
            return

        taufac = properties.Property().get_property("geninp.tau_mass_factor",float)
        ndiv = properties.Property().get_property("geninp.tau_ndiv",float)
        if taufac is None or ndiv is None:
            return

        thres=properties.Property().get_property('geninp.neb.thres',float)
        if thres is None:
            thres = 0.001
        spring=properties.Property().get_property('geninp.neb.spring',float)
        if spring is None:
            spring  = 0.3

        method='quench'
        dt = config.get_dt_estimate(self.atomic_coordinates,taufac,ndiv)*0.5

        target_input.get_root_entry().get_entry('control').add_entry(\
                     phase.PrimitiveEntry('multiple_replica_mode','on'))  
        target_input.get_root_entry().get_entry('control').add_entry(\
                     phase.PrimitiveEntry('multiple_replica_max_iteration',5000))
        
        multiblock = phase.BlockEntry('multiple_replica')
        target_input.get_root_entry().add_entry(multiblock)

        accuracyblock = phase.BlockEntry('accuracy')
        multiblock.add_entry(accuracyblock)
        accuracyblock.add_entry(phase.PrimitiveEntry('dt',dt))
        accuracyblock.add_entry(phase.PrimitiveEntry('neb_time_integral',method))
        accuracyblock.add_entry(phase.PrimitiveEntry('neb_convergence_condition',3))
        accuracyblock.add_entry(phase.PrimitiveEntry('neb_convergence_threshold',thres))
        
        constraintblock = phase.BlockEntry('constraint')        
        multiblock.add_entry(constraintblock)
        constraintblock.add_entry(phase.PrimitiveEntry('sp_k_init',spring))

        structureblock = phase.BlockEntry('structure')
        multiblock.add_entry(structureblock)
        structureblock.add_entry(phase.PrimitiveEntry('number_of_replicas',self.nreplica))
        structureblock.add_entry(phase.PrimitiveEntry('endpoint_images','directin'))

        replicablock = phase.BlockEntry('replicas')
        structureblock.add_entry(replicablock)
        replicatable = phase.TableEntry()
        replicablock.add_entry(replicatable)

        for i in range(self.nreplica):
            row={'replica_number':str(i+1),'howtogive_coordinates':'proportional'}
            replicatable.add_row_data(row)

        atom0 = phase.BlockEntry('atom_list_end0')
        structureblock.add_entry(atom0)
        atom0table = phase.TableEntry()
        atom0.add_entry(atom0table)
        self.fill_neb_atom_list(self.atomic_coordinate_end0,atom0table)

        atom1 = phase.BlockEntry('atom_list_end1')
        structureblock.add_entry(atom1)
        atom1table = phase.TableEntry()
        atom1.add_entry(atom1table)
        self.fill_neb_atom_list(self.atomic_coordinate_end1,atom1table)

    def fill_neb_atom_list(self,coordinates,atmtable):
        if coordinates is None:
            logger.warn('invalid atomic coordinates')
            return

        for coord in coordinates:
            elem = coord.get_element_name()
            pos = coord.get_pos(mode=config.FRACTIONAL)
            row = {'element':elem,'rx':pos[0],'ry':pos[1],'rz':pos[2]}
            atmtable.add_row_data(row)
            
    def run(self):
        self.generate_input()

if __name__ == "__main__":
    piou.run_main(InputGenerator)

