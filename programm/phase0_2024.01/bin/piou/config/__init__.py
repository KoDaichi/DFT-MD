''' module which defines constants, functions and classes for 
    manipulating atomic configuration data '''

import logging
import os
import math

from piou.data import constants
from piou.data import elements

from piou.util import pyutil
from piou.util import vector_op
from piou.data import properties

logger = logging.getLogger('piou.config')

FRACTIONAL = 0
CARTESIAN  = 1

BOHR     = 0
ANGSTROM = 1
NM       = 2

def resolve_frame_no(no_specified,numframe):
    if no_specified is not None and no_specified>=0 and no_specified<numframe-1:
        return no_specified
    return numframe-1

def pack(pos):
    ''' pack the specified position within the unit cell. fractional coordinates are assumed. '''
    ret = [0.0,0.0,0.0]
    for i in range(len(pos)):
        ret[i] = pos[i]-math.floor(pos[i])
    return ret

def get_length_unit_conversion_factor(from_unit,to_unit,invert=False):
    ''' get the length conversion factor. from_unit and to_unit should be 
    one of config.BOHR, config.ANGSTROM or config.NM. set invert to True
    if you want the inverted value. '''
    if from_unit is None or to_unit is None or from_unit==to_unit:
        return 1.0

    convfactor = 1.0
    if to_unit == BOHR and from_unit==ANGSTROM:
        convfactor = constants.Angstrom_2_Bohr
    elif to_unit == BOHR and from_unit == NM:
        convfactor = constants.nm_2_Angstrom * constants.Angstrom_2_Bohr
    elif to_unit == ANGSTROM and from_unit == BOHR:
        convfactor = constants.Bohr_2_Angstrom
    elif to_unit == ANGSTROM and from_unit == NM:
        convfactor = constants.nm_2_Angstrom
    elif to_unit == NM and from_unit == ANGSTROM:
        convfactor = constants.Angstrom_2_nm
    elif to_unit == NM and from_unit == BOHR:
        convfactor = constants.Bohr_2_Angstrom * constants.Angstrom_2_nm

    if invert:
        convfactor = 1.0/convfactor

    return convfactor

def get_reference_kmesh(atomic_conf,dk):
    ''' get the reference k-mesh. will return None if the unitcell is not set properly. '''
    if not atomic_conf.is_cell_ready():
        return (None,None,None)

    unitcell = atomic_conf.get_unit_cell()
    reclat = unitcell.get_reciprocal_lattice(to_unit=BOHR)
    if reclat is None:
        return (None,None,None)

    (b1,b2,b3) = (vector_op.norm(reclat[0]),\
                  vector_op.norm(reclat[1]),\
                  vector_op.norm(reclat[2]))
    (rb1,rb2,rb3) = (int(round(b1/dk)),\
                     int(round(b2/dk)),\
                     int(round(b3/dk)))
    return (rb1,rb2,rb3)

def get_dt_estimate(atomic_conf,taufac,ndiv):
    ''' estimate dt from the lightest atom defined, and return it.
        will return None if any errors are encountered. '''
    if ndiv == 0:
        return None
    tau = get_tau(atomic_conf,taufac)
    if tau==None:
        return None
    return tau/ndiv

def get_qmass_estimate(atomic_conf,taufac,nat,target_temperature):
    ''' estimate the mass of the Nose thermostat, and return it in atomic units. '''
    tau = get_tau(atomic_conf,taufac)
    return 2.0*3.0*nat*constants.kbt*target_temperature*(0.5*tau/math.pi)**2

def get_tau(atomic_conf,taufac):
    ''' get the estimate for the characteristic time of vibration for the specified 
    atomic configuration in atomic units '''
    masslist = []
    for atom in atomic_conf:
        if not atom.is_mobile():
            continue
        try:
            masslist.append(atomic_conf.\
            element_info_dict[atom.get_element_name()]['mass'])
        except:
            pass
        
    if len(masslist)==0:
        return None

    masslist.sort()
    logger.debug("mass of the lightest atom : "+str(masslist[0]))
    logger.debug("mass of the heaviest atom : "+str(masslist[len(masslist)-1]))
    return taufac * math.sqrt(masslist[0])

def get_fract_coord(atom_config, cartx, carty, cartz):
    unit_cell = atom_config.get_unit_cell()
    if unit_cell is None:
        return None
    invmat = unit_cell.get_inverse_lattice_vector()
    if invmat is None:
        return
    rx = invmat[0][0]*cartx+invmat[0][1]*carty+invmat[0][2]*cartz
    ry = invmat[1][0]*cartx+invmat[1][1]*carty+invmat[1][2]*cartz
    rz = invmat[2][0]*cartx+invmat[2][1]*carty+invmat[2][2]*cartz
    return (rx, ry, rz)

class AtomConfigGenerator:
    ''' interface class which converts itself to an instance of the AtomConfig class '''

    def __init__(self):
        self.subtype = None
        self.frame = None
        self.valid=True

    def get_name(self):
        ''' return the 'name' for this atomic coniguration generator '''
        return 'AtomConfigGeneator'

    def get_all_subtypes(self,mode='w'):
        ''' return the possible subtypes as a list. if subtype support is unavailable, will return None '''
        return None

    def get_defaultname(self,arg=None):
        ''' return the default name '''
        return None

    def supports_import(self):
        ''' check if this generator supports import '''
        return True

    def supports_export(self):
        ''' check if this generator supports export '''
        return True

    def check_support(self):
        ''' implement this method if this generator has some kind of dependency. '''
        return True

    def select_subtype(self,mode=None):
        ''' if 'subtype' should be enabled, implement this method so that the subtype can be specified by the user. 
        'mode' will be either 'r' or 'w' '''
        return None

    def set_subtype(self,subtype):
        ''' set the 'subtype' for this atomic configuration '''
        self.subtype = subtype

    def get_subtype(self):
        ''' get the 'subtype' for this atomic configuration '''
        return self.subtype

    def get_atomic_configuration(self,nframe=-1):
        ''' get the atomic configuration associated with this data structure.
            returns the last frame if a negative value is specified for 
            the 'nframe' argument '''
        #if self.frame is None:
        self.frame = self._gen_atomic_configuration()

        if self.frame is None:
            return None

        index = nframe
        if index<0:
            index = len(self.frame)-1
        if index<0: 
            return None
        return self.frame[index]

    def get_num_frames(self):
        ''' get the number of frames defined '''
        if self.frame is None:
            self.frame = self._gen_atomic_configuration()
        return len(self.frame)

    def supports_frame(self):
        ''' override and return True if this coord-generator creates frame data, or in other words,
            get_num_frames() can be greater than 1.'''
        return False

    def supports_volumetric_data(self):
        ''' override and return True if this coord-generator supports volumetric data, ie Gaussian Cube format '''
        return False

    def get_frame(self):
        ''' get reference to the list of all the atomic configuration defined in this class '''
        if self.frame is None:
            self.frame = self._gen_atomic_configuration()
        return self.frame

    def _gen_atomic_configuration(self):
        ''' subclass must implement this method. should return the list of an instance of 
            piou.data.config.AtomConfig'''
        raise NotImplementedError()

    def reset_atomic_configuration(self):
        self.frame = None

    def export_atomic_configuration(self,atomic_coordinates,to_file=None,frame_no=None,all_frames=False):
        ''' generate an instance of the current class from the atomic configuration data
            specified by the arg. subclass must implement this method. '''
        raise NotImplementedError()

    def is_valid(self):
        ''' return false if this 'inputgenerator' is invalid for some reason. '''
        return self.valid

    def set_is_valid(self,val):
        ''' set whether this atomic configuration is valid '''
        self.valid = val

class AtomConfig:
    ''' the class which holds atomic configuration data 
        since this class implements the iterator protocal, 
        it's possible to do something like:
        for atom in atom_config:
            ... do something to atom ...  
        where atom_config is an instance of this class, and atom is
        an instance of Atom '''
       


    def __init__(self,coordinate_system=FRACTIONAL, length_unit=BOHR, \
                 unit_cell=None,a_vector=None,b_vector=None,c_vector=None,\
                 a=None,b=None,c=None,alpha=None,beta=None,gamma=None):
        ''' you can specify the 'coordinate_system', length_unit and the
            unit cell from this constructor. '''
        self.element_info_dict={}
        self.coordinate_system = FRACTIONAL
        self.length_unit = BOHR
        self.atom_list = []
        self.unit_cell = None
        self.volumetric_data = None
        self.total_energy=None
        self.index=-1
        self.coordinate_system = coordinate_system
        self.length_unit = length_unit
        self.set_unit_cell(cell=unit_cell,a_vector=a_vector,b_vector=b_vector,c_vector=c_vector,\
                           a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,msg=False)
        if self.unit_cell is not None and unit_cell is None:
            self.unit_cell.set_length_unit(length_unit)

    def pack(self,abc=None):
        ''' pack all atoms into the unit cell. '''
        for atom in self:
            p = pack(atom.get_pos(mode=FRACTIONAL))
            if abc is None:
                pp = p
            else:
                pp = atom.get_pos(mode=FRACTIONAL)
                for i in range(3):
                    if abc[i]:
                        pp[i] = p[i]
            atom.set_pos(pp,mode=FRACTIONAL)

    def set_coordinate_system(self,coordsys):
        ''' set the coordinate system for this atomic configuration '''
        self.coordinate_system = coordsys

    def get_coordinate_system(self):
        ''' get the coordinate system for this atomic configuration. will return either FRACTIONAL or CARTESIAN'''
        return self.coordinate_system

    def add_volumetric_data(self,voldata):
        ''' set the specified volumetric data. the volumetric data must be an instance of VolumetricData.'''
        if not isinstance(voldata,VolumetricData):
            logger.error('volumetric data must be an instance of VolumetricData')
            return
        if self.volumetric_data is None:
            self.volumetric_data=[]
        self.volumetric_data.append(voldata)

    def remove_volumetric_data(self,voldata=None,index=0):
        ''' remove the specified volumetric from list '''
        if self.volumetric_data is None:
            return
        if voldata is not None and voldata in self.volumetric_data:
            self.volumetric_data.remove(voldata)
        else:
            if index>=0 and index<len(self.volumetric_data):
                self.volumetric_data.pop(index)

    def get_volumetric_data(self):
        ''' get the list of volumetric data associated with this atomic coordinates.
        will return None if no volumetric data are associated to this atomic coordinates'''
        return self.volumetric_data

    def get_natm2(self):
        ''' obtain the number of atoms with the weight parameter taken into account '''
        nat = 0
        for atom in self:
            aa = 1
            if atom.get_weight()==2:
                aa=2
            nat += aa
        return nat

    def update_cart(self):
        natm = self.get_num_atom()
        for i in range(natm):
            atom = self.get_atom_at(i)
            pos = atom.get_pos(mode=FRACTIONAL)
            atom.set_pos(pos=[pos[0],pos[1],pos[2]],mode=FRACTIONAL)

    def update_fract(self):
        natm = self.get_num_atom()
        for i in range(natm):
            atom = self.get_atom_at(i)
            pos = atom.get_pos(mode=CARTESIAN)
            atom.set_pos(pos=[pos[0],pos[1],pos[2]],mode=CARTESIAN)

    def to_super(self,na,nb,nc):
        natmorg = self.get_num_atom()
        orgcell = self.get_unit_cell().get_lattice_vector()
        uni = self.get_unit_cell().get_length_unit()
        avec=[orgcell[0][0]*na,orgcell[0][1]*na,orgcell[0][2]*na]
        bvec=[orgcell[1][0]*nb,orgcell[1][1]*nb,orgcell[1][2]*nb]
        cvec=[orgcell[2][0]*nc,orgcell[2][1]*nc,orgcell[2][2]*nc]
        self.set_unit_cell(a_vector=avec,b_vector=bvec,c_vector=cvec)
        self.get_unit_cell().set_length_unit(uni)
        for i in range(natmorg):
            atom = self.get_atom_at(i)
            pos = atom.get_pos(mode=CARTESIAN)
            atom.set_pos(pos=[pos[0],pos[1],pos[2]],mode=CARTESIAN)

        for i in range(na):
            for j in range(nb):
                for k in range(nc):
                    if i==0 and j==0 and k==0:
                        continue
                    trans=[0,0,0]
                    for l in range(3):
                        trans[l] = i*orgcell[0][l]+j*orgcell[1][l]+k*orgcell[2][l]
                    for l in range(natmorg):
                        atom = self.get_atom_at(l)
                        cpos = atom.get_pos(CARTESIAN)
                        newpos = [cpos[0]+trans[0],cpos[1]+trans[1],cpos[2]+trans[2]]
                        newatom = Atom(self,element=atom.get_element_name(),rx=newpos[0],ry=newpos[1],rz=newpos[2], \
                                  mobile=atom.is_mobile(),weight=atom.get_weight(),mode=CARTESIAN)
                        self.add_atom(atom=newatom)

    def get_center_of_mass(self,mode=FRACTIONAL,to_unit=BOHR):
        ''' get the center of mass for this atomic configuration '''
        M=0.0
        pospos=[0.0,0.0,0.0]
        for atom in self:
            pos=atom.get_pos(mode=mode,to_unit=to_unit)
            mass = atom.get_mass()
            M += mass
            for i in range(len(pos)):
                pospos[i] += pos[i]*mass
        pospos=[pospos[0]/M,pospos[1]/M,pospos[2]/M]

        return pospos

    def shift(self,dshift,mode=FRACTIONAL,unit=BOHR):
        ''' shift all the coordinates of the current atomic configuration '''
        if len(dshift)<3:
            logger.error("dshift must have atleast three elements")
            return

        for atom in self:
            pos = atom.get_pos(mode=mode,to_unit=unit)
            for i in range(len(pos)):
                pos[i] += dshift[i]
            atom.set_pos(pos,mode=mode,theunit=unit)

    def get_bond_angle_list(self,cutoff_length_sq=None):
        ''' get a list of bond angle. each element in the list is a four-value tuple,
            where the first element is the bond angle, while remaining elements are the reference to 
            the corresponding atoms.  will return None if a valid unit cell is not defined. '''
        unitcell = self.get_unit_cell()
        if unitcell is None:
            return None

        [a_vector,b_vector,c_vector] = unitcell.get_lattice_vector(to_unit=BOHR)
        if a_vector is None or b_vector is None or c_vector is None:
            return None

        ret = []
        for i in range(len(self.atom_list)):
            for j in range(len(self.atom_list)):
                if i>=j:
                    continue
                for k in range(len(self.atom_list)):
                    if j>=k:
                        continue
                    atomi=self.atom_list[i]
                    atomj=self.atom_list[j]
                    atomk=self.atom_list[k]
                    posi=atomi.get_pos(mode=FRACTIONAL,to_unit=BOHR)
                    posj=atomj.get_pos(mode=FRACTIONAL,to_unit=BOHR)
                    posk=atomk.get_pos(mode=FRACTIONAL,to_unit=BOHR)
                    diffij = [posi[0]-posj[0],posi[1]-posj[1],posi[2]-posj[2]]
                    diffik = [posi[0]-posk[0],posi[1]-posk[1],posi[2]-posk[2]]
                    for l in range(3):
                        diffij[l] -= math.floor(diffij[l])
                        diffik[l] -= math.floor(diffik[l])
                    for l in range(3):
                        if math.fabs(diffij[l])>0.5:
                            diffij[l] -= 1
                        elif diffij[l]<-0.5:
                            diffij[l] += 1
                        if math.fabs(diffik[l])>0.5:
                            diffik[l] -= 1
                        elif diffik[l]<-0.5:
                            diffik[l] += 1
                    xij = a_vector[0]*diffij[0]+b_vector[0]*diffij[1]+c_vector[0]*diffij[2]
                    yij = a_vector[1]*diffij[0]+b_vector[1]*diffij[1]+c_vector[1]*diffij[2]
                    zij = a_vector[2]*diffij[0]+b_vector[2]*diffij[1]+c_vector[2]*diffij[2]

                    xik = a_vector[0]*diffik[0]+b_vector[0]*diffik[1]+c_vector[0]*diffik[2]
                    yik = a_vector[1]*diffik[0]+b_vector[1]*diffik[1]+c_vector[1]*diffik[2]
                    zik = a_vector[2]*diffik[0]+b_vector[2]*diffik[1]+c_vector[2]*diffik[2]

                    r2ij = xij**2+yij**2+zij**2
                    r2ik = xik**2+yik**2+zik**2

                    if r2ij<cutoff_length_sq and r2ik<cutoff_length_sq:
                        rij=math.sqrt(r2ij)
                        rik=math.sqrt(r2ik)
                        dp = vector_op.dot_product((xij,yij,zij),(xik,yik,zik))
                        cosin = dp/(rij*rik)
                        if cosin>1:
                            cosin=1
                        elif cosin<-1:
                            cosin=-1
                        ret.append((math.acos(cosin),atomi,atomj,atomk))
        return ret
        
    def get_coordination_number(self,atom,cutoff_length_sq,alist=None,pbc=True):
        unitcell = self.get_unit_cell()
        if unitcell is None:
            return None

        [a_vector,b_vector,c_vector] = unitcell.get_lattice_vector(to_unit=BOHR)
        if a_vector is None or b_vector is None or c_vector is None:
            return None

        thealist = self.atom_list
        if alist is not None:
            thealist = alist
        ret = 0
        for atm in thealist:
            if atm==atom:
                continue
            if atm.is_dangling() or atom.is_dangling():
                continue
            posi = atom.get_pos(mode=FRACTIONAL,to_unit=BOHR)
            posj = atm.get_pos(mode=FRACTIONAL,to_unit=BOHR)
            posdiff = [posi[0]-posj[0],posi[1]-posj[1],posi[2]-posj[2]]
            if pbc:
                for k in range(3):
                    posdiff[k] -= math.floor(posdiff[k])
                    for k in range(3):
                        if math.fabs(posdiff[k])>0.5:
                            posdiff[k] -= 1
                        elif posdiff[k]<-0.5:
                            posdiff[k] += 1

            x = a_vector[0]*posdiff[0]+b_vector[0]*posdiff[1]+c_vector[0]*posdiff[2]
            y = a_vector[1]*posdiff[0]+b_vector[1]*posdiff[1]+c_vector[1]*posdiff[2]
            z = a_vector[2]*posdiff[0]+b_vector[2]*posdiff[1]+c_vector[2]*posdiff[2]

            r2 = x**2+y**2+z**2
            if r2<cutoff_length_sq:
                ret += 1
        return ret

    covrads = None
    def get_bondlenth_list(self,cutoff_factor=1.0,cutoff_length_sq=None,pbc=True,tlist=None,register_neighbor=False,alist=None,alist2=None):
        ''' get a list of bond length. each element in the list is a three-value tuple,
            where the first element is the bond length, while the second and third
            elements are the reference to the corresponding atoms.  will return None
            if a valid unit cell is not defined. '''
        unitcell = self.get_unit_cell()
        if unitcell is None:
            return None

        [a_vector,b_vector,c_vector] = unitcell.get_lattice_vector(to_unit=BOHR)
        if a_vector is None or b_vector is None or c_vector is None:
            return None

        build_cutoff = cutoff_length_sq is None
        if self.covrads is None:
            self.covrads = {}
            for atom in self:
                covrad = elements.ElementInfo().get_element_attribute(\
                         atom.get_element_name(),'covalent_radius')
                if covrad is None:
                    continue
                self.covrads[atom.get_element_name()] = covrad

        target_atoms = self.atom_list
        target_atoms2 = self.atom_list
        if alist is not None:
            target_atoms = alist
        if alist2 is not None:
            target_atoms2 = alist2
        ret = []
        for i in range(len(target_atoms)):
            for j in range(len(target_atoms2)):
                if i>=j and alist is None:
                    continue

                if tlist is not None and not (i in tlist and j in tlist):
                    continue

                atomi = target_atoms[i]
                atomj = target_atoms2[j]

                if not (atomi.get_element_name() in self.covrads or \
                        atomj.get_element_name() in self.covrads):
                    continue

                posi = atomi.get_pos(mode=FRACTIONAL,to_unit=BOHR)
                posj = atomj.get_pos(mode=FRACTIONAL,to_unit=BOHR)
                posdiff = [posi[0]-posj[0],posi[1]-posj[1],posi[2]-posj[2]]

                if build_cutoff:
                    covi = self.covrads[atomi.get_element_name()]
                    covj = self.covrads[atomj.get_element_name()]
                    cutoff_length_sq = ((covi+covj)*cutoff_factor)**2

                if pbc:
                    for k in range(3):
                        posdiff[k] -= math.floor(posdiff[k])
                        for k in range(3):
                            if math.fabs(posdiff[k])>0.5:
                                posdiff[k] -= 1
                            elif posdiff[k]<-0.5:
                                posdiff[k] += 1

                x = a_vector[0]*posdiff[0]+b_vector[0]*posdiff[1]+c_vector[0]*posdiff[2]
                y = a_vector[1]*posdiff[0]+b_vector[1]*posdiff[1]+c_vector[1]*posdiff[2]
                z = a_vector[2]*posdiff[0]+b_vector[2]*posdiff[1]+c_vector[2]*posdiff[2]

                r2 = x**2+y**2+z**2
                if r2<cutoff_length_sq:
                    bl = math.sqrt(r2)
                    if register_neighbor:
                        atomi.add_neighbor(atomj,bl)
                        if alist is None:
                            atomj.add_neighbor(atomi,bl)
                    ret.append( (bl,atomi,atomj) )
        return ret

    def get_length_unit(self):
        ''' get the length unit associated with this atomic configuration '''
        return self.length_unit

    def set_length_unit(self,un):
        ''' set the length unit associted with this atomic configuration '''
        self.length_unit = un

    def add_atom(self,element=None,rx=None,ry=None,rz=None,fx=None,fy=None,fz=None,mobile=False,weight=1,attr=None,atom=None,idi=-1,coordinate_system=None):
        ''' add atom with the specified attributes to list along with the reference to myself. '''
        if atom is not None:
            self.atom_list.append(atom)
            atom.atom_config=self
            return
        sys = self.coordinate_system
        if coordinate_system is not None:
            sys = coordinate_system
        atom = Atom(self,element,rx,ry,rz,mobile,weight,sys)
        atom.fx = fx
        atom.fy = fy
        atom.fz = fz
        atom.set_id(idi)
        self.element_info_dict[element] = elements.ElementInfo().get_element_info(element)
        self.atom_list.append(atom)
        if attr is not None:
            for k,v in attr.items():
                atom.set_attribute(k,v)

    def remove_atom(self,atom):
        ''' remove the specified atom from list (if it exists) '''
        #if atom in self.atom_list:
        #    self.atom_list.remove(atom)
        try:
            self.atom_list.remove(atom)
        except ValueError:
            pass

    def get_unit_cell(self):
        ''' get the instance of the unit cell '''
        return self.unit_cell

    def bond_list(self,rcut=None,pbc=True,register_neighbor=False):
        neihborind=(-1,0,1)
        self.generate_cell_index(rcut)
        self.generate_atom_cell_index()
        rc2 = rcut*rcut
        [a_vector,b_vector,c_vector] = self.get_unit_cell().get_lattice_vector(to_unit=BOHR)
        ret = []
        natm=self.get_num_atom()
        if self.covrads is None:
            self.covrads = {}
            for atom in self:
                covrad = elements.ElementInfo().get_element_attribute(\
                         atom.get_element_name(),'covalent_radius')
                if covrad is None:
                    continue
                self.covrads[atom.get_element_name()] = covrad

        for i0 in range(self.n0):
            for i1 in range(self.n1):
                for i2 in range(self.n2):
                    atom_in_cell = self.atom_cell_index[i0][i1][i2]
                    for atom in atom_in_cell:
                        aind=atom.index()
                        for in0 in neihborind:
                            nn0 = i0+in0
                            if nn0>=self.n0:
                                nn0=0
                            elif nn0<0:
                                nn0=self.n0-1
                            for in1 in neihborind:
                                nn1 = i1+in1
                                if nn1>=self.n1:
                                    nn1=0
                                elif nn1<0:
                                    nn1=self.n1-1
                                for in2 in neihborind:
                                    nn2 = i2+in2
                                    if nn2>=self.n2:
                                        nn2=0
                                    elif nn2<0:
                                        nn2=self.n2-1
                                    atom_in_cell2 = self.atom_cell_index[nn0][nn1][nn2]
                                    for atom2 in atom_in_cell2:
                                        aind2=atom2.index()
                                        if aind==aind2 and in0==0 and in1==0 and in2==0:
                                            continue
                                        if aind2>aind:
                                            continue
                                        if rcut is None:
                                            covi = self.covrads[atomi.get_element_name()]
                                            covj = self.covrads[atomj.get_element_name()]
                                            rc2 = ((covi+covj)*cutoff_factor)**2
                                        posi = atom.get_pos(mode=FRACTIONAL,to_unit=BOHR)
                                        posj = atom2.get_pos(mode=FRACTIONAL,to_unit=BOHR)
                                        posdiff = [posi[0]-posj[0],posi[1]-posj[1],posi[2]-posj[2]]
                                        if pbc:
                                            for k in range(3):
                                                posdiff[k] -= math.floor(posdiff[k])
                                                for k in range(3):
                                                    if posdiff[k]>0.5:
                                                        posdiff[k] -= 1
                                                    elif posdiff[k]<-0.5:
                                                        posdiff[k] += 1
                                        x = a_vector[0]*posdiff[0]+b_vector[0]*posdiff[1]+c_vector[0]*posdiff[2]
                                        y = a_vector[1]*posdiff[0]+b_vector[1]*posdiff[1]+c_vector[1]*posdiff[2]
                                        z = a_vector[2]*posdiff[0]+b_vector[2]*posdiff[1]+c_vector[2]*posdiff[2]
                                        r2 = x**2+y**2+z**2
                                        if r2<rc2:
                                            bl = math.sqrt(r2)
                                            if register_neighbor:
                                                atom.add_neighbor(atom2,bl)
                                                atom2.add_neighbor(atom,bl)
                                            ret.append( (bl,atom,atom2) )
        return ret

    atom_cell_index={}
    def generate_atom_cell_index(self):
        self.atom_cell_index={}
        for i0 in range(self.n0):
            self.atom_cell_index[i0] = {}
            for i1 in range(self.n1):
                self.atom_cell_index[i0][i1] = {}
                for i2 in range(self.n2):
                    self.atom_cell_index[i0][i1][i2] = []
        for atom in self:
            [i0,i1,i2] = self.where_am_i(atom)
            self.atom_cell_index[i0][i1][i2].append(atom)
            
    def where_am_i(self,atom):
        fpos = atom.get_pos(mode=FRACTIONAL)
        fposp = pack(fpos)
        nn0 = int(fposp[0]*self.n0)
        nn1 = int(fposp[1]*self.n1)
        nn2 = int(fposp[2]*self.n2)
        if nn0==self.n0:
            nn0 = self.n0-1
        if nn1==self.n1:
            nn1 = self.n1-1
        if nn2==self.n2:
            nn2 = self.n2-1
        if nn0<0:
            nn0 = 0
        if nn1<0:
            nn1 = 0
        if nn2<0:
            nn2 = 0
        return [nn0,nn1,nn2]

    n0=1
    n1=1
    n2=1
    def generate_cell_index(self,rcut):
        [a,b,c,alph,beta,gamma] = self.get_unit_cell().get_lattice_constant(to_unit=BOHR)
        self.n0 = int(a/rcut)
        self.n1 = int(b/rcut)
        self.n2 = int(c/rcut)
    
    def is_cell_ready(self):
        ''' determine whether the unit cell is ready or not '''
        if self.unit_cell is None:
            return False
        return self.unit_cell.ready()

    def set_unit_cell(self,cell=None,\
                      a_vector=None,b_vector=None,c_vector=None,\
                      a=None,b=None,c=None,alpha=None,beta=None,gamma=None,msg=True):
        ''' set the unit cell associated with this atomic configuration.
            This can be done by one of three ways :
            1. instantiate the UnitCell class and specify it  
            2. specify the three lattice vectors. each vector is a list holding exactly three floats.
            3. specify the lattice constants directly 
            if msg is True, then an error message will be printed on unit cell generation failure. '''
        if cell is not None:
            self.unit_cell = cell
        elif a_vector is not None and b_vector is not None and c_vector is not None:
            self.unit_cell = UnitCell(a_vector=a_vector,b_vector=b_vector,c_vector=c_vector)
        elif a is not None and b is not None and c is not None and \
             alpha is not None and beta is not None and gamma is not None:
            self.unit_cell = UnitCell(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)
        else:
            if msg:
                logger.error("you must specify one of : "+\
                             " 1. an instance of UnitCell"+\
                             " 2. the 3 lattice vectors"+\
                             " 3. the 6 lattice constants")

    def get_atom_list(self,elems=None):
        ''' get the list of atoms. if elems is specified, it will be used as an element filter. '''
        if elems is None:
            return self.atom_list
        else:
            ret=[]
            for elem in elems:
                for atom in self:
                    if atom.get_element_name()==elem:
                        ret.append(atom)
            return ret

    def get_atom_at(self,index):
        ''' get the atom at the specified index. will return None if index is out of range.'''
        if index>=len(self) or index<0:
            return None
        return self.atom_list[index]
    
    def get_num_atom(self):
        ''' get the number of atoms included in the current atomic configuration '''
        return len(self.atom_list)

    def set_total_energy(self,en):
        ''' set the total energy associated with this atomic configuration. '''
        self.total_energy = en

    def get_total_energy(self):
        ''' obtain the total energy associated with this atomic configuration.
            if the total energy is unset, will return None. '''
        return self.total_energy

    def __str__(self):
        ret = "\n"+"status of the atomic configuration\n"
        ret += "coordinate system : "
        if self.coordinate_system == FRACTIONAL:
            ret += "fractional"
        else:
            ret += "cartesian"
        ret += "\n"

        ret += "length unit : "
        if self.length_unit == BOHR:
            ret += "bohr"
        elif self.length_unit == ANGSTROM:
            ret += "angstrom"
        else:
            ret += "nm"
        ret += "\n"

        if self.is_cell_ready():
            ret += str(self.unit_cell)
        ret += "atoms\n"
        ret += "  number of atoms defined : "+str(len(self))+"\n"
        for atom in self:
            ret += "  "+str(atom)+"\n"
        return ret

    def __iter__(self):
        ''' implementation of the iterator protocal '''
        return self

    def __len__(self):
        ''' returns the number of atoms '''
        return len(self.atom_list)

    def next(self):
        ''' implementation of the iterator protocal '''
        if self.index==len(self.atom_list)-1:
            self.index = -1
            raise StopIteration
        self.index += 1 
        return self.atom_list[self.index]

    def __next__(self):
        ''' implementation of the iterator protocal '''
        if self.index==len(self.atom_list)-1:
            self.index = -1
            raise StopIteration
        self.index += 1 
        return self.atom_list[self.index]

    def has_inversion_symmetry(self):
        if not self.is_cell_ready():
            return False

        eps = constants.very_small_value*100
        flg=[]
        for atom in self.atom_list:
            flg.append(False)

        for i in range(len(self.atom_list)):
            atom = self.atom_list[i]
            if flg[i]:
                continue
            if atom.get_weight()==2:
                flg[i]=True
                continue
            pos = pack(atom.get_pos(mode=FRACTIONAL))
            invpos=pack([-pos[0],-pos[1],-pos[2]])
            for j in range(len(self.atom_list)):
                if i>j or flg[j]:
                    continue
                atomj=self.atom_list[j]
                posj=pack(atomj.get_pos(mode=FRACTIONAL))
                diff=vector_op.diff(invpos,posj)
                if math.fabs(diff[0])<eps and math.fabs(diff[1])<eps and\
                    math.fabs(diff[2])<eps and\
                    atom.get_element_name()==atomj.get_element_name():
                    flg[i]=True
                    flg[j]=True
                    break
        return not (False in flg)

    def get_elem_list(self):
        ''' obtain a list of elements '''
        ret = []
        for atom in self:
            if not atom.get_element_name() in ret:
                ret.append(atom.get_element_name())
        return ret

    def get_elem_typ_map(self):
        ''' obtain a dictionary which maps the element name to an integer ID 
            ie map['Si']=1, map['O']=2, ... '''
        ret = {}
        the_id=0
        for atom in self:
            if not atom.get_element_name() in ret:
                the_id+=1
                ret[atom.get_element_name()] = the_id
        return ret

class Atom:
    ''' class which holds attributes for each atom '''

    def add_neighbor(self,atom,dist):
        ''' add neighbor atoms. 
            arg : atom the atom to be regarded as neighbor. 
            arg : dist the distance between myself and atom. '''
        if self.neighbor is None:
            self.neighbor = []
        self.neighbor.append((atom,dist))

    def set_id(self,idi):
        self.idi = idi 
    
    def get_neighbor(self):
        ''' return a list of neighbors for this atom. 
            element of the neighbor list is a tuple of the form (atom obj, distance)'''
        if self.neighbor is None:
            self.neighbor = []
        return self.neighbor

    def __init__(self,atom_config,element,rx,ry,rz,mobile=False,weight=1,mode=FRACTIONAL,mass=None):
        ''' specify the element name and the coordinates.  '''
        self.atom_config = atom_config
        self.element=''
        self.rx=0.0
        self.ry=0.0
        self.rz=0.0
        self.fx=0.0
        self.fy=0.0
        self.fz=0.0
        self.rx_cart=0.0
        self.ry_cart=0.0
        self.rz_cart=0.0
        self.mobile=False
        self.weight = 1

        self.mass = 50000
        self.neighbor = None
        self.valid_atm = False
        self.idi = -1
        self.atom_attr= {
         'num_layer':None
        ,'aldos':None
        ,'thermo_group':None
        ,'proj_group':None
    
        ,'forcex':None
        ,'forcey':None
        ,'forcez':None

        ,'velocx':None
        ,'velocy':None
        ,'velocz':None
        ,'id':None
        }

        self.atom_attr_type= {
         'num_layer':'int'
        ,'aldos':'int'
        ,'thermo_group':'int'
        ,'proj_group':'int'

        ,'forcex':'float'
        ,'forcey':'float'
        ,'forcez':'float'

        ,'velocx':'float'
        ,'velocy':'float'
        ,'velocz':'float'
        ,'id':'int'
        }
        try:
            rx = float(str(rx).lower().replace('d','e'))
            ry = float(str(ry).lower().replace('d','e'))
            rz = float(str(rz).lower().replace('d','e'))
        except ValueError:
            logger.error("the coordinates are non-float : "+str(rx)+', '+str(ry)+' '+str(rz))
            return

        self.set_element_name(element)

        self.set_pos(pos=[rx,ry,rz],mode=mode)
    
        self.valid_atm = True
        if mobile:
            self.mobile = True
        else:
            self.mobile = False

        if weight==2:
            self.weight = 2
        else:
            self.weight = 1

        if mass is not None:
            self.mass =mass

    def is_valid(self):
        return self.valid_atm

    def get_mass(self):
        ''' get the mass associated with this atom '''
        return self.mass

    def index_str(self):
        ind = self.index()
        return str(ind).rjust(5)

    def index(self):
        ''' return the index of the current atom in the parent list '''
        return self.atom_config.atom_list.index(self)
    
    def is_dangling(self):
        return self.atom_config is None

    def deregister(self):
        self.atom_config = None

    def detach(self):
        ''' remove the current atom from configuration '''
        if self.atom_config is not None:
            self.atom_config.remove_atom(self)
            self.atom_config=None

    def _gen_cart_pos(self):
        unit_cell = self.atom_config.get_unit_cell()
        if unit_cell is None:
            return

        [a_vector,b_vector,c_vector] = unit_cell.get_lattice_vector()
        if a_vector is None or b_vector is None or c_vector is None:
            return

        factor = 1.0
        if unit_cell.get_length_unit()!=self.atom_config.length_unit:
            factor = get_length_unit_conversion_factor(\
            from_unit=unit_cell.get_length_unit(),to_unit=self.atom_config.length_unit,invert=False)

        self.rx_cart = \
        (self.rx * a_vector[0] + self.ry * b_vector[0] + self.rz * c_vector[0])*factor

        self.ry_cart = \
        (self.rx * a_vector[1] + self.ry * b_vector[1] + self.rz * c_vector[1])*factor

        self.rz_cart = \
        (self.rx * a_vector[2] + self.ry * b_vector[2] + self.rz * c_vector[2])*factor

        eps = constants.very_small_value
        if math.fabs(self.rx_cart)<eps:
            self.rx_cart=0
        if math.fabs(self.ry_cart)<eps:
            self.ry_cart=0
        if math.fabs(self.rz_cart)<eps:
            self.rz_cart=0

    def _gen_fract_pos(self):
        unit_cell = self.atom_config.get_unit_cell()
        if unit_cell is None:
            return

        factor = 1.0
        if unit_cell.get_length_unit()!=self.atom_config.length_unit:
            factor = get_length_unit_conversion_factor(\
            from_unit=unit_cell.get_length_unit(),to_unit=self.atom_config.length_unit,invert=True)

        invmat = unit_cell.get_inverse_lattice_vector()
        if invmat is None:
            return

        self.rx = \
        (invmat[0][0]*self.rx_cart+invmat[0][1]*self.ry_cart+invmat[0][2]*self.rz_cart)*factor

        self.ry = \
        (invmat[1][0]*self.rx_cart+invmat[1][1]*self.ry_cart+invmat[1][2]*self.rz_cart)*factor

        self.rz = \
        (invmat[2][0]*self.rx_cart+invmat[2][1]*self.ry_cart+invmat[2][2]*self.rz_cart)*factor

        eps = constants.very_small_value
        if math.fabs(self.rx)<eps:
            self.rx=0
        if math.fabs(self.ry)<eps:
            self.ry=0
        if math.fabs(self.rz)<eps:
            self.rz=0

    def set_attribute(self,attrname,attrval):
        ''' set the specified attribute to the specified value. 
            Will do nothing if the specified attribute is undefined. '''
        if attrname in self.atom_attr:
            typ = self.atom_attr_type[attrname]
            try:
                if typ=='string':
                    self.atom_attr[attrname] = str(attrval)
                elif typ=='float':
                    self.atom_attr[attrname] = float(str(attrval).lower().replace('d','e'))
                elif typ=='int':
                    self.atom_attr[attrname] = int(attrval)
                elif typ=='bool':
                    self.atom_attr[attrname] = pyutil.parse_bool(attrval)
                else:
                    logger.warn("could not find type : "+typ+". assuming 'string'")
                    self.atom_attr[attrname] = str(attrval)
            except ValueError:
                logger.error("'type' is not consistent : ["+str(attrname)+" = "+str(attrval)+"]")
        else:
            logger.error("could not find attribute : "+attrname)

    def set_element_name(self,elem):
        ''' set the element name. '''
        #if not elements.ElementInfo().element_exists(elem):
        #    logger.warn("unknown element : "+str(elem))
        self.element = elem

    def get_element_name(self):
        ''' get the element name of this atom '''
        return self.element

    def set_pos(self,pos,mode=FRACTIONAL,theunit=None):
        ''' set the position for this atom. if len(pos)!=0 or any one of 
            the coordinates cannot be transformed to float, will do nothing. 
            'theunit' makes sense only for mode=CARTESIAN '''
        eps = constants.very_small_value

        if len(pos)!=3:
            return
        if mode==FRACTIONAL:
            self.rx = pos[0]
            self.ry = pos[1]
            self.rz = pos[2]
            if math.fabs(self.rx)<eps:
                self.rx=0
            if math.fabs(self.ry)<eps:
                self.ry=0
            if math.fabs(self.rz)<eps:
                self.rz=0
            self._gen_cart_pos()
        else:
            unit = self.atom_config.length_unit
            if theunit is not None and theunit != unit:
                unit=theunit
            factor = get_length_unit_conversion_factor(unit,self.atom_config.length_unit)
            self.rx_cart = pos[0]*factor
            self.ry_cart = pos[1]*factor
            self.rz_cart = pos[2]*factor
            if math.fabs(self.rx_cart)<eps:
                self.rx_cart=0
            if math.fabs(self.ry_cart)<eps:
                self.ry_cart=0
            if math.fabs(self.rz_cart)<eps:
                self.rz_cart=0
            self._gen_fract_pos()

    def get_pos(self,mode=FRACTIONAL,to_unit=None):
        ''' get the position of this atom as a python list. 
            to_unit will be ignored if mode is equal to FRACTIONAL'''
        if mode==FRACTIONAL:
            return [self.rx, self.ry, self.rz]
        else:
            if to_unit==None or to_unit == self.atom_config.length_unit:
                return [self.rx_cart,self.ry_cart,self.rz_cart]
            factor = get_length_unit_conversion_factor(self.atom_config.length_unit,to_unit)
            return [self.rx_cart*factor, self.ry_cart*factor, self.rz_cart*factor]

    def set_mobile(self,mob):
        ''' setter for the mobile attribute. '''
        self.mobile = mob

    def is_mobile(self):
        ''' returns True if this atom is 'mobile', False otherwise '''
        return self.mobile

    def get_weight(self):
        ''' get the 'weight' attribute associated with this atom. '''
        return self.weight

    def set_weight(self,wei):
        ''' set the 'weight' attribute associated with this atom. '''
        if wei==1 or wei==2:
            self.weight=wei
        else:
            logger.error("'weight' must either be 1 or 2.")

    def get_attribute(self,attrname):
        ''' get the specified attribute associated with 
            this atom by a two-value tuple. 
            the first element is the value, while the
            second element is the associated type. 
            Will returns None if no such attr exist, 
            or if the corresponding attribute is None.'''
        if not attrname in list(self.atom_attr.keys()):
            return None
        elif self.atom_attr[attrname] is None:
            return None
        return (self.atom_attr[attrname],self.atom_attr_type[attrname])

    def __str__(self):
        ret  = str("element : "+self.element+",").ljust(14)
        ret += str(" pos (fract) : ("+str(self.rx)+\
                         "  "+str(self.ry)+\
                         "  "+str(self.rz)+")").ljust(58)
        ret += str(" pos (cart) : ("+str(self.rx_cart)+\
                         "  "+str(self.ry_cart)+\
                         "  "+str(self.rz_cart)+")").ljust(62)
        ret += str(", mobile : "+str(self.mobile)+",").ljust(18)
        ret += " weight : "+str(self.weight)
        return ret

class UnitCell:
    ''' the class which represents the unit cell '''

    def __ne__(self,other):
        return not self.__eq__(other)

    def __eq__(self,other):
        ''' will return True if all the elements of the lattice vector is smaller than 
        constants.very_small_value.'''
        if not isinstance(other,UnitCell):
            return False
        (avec0,bvec0,cvec0) = self.get_lattice_vector()
        if avec0 is None or bvec0 is None or cvec0 is None:
            return False
        (avec1,bvec1,cvec1) = other.get_lattice_vector()
        if avec1 is None or bvec1 is None or cvec1 is None:
            return False
        eps = constants.very_small_value
        for ia in range(3):
            if math.fabs(avec0[ia]-avec1[ia])>eps:
                return False
        for ib in range(3):
            if math.fabs(bvec0[ib]-bvec1[ib])>eps:
                return False
        for ic in range(3):
            if math.fabs(cvec0[ic]-cvec1[ic])>eps:
                return False
        return True

    def __init__(self,a_vector=None,b_vector=None,c_vector=None,\
                 a=None,b=None,c=None,alpha=None,beta=None,gamma=None):
        ''' initialize the unit cell by either the cell vector or the lattice constants '''
        self.inver=None
        self.reclat = None
        self.a_vector=[0,0,0]
        self.b_vector=[0,0,0]
        self.c_vector=[0,0,0]
        self.a=0
        self.b=0
        self.c=0
        self.alpha=0
        self.beta=0
        self.gamma=0
        self.length_unit = BOHR
        self.set_cell(a_vector=a_vector,b_vector=b_vector,c_vector=c_vector,\
                      a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)

    def ready(self):
        ''' determine whether this unit cell is 'ready' or not '''
        [a,b,c] = self.get_lattice_vector()
        if a is None or b is None or c is None:
            return False

        latconst = self.get_lattice_constant()
        if len(latconst)<6:
            return False

        v = self.get_volume()
        if v < 1e-10:
            return False
        return True

    def set_length_unit(self,length_unit):
        ''' set the length unit for the unit cell '''
        self.length_unit = length_unit

    def get_length_unit(self):
        ''' obtain the length unit for this unit cell'''
        return self.length_unit

    def set_cell(self,a_vector=None,b_vector=None,c_vector=None,\
                 a=None,b=None,c=None,alpha=None,beta=None,gamma=None):
        ''' set the unit cell, either by cell vector (a three-value list of floats)
            or by the lattice constant. if a, b and c vectors are all non-None
            the former is assumed, while the latter will be assumed if a,b,c,alpha,gamma,beta
            are all non-None. will do nothing if either do not apply. '''
        if a_vector is not None and b_vector is not None and c_vector is not None:
            if len(a_vector)<3 or len(b_vector)<3 or len(c_vector)<3:
                logger.error('cell vector must have at least three elements')
                return
            try:
                for i in range(3):
                    a_vector[i] = float(str(a_vector[i]).lower().replace('d','e'))
                    b_vector[i] = float(str(b_vector[i]).lower().replace('d','e'))
                    c_vector[i] = float(str(c_vector[i]).lower().replace('d','e'))
            except ValueError:
                logger.error('cell vector must be a float')
                return
            self.a_vector = a_vector
            self.b_vector = b_vector
            self.c_vector = c_vector
            self._gen_lattice_constant()
            self._gen_inverse()
            self._gen_reciprocal_lattice()
        elif a is not None and b is not None and c is not None and \
             alpha is not None and beta is not None and gamma is not None:
            try:
                a = float(str(a).lower().replace('d','e'))
                b = float(str(b).lower().replace('d','e'))
                c = float(str(c).lower().replace('d','e'))
                alpha  = float(str(alpha).lower().replace('d','e'))
                beta   = float(str(beta).lower().replace('d','e'))
                gamma  = float(str(gamma).lower().replace('d','e'))
            except ValueError:
                logger.error('lattice constants must be a float')
                return
            self.a = a
            self.b = b
            self.c = c
            self.alpha = alpha
            self.beta = beta
            self.gamma = gamma
            self._gen_lattice_vector()
            self._gen_inverse()
            self._gen_reciprocal_lattice()

    def set_lattice_vector(self,avec,bvec,cvec):
        ''' set the lattice vector. '''
        self.a_vector = avec
        self.b_vector = bvec
        self.c_vector = cvec
        self._gen_lattice_constant()
        self._gen_inverse()
        self._gen_reciprocal_lattice()

    def set_lattice_constant(self,a,b,c,alpha,beta,gamma):
        ''' set the lattice constant '''
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self._gen_lattice_vector()
        self._gen_inverse()
        self._gen_reciprocal_lattice()

    def _gen_inverse(self):
        mat = []
        for i in range(3):
            mat.append([self.a_vector[i],self.b_vector[i],self.c_vector[i]])
        self.inver = vector_op.get_inverse_matrix(mat)

    def _gen_lattice_constant(self):
        if self.a_vector is None or self.b_vector is None or self.c_vector is None:
            return

        self.a = vector_op.norm(self.a_vector)
        self.b = vector_op.norm(self.b_vector)
        self.c = vector_op.norm(self.c_vector)
        if self.a==0 or self.b==0 or self.c==0:
            logger.error("failed to calculate the length of cell.")
            return

        cosbeta = vector_op.dot_product(self.c_vector,self.a_vector)/(self.a*self.c)
        cosalph = vector_op.dot_product(self.c_vector,self.b_vector)/(self.c*self.b)
        cosgamm = vector_op.dot_product(self.a_vector,self.b_vector)/(self.a*self.b)

        self.beta  = math.acos(cosbeta) * constants.rad_2_angle
        self.alpha = math.acos(cosalph) * constants.rad_2_angle
        self.gamma = math.acos(cosgamm) * constants.rad_2_angle

    def _gen_lattice_vector(self):
        if self.a is None or self.b is None or self.c is None or \
           self.alpha is None or self.beta is None  or self.gamma is None:
            return

        eps = constants.very_small_value

        self.a_vector = [0.0,0.0,0.0]
        self.b_vector = [0.0,0.0,0.0]
        self.c_vector = [0.0,0.0,0.0]

        alpha = self.alpha*constants.angle_2_rad
        beta  = self.beta*constants.angle_2_rad
        gamma = self.gamma*constants.angle_2_rad

        self.a_vector[0] = 1.0
        self.b_vector[0] = math.cos(gamma)
        self.c_vector[0] = math.cos(beta)

        self.b_vector[1] = math.sin(gamma)
        temp = (math.cos(alpha) - math.cos(beta) * math.cos(gamma))/math.sin(gamma)
        self.c_vector[1] = temp
        self.c_vector[2] = math.sqrt(math.sin(self.beta*constants.angle_2_rad)**2 - \
                                     temp * temp)
        for i in range(3):
            self.a_vector[i] *= self.a
            self.b_vector[i] *= self.b
            self.c_vector[i] *= self.c
            if math.fabs(self.a_vector[i])<eps:
                self.a_vector[i] = 0
            if math.fabs(self.b_vector[i])<eps:
                self.b_vector[i] = 0
            if math.fabs(self.c_vector[i])<eps:
                self.c_vector[i] = 0

    def get_inverse_lattice_vector(self,to_unit=None):
        ''' get the inverse of the lattice vector (remark: not the reciprocal lattice) '''
        if to_unit is None or to_unit == self.length_unit:
            return self.inver
        factor = get_length_unit_conversion_factor(self.length_unit,to_unit,invert=True)
        ret = []
        for i in range(3):
            row = []
            for j in range(3):
                row.append(self.inver[i][j]*factor)
            ret.append(row)
        return ret

    def get_lattice_vector(self,to_unit=None):
        ''' get the lattice vector expressed by a 3x3 list of floats '''
        if to_unit is None or to_unit==self.length_unit:
            return [self.a_vector,self.b_vector,self.c_vector]
        factor = get_length_unit_conversion_factor(self.length_unit,to_unit)
        orig = [self.a_vector,self.b_vector,self.c_vector]
        ret = []
        for i in range(3):
            row = []
            for j in range(3):
                row.append(orig[i][j]*factor)
            ret.append(row)
        return ret

    def get_lattice_constant(self,to_unit=None):
        ''' get the lattice constant. the ordering of the list is :
            [a,b,c,alpha,beta,gamma]'''
        if to_unit is None or to_unit == self.length_unit:
            return [self.a,self.b,self.c,self.alpha,self.beta,self.gamma]
        factor = get_length_unit_conversion_factor(self.length_unit,to_unit)
        return [self.a*factor,self.b*factor,self.c*factor,self.alpha,self.beta,self.gamma]

    def get_reciprocal_lattice(self,to_unit=None):
        ''' get the reciprocal lattice expressed by a 3x3 list of floats. '''
        if to_unit is None or to_unit==self.length_unit:
            return self.reclat
        if self.reclat is None:
            return None

        factor = get_length_unit_conversion_factor(self.length_unit,to_unit,invert=True)
        ret = []
        for i in range(3):
            row = []
            for j in range(3):
                row.append(self.reclat[i][j]*factor)
            ret.append(row)
        return ret

    def _gen_reciprocal_lattice(self):
        vol = self.get_volume()
        if vol is None:
            logger.error("failed to calculate the volume of cell.")
            return None

        factor = 2.0 * math.pi/vol
        b_cross_c = vector_op.cross_product(self.b_vector,self.c_vector,factor)
        c_cross_a = vector_op.cross_product(self.c_vector,self.a_vector,factor)
        a_cross_b = vector_op.cross_product(self.a_vector,self.b_vector,factor)

        self.reclat = [b_cross_c,c_cross_a,a_cross_b]

    def get_volume(self,to_unit=None):
        ''' get the volume of the unit cell. Will return None if any errors are encountered. '''
        if self.a_vector is None or self.b_vector is None or self.c_vector is None:
            if a is None or b is None or c is None or \
               alpha is None or beta is None or gamma is None:
                logger.error("invalid cell specification. cannot continue...")
                return None 
            self._gen_lattice_vector()
        
        factor = get_length_unit_conversion_factor(self.length_unit,to_unit)
        avec = [self.a_vector[0]*factor,self.a_vector[1]*factor,self.a_vector[2]*factor]
        bvec = [self.b_vector[0]*factor,self.b_vector[1]*factor,self.b_vector[2]*factor]
        cvec = [self.c_vector[0]*factor,self.c_vector[1]*factor,self.c_vector[2]*factor]

        b_cross_c=vector_op.cross_product(bvec,cvec)
        if b_cross_c is not None:
            return math.fabs(vector_op.dot_product(avec,b_cross_c))
        else:
            return None

    def __str__(self):
        ret = "unit cell\n"

        ret += "  length unit : "
        if self.get_length_unit()==BOHR:
            ret += "bohr"
        elif self.get_length_unit()==ANGSTROM:
            ret += "angstrom"
        elif self.get_length_unit()==NM:
            ret += "nm"
        ret += "\n"
       
        if self.a_vector is not None and self.b_vector is not None and self.c_vector is not None:
            ret +=  "  cell vector\n"
            ret +=  "    a-vector : ("\
            +str(self.a_vector[0])+" "+str(self.a_vector[1])+" "+str(self.a_vector[2])+")\n"
            ret +=  "    b-vector : ("\
            +str(self.b_vector[0])+" "+str(self.b_vector[1])+" "+str(self.b_vector[2])+")\n"
            ret +=  "    c-vector : ("\
            +str(self.c_vector[0])+" "+str(self.c_vector[1])+" "+str(self.c_vector[2])+")\n"
        if self.a is not None and self.b is not None and self.c is not None and \
           self.alpha is not None and self.beta is not None and self.gamma is not None:
            ret +=  "  lattice constant\n"
            ret +=  "    a     : "+str(self.a)+"\n"
            ret +=  "    b     : "+str(self.b)+"\n"
            ret +=  "    c     : "+str(self.c)+"\n"
            ret +=  "    alpha : "+str(self.alpha)+"\n"
            ret +=  "    beta  : "+str(self.beta)+"\n"
            ret +=  "    gamma : "+str(self.gamma)+"\n"
        reclat = self.get_reciprocal_lattice()
        if reclat is not None:
            ret += "  reciprocal lattice\n"
            ret += "    b1  : ("+str(reclat[0][0])+" "+str(reclat[0][1])+" "+str(reclat[0][2])+")\n"
            ret += "    b2  : ("+str(reclat[1][0])+" "+str(reclat[1][1])+" "+str(reclat[1][2])+")\n"
            ret += "    b3  : ("+str(reclat[2][0])+" "+str(reclat[2][1])+" "+str(reclat[2][2])+")\n"
            ret += "   |b1| : "+str(vector_op.norm(reclat[0]))+"\n"
            ret += "   |b2| : "+str(vector_op.norm(reclat[1]))+"\n"
            ret += "   |b3| : "+str(vector_op.norm(reclat[2]))+"\n"
        vol = self.get_volume()
        if vol is not None:
            ret += "  cell volume : "+str(vol)+"\n"

        return ret
 
class VolumetricData:
    ''' the class which holds volumetric data, ie, charge density. '''
    def __init__(self,unitcell=None,\
        a=None,b=None,c=None,alpha=None,beta=None,gamma=None,\
        a_vector=None,b_vector=None,c_vector=None,
        n1=None,n2=None,n3=None,vdata=None):
        self.unitcell=None
        self.ndiv=[]
        self.div=[]
        self.vdata=[]
        self.origin=[0.0,0.0,0.0]
        if unitcell is not None:
            self.unitcell = unitcell
        else:
            self.unitcell = UnitCell(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
            a_vector=a_vector,b_vector=b_vector,c_vector=c_vector)
        if n1 is not None and n2 is not None and n3 is not None:
            self.ndiv.append(n1)
            self.ndiv.append(n2)
            self.ndiv.append(n3)
            self._init_div()
        if vdata is not None:
            self.vdata = vdata
        else:
            self._init_vdata()

    def set_origin(self,origin):
        ''' set the amount of the 'origin' for this volumetric data '''
        self.origin = origin

    def get_origin(self):
        ''' get the amount of the 'origin' for this volumetric data '''
        return self.origin

    def _init_vdata(self):
        if len(self.ndiv)==3 and len(self.div)==3:
            self.vdata = []
            ndata = self.ndiv[0]*self.ndiv[1]*self.ndiv[2]
            for i in range(ndata):
                self.vdata.append(0.0)

    def _init_div(self):
        if len(self.ndiv)==3 and self.unitcell is not None and self.unitcell.ready():
            (avec,bvec,cvec) = self.unitcell.get_lattice_vector()
            self.div=[]

            firstvec=[]
            secondvec=[]
            thirdvec=[]
            for i in range(3):
                firstvec.append(  avec[i]/float(self.ndiv[i]))
                secondvec.append( bvec[i]/float(self.ndiv[i]))
                thirdvec.append(  cvec[i]/float(self.ndiv[i]))
            self.div.append(firstvec)
            self.div.append(secondvec)
            self.div.append(thirdvec)

    def has_enough_data(self):
        ''' check whether the volumetric data specified has enough data in it'''
        if len(self.ndiv)<3:
            return False
        if self.vdata is None:
            return False
        if len(self.div)<3:
            self._init_div()
        ndata = self.ndiv[0]*self.ndiv[1]*self.ndiv[2]
        return len(self.vdata)==ndata

    def set_vdata(self,vdata):
        ''' set the list of volumetric data '''
        self.vdata = vdata

    def set_vdata_at(self,n1,n2,n3,vdatum):
        ''' set the data at the specified position. note that 'vdatum' should be a float'''
        if not has_enough_data():
            logger.error('the volumetric data is not ready!')
            return
        self.vdata[self.get_1d_index(n1,n2,n3)] = vdatum

    def set_ndiv(self,n1,n2,n3):
        ''' set the division number of each cell '''
        ndiv=[]
        ndiv.append(n1)
        ndiv.append(n2)
        ndiv.append(n3)
        self._init_div()
        self._init_vdata()

    def get_ndiv(self):
        ''' get the division numbers of each cell '''
        return self.ndiv

    def get_div(self,to_unit=None):
        ''' get the length of each voxel. note that this is calculated from the unit cell and ndiv. 
        the length unit is the same as that for the associated unit cell. '''
        if to_unit is None or to_unit == self.unitcell.get_length_unit():
            return self.div
        factor = get_length_unit_conversion_factor(\
        from_unit=self.unitcell.get_length_unit(),to_unit=to_unit)
        tmp=[]
        for t in self.div:
            tmptmp=[]
            for tt in t:
                tmptmp.append(factor*tt)
            tmp.append(tmptmp)
        return tmp

    def set_unit_cell(self,unitcell=None,a=None,b=None,c=None,alpha=None,beta=None,gamma=None,\
        a_vector=None,b_vector=None,c_vector=None):
        ''' set the unit cell '''
        if unitcell is not None:
            self.unitcell = unitcell
        else:
            self.unitcell = UnitCell(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
            a_vector=a_vector,b_vector=b_vector,c_vector=c_vector)
        self._init_div()
        self._init_vdata()

    def get_unit_cell(self):
        ''' get the unitcell instance associated with this volumetric data '''
        return self.unitcell

    def get_vdata(self):
        ''' get the reference to the list of volumetric data '''
        return self.vdata

    def get_vdata_at(self,n1,n2,n3):
        ''' get the volumetric data at the specified position. 
        will return None if any errors are encountered. '''
        if (len(self.ndiv))<3:
            return None
        if n1>self.ndiv[0] or n2>self.ndiv[1] or n3>self.ndiv[2]:
            return None
        if not self.has_enough_data():
            return None
        return self.vdata[self.get_1d_index(n1,n2,n3)]

    def get_1d_index(self,n1,n2,n3):
        ''' convert the 3 indeces to an 1-d index.
        the index will be calculated as: (n1-1)*N2*N3+(n2-1)*N3+n3,
        which means that the order of the change for the
        3-d indeces are n3->n2->n1 (same as the Gaussian Cube format) '''
        if len(self.ndiv)<3:
            return None
        nn2=self.ndiv[1]
        nn3=self.ndiv[2]
        return int(n1*nn2*nn3+n2*nn3+n3)

    def __str__(self):
        strresult='status of the volumetric data\n'
        if self.unitcell is not None and self.unitcell.ready():
            strresult+= 'associated unit cell : '+str(self.unitcell)+'\n'
        else:
            strresult+='the unit cell is not ready\n'

        if len(self.ndiv)>=3:
            strresult += 'division for each dimension : '+\
            str(self.ndiv[0])+' '+str(self.ndiv[1])+' '+str(self.ndiv[2])+'\n'
        else:
            strresult += 'division for each dimension is undefined.\n'

        if len(self.div)>=3:
            d0 = self.div[0]
            d1 = self.div[1]
            d2 = self.div[2]
            strresult += 'vector for the first  dimension of the voxel : '+str(d0[0])+' '+str(d0[1])+' '+str(d0[2])+'\n'
            strresult += 'vector for the second dimension of the voxel : '+str(d1[0])+' '+str(d1[1])+' '+str(d1[2])+'\n'
            strresult += 'vector for the third  dimension of the voxel : '+str(d2[0])+' '+str(d2[1])+' '+str(d2[2])+'\n'
        else:
            strresult+='the vector for each dimension of the voxel is undefined.\n'

        return strresult

