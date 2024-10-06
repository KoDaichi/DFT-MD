
from piou.util import phase_util
from piou.util import vector_op
from piou.data import properties

import piou.data.elements
import piou.data.constants

from piou import config

import logging
import os
import math

logger = logging.getLogger("piou.config.functions")

def check_cpumax(args):
    fact = properties.Property().get_property("inpcheck.cpumax_warn",float)
    if fact is None:
        return True

    cpum = args['currelem']
    if cpum is None:
        return True
    c = config.phase.get_value_in_default_units(cpum)
    if c <fact:
        logger.warn("[control.cpumax] = "+str(cpum.get_value())+ " sec ... are you sure?")
        args['input'].increment_nwarn()
    return True

def check_tetrahedral_dos(args):
    logger.debug("checking if tetrahedral dos is propertly configured...")
    root = args['input'].get_root_entry()
    kmethod = root.get_entry("accuracy.ksampling.method")
    bok=True

    if kmethod is None:
        logger.warn("[accuracy.ksampling.method] must be defined and its value must be 'mesh' "+\
                    "in order to calculate the dos by the tetrahedron method")
        args['input'].increment_nwarn()
    elif kmethod.get_value() !="mesh":
            bok=False
            logger.warn("[accuracy.ksampling.method] must be 'mesh'"+\
                        " in order to calculate the dos by the tetrahedron method")
            args['input'].increment_nwarn()
    else:
        bok=True

    smear = root.get_entry("accuracy.smearing.method")
    if smear==None:
        logger.warn(\
       "[accuracy.smearing.method] must be defined and its value must be 'tetrahedral' "+\
       "in order to calculate the dos by the tetrahedral method")
        args['input'].increment_nwarn()
        bok=False
    else:
        comp = smear.get_value().strip().lower()
        if not (comp == "tetrahedral" or comp == "tetrahedron" or \
                comp == "improved_tetrahedron" or comp == "t"):
            logger.warn( \
            "[accuracy.smearing.method] must be set to 'tetrahedral' in order to calculate "+\
            "the dos by the tetrahedral method")
            args['input'].increment_nwarn()
            bok=False

    if not bok:
        logger.warn("tetrahedron method is not configured propertly; the gaussain broadening method will be used")
    else:
        logger.debug("tetrahedron method seems to be configured properly.")

    return bok

def check_cutoff_cd(args):
    cutoff_cd = args['currelem']
    cutoff_wf = args['input'].get_root_entry().get_entry('accuracy.cutoff_wf')
    wf = config.phase.get_value_in_default_units(cutoff_wf)
    bo=True
    if wf is None:
        logger.error("invalid specification : "+str(cutoff_wf))
        args['input'].increment_nerr()
        bo=False
    cd = config.phase.get_value_in_default_units(cutoff_cd)
    if cd is None:
        logger.error("invalid specification : "+str(cutoff_cd))
        args['input'].increment_nerr()
        bo=False
    if not bo:
        return False
   
    ok = cd >= 4.0 * wf

    if not ok:
        logger.error("")
        logger.error("the cutoff energy for the charge density : "+\
                      str(cutoff_cd).strip()+" [hartree] is not enough.")
        args['input'].increment_nerr()
    return ok

def check_consistency_of_wavefunction_solver(args):
    ''' check consistency of solver. for instance, 'prec' must be 'off' for the davidson solver. '''
    ok = True
    wfsolvers = args['currelem']
    rmmblock = args['input'].get_root_entry().get_entry('wavefunction_solver.rmm')
    if rmmblock == None:
        rmmblock = args['input'].get_root_entry().get_entry('wf_solver.rmm')
    nsolver = wfsolvers.get_num_data()
    minus_found=False
    maxcmix=1
    maxid=0
    for i in range(nsolver):
        if minus_found:
            logger.warn("found a solver defined after a solver whose 'till_n' attribute is smaller than 0.")
            logger.warn("this is not an error, but mind that such solver will never be used.")
            logger.warn('')
            args['input'].increment_nwarn()

        tilln = wfsolvers.get_data(i,'till_n')
        ii = i+1
        if wfsolvers.get_data(i,'id') is not None:
            ii = wfsolvers.get_data(i,'id')
        elif wfsolvers.get_data(i,'no') is not None:
            ii = wfsolvers.get_data(i,'no')
        if ii==1 and tilln is not None and tilln==1:
            mat = args['input'].get_root_entry().get_entry('accuracy.initial_wavefunctions')
            if mat is not None and mat.get_value() == 'matrix_diagon':
                logger.warn(\
                "till_n for the first solver should be greater than 1 when using 'matrix_diagon' "+\
                "for the initial WF generation")
                logger.warn('')
                args['input'].increment_nwarn()

        sol = wfsolvers.get_data(i,'sol')
        if sol!=None and sol.lower()=='davidson':
            prec = wfsolvers.get_data(i,'prec')
            if prec:
                logger.error("prec must be 'off' when using the davidson solver")
                args['input'].increment_nerr()
                ok=False
            submat = wfsolvers.get_data(i,'submat')
            if submat:
                logger.error("submat must be 'off' when using the davidson solver")
                args['input'].increment_nerr()
                ok=False
            if not ok:
                logger.error("")

        #if sol is not None and sol.lower().startswith('rmm'):
        #    if (rmmblock is None or rmmblock.get_entry('edelta_change_to_rmm') is None):
        #        logger.warn("'edelta_change_to_rmm' variable not found.")
        #        logger.warn("it is recommended to specify this variable when using the rmm solver.")
        #        logger.warn('')
        #        args['input'].increment_nwarn()
        #    submat = wfsolvers.get_data(i,'submat')
        #    if not submat:
        #        logger.warn("subspace rotation is 'off' ")
        #        logger.warn("it is strongly recommended to use subspace rotation when using the rmm solver.")
        #        logger.warn('')

        if tilln is not None and isinstance(tilln,int) and tilln<0:
            minus_found = True

        cmix=wfsolvers.get_data(i,'cmix')
        if cmix is not None and cmix!='*' and cmix>maxcmix:
            maxcmix = cmix
            maxid=i

    cmixtable = args['input'].get_root_entry().get_entry('charge_mixing.mixing_methods.table')
    if cmixtable is not None:
        n = cmixtable.get_num_data()
        nn = n
        for i in range(n):
            nnn = cmixtable.get_data(i,'no')
            if not isinstance(nnn,int):
                nnn = i
            if nnn is not None and nnn>nn:
                nn = nnn
            nnn = cmixtable.get_data(i,'id')
            if nnn is not None and nnn>nn:
                nn = nnn
        if nn<maxcmix:
            logger.error("'cmix' for solver no. "+str(maxid+1) +' is invalid : '+str(maxcmix))
            args['input'].increment_nerr()
            ok=False
    elif cmixtable is None and maxcmix >1:
        logger.error('max. cmix ('+str(maxcmix)+') is too large' )
        args['input'].increment_nerr()
        ok=False
            
    return ok

def check_charge_mixing(args):
    cmix_table = args['currelem']
    ncmix = cmix_table.get_num_data()
    ok=True
    for i in range(ncmix):
        nbmix = cmix_table.get_data(i,'nbmix')
        if nbmix is None:
            logger.error('nbmix is undefined')
            return False
        try:
            nbmix = int(nbmix)
        except ValueError:
            logger.error('invalid nbmix parameter '+str(nbmix))
            return False
        if nbmix is not None and nbmix <3:
            logger.warn("'nbmix' for charge mixing scheme no. "+str(i+1)+' is too small : '+str(nbmix))
            args['input'].increment_nwarn()
            ok=False
    return ok

def check_xctype(args):
    xcfrompp = phase_util.get_xctype_from_pp(args['input'])
    if xcfrompp is None:
        logger.error("XC type defined in the pp file is invalid")
        args['input'].increment_nerr()
        return False

    xctype = args['input'].get_root_entry().get_entry('accuracy.xctype')

    if xctype is None or xctype.get_value().lower().startswith('vdwdf'):
        return True

    if xctype.get_value().lower() != xcfrompp.lower():
        if xctype.get_value().lower()=='ggapbex' and xcfrompp.lower() == 'ggapbe':
            return True
        if 'pbe' in xctype.get_value().lower() and xcfrompp.lower() == 'ggapbe':
            return True
        logger.error("the XC type defined in F_INP ["+xctype.get_value() +\
                     "] is inconsistent with that defined in the pp file ["+xcfrompp+"]")
        args['input'].increment_nerr()
        return False

    return True

def check_ppfiles(args):
    inp = args['input']
    ret = phase_util.ppfiles_exist(inp)
    if not ret:
        inp.increment_nerr()
    return ret
    
def check_num_bands(args):
    fact = properties.Property().get_property("inpcheck.band_tolerance_factor",float)
    if fact is None:
        return True

    nbdefined = args['input'].get_root_entry().get_entry('accuracy.num_bands')
    if nbdefined is None:
        logger.error('number of bands is not defined')
        args['input'].increment_nerr()
        return False

    nv = phase_util.get_num_valence_electrons(args['input'])
    if nv is None:
        logger.error("failed to get the number of valence electrons")
        args['input'].increment_nerr()
        return False

    nbmin = nv * 0.5

    if nbdefined.get_value() < nbmin:
        logger.error("the number of bands is insufficient. defined : "+\
                     str(nbdefined.get_value())+", minimum : "+str(int(nbmin))+\
                    ", recommended : "+str(int(nbmin*fact)))
        args['input'].increment_nerr()
        return False

    if int(nbdefined.get_value()) < int(nbmin*fact):
        logger.warn("the number of bands may be insufficient. defined : "+\
                     str(nbdefined.get_value())+", minimum : "+str(int(nbmin))+\
                    ", recommended : "+str(int(nbmin*fact)))
        args['input'].increment_nwarn()

    return True

def all_elements_are_defined(args):
    ok=True

    input = args['input']
    elem_table = args['currelem']
    atom_table = args['input'].get_root_entry().get_entry('structure.atom_list.atoms.table')

    if atom_table is None:
        logger.error('atom table is undefined.')
        args['input'].increment_nerr()
        return False

    nat = atom_table.get_num_data()
    elem_dict={}
    for i in range(nat):
        elem = atom_table.get_data(i,'element')
        if not elem in elem_dict:
            elem_dict[elem] = False

    nelem = elem_table.get_num_data()
    for i in range(nelem):
        elemname = elem_table.get_data(i,'element')
        if elemname is None:
            logger.error("element name undefined")
            args['input'].increment_nerr()
            return False

        if elemname in list(elem_dict.keys()):
            elem_dict[elemname] = True

    for k,v in elem_dict.items():
        if not v and k is not None:
            logger.error("element : '"+k+"' undefined in ["+elem_table.get_long_name()+"]")
            args['input'].increment_nerr()
            ok=False

    return ok

def projector_exists(args):
    proj0 = args['input'].get_root_entry().get_entry('accuracy.hubbard.projectors.table')
    proj1 = args['input'].get_root_entry().get_entry('accuracy.projector_list.projectors.table')
    if proj0 is None or proj1 is None:
        return True
    n0 = proj0.get_num_data()
    n1 = proj1.get_num_data()
    ok=True
    for i in range(n0):
        ii = proj0.get_data(i,'no')
        if ii is None:
            continue
        found = False
        for j in range(n1):
            jj = proj1.get_data(j,'group')
            if jj is None:
                continue
            if jj==ii:
                found=True
                break
        if not found:
            logger.error("projector : "+str(ii)+" specified in accuracy.hubbard.projectors.table "+\
                         "is undefined under accuracy.hubbard.projector_list.projectors.table.")
            args['input'].increment_nerr()
            ok=False
    return ok

def valid_proj_group_assigned(args):
    projs = phase_util.get_atom_attributes(args['input'],'proj_group')
    projtable = args['input'].get_root_entry().get_entry('accuracy.projector_list.projectors.table')
    group = phase_util.get_attributes(projtable,'group')
    if projs is None or group is None:
        return True
    ok = True
    for i in range(len(projs)):
        proj = projs[i]
        if not proj in group:
            if not (proj=='*' or proj==0):
                logger.error("undefined proj_group : "+str(proj)+" assigned to atom "+str(i+1))
                args['input'].increment_nerr()
                ok=False
    return ok

def check_proj_group(args):
    projtable = args['input'].get_root_entry().get_entry('accuracy.projector_list.projectors.table')
    if projtable is None:
        logger.error("projectors undefined.")
        return False
    group = phase_util.get_attributes(projtable,'group')
    l = phase_util.get_attributes(projtable,'l')
    if group is None or l is None:
        return True
            
    for i in range(len(l)):
        for j in range(len(l)):
            if j>=i:
                continue
            if l[i]==l[j] and group[i]==group[j]:
                logger.error("the 'l' attribute for projector "+str(i+1)+" and "+str(j+1)+" are the same "+\
                "and belong to the same group.")
                args['input'].increment_nerr()
                return False
    return True

def check_num_layer(args):
    naxis = 1
    axi = args['input'].get_root_entry().get_entry('postprocessing.ldos.layerdos.normal_axis')
    if axi is not None:
        naxis = axi.get_value()

    nlayer = phase_util.get_atom_attributes(args['input'],'num_layer')
    if nlayer is None:
        return False

    xyz = 'rx'
    if naxis==2:
        xyz = 'ry'
    if naxis==3:
        xyz = 'rz'
    coord = phase_util.get_atom_attributes(args['input'],xyz)
    if coord is None:
        return False

    nlayer_coord={}
    for layer in nlayer:
        nlayer_coord[layer]=[]
    for i in range(len(nlayer)):
        nlayer_coord[nlayer[i]].append(coord[i])

    min = {}
    max = {}
    for k,v in nlayer_coord.items():
        minv=10000000
        maxv=-10000000
        for vv in v:
            if vv<minv:
                minv=vv
            if vv>maxv:
                maxv=vv
        min[k] = minv
        max[k] = maxv

    b=True
    for k,v in nlayer_coord.items():
        for k0,v0 in nlayer_coord.items():
            if k<=k0:
                continue
            if not (min[k0]>max[k] or max[k0] <min[k]):
                logger.error("inconsistent num_layer attribute : "+str(k0)+", "+str(k))
                args['input'].increment_nerr()
                b=False
    return b

def check_mass(args):
    fact = properties.Property().get_property("inpcheck.mass_tolerance_factor",float)
    if fact is None:
        return True

    elemtable = args['input'].get_root_entry().get_entry('structure.element_list.table')
    elemblock = args['input'].get_root_entry().get_entry('structure.element_list')
    if elemtable is None:
        return True
    uni = elemblock.get_default_unit('mass')
    nelem = elemtable.get_num_data()
    for i in range(nelem):
        elemname = elemtable.get_data(i,'element')
        if elemname is None:
            logger.error("element name undefined.")
            return True
        mass = elemtable.get_data(i,'mass')
        if mass is None:
            return True

        if mass is None or mass == '*':
            mass = 51577.50

        if uni == 'atomic_mass':
            mass *= piou.data.elements.ElementInfo().to_au_mass

        compmass = piou.data.elements.ElementInfo().get_element_attribute(elemname,"mass")
        if compmass is None:
            logger.error("invalid element name : "+elemname)
            args['input'].increment_nerr()
            return False
        if mass<compmass/fact:
            logger.warn("the mass "+str(mass)+" [au] of element ["+elemname+"] might be too light. "+\
                        "reference mass : "+str(compmass)+" [au].")
            args['input'].increment_nwarn()
        elif mass>compmass*fact:
            logger.warn("the mass "+str(mass)+" [au] of element ["+elemname+"] might be too heavy. "+\
                        "reference mass : "+str(compmass)+" [au].")
            args['input'].increment_nwarn()

    return True

def check_qmass(args):
    taufac = properties.Property().get_property("inpcheck.tau_mass_factor",float)
    if taufac is None:
        return True
    lmass = properties.Property().get_property("inpcheck.lightest_mass",float)
    if lmass is None:
        return True

    thermo_group = phase_util.get_atom_attributes(args['input'],'thermo_group')
    if thermo_group is None:
        logger.warn("the 'thermo_group' attribute is unassinged")
        args['input'].increment_nwarn()
        return True

    thermostats = args['currelem']
    mass = phase_util.get_mass_list(args['input'])
    wei = phase_util.get_atom_attributes(args['input'],'weight')
    nthermo = thermostats.get_num_data()
    for i in range(nthermo):
        id = thermostats.get_data(i,'id')
        if id is None:
            id = thermostats.get_data(i,'no')
        if id is None:
            id = i+1

        natm = len(thermo_group)
        mass_list_now = []
        nat_now = 0
        for j in range(natm):
            if id == thermo_group[j]:
                if wei is not None and wei==2:
                    nat_now += 2
                else:
                    nat_now += 1
                mass_list_now.append(mass[j])
        if len(mass_list_now) == 0:
            logger.error("thermo_group assignment is not valid")
            args['input'].increment_nerr()
            return False
        mass_list_now.sort()
        mass0 = mass_list_now[0]
        logger.debug("number of atoms associated with thermostat "+str(id)+" : "+str(nat_now))
        logger.debug("mass of the lightest element associated with this thermostat: "\
                      +str(mass0))
        logger.debug("mass of the heaviest element associated with this thermostat: "\
                      +str(mass_list_now[len(mass_list_now)-1]))
 
        temperature = thermostats.get_corresponding_data('id',id,'temp')
        if temperature is None:
            temperature = thermostats.get_corresponding_data('no',id,'temp')
        if temperature is None:
            temperature = thermostats.get_data(i,'temp')

        qmass = thermostats.get_corresponding_data('id',id,'qmass')
        if qmass is None:
            qmass = thermostats.get_corresponding_data('no',id,'qmass')
        if qmass is None:
            qmass = thermostats.get_data(i,'qmass')

        if temperature is None or qmass is None:
            logger.error("thermostat "+str(id)+" is invalid.")
            args['input'].increment_nerr()
            continue
        logger.debug("temperature and qmass associate with thermostat "+str(id)+" : "+\
                      str(temperature)+", "+str(qmass))
        
        tau = taufac * math.sqrt(mass0)
        logger.debug("estimated tau : "+str(tau))
#        qmassref = 2.0*3.0*nat_now*piou.data.constants.kbt*temperature*(0.5*tau/math.pi)**2
        qmassref = config.get_qmass_estimate(args['input'].get_atomic_configuration(),taufac,nat_now,temperature)
        logger.debug("reference value for qmass : "+str(qmassref))
        if qmassref>qmass:
            logger.warn("qmass "+str(qmass)+" [au] for thermostat "+str(id)+" may be too light. "+\
                        "reference value : "+str(qmassref)+" [au].")
            args['input'].increment_nwarn()
    return True

def check_dt(args):
    taufac = properties.Property().get_property("inpcheck.tau_mass_factor",float)
    if taufac is None:
        return True

    ndiv = properties.Property().get_property("inpcheck.tau_ndiv",float)
    if ndiv is None:
        return True

    dt = args['input'].get_root_entry().get_entry('structure_evolution.dt')
    if dt is None:
        return True
    fdt = config.phase.get_value_in_default_units(dt)
    aconf = args['input'].get_atomic_configuration()
    if aconf is None:
        return True
    dtref = config.get_dt_estimate(aconf,taufac,ndiv)
    if dtref is None:
        return True
    logger.debug("reference dt : "+str(dtref))
    logger.debug("defined dt   : "+str(fdt))

    if fdt>dtref*1.1:
        logger.warn("structure_evolution.dt may be too large. defined value : "+str(fdt)+\
                    " [au], reference value : "+str(dtref)+" [au].")
        args['input'].increment_nwarn()

    return True

def check_scf_convergence_criteria(args):
    root = args['input'].get_root_entry()
    deinp = config.phase.get_value_in_default_units(\
            root.get_entry("accuracy.scf_convergence.delta_total_energy"))

    if deinp is None:
        deinp = 1e-10

    de = properties.Property().get_property("inpcheck.scf_criteria.scf",float)

    mobile = phase_util.get_atom_attributes(args['input'],'mobile')
    strwarn = "[accuracy.scf_convergence.delta_total_energy] might be too large"

    if mobile is None:
        if de is not None and deinp>de:
            logger.warn(strwarn+" for a static calculation.")
            logger.warn("specified : "+str(deinp)+" [hartree], reference : "+str(de)+" [hartree].")
            args['input'].increment_nwarn()
        return True
    smob = 0
    if mobile is not None:
        for m in mobile:
            if m is None:
                logger.error('invalid value for the mobile attribute : '+str(m))
                return False
            try:
                smob += int(m)
            except ValueError:
                logger.error('invalid value for the mobile attribute : '+str(m))
                return False
    elif mobile is not None and sum(mobile)==0:
        if de is not None and deinp>de:
            logger.warn(strwarn+" for a static calculation.")
            logger.warn("specified : "+str(deinp)+" [hartree], reference : "+str(de)+" [hartree].")
            args['input'].increment_nwarn()
        return True
        
    swphonon = root.get_entry("phonon.sw_phonon")
    swcalfor = root.get_entry("phonon.sw_calc_force")
    if swphonon is not None and swcalfor is not None and \
       swphonon.get_value() and swcalfor.get_value():
        de = properties.Property().get_property("inpcheck.scf_criteria.phonon",float)
        if de is not None and deinp>de:
            logger.warn(strwarn+" for a phonon calculation.")
            logger.warn("specified : "+str(deinp)+" [hartree], reference : "+str(de)+" [hartree].")
            args['input'].increment_nwarn()
        return True

    strevlm = root.get_entry("structure_evolution.method")
    if strevlm is not None and phase_util.is_md(strevlm.get_value()):
        de = properties.Property().get_property("inpcheck.scf_criteria.md",float)
        if de is not None and deinp>de:
            logger.warn(strwarn+" for a molecular-dynamics simulation.")
            logger.warn("specified : "+str(deinp)+" [hartree], reference : "+str(de)+" [hartree].")
            args['input'].increment_nwarn()
        return True

    de = properties.Property().get_property("inpcheck.scf_criteria.stropt",float)
    if de is not None and deinp>de:
        logger.warn(strwarn+" for structural optimization.")
        logger.warn("specified : "+str(deinp)+" [hartree], reference : "+str(de)+" [hartree].")
        args['input'].increment_nwarn()
    return True

def check_mobile(args):
    mobile = phase_util.get_atom_attributes(args['input'],'mobile')
    method = args['input'].get_root_entry().get_entry('structure_evolution.method').get_value()
    if mobile is None:
        logger.error('mobile is None')
        return False
    if None in mobile: 
        return False
    if mobile is None or sum(mobile)==0:
        if phase_util.is_md(method):
            logger.warn("the 'mobile' attribute for all atoms are set to 'off', "+\
                        "which is unlikely for a molecular-dynamics simulation.")
            args['input'].increment_nwarn()
        else:
            logger.warn("the 'mobile' attribute for all atoms are set to 'off', "+\
                        "which is unlikely for structural optimization.")
            args['input'].increment_nwarn()
    return True

def check_k_mesh(args):
    method = args['input'].get_root_entry().get_entry('accuracy.ksampling.method')
    kmesh_dk = properties.Property().get_property("inpcheck.kmesh_dk",float)
    if kmesh_dk is None:
        logger.error("invalid property : inpcheck.kmesh_dk")
        return True

    if method is not None and (method.get_value().lower() != 'monk' and\
                               method.get_value().lower() != 'gamma' and\
                               method.get_value().lower() != 'mesh' ):
        return True

    bgamma=False
    if method is not None and method.get_value().lower() == 'gamma':
        bgamma=True
    
    f=True
    n1 = args['input'].get_root_entry().get_entry('accuracy.ksampling.mesh.nx')
    if n1 is None:
        f = False
        n1 = args['input'].get_root_entry().get_entry('accuracy.ksampling.mp_index.n1')
    if n1 is None:
        return True

    if f:
        n2 = args['input'].get_root_entry().get_entry('accuracy.ksampling.mesh.ny')
        n3 = args['input'].get_root_entry().get_entry('accuracy.ksampling.mesh.nz')
    else:
        n2 = args['input'].get_root_entry().get_entry('accuracy.ksampling.mp_index.n2')
        n3 = args['input'].get_root_entry().get_entry('accuracy.ksampling.mp_index.n3')
    if n2 is None or n3 is None and not bgamma:
        return True

    if not bgamma:
        in1=n1.get_value()
        in2=n2.get_value()
        in3=n3.get_value()
        if in1==0 or in2==0 or in3==0:
            logger.error("invalid kmesh")
            return True
    else:
        in1=1
        in2=1
        in3=1

    aconfig = args['input'].get_atomic_configuration()
    if aconfig is None:
        return True

    (rb1,rb2,rb3) = config.get_reference_kmesh(aconfig,kmesh_dk)
    if rb1 is None or rb2 is None or rb3 is None:
        return True

    try:
        in1 = int(in1)
        in2 = int(in2)
        in3 = int(in3)
    except ValueError:
        logger.error('non-int kmesh : '+str(in1)+', '+str(in2)+', '+str(in3))
        return False
    logger.debug("rb1, rb2, rb3 : "+str(rb1)+", "+str(rb2)+", "+str(rb3))

    if rb1>in1:
        logger.warn("kmesh for the first reciprocal lattice  : ["+str(in1)+"] may not be enough. "+\
                    "reference value : ["+str(rb1)+"].")
        args['input'].increment_nwarn()

    if rb2>in2:
        logger.warn("kmesh for the second reciprocal lattice : ["+str(in2)+"] may not be enough. "+\
                    "reference value : ["+str(rb2)+"].")
        args['input'].increment_nwarn()

    if rb3>in3:
        logger.warn("kmesh for the third reciprocal lattice  : ["+str(in3)+"] may not be enough. "+\
                    "reference value : ["+str(rb3)+"].")
        args['input'].increment_nwarn()

    return True

def check_bond_length(args):
    conf = args['input'].get_atomic_configuration()
    if conf is None:
        return True
    fact = properties.Property().get_property("inpcheck.bond_length_factor",float)
    if fact is None:
        return True
    for at in conf:
        if not at.is_valid():
            return True
    bl = conf.get_bondlenth_list(cutoff_factor=fact)
    if bl is None or len(bl)==0:
        return True
    for b in bl:
        logger.warn("bond length between atom "+str(b[1].index()+1)+" ("+b[1].get_element_name()+\
                    ") and "+str(b[2].index()+1)+" ("+b[2].get_element_name()+") : "+\
                    str(str(b[0])+" [bohr] ").rjust(25) +"may be too short.")
        args['input'].increment_nwarn()

    return True

def check_fract(args):
    conf = args['input'].get_atomic_configuration()
    if conf is None:
        return True

    fact = properties.Property().get_property("inpcheck.max_fract_coord",float)
    if fact is None:
        return True

    found=False
    for atom in conf:
        pos = atom.get_pos(mode=config.FRACTIONAL)
        valid = True
        if math.fabs(pos[0]) > fact or math.fabs(pos[1]) > fact or math.fabs(pos[2]) > fact:
            valid = False
        if not valid:
            logger.warn("fractional coordinate ("+str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+\
                        ") of atom "+str(atom.index()+1)+" may not be valid.")
            found=True
            
    if found:
        args['input'].increment_nwarn()

    return True

def check_inversion(args):
    conf = args['input'].get_atomic_configuration()
    if conf.has_inversion_symmetry():
        inv=args['input'].get_root_entry().get_entry('structure.symmetry.sw_inversion')
        if inv is None or not inv.get_value():
            logger.warn(\
            "the system has inversion symmetry, but [structure.symmetry.sw_inversion] is not 'on'")
            args['input'].increment_nwarn()
    else:
        inv=args['input'].get_root_entry().get_entry('structure.symmetry.sw_inversion')
        if inv is not None and inv.get_value():
            logger.error(\
            "the system does not have inversion symmetry, but [structure.symmetry.sw_inversion] is 'on'")
            args['input'].increment_nerr()
            return False

    return True

