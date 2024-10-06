'''define constants for the package 'piou (PHASE I/O utilities) ' '''
import math

Hartree_2_eV = 27.21139615

eV_2_Hartree = 1/Hartree_2_eV
cminv_2_eV      = 1.233984185e-4
cminv_2_Hartree = cminv_2_eV*eV_2_Hartree
cminv_2_THz     = 1.0/33.0
eV_2_cminv      = 1/cminv_2_eV
rad_2_angle     = 180.0/math.pi
angle_2_rad     = math.pi/180.0

Bohr_2_Angstrom = 0.5291772480

Angstrom_2_Bohr = 1/Bohr_2_Angstrom

Angstrom_2_nm   = 0.1
nm_2_Angstrom   = 10.0

very_small_value = 1e-12

name                 = 'name'
has_table            = 'has_table'
choice               = 'choice'
description          = 'description'
num_elements         = 'num_elements'
columns              = 'columns'
entry_type           = 'entry_type'
entry_type_block     = 'block'
entry_type_primitive = 'primitive'
entry_type_table     = 'table'
entry_type_vector    = 'vector'

val_type             = 'val_type'
val_type_float       = 'float'
val_type_int         = 'int'
val_type_string      = 'string'
val_type_bool        = 'bool'
val_type_fract       = 'fract'

val_range            = 'val_range'

unit_type            = 'unit_type'
unit_type_none       = 'none'
unit_type_length     = 'length'
unit_type_energy     = 'energy'
unit_type_time       = 'time'
unit_type_longtime   = 'longtime'
unit_type_velocity   = 'velocity'
unit_type_force      = 'force'
unit_type_pressure   = 'pressure'
unit_type_mass       = 'mass'
unit_type_angle      = 'angle'
unit_type_temperature= 'temperature'

unit_spec={
unit_type_none:[''],
unit_type_length:['bohr','angstrom','nm'],
unit_type_energy:['hartree','ev','rydberg'],
unit_type_time:['au_time','fs','ps','ns'],
unit_type_longtime:['sec','min','hour','day'],
unit_type_velocity:['bohr/au_time','bohr/fs','angstrom/fs','angstrom/au_time','nm/fs','nm/au_time'],
unit_type_force:['hartree/bohr','hartree/angstrom','hartree/nm','ev/angstrom','ev/bohr','ev/nm','rydberg/bohr','rydberg/angstrom','rydberg/nm'],
unit_type_pressure:['hartree/bohr3','hartree/angstrom3','hartree/nm3','ev/angstrom3','ev/bohr3','ev/nm3','rydberg/angstrom3','rydberg/bohr3','rydberg/nm3'],
unit_type_mass:['au_mass','atomic_mass'],
unit_type_angle:['degree','radian'],
unit_type_temperature:['k']
}

unit_conversion = {
 '':1.0

,'bohr':1.0
,'angstrom':0.5291772480
,'nm':0.05291772480

,'hartree':1.0
,'rydberg':2.0
,'ev':27.211396150

,'au_time':1.0
,'fs':2.418884327e-2
,'ps':2.418884327e-5
,'ns':2.418884327e-8

,'sec':1.0
,'min':1.0/60.0
,'hour':1.0/3600.0
,'day':1.0/86400.0

,'au_mass':9.1093897e-31
,'atomic_mass':1.66053e-27

}

au_mass=9.1093897e-31
atomic_mass=1.66053e-27

kbt=3.16676087e-6 # kbt in hartree 

default_units = {
unit_type_none:'',
unit_type_length:'bohr',
unit_type_energy:'hartree',
unit_type_time:'au_time',
unit_type_longtime:'sec',
unit_type_velocity:'bohr/au_time',
unit_type_force:'hartree/bohr',
unit_type_pressure:'hartree/bohr3',
unit_type_mass:'au_mass',
unit_type_angle:'degree',
unit_type_temperature:'k'
}

def get_unit_factor(from_unit, to_unit):
    fr = unit_conversion[from_unit.strip().lower()]
    to = unit_conversion[to_unit.strip().lower()]
    fac=1
    if fr!=None and to!=None:
        fac = to/fr
    return fac

condition  = 'condition'
importance = 'importance'
importance_required='required'
importance_recommended = 'recommended'
importance_optional = 'optional'

begin_tag='[begin]'
end_tag = '[end]'

NOT_READY=1

BLOCK = 'block'
TABLE = 'table'
PRIMITIVE = 'primitive'
VECTOR = 'vector'

import logging

def load_spec_dict(filename):
    ret={}
    try:
        f=open(filename,'r')
        nam=''
        tmpdict={}
        begin=False
        for line in f:
            line=line.strip()
            vals = line.split('=')
            if line.startswith(begin_tag):
                begin=True
            elif line.startswith(end_tag):
                tmpdict={}
            elif len(vals)>=2 and begin:
                if vals[0]==name:
                    ret[vals[1]] = tmpdict
                else:
                    tmpdict[vals[0]]=vals[1]
        f.close()
        return ret
    except IOError:
        logger.error('failed read of: '+filename)
        return None
