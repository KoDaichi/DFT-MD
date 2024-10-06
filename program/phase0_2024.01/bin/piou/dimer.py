#!/usr/bin/env python
from datetime import datetime
from piou import config
from piou.config import phase
from piou.data import elements
import optparse
import math
import os
import shutil
import sys

dalton = 1822.877332827
Hartree_2_eV = 27.21139615
eV_2_Hartree = 1/Hartree_2_eV
cminv_2_eV      = 1.233984185e-4
cminv_2_Hartree = cminv_2_eV*eV_2_Hartree

Bohr_2_Angstrom = 0.5291772480

def get_energy():
    nfefn = open("nfefn.data")
    energy = 0.0
    for line in nfefn:
        words = line.split()
        try:
            energy = float(words[2])
        except ValueError:
            pass
    return energy

def run_phase(phase,np,ne,nk):
    comm = "mpirun -n "+str(np)+" "+phase+" ne="+str(ne)+" nk="+str(nk)
    print "executing : "+comm
    sys.stdout.flush()
    os.system(comm)

def morse(best_bl=None,best_ene=None):
    if best_bl is None or best_ene is None:
        best_bl = -1
        best_ene = +1000
        blfile = open('bl_energy.data')
        for line in blfile:
            words = line.strip().split()
            bl = float(words[0])
            ene = float(words[1])
            if ene<best_ene:
                best_ene = ene
                best_bl = bl
    print ""
    print "perform fit to the morse potential"
    os.system("rm -f fit.log")
    gnstr  = "f(x) = a*(1-exp(-b*(x-c)))**2+d\n"
    gnstr += "a=0.1\n"
    gnstr += "b=0.1\n"
    gnstr += "c="+str(best_bl)+"\n"
    gnstr += "d="+str(best_ene)+"\n"
    gnstr += "fit f(x) \'bl_energy.data\' using 1:2 via a,b,c,d\n"
    gnstr += "plot \'bl_energy.data' u 1:($2*27.2114) not w p ps 2 pt 6\n"
    gnstr += "replot f(x)*27.2114 not w l lt 1\n"
    gnstr += "set ylabel \'energy (ev)\'\n"
    gnstr += "set xlabel \'bond length (\\305)\'\n"
    gnstr += "set encoding iso\n"
    gnstr += "set term post eps color solid enhanced 18\n"
    gnstr += "set output \'ebl.eps\' \n"
    gnstr += "replot"

    gnfile = open("fit.gn","w")
    gnfile.write(gnstr+"\n")
    gnfile.close()

    print "begin fit..."

    os.system("gnuplot < fit.gn")

    fitlog = open('fit.log')
    in_final = False
    keys=["a","b","c","d"]
    params={}
    for line in fitlog:
        if line.startswith("Final set of parameters"):
            in_final=True
        if not in_final:
            continue
        if line.find("=")<0:
            continue
        words = line.split()
        for i in range(len(keys)):
            if line.startswith(keys[i]):
                params[keys[i]] = float(words[2])
    print "...done."
    for i in range(len(keys)):
        print keys[i]+" = "+str(params[keys[i]])
    print ""
    best_bl = params["c"]
    spring_const = math.pow(params["b"],2)*2.0*params["a"]
    mdenom=0
    mnum=1
    for ename in elem_nam:
        m = elements.ElementInfo().get_element_attribute(ename,'mass')
        mdenom += m
        mnum *= m
    mu = mnum/mdenom

    vib = 0.5*math.sqrt(spring_const/mu)/cminv_2_Hartree
    print "optimized bond length => "+str(best_bl)+" Angstrom"
    print "vibrational frequency => "+str(vib)+" cm^{-1}"
    print ""
    return best_bl

def dimer_main():
    parser = optparse.OptionParser(usage="\n%prog [OPTION] ",version="%prog 0.5")
    parser.add_option("-n","--np",dest="np",default=1)
    parser.add_option("-e","--ne",dest="ne",default=1)
    parser.add_option("-k","--nk",dest="nk",default=1)
    parser.add_option("-p","--phase",dest="phase",default="/s0/home2/jkoga/phase_dev/open66_wk/phase")
    parser.add_option("-d","--dir",dest="dir",default=".")
    parser.add_option("","--dryrun",dest="dryrun",action="store_true",default=False)
    parser.add_option("","--start",dest="start",type=float,default=0.8)
    parser.add_option("","--end",dest="end",type=float,default=1.2)
    parser.add_option("","--npoints",dest="npoints",type=float,default=15)

    (options, args) = parser.parse_args()

    os.chdir(options.dir)

    print "calculating the optimal bond-length..."
    shutil.copyfile("nfinp.data","nfinp.data.org")
    nfinp = phase.Input()
    conf = nfinp.get_atomic_configuration()
    pos0 = conf.get_atom_at(0).get_pos(config.CARTESIAN)
    pos1 = conf.get_atom_at(1).get_pos(config.CARTESIAN)
    blen = 0.0
    for i in range(3):
        blen += math.pow(pos0[i]-pos1[i],2)
    blen = math.sqrt(blen)
    blen0 = options.start*blen
    blen1 = options.end*blen
    dlen = (blen1-blen0)/float(options.npoints-1)
    if not options.dryrun:
        bl_ene = open("bl_energy.data","w")
    best_bl = -1
    best_ene = +10000000
    elem_nam=[]
    for i in range(options.npoints):
        blen = blen0+dlen*i
        print ""
        print "calculating the total energy for bond-length : "+str(blen)+" A"
        atmtable = nfinp.get_root_entry().get_entry("structure.atom_list.atoms.table")
        atom0 = atmtable.get_row_data(0)
        atom1 = atmtable.get_row_data(1)
        if(i==0):
            elem_nam.append(atom0['element'])
            elem_nam.append(atom1['element'])
        atom0['rx'] = 0.0
        atom0['ry'] = 0.0
        atom0['rz'] = blen*0.5
        atom1['rx'] = 0.0
        atom1['ry'] = 0.0
        atom1['rz'] = -blen*0.5
        if 'mobile' in atom0.keys():
            atom0['mobile'] = False
        if 'mobile' in atom1.keys():
            atom1['mobile'] = False
        atmtable.reset_data()
        atmtable.add_row_data(atom0)
        atmtable.add_row_data(atom1)
        nfinp.save()
        if not options.dryrun:
            run_phase(options.phase,options.np,options.ne,options.nk)
        ene = get_energy()
        if ene<best_ene:
            best_ene = ene
            best_bl = blen
        if not options.dryrun:
            bl_ene.write(str(blen)+" "+str(ene)+"\n")
            print "calculated energy : "+str(ene)+" hartree"
        if not options.dryrun:
            bl_ene.flush()
    if not options.dryrun:
        bl_ene.close()
    print "...done"

    best_bl = morse(best_bl,best_ene)

    print "calculating the energy at the optimized bond-length"
    atmtable = nfinp.get_root_entry().get_entry("structure.atom_list.atoms.table")
    atom0 = atmtable.get_row_data(0)
    atom1 = atmtable.get_row_data(1)
    atom0['rx'] = 0.0
    atom0['ry'] = 0.0
    atom0['rz'] = best_bl*0.5
    atom1['rx'] = 0.0
    atom1['ry'] = 0.0
    atom1['rz'] = -best_bl*0.5
    atmtable.reset_data()
    atmtable.add_row_data(atom0)
    atmtable.add_row_data(atom1)
    nfinp.save()
    run_phase(options.phase,options.np,options.ne,options.nk)
    best_ene = get_energy()
    print "optimized energy : "+str(best_ene)+" hartree"
    print "...done"

    shutil.copyfile("nfinp.data.org","nfinp.data")

    try:
        os.makedirs("atom")
    except OSError:
        pass
    shutil.copyfile("file_names.data","atom/file_names.data")
    shutil.copyfile("nfinp.data","atom/nfinp.data")
    os.chdir("atom")
    fnamesdata = open("file_names.data")
    newfnames = ""
    for fname in fnamesdata:
        if fname.strip().startswith("/"):
            continue
        if fname.startswith("F_POT"):
            words = fname.split("=")
            path = words[1].strip()
            newpp = words[0]+"= \'../"+path[1:-1]+"\'"
            newfnames += newpp+"\n"
        else:
            newfnames += fname
    newfnames += "/"
    fnamesdata.close()

    fnamesdata = open("file_names.data","w")
    fnamesdata.write(newfnames)
    fnamesdata.close()

    fnames = open('file_names.data')
    for fname in fnames:
        ff = fname.split("=")
        if ff[0].strip().startswith("F_POT"):
            pot = ff[1].strip()[1:-1]
    fnames.close()
    print "potential : "+pot

    fpot = open(pot)
    reading_elec = False
    reading_pp = False
    line_count=0
    num_elec_conf=0
    nup = 0
    ndown = 0
    electron_config={}
    pp_config = {}
    norb = -1

    for line in fpot:
        ll = line.strip()
        words = ll.split()
        if len(words)<2:
            continue
        if ll.startswith("electron_config"):
            line_count=0
            reading_elec = True
            num_elec_conf = int(words[1].strip())
            continue
        if ll.startswith("pseudo_potential"):
            reading_pp = True
        if not reading_elec and not reading_pp:
            continue

        if reading_elec:
            line_count+=1
            ne = float(words[1])
            electron_config[words[0].strip()] = ne
            if line_count==num_elec_conf:
                reading_elec = False
        if reading_pp:
            if words[0].strip() == "orbitals":
                norb = int(words[1])
            if norb<0:
                continue 
            if words[0].strip() in electron_config.keys():
                pp_config[words[0].strip()] = electron_config[words[0].strip()]

    fpot.close()

    for key in pp_config.keys():
        nv = pp_config[key]
        norb = -1
        if key.find("s")>0:
            norb = 1
        if key.find("p")>0:
            norb = 3
        if key.find("d")>0:
            norb = 5
        if key.find("f")>0:
            norb = 7

        nvh = nv*0.5
        nu=0.0
        nd=0.0
        if nv<norb:
            nu   = float(nv)
        else:
            nu = float(norb)
            nd = float(nv-norb)
        print 'number of up   electrons in orbital '+key+" = "+str(nu)
        print 'number of down electrons in orbital '+key+" = "+str(nd)
        nup += nu
        ndown += nd

    print 'nup, ndown : '+str(nup)+", "+str(ndown)
    zeta = (nup-ndown)/(nup+ndown)
    nb = (nup+ndown)*1.2
    print '=> zeta : '+str((nup-ndown)/(nup+ndown))

    nfinp = phase.Input()

    root = nfinp.get_root_entry()
    elem_table = root.get_entry("structure.element_list.table")

    if not 'zeta' in elem_table.get_identifiers():
        elem_table.add_identifier('zeta')

    for i in range(elem_table.get_num_data()):
        elem_table.set_data(i,'zeta',zeta)

    cmix = root.get_entry("charge_mixing.mixing_methods.table")
    cmix.set_data(0,'rmxs',0.1)
    cmix.set_data(0,'rmxe',0.1)

    smethod = root.get_entry("accuracy.ksampling.method")
    if smethod is None:
        smethod = phase.PrimitiveEntry("method")
        root.get_entry("accuracy.ksampling").add_entry(smethod)
    smethod.set_value("gamma")

    symmblock = root.get_entry("structure.symmetry")
    latsys = None
    lattice_system = ""
    if symmblock is not None:
        tsp = symmblock.get_entry("tspace")
        if tsp is not None:
            latsys = tsp.get_entry("lattice_system")
            if latsys is not None:
                lattice_system = latsys.get_value()
                tsp.remove_entry(latsys)

    mag = root.get_entry("structure.magnetic_state")
    if mag is None:
        mag = phase.PrimitiveEntry("magnetic_state")
        root.get_entry("structure").add_entry(mag)
    mag.set_value("ferro")

    iwav = root.get_entry("accuracy.initial_wavefunctions")
    if iwav is None:
        iwav = phase.PrimitiveEntry("initial_wavefunctions")
        root.get_entry("accuracy").add_entry(iwav)
    iwav.set_value("atomic_orbitals")

    ichg = root.get_entry("accuracy.initial_charge_density")
    if ichg is None:
        ichg = phase.PrimitiveEntry("initial_charge_density")
        root.get_entry("accuracy").add_entry(ichg)
    ichg.set_value("atomic_charge_density")

    root.get_entry("accuracy.num_bands").set_value(int(nb+1))

    atmtable = root.get_entry("structure.atom_list.atoms.table")
    nat_per_cell = atmtable.get_num_data()
    atom0 = atmtable.get_row_data(0)
    atom0['rx'] = 0.0
    atom0['ry'] = 0.0
    atom0['rz'] = 0.0
    atmtable.reset_data()
    atmtable.add_row_data(atom0)

    root.get_entry("structure").add_entry(phase.BlockEntry("ferromagnetic_state"))
    root.get_entry("structure.ferromagnetic_state").add_entry(phase.PrimitiveEntry("sw_fix_total_spin"))
    root.get_entry("structure.ferromagnetic_state.sw_fix_total_spin").set_value("on")
    root.get_entry("structure.ferromagnetic_state").add_entry(phase.PrimitiveEntry("total_spin"))
    root.get_entry("structure.ferromagnetic_state.total_spin").set_value(nup-ndown)
    root.get_entry("structure.ferromagnetic_state").add_entry(phase.PrimitiveEntry("spin_fix_period"))
    root.get_entry("structure.ferromagnetic_state.spin_fix_period").set_value(15)

    nfinp.save()

    print "performing atom calculations..."
    if not options.dryrun:
        run_phase(options.phase,options.np,options.ne,options.nk)
    atomenergy = get_energy()
    print "...done"
    print "energy of the atom : "+str(atomenergy)+" hartree"
    print ""

    print "bond-length           : "+str(best_bl)+" Angstrom"
    print "cohesive energy       : "+str((atomenergy*2-best_ene)*Hartree_2_eV)+" eV/dimer"
    print "vibrational frequency : "+str(vib)+" cm^{-1}"

    os.chdir("..")

if __name__=='main':
    dimer_main()
