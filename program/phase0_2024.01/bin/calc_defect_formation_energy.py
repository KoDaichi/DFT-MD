#!/usr/bin/env python3

import sys, math;
import argparse
import subprocess
import shlex

class calc_formation_energy:
    def __init__( self, fname_in, fname_out, de, emin, emax, vmin, vmax, image_format ):
        self.Hartree = 27.2116;
        self.vbm = 0.0;     self.ene_host = 0.0;
        self.emin = emin;   self.emax = emax;    self.de = de;
        self.vmin = vmin;   self.vmax = vmax;
        self.banad_gap = -100;

        self.num_elem_defect = 0;     self.elem_defect = [];
        self.num_atom_defect = [];    self.chem_pot_defect = [];

        self.num_charge_state = 0;    self.q_charge_state = [];
        self.ene_charge_state = [];   self.corr_charge_state = [];

        in_f = open( fname_in, "r" );
        self.read_input( in_f );
        in_f.close();

        fname_out_ene_all = fname_out +'.qdep'
        fname_out_ene_min = fname_out +'.min'
        out_f1 = open( fname_out_ene_all, "w" );
        out_f2 = open( fname_out_ene_min, "w" );
        self.write_formation_energy( out_f1, out_f2 );
        out_f1.close();
        self.write_charge_transition_level( out_f2 );
        out_f2.close();

        fname_out_gnu = fname_out +'.gnu'
        out_f = open( fname_out_gnu, "w" );
        self.write_gnuplot_file( out_f, fname_out_ene_all, fname_out_ene_min,
                                 fname_out, image_format );
        out_f.close();

        cmd = 'gnuplot '+fname_out_gnu
        tokens = shlex.split(cmd)
        subprocess.run(tokens)

    def write_charge_transition_level( self, out_f ):
        out_f.write('\n')
        out_f.write('{:s}\n'.format("# Charge   Transtion level [eV]"))
        for j in range( self.num_charge_state -1 ):
            e_form_1 = self.ene_charge_state[j] +self.corr_charge_state[j]
            e_form_2 = self.ene_charge_state[j+1] +self.corr_charge_state[j+1]
            dq = self.q_charge_state[j+1] -self.q_charge_state[j];
            c1 = ( e_form_1 -e_form_2 )/ float(dq) -self.vbm;

            out_f.write('{:s}{:2d}{:s}'.format("#",self.q_charge_state[j+1], "/"))
            out_f.write('{:2d}'.format(self.q_charge_state[j]) )
            out_f.write('{:20.5f} \n'.format(c1))

    def write_gnuplot_file( self, out_f, ene_data_file, ene_min_file, fname_out,
                            image_format ):
        out_f.write('{:s}\n'.format("# Formation energy"));
        out_f.write('{:s}\n'.format("se gr"));
        out_f.write('{:s}\n'.format("set key reverse opaque box"));
        out_f.write('{:s}\n'.format("set st da line"));
        if self.emax < 3.0:
            out_f.write('{:s}\n'.format("set format x '%3.1f'"))
#        out_f.write('{:s}\n'.format("set format y '%3.1f'"))
        out_f.write('{:s}\n'.format("set xlabel 'Fermi energy (eV)'"));
        out_f.write('{:s}\n'.format("set ylabel 'Formation energy (eV)'"));
        out_f.write('{:s}{:8.3f}{:s}{:8.3f}{:s}\n'
                    .format("set xr [",self.emin,":",self.emax,"]"));
        out_f.write('{:s}{:8.3f}{:s}{:8.3f}{:s}\n'
                    .format("set yr [",self.vmin,":",self.vmax,"]"));
        out_f.write('{:s}{:8.3f}{:s}{:s}{:8.3f}{:s}\n'
                    .format("set arrow 1 from ", 0.0, ", graph 0", 
                            " to ", 0.0, ", graph 1 nohead") )
        out_f.write('{:s}{:8.3f}{:s}{:s}{:8.3f}{:s}\n'
                    .format("set arrow 2 from ", self.band_gap, ", graph 0", 
                            " to ", self.band_gap, ", graph 1 nohead") )

        if image_format == "png":
            fname_png = fname_out + ".png"
            out_f.write('{:s}\n'.format("set term pngcairo enhanced"))
            out_f.write('{:s}{:s}{:s}\n'.format("set out '",fname_png,"'"))
        elif image_format == "eps":
            fname_eps = fname_out + ".eps"
            out_f.write('{:s}\n'.format("set term postscript color enhanced"))
            out_f.write('{:s}{:s}{:s}\n'.format("set out '",fname_eps,"'"))

        out_f.write('{:s}\n'.format("plot \\"))

        for j in range( self.num_charge_state ):
            out_f.write('  {:s}{:s}{:s}'.format("'",ene_data_file,"'"))
            out_f.write('{:s}{:d}{:s}'.format(" us 1:",j+2," ") )
            out_f.write('{:s}{:2d}{:s}'.
                        format("title 'q=", self.q_charge_state[j]," '"));
            if j <= self.num_charge_state -1:
                out_f.write('{:s}\n'.format(", \\"));

        out_f.write('{:s}{:s}{:s}'.format("  '",ene_min_file,"'"))
#        out_f.write('{:s}'.format(" us 1:2 title 'minimum'"));
        out_f.write('{:s}'.format(" us 1:2 title 'min.  '"));
        out_f.write('{:s}\n'.format(" lw 3"));
        
    def write_formation_energy( self, out_f1, out_f2 ):
        out_f1.write('{:s}\n'.format("# Formation energy"));
        out_f1.write('{:s}    '.format("# Ef (eV)"));
        for j in range( self.num_charge_state ):
            out_f1.write('{:s}{:d}       '.format("q=", self.q_charge_state[j]));
        out_f1.write('\n')

        out_f2.write('{:s}\n'.format("# Formation energy"));
        out_f2.write('{:s}    '.format("# Ef (eV)"));
        out_f2.write('{:s}\n'.format("min"))

        e1 = ( self.emax -self.emin ) /self.de
        num_ene_points = int( e1 )+1
        for i in range (num_ene_points):
            ene = self.emin +self.de *i
            out_f1.write('{:10.5f}'.format(ene));
            out_f2.write('{:10.5f}'.format(ene));

            e_form_min = 1.0E10
            for j in range( self.num_charge_state ):
                e_form = self.ene_charge_state[j] -self.ene_host
                for k in range( self.num_elem_defect ):
                    e_form = e_form \
                             -self.num_atom_defect[k] *self.chem_pot_defect[k]

                e_form = e_form +self.corr_charge_state[j] \
                         +self.q_charge_state[j]*( ene +self.vbm )
                if e_form < e_form_min :
                    e_form_min = e_form

                out_f1.write('{:10.5f}'.format(e_form));

            out_f1.write('\n');

            out_f2.write('{:10.5f}'.format(e_form_min));
            out_f2.write('\n');

    def read_input( self, in_f ):
        num_chem_pot_data = 0;   chem_pot_elem = [];    chem_pot_val = []

        while True:
            line = in_f.readline()
            if not line: break;
                
            if line.find("&VBM") == 0 :
                line = in_f.readline()
                self.vbm = float(line)

            if line.find("&band_gap") == 0 :
                line = in_f.readline()
                self.band_gap = float(line)

            if line.find("&Defects") == 0 :
                while True:
                    line = in_f.readline()
                    data = line.split();   num_data = len(data)
                    if num_data == 0: break

                    self.elem_defect.append( data[0] )
                    self.num_atom_defect.append( int(data[1]) )
                    self.num_elem_defect = self.num_elem_defect +1

            if line.find("&Chemical_po") == 0 :
                while True:
                    line = in_f.readline()
                    data = line.split();  num_data = len(data)
                    if num_data == 0:  break

                    chem_pot_elem.append( data[0] )
                    chem_pot_val.append( float(data[1]) *self.Hartree )
                    num_chem_pot_data = num_chem_pot_data +1

            if line.find("&Host_supercell_energy") == 0 :
                line = in_f.readline()
                self.ene_host = float(line) *self.Hartree;

            if line.find("&Defective_supercell") == 0 :
                while True:
                    line = in_f.readline()
                    data = line.split();   num_data = len(data)
                    if num_data == 0:  break

                    self.q_charge_state.append( int(data[0]) )
                    self.ene_charge_state.append( float(data[1]) *self.Hartree )
                    self.num_charge_state = self.num_charge_state +1

            if line.find("&Correction") == 0 :
                while True:
                    line = in_f.readline()
                    data = line.split(); num_data = len(data)
                    if num_data == 0: break

                    self.corr_charge_state.append( float(data[1]) );

        for i in range( self.num_elem_defect ):
            for j in range( num_chem_pot_data ):
                if chem_pot_elem[j] == self.elem_defect[i]:
                    self.chem_pot_defect.append( chem_pot_val[j] )

###################################### Main ########
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Caulculte formation energy of a given defect');

    parser.add_argument( 'input', help='input file');
    parser.add_argument( '-o', '--outfile',
                         help='output file name (default:result)',
                         default='result_formation_energy');
    parser.add_argument( '--de', 
                         help='divisoin of energy scale (default:0.01)',
                         type=float, default='0.01');
    parser.add_argument( '--emin', 
                         help='minimum of Fermi energy (default:-1.0)',
                         type=float, default='-1.0');
    parser.add_argument( '--emax', 
                         help='maximum of Fermi energy (default:6.0)',
                         type=float, default='6.0');
    parser.add_argument( '--vmin', 
                         help='minimum of Formation energy (default:-5.0)',
                         type=float, default='-5.0');
    parser.add_argument( '--vmax', 
                         help='maximum of Formation energy (default:5.0)',
                         type=float, default='5.0');
    parser.add_argument( '--image_format', 
                         help='image format of figure (default:png)',
                         default='png');

    args = parser.parse_args();

    calc_formation_energy( args.input, args.outfile, args.de, args.emin, args.emax,
                           args.vmin,  args.vmax,
                           args.image_format );


