#!/usr/bin/env python3
#
##################################################
#
# phonon_band_atom_proj.py
#   version : 1.02
#
#   contact address :  Phase System Consortium
#
##################################################

import sys, math, re;
import argparse, subprocess, shlex

class atom_projected_phonon_band:
    def __init__( self, phonon_file, qpt_file, 
                  atom_id, element, key, z_range, mode_sym, 
                  neglect_mass, disp_squared, proj_qdir, unfolding,
                  e_range, e_inc, 
                  unit, plot_style, circle_scale,
                  cb_range, fig_format, outfile, ref_file,
                  ndiv_erange_map, broadening_width_map, threshold ):

        self.natom = 0;        self.nmodec = 0;          self.nqpt = 0;
        self.atom_mass = [];   self.atom_name = [];      self.atom_key = [];
        self.coord_x = [];     self.coord_y = [];        self.coord_z = [];

        self.num_special_qpt = 0
        self.special_qpt_id = [];   self.special_qpt_name = []
        self.dq_list = [];          self.dq_max = 0.0

        ph_in = open( phonon_file, "r" );
        self.read_header_mode_data( ph_in )

        self.qcart = [ [0.0 for j in range(3) ] for i in range(self.nqpt) ]
        self.qfrac = [ [0.0 for j in range(3) ] for i in range(self.nqpt) ]
        self.freq = [ [0.0 for j in range(self.nmode)] for i in range(self.nqpt)]

        self.mode_sym = [ ["" for j in range(self.nmode)] for i in range(self.nqpt)]
        if unfolding == 1:
            self.unfolding_weight = [ [0.0 for j in range(self.nmode)] \
                                      for i in range(self.nqpt)]
        if proj_qdir == 1:
            self.longitudinal_weight = [ [0.0 for j in range(self.nmode)] \
                                         for i in range(self.nqpt)]

        self.weight = [ [0.0 for j in range(self.nmode)] for i in range(self.nqpt)]
        self.norm   = [ [0.0 for j in range(self.nmode)] for i in range(self.nqpt)]

        self.factor = [ 0.0 for i in range(self.natom) ]
        
        self.re_eigenvec = [ [ [ [ 0.0 for m in range(3)] \
                                for k in range(self.natom) ] \
                               for j in range(self.nmode)] \
                             for i in range(self.nqpt)]
        self.im_eigenvec = [ [ [ [ 0.0 for m in range(3)] \
                                for k in range(self.natom) ] \
                               for j in range(self.nmode)] \
                             for i in range(self.nqpt)]

        self.qdir = [ [0.0 for j in range(3) ] for i in range(self.nqpt) ]

        self.set_factor( atom_id, element, key, z_range )

        self.read_body_mode_data( ph_in, unit, unfolding )
        self.set_qdir()

        if proj_qdir == 1:
            self.longitudinal_weight = [ [0.0 for j in range(self.nmode)] \
                                         for i in range(self.nqpt)]
            self.set_longitudinal_weight()

        self.set_weight_mode_data( neglect_mass, disp_squared, proj_qdir, \
                                   unfolding, mode_sym )

        self.set_dq_list()
        self.read_special_qpt( qpt_file )

        if e_range == "none":
            emin = 1.0E10;  emax = -1.0E10
            for i in range(self.nqpt):
                for j in range(self.nmode):
                    emin = min( emin, self.freq[i][j] )
                    emax = max( emax, self.freq[i][j] )
            emax = emax *1.02
        else:
            emin = float( e_range[0] );  emax = float( e_range[1] )

        if cb_range == "none":
            cbmin = 1.0E10;  cbmax = -1.0E10
            for i in range(self.nqpt):
                for j in range(self.nmode):
                    cbmin = min( cbmin, self.weight[i][j] )
                    cbmax = max( cbmax, self.weight[i][j] )
            cbmax = cbmax *1.02
        else:
            cbmin = float( cb_range[0] );  cbmax = float( cb_range[1] )

        out_data = outfile +".dat";    out_gnu  = outfile +".gnu"

        if ref_file == "none":
            mapping_mode = 0;            out_map = ""
        else:
            mapping_mode = 1;            out_map = 'reference.map'
            ref_in = open( ref_file, "r" );
            self.nqpt_ref = 0
            self.read_header_mode_ref( ref_in )
            self.qcart_ref = [ [0.0 for j in range(3) ] for i in range(self.nqpt_ref) ]
            self.qfrac_ref = [ [0.0 for j in range(3) ] for i in range(self.nqpt_ref) ]
            self.freq_ref = [ [0.0 for j in range(self.nmode_ref)] 
                              for i in range(self.nqpt_ref)]
            self.read_body_mode_ref( ref_in, unit )
            self.set_ref_map( unit, emin, emax, out_map,
                              ndiv_erange_map, broadening_width_map )

        self.write_weight_data( out_data, unit, mapping_mode, threshold )
        self.write_gnu_file( outfile, out_gnu, out_data, emin, emax, e_inc, 
                             unit, plot_style, circle_scale, 
                             cb_range, cbmin, cbmax, fig_format,
                             mapping_mode, out_map )
        self.run_gnuplot( out_gnu )

    def set_ref_map( self, unit, emin, emax, out_map, 
                     ndiv_erange_map, broadening_width_map ):
        qnow = [0,0]; qold = [-99,-99]
        count = 0;  dq = 0.0

        dq_list_2d = [];

        for i in range( self.nqpt_ref ):
            qnow[0] = self.qcart_ref[i][0]
            qnow[1] = self.qcart_ref[i][1]
            if qnow[0] != qold[0] or qnow[1] != qold[1]:
                count = count +1
                c1 = ( qnow[1] -qold[1] )**2 +( qnow[0] -qold[0] )**2
                c1 = math.sqrt( c1 )
                if i == 0:
                    dq_list_2d.append(0.0)
                else:
                    dq = dq +c1
                    dq_list_2d.append(dq)

            qold[0] = qnow[0];  qold[1] = qnow[1]
        
        nqpt_ref_2d = count

        if ndiv_erange_map == "none":
            de = ( emax -emin ) /1000.0;
            ne_div = int( ( emax -emin )/de ) +1
        else:
            ne_div = int(ndiv_erange_map)
            de = ( emax -emin ) /ne_div

        ovp = broadening_width_map
        
        itmp = 0
        spectr = [ [0.0 for j in range(ne_div)] for i in range(nqpt_ref_2d)]
        nn = int( self.nqpt_ref / nqpt_ref_2d )

        for i in range( nqpt_ref_2d ):
            emin_wk = [  1.0E10 for j in range(self.nmode_ref) ]
            emax_wk = [ -1.0E10 for j in range(self.nmode_ref) ]

            for k in range( nn ):
                for j in range( self.nmode_ref ):
                    itmp = nn *i +k
                    c1 = min( emin_wk[j], self.freq_ref[itmp][j] )
                    c2 = max( emax_wk[j], self.freq_ref[itmp][j] )
                    emin_wk[j] = c1;  emax_wk[j] = c2

            for n in range(ne_div):
                count = 0
                for j in range( self.nmode_ref ):
                    e1 = emin +de *n
                    if e1 >= emin_wk[j]-de*ovp and e1 < emax_wk[j]+de*ovp:
                        count = count +1
                if count > 0:
                    spectr[i][n] = 1

        out_f = open( out_map, "w" );
        if unit == "cm-1":
            out_f.write("{:s}\n".format("#      dq           freq[cm-1]       weight"))
        elif unit == "THz":
            out_f.write("{:s}\n".format("#      dq           freq[THz]        weight"))
        elif unit == "meV":
            out_f.write("{:s}\n".format("#      dq           freq[meV]        weight"))

        for i in range(nqpt_ref_2d):
            dq = dq_list_2d[i]
            for j in range(ne_div):
                e1 = emin +de *j
                c1 = spectr[i][j]
                out_f.write("{:15.8f}{:15.8f}{:15.8f}\n".format(dq,e1,c1) )
            out_f.write("\n")
        out_f.close()
        
    def set_factor( self, atom_id, element, key, z_range ):
        num_flag = 0
        count = [ 0 for i in range(self.natom) ]

        if atom_id != 'none':
            num_flag = num_flag +1

            for i in range(len(atom_id)):
                itmp = atom_id[i].split('-');   inum = len(itmp)
                if inum == 2:
                    atom_s = int(itmp[0])-1;  atom_e = int(itmp[1])
                    for j in range(atom_s,atom_e):
                        if j < 0 or j > self.natom -1:
                            pass
                        else:
                            count[j] = count[j] +1
                else:
                    atom_s = int(itmp[0]) -1;    j = atom_s
                    if j < 0 or j > self.natom -1: 
                        pass
                    else:
                        count[j] = count[j] +1

        if element != 'none':
            num_flag = num_flag +1
            for i in range(len(element)):
                for k in range( self.natom ):
                    if self.atom_name[k] == element[i]:
                        count[k] = count[k] +1

        if key != 'none':
            num_flag = num_flag +1

            for i in range(len(key)):
                itmp = key[i].split('-');   inum = len(itmp)
                if inum == 2:
                    key_s = int(itmp[0]);  key_e = int(itmp[1]) +1
                    for j in range(key_s,key_e):
                        for k in range( self.natom ):
                            if self.atom_key[k] == j:
                                count[k] = count[k] +1
                else:
                    key_s = int(itmp[0]);    j = key_s
                    for k in range( self.natom ):
                        if self.atom_key[k] == j:
                            count[k] = count[k] +1

        if z_range != 'none':
            num_flag = num_flag +1
            for i in range(self.natom):
                zmin = float(z_range[0]); zmax = float(z_range[1])
                if self.coord_z[i] >= zmin and self.coord_z[i] <= zmax:
                    count[i] = count[i] +1

        for i in range(self.natom):
            if count[i] == num_flag:
                self.factor[i] = 1.0

    def set_dq_list( self ):
        for i in range( self.nqpt ):
            if i == 0:
                dq = 0.0
            else:
                dx = self.qcart[i][0] -self.qcart[i-1][0]
                dy = self.qcart[i][1] -self.qcart[i-1][1]
                dz = self.qcart[i][2] -self.qcart[i-1][2]
                dq = dq +math.sqrt( dx**2 +dy**2 +dz**2 )
            self.dq_list.append( dq )

        self.dq_max = dq

    def read_special_qpt( self, qpt_file ):
        qpt_in = open( qpt_file, "r" );
        for i in range(4):
            line = qpt_in.readline()

        last = 0
        while True:
            line = qpt_in.readline()
            if not line:
                eof_flag = 1;  break;
                               
            data = line.split()
            if data[4] == "#":
                found = 0
                qx = float( data[0] ) /float( data[3] )
                qy = float( data[1] ) /float( data[3] )
                qz = float( data[2] ) /float( data[3] )

                data = re.split('[#\n]',line);
                for i in range( last,self.nqpt ):
                    dx = qx -self.qfrac[i][0]
                    dy = qy -self.qfrac[i][1]
                    dz = qz -self.qfrac[i][2]
                    d2 = dx**2 +dy**2 +dz**2
                    if d2 < 1.0E-6:
                        found = 1
                        self.special_qpt_id.append( i )
                        self.special_qpt_name.append( data[1] )
                        self.num_special_qpt = self.num_special_qpt +1
                        last = i;  break

    def run_gnuplot( self, out_gnu ):
        cmd = 'gnuplot '+out_gnu
        tokens = shlex.split(cmd)
        subprocess.run(tokens)

    def write_gnu_file( self, outfile, out_gnu, out_data, emin, emax, e_inc, 
                        unit, plot_style, circle_scale, 
                        cb_range, cbmin, cbmax, fig_format,
                        mapping_mode, out_map ):
        out_f = open( out_gnu, "w" );
        out_f.write("set nokey \n");

        if unit == "cm-1":
            out_f.write("set ylabel 'Frequency (cm^{-1})' \n")
        elif unit == "THz":
            out_f.write("set ylabel 'Frequency (THz)' \n")
        elif unit == "meV":
            out_f.write("set ylabel 'Frequency (meV)' \n")

        fname = outfile +'.' +fig_format
        if fig_format == "eps":
            out_f.write("set terminal postscript eps enhanced color solid \n")
            out_f.write(f"set out '{fname}' \n")

        elif fig_format == "png":
#            out_f.write("set terminal pngcairo enhanced size 320, 480 \n")
            out_f.write("set terminal pngcairo enhanced \n")
            out_f.write(f"set out '{fname}' \n")
        else:
            pass

        if mapping_mode == 0:
            out_f.write("set palette rgb 22,13,-31 \n");
            out_f.write("set style line 100 lc palette \n")
        else:
            out_f.write("set palette gray \n");
            out_f.write("set style line 100 lc palette \n")

#        out_f.write(f"set xr [0.0:{self.dq_max}]\n" )
        out_f.write("{:s}{:8.5f}{:s}{:8.5f}{:s}\n".format(
            "set xr [",0.0,":",self.dq_max,"]") )

        out_f.write(f"set yr [{emin}:{emax}]\n" )
        if e_inc != "none":
            e1 = float( e_inc[0] )
            out_f.write(f"set ytics {e1}\n" )

        if cb_range != "none":
            out_f.write(f"set cbrange [{cbmin}:{cbmax}]\n" )

        if mapping_mode == 0:
            out_f.write("set gr \n");
            for i in range(self.num_special_qpt):
                j = self.special_qpt_id[i]
                out_f.write(f"set arrow from {self.dq_list[j]},{emin} ")
                out_f.write(f"to {self.dq_list[j]},{emax} nohead lt -1 \n")
        else:
            out_f.write("set gr front \n");
            out_f.write("set tics front \n")
            for i in range(self.num_special_qpt):
                j = self.special_qpt_id[i]
                out_f.write(f"set arrow front from {self.dq_list[j]},{emin} ")
                out_f.write(f"to {self.dq_list[j]},{emax} nohead lt -1 \n")

        out_f.write("{:s}".format("set xtics (") )
        for i in range(self.num_special_qpt):
            j = self.special_qpt_id[i]
            out_f.write("{:s}{:s}{:s}{:8.5f}".format(
                '"', self.special_qpt_name[i], '"' ,self.dq_list[j] ) )

            if i < self.num_special_qpt-1:
                out_f.write("{:s}".format(",") )
            else:
                out_f.write("{:s}\n".format(")") )

        if fig_format == "eps":
            transp = 0.8
        else:
            transp = 0.6


        if mapping_mode == 0:
            if plot_style == 1:
                ctmp = 0.01
                out_f.write(f"plot '{out_data}' u 1:2:({circle_scale}*{ctmp}*$3) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write( " lc rgb 'web-blue', \\\n")
            elif plot_style == 2:
                ctmp = 0.005
                out_f.write(f"plot '{out_data}' u 1:2:({circle_scale}*{ctmp}):($3) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write(" ls 100, \\\n")
            elif plot_style == 3:
                ctmp = 0.5 /( cbmax -cbmin ) 
                ctmp = 0.01
                out_f.write(f"plot '{out_data}' u 1:2:({circle_scale}*{ctmp}*$3):($3) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write(" ls 100, \\\n")

            out_f.write(f"     '{out_data}' u 1:2 w l lw 0.2 lc 'black' \n")

        else:
            out_f.write(f"set cbrange [0:1] \n" )
            out_f.write(f"unset colorbox \n" )
            out_f.write(f"plot '{out_map}' u 1:2:(1-0.2*$3) w image, \\\n")
            out_f.write(f"     '{out_data}' u 1:2 w l lt 7 lw 2 lc rgb 'web-blue'\\\n")

        out_f.close()

    def set_qdir( self ):
        vec = [ 0.0 for j in range(3) ]

        for iq in range(self.nqpt-1):
            for k in range(3):
                vec[k] = self.qcart[iq+1][k] -self.qcart[iq][k]
            c1 = math.sqrt( vec[0]**2 +vec[1]**2 +vec[2]**2 )
            for k in range(3):
                self.qdir[iq][k] = vec[k] /c1

        for k in range(3):
            self.qdir[self.nqpt-1][k] = self.qdir[self.nqpt-2][k]

        for iq in range(self.nqpt):
            for k in range(3):
                vec[k] = self.qcart[iq][k]
            c1 = math.sqrt( vec[0]**2 +vec[1]**2 +vec[2]**2 )

            if c1 > 1E-8:
                for k in range(3):
                    self.qdir[iq][k] = vec[k] /c1
                

    def write_weight_data( self, outfile, unit, mapping_mode, threshold ):
        out_f = open( outfile, "w" );

        if unit == "cm-1":
            out_f.write("{:s}\n".format("#      dq           freq[cm-1]       weight"))
        elif unit == "THz":
            out_f.write("{:s}\n".format("#      dq           freq[THz]        weight"))
        elif unit == "meV":
            out_f.write("{:s}\n".format("#      dq           freq[meV]        weight"))

        if mapping_mode == 0:
            for j in range( self.nmode ):
                for i in range( self.nqpt ):
                    if self.norm[i][j] == 0.0:
                        c1 = 0.0
                    else:
                        c1 = self.weight[i][j] /self.norm[i][j]

                    out_f.write("{:15.8f}{:15.8f}{:15.8f}\n".format( 
                        self.dq_list[i], self.freq[i][j], c1 ) )
                out_f.write("\n")

        else:
            for j in range( self.nmode ):
                flag = 0;  flag_old = 0
                for i in range( self.nqpt ):
                    if self.norm[i][j] == 0.0:
                        c1 = 0.0
                    else:
                        c1 = self.weight[i][j] /self.norm[i][j]

                    if c1 > threshold:
                        flag = 1
                        if flag_old == 0: out_f.write("\n")
                        out_f.write("{:15.8f}{:15.8f}{:15.8f}\n".format( 
                            self.dq_list[i], self.freq[i][j], c1 ) )
                    else:
                        flag = 0

                    flag_old = flag
                        
        out_f.close()

    def read_header_mode_ref( self, ref_in ):
        eof_flag = 0
        while True:
            line = ref_in.readline()
            if not line:
                eof_flag = 1;  break;

            if line.find('Nqvec') > 0:
                data = line.split();
                ntarget = len(data)-1;
                self.nmode_ref = int( data[1] )
                self.nqpt_ref = int( data[ntarget] )
                break;

    def read_body_mode_ref( self, ref_in, unit ):
        eof_flag = 0;    mode_flag = 0

        conv_from_cm_to_THz = 0.02998E0;    conv_from_cm_to_meV = 0.12398E0

        while True:
            line = ref_in.readline()
            if not line:
                eof_flag = 1;  break;

            if line.find('iq=') > 0:
                mode_flag = 0;
                data = line.split()
                qpt_id = int( data[1] ) -1

                data = re.split('[(,)]',line);
                self.qcart_ref[qpt_id][0] = float(data[5])
                self.qcart_ref[qpt_id][1] = float(data[6])
                self.qcart_ref[qpt_id][2] = float(data[7])

            if mode_flag == 1:
                if line.find('hbar') > 0:
                    data = line.split()
                    ntarget = len(data)-2;
                    if unit == "cm-1":
                        self.freq_ref[qpt_id][mode_id] = float( data[ntarget] )
                    elif unit == "THz":
                        self.freq_ref[qpt_id][mode_id] = float( data[ntarget] ) \
                                                         *conv_from_cm_to_THz
                    elif unit == "meV":
                        self.freq_ref[qpt_id][mode_id] = float( data[ntarget] ) \
                                                         *conv_from_cm_to_meV

            if line.find('n=') > 0:
                data = line.split()
                mode_id = int( data[1] ) -1
                mode_flag = 1;

    def read_body_mode_data( self, ph_in, unit, unfolding ):
        eof_flag = 0;    mode_flag = 0

        conv_from_cm_to_THz = 0.02998E0;    conv_from_cm_to_meV = 0.12398E0

        while True:
            line = ph_in.readline()
            if not line:
                eof_flag = 1;  break;

            if line.find('iq=') > 0:
                mode_flag = 0;
                data = line.split()
                qpt_id = int( data[1] ) -1

                data = re.split('[(,)]',line);
                self.qfrac[qpt_id][0] = float(data[1])
                self.qfrac[qpt_id][1] = float(data[2])
                self.qfrac[qpt_id][2] = float(data[3])
                self.qcart[qpt_id][0] = float(data[5])
                self.qcart[qpt_id][1] = float(data[6])
                self.qcart[qpt_id][2] = float(data[7])

            if mode_flag == 1:
                if line.find('hbar') > 0:
                    data = line.split()
                    ntarget = len(data)-2;
                    if unit == "cm-1":
                        self.freq[qpt_id][mode_id] = float( data[ntarget] )
                    elif unit == "THz":
                        self.freq[qpt_id][mode_id] = float( data[ntarget] ) \
                                                     *conv_from_cm_to_THz
                    elif unit == "meV":
                        self.freq[qpt_id][mode_id] = float( data[ntarget] ) \
                                                     *conv_from_cm_to_meV

                elif line.find('n=') <= 0:
                    data = line.split()
                    if len(data) == 4:
                        isum = isum +1
                        atom_id = int( data[0] ) -1
                        if isum <= self.natom:
                            for m in range(3):
                                self.re_eigenvec[qpt_id][mode_id][atom_id][m] \
                                = float( data[m+1] )
                        else:
                            for m in range(3):
                                self.im_eigenvec[qpt_id][mode_id][atom_id][m] \
                                = float( data[m+1] )

            if line.find('n=') > 0:
                data = line.split()
                mode_id = int( data[1] ) -1
                self.mode_sym[qpt_id][mode_id] = data[2]
                if unfolding == 1: 
                    self.unfolding_weight[qpt_id][mode_id] = float(data[5])

                isum = 0;   mode_flag = 1;

    def set_longitudinal_weight( self ):
        for iq in range(self.nqpt):
            for n in range(self.nmode):
                ctmp1 = 0.0;    ctmp2 = 0.0
                for i in range(self.natom):
                    fac1 = self.atom_mass[i]
                    c1 = 0.0;  c2 = 0.0
                    for k in range(3):
                        c1 = c1 +self.re_eigenvec[iq][n][i][k] \
                             *self.qdir[iq][k]
                        c2 = c2 +self.im_eigenvec[iq][n][i][k] \
                             *self.qdir[iq][k]
                    ctmp1 = ctmp1 +( c1**2 +c2**2 ) /fac1

                    c1 = 0.0;  c2 = 0.0
                    for k in range(3):
                        c1 = c1 +self.re_eigenvec[iq][n][i][k]**2
                        c2 = c2 +self.im_eigenvec[iq][n][i][k]**2
                    ctmp2 = ctmp2 +( c1 +c2 ) /fac1

                self.longitudinal_weight[iq][n] = math.sqrt( ctmp1 /ctmp2 )

    def set_weight_mode_data( self, neglect_mass, disp_squared, proj_qdir, unfolding,
                              mode_sym ):
        vec = [0.0 for k in range(3)]

        for iq in range(self.nqpt):
            for n in range(self.nmode):
                if mode_sym != 'none':
                    if ( self.mode_sym[iq][n] in mode_sym ):
                        pass
                    else:
                        continue
                else:
                    pass
    
                for i in range(self.natom):
                    dr2 = 0.0
                    for k in range(3):
                        dr2 = dr2 +self.re_eigenvec[iq][n][i][k]**2 \
                              +self.im_eigenvec[iq][n][i][k]**2
                        
                    if neglect_mass == 0:
                        dr2 = dr2 /self.atom_mass[i]
                    
                    if disp_squared == 1:
                        d2 = dr2
                    else:
                        d2 = math.sqrt(dr2)

                    c2 = self.factor[i] *d2
                        
                    self.weight[iq][n] = self.weight[iq][n] +c2
                    self.norm[iq][n]   = self.norm[iq][n] +d2
                
                if unfolding == 1:
                    self.weight[iq][n] = self.weight[iq][n] \
                                         *self.unfolding_weight[iq][n]
                if proj_qdir == 1:
                    self.weight[iq][n] = self.weight[iq][n] \
                                         *self.longitudinal_weight[iq][n]
                    
    def read_header_mode_data( self, ph_in ):
        eof_flag = 0;     atom_flag = 0;    mode_flag = 0
        count = 0;
        conv_from_aumass_to_atommass = 1822.877332825052;

        while True:
            line = ph_in.readline()
            if not line:
                eof_flag = 1;  break;

            if atom_flag == 1:
                data = line.split();
                ndata = len(data);
                self.coord_x.append( float(data[1]) )
                self.coord_y.append( float(data[2]) )
                self.coord_z.append( float(data[3]) )
                self.atom_mass.append( float(data[4]) /conv_from_aumass_to_atommass );
                self.atom_name.append( data[5] );
                if ndata == 7: self.atom_key.append( int(data[6]) )

                count = count +1;
                if count == self.natom:     
                    atom_flag = 0;  count = 0

            if line.find('Natom') > 0 and line.find('Nqvec') <= 0:
                data = line.split();
                ntarget = len(data)-1;
                self.natom = int( data[ntarget] )
                atom_flag = 1

            if line.find('Nqvec') > 0:
                data = line.split();
                ntarget = len(data)-1;
                self.nmode = int( data[1] )
                self.nqpt = int( data[ntarget] )
                break;
                

###################################### Main ########
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='phonon dispersion with weight')

    parser.add_argument( 'phonon_file', help='mode.data' )
    parser.add_argument( 'qpt_file', help='q-point generation control file' )
    parser.add_argument( '--atom_id', help='list of atom ids',
                         default='none', nargs='*' )
    parser.add_argument( '--element', help='list of elements',
                         default='none', nargs='*' )
    parser.add_argument( '--key', help='list of atom keys',
                         default='none', nargs='*' )
    parser.add_argument( '--z_range', help='height range (min and max)',
                         default='none', nargs=2 )
    parser.add_argument( '--mode_sym', help='mode symmetry',
                         default='none', nargs='*' )

    parser.add_argument( '--neglect_mass',
                         help='atom mass is assumed to be 1.00',
                         action='store_true')
    parser.add_argument( '--disp_squared',
                         help='calc square of atomic displacement',
                         action='store_true')
    parser.add_argument( '--proj_qdir',
                         help='projection along q vector',
                         action='store_true')
    parser.add_argument( '--unfolding',
                         help='plot with phonon band-unfolding weight',
                         action='store_true')

    parser.add_argument( '--e_range', help='energy range (min and max)',
                         default='none', nargs=2 )
    parser.add_argument( '--e_inc', help='energy increment',
                         default='none', nargs=1 )

    parser.add_argument( '--unit', 
                         help='unit (meV, THz, cm-1) (default:cm-1)',
                         default='cm-1' )
    parser.add_argument( '--plot_style', 
                         help='plot style (1,2,3) (default:1)',
                         type=int, default=1 )
    parser.add_argument( '--circle_scale', 
                         help='circle scale for plotting (default:1.0)',
                         type=float, default=1.0 )

    parser.add_argument( '--cb_range', help='color range (min and max)',
                         default='none', nargs=2 )

    parser.add_argument( '--fig_format', 
                         help='figure format (eps,png) (default:eps)',
                         default='eps' )
    parser.add_argument( '--out_file', 
                         help='output filename (default:atom_projected_phonon_band)',
                         default='atom_projected_phonon_band' )

    parser.add_argument( '--ref_file', 
                         help='phonon file for bulk',
                         default='none' )
    parser.add_argument( '--ndiv_erange_map', 
                         help='divison of energy range for mapping bulk phonon',
                         default='none' )
    parser.add_argument( '--broadening_width_map', 
                         help='broadening width for mapping bulk phonon',
                         type=int, default=5 )
    parser.add_argument( '--threshold', 
                         help='threshold for plotting when ref_file is used (default:0.01)',
                         type=float, default=0.01 )

    args = parser.parse_args()

    atom_projected_phonon_band( args.phonon_file, args.qpt_file, 
                                args.atom_id, 
                                args.element, args.key, args.z_range, args.mode_sym,
                                args.neglect_mass, args.disp_squared,
                                args.proj_qdir, args.unfolding,
                                args.e_range, args.e_inc, 
                                args.unit, args.plot_style, args.circle_scale,
                                args.cb_range, args.fig_format, args.out_file,
                                args.ref_file, args.ndiv_erange_map, 
                                args.broadening_width_map,
                                args.threshold )
