#!/usr/bin/env python3
#
##################################################
#
# band_local_decomp.py
#   version : 1.02
#
#   contact address :  Phase System Consortium
#
##################################################

import sys, math, re;
import argparse, subprocess, shlex

class plot_lband:
    def __init__( self, band_file, kpt_file, lband_weight_file,
                  atom_id, layer_id,
                  e_range, e_inc, unit,
                  plot_style, circle_scale, cb_range,
                  line_width, window_scale, vertical,
                  fig_format, out_file, 
                  with_fermi, no_dispersion, threshold ):

        self.nkpt = 0;       self.nband = 0;     self.nspin = 0;
        self.fermi = 0.0;
        self.natom = 0;      self.nlayer = 0;
        self.ndim_magmom = 1;

        self.recip_vec1 = [ 0.0 for j in range(3) ]
        self.recip_vec2 = [ 0.0 for j in range(3) ]
        self.recip_vec3 = [ 0.0 for j in range(3) ]
        self.num_special_kpt = 0
        self.special_kpt_id = [];   self.special_kpt_name = []

        band_in = open( band_file, "r" );
        self.read_header_band_energy( band_in )

        self.kfrac = [ [0.0 for j in range(3) ] for i in range(self.nkpt) ]
        self.kcart = [ [0.0 for j in range(3) ] for i in range(self.nkpt) ]
        self.bands = [ [0.0 for j in range(self.nband) ] \
                       for i in range(self.nkpt) ]

        self.read_body_band_energy( band_in )
        band_in.close()

        self.dk_list = [];       self.dk_max = 0.0
       
        self.read_special_kpt( kpt_file )
        self.set_kpt_cart()
        self.set_dk_list()
        self.conv_unit( unit )

        lband_weight_in = open( lband_weight_file, "r" );
        self.read_header_lband_weight( lband_weight_in )

        if self.natom > 0:
            self.weight_atom = [ [ [ 0.0 for k in range(self.natom) ] \
                                   for j in range(self.nband) ] \
                                 for i in range(self.nkpt) ] 
        if self.nlayer > 0:
            self.weight_layer = [ [ [ 0.0 for k in range(self.nlayer) ] \
                                    for j in range(self.nband) ] \
                                  for i in range(self.nkpt) ] 

        self.read_body_lband_weight( lband_weight_in )
        lband_weight_in.close()

        if self.natom > 0:
            self.factor_atom  = [ 0 for i in range(self.natom) ]
        if self.nlayer > 0:
            self.factor_layer = [ 0 for i in range(self.nlayer) ]

        self.weight_plot = [ [0.0 for j in range(self.nband)] for i in range(self.nkpt)]
        self.set_factor( atom_id, layer_id )

        out_data = out_file +".dat";    out_gnu  = out_file +".gnu"

        e_origin = self.fermi
        if e_range == "none":
            emin = 1.0E10;  emax = -1.0E10
            for i in range(self.nkpt):
                for j in range(self.nband):
                    emin = min( emin, self.bands[i][j] )
                    emax = max( emax, self.bands[i][j] )
            emin = emin -e_origin
            emax = emax -e_origin
        else:
            emin = float( e_range[0] );  emax = float( e_range[1] )
            
        self.write_weight_data( out_data, unit, e_origin, atom_id, layer_id )

        self.write_gnu_file( out_file, out_gnu, out_data,
                             emin, emax, e_inc, unit, plot_style, circle_scale,
                             cb_range, line_width, window_scale, vertical,
                             with_fermi, no_dispersion, threshold, fig_format )
        self.run_gnuplot( out_gnu )

    def run_gnuplot( self, out_gnu ):
        cmd = 'gnuplot '+out_gnu
        tokens = shlex.split(cmd)
        subprocess.run(tokens)

    def write_gnu_file( self, outfile, out_gnu, out_data,
                        emin, emax, e_inc, unit, plot_style, circle_scale,
                        cb_range, line_width, window_scale, vertical,
                        with_fermi, no_dispersion, threshold, fig_format ):

        out_f = open( out_gnu, "w" );
        out_f.write("set nokey \n");

        if unit == "eV":
            out_f.write("set ylabel 'Energy (eV)' \n")
        elif unit == "Ry":
            out_f.write("set ylabel 'Energy (Ry)' \n")
        elif unit == "Ha":
            out_f.write("set ylabel 'Energy (Ha)' \n")

        out_f.write(f"circle_scale={circle_scale} \n")
        out_f.write(f"threshold={threshold} \n")

        fname = outfile +'.' +fig_format

        if window_scale == "none":
            wscale_x = 1.0;  wscale_y = 1.0
        else:
            wscale_x = float( window_scale[0]);  wscale_y = float( window_scale[1] )

        if fig_format == "eps":
            w1 = 5;  w2 = 3.5
            out_f.write(f"wsx={wscale_x} \n")
            out_f.write(f"wsy={wscale_y} \n")
            out_f.write(f"set terminal postscript eps enhanced color solid size")
            out_f.write(f" {w1}*wsx, {w2}*wsy \n")
#            out_f.write(f"set terminal postscript eps enhanced color solid \n")
            out_f.write(f"set out '{fname}' \n")

        elif fig_format == "png":
            w1 = 640;  w2 = 480
            out_f.write(f"wsx={wscale_x} \n")
            out_f.write(f"wsy={wscale_y} \n")
            out_f.write(f"set terminal pngcairo enhanced size {w1}*wsx, {w2}*wsy \n")
#            out_f.write("set terminal pngcairo enhanced \n")
            out_f.write(f"set out '{fname}' \n")
        else:
            pass

        out_f.write("set palette rgb 22,13,-31 \n");
        out_f.write("set style line 100 lc palette \n")

        out_f.write(f"emin={emin}\n" )
        out_f.write(f"emax={emax}\n" )

        out_f.write("{:s}{:8.5f}{:s}{:8.5f}{:s}\n".format(
            "set xr [",0.0,":",self.dk_max,"]") )

        out_f.write(f"set yr [emin:emax]\n" )
        if e_inc != "none":
            e1 = float( e_inc[0] )
            out_f.write(f"set ytics {e1}\n" )

        out_f.write("set gr front \n");

        if with_fermi == True:
            out_f.write(f"set arrow front from 0.0,0.0 ")
            out_f.write(f"to {self.dk_max},{0.0} nohead lt -1 lw 0.4 \n")

        for i in range(self.num_special_kpt):
            j = self.special_kpt_id[i]
            out_f.write(f"set arrow front from {self.dk_list[j]},emin ")
            out_f.write(f"to {self.dk_list[j]},emax nohead lt -1 \n")

        out_f.write("{:s}".format("set xtics (") )
        for i in range(self.num_special_kpt):
            j = self.special_kpt_id[i]
            out_f.write("{:s}{:s}{:s}{:8.5f}".format(
                '"', self.special_kpt_name[i], '"' ,self.dk_list[j] ) )

            if i < self.num_special_kpt-1:
                out_f.write("{:s}".format(",") )
            else:
                out_f.write("{:s}\n".format(")") )

        if fig_format == "eps":
            transp = 0.8
        else:
            transp = 0.6

        lw = line_width
#        out_f.write(f"set size {window_width},1.0 \n")

        if vertical == False:
            lay_x = 1;    lay_y = 2
        else:
            lay_x = 2;    lay_y = 1

        if cb_range == "none":
            cbmin = 0.0;  cbmax = 0.5
        else:
            cbmin = float( cb_range[0] );  cbmax = float( cb_range[1] )

        if self.nspin == 1:
            if plot_style == 1:
                ctmp = 0.01
                out_f.write(f"plot '{out_data}' u 1:2:({circle_scale}*{ctmp}*sqrt($3**2)) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write( " lc rgb 'web-blue'")
            elif plot_style == 2:
                ctmp = 0.005
#                out_f.write(f"plot '{out_data}' u 1:2:({circle_scale}*{ctmp}):($3) \\\n")
                out_f.write(f"plot '{out_data}' u 1:($3>=threshold ? $2:1/0):(circle_scale*{ctmp}):($3) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write(" ls 100")
            elif plot_style == 3:
                ctmp = 0.5 /( cbmax -cbmin )
                ctmp = 0.01
#                out_f.write(f"plot '{out_data}' u 1:2:({circle_scale}*{ctmp}*$3):($3) \\\n")
                out_f.write(f"plot '{out_data}' u 1:($3>=threshold ? $2:1/0):(circle_scale*{ctmp}*(sqrt(($3-threshold)**2))):($3) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write(" ls 100")

            if no_dispersion == False:
                out_f.write( ", \\\n")
                out_f.write(f"     '{out_data}' u 1:2 w l lw 0.2 lc 'black' \n")
            else:
                out_f.write( "\n")

        elif self.nspin == 2:
            out_f.write(f"set multiplot layout {lay_x},{lay_y} \n")

            if plot_style == 1:
                ctmp = 0.01
                out_f.write(f"plot '{out_data}' u 1:2:({circle_scale}*{ctmp}*sqrt($3**2)) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write( " lc rgb 'web-blue'")

            elif plot_style == 2:
                ctmp = 0.01
                out_f.write(f"plot '{out_data}' u 1:($3>=threshold ? $2:1/0):(circle_scale*{ctmp}):($3) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write(" ls 100")

            elif plot_style == 3:
                ctmp = 0.5 /( cbmax -cbmin )
                ctmp = 0.01
                out_f.write(f"plot '{out_data}' u 1:($3>=threshold ? $2:1/0):(circle_scale*{ctmp}*(sqrt(($3-threshold)**2))):($3) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write(" ls 100")

            if no_dispersion == False:
                out_f.write(", \\\n")
                out_f.write(f"     '{out_data}' u 1:2 w l lw 0.2 lc 'black' \n")
            else:
                out_f.write("\n")
                
            if plot_style == 1:
                ctmp = 0.01
                out_f.write(f"plot '{out_data}' u 1:4:({circle_scale}*{ctmp}*sqrt($5**2)) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write( " lc rgb 'web-blue'")

            elif plot_style == 2:
                ctmp = 0.01
                out_f.write(f"plot '{out_data}' u 1:($5>=threshold ? $4:1/0):(circle_scale*{ctmp}):($5) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write(" ls 100")

            elif plot_style == 3:
                ctmp = 0.5 /( cbmax -cbmin )
                ctmp = 0.01
                out_f.write(f"plot '{out_data}' u 1:($5>=threshold ? $4:1/0):(circle_scale*{ctmp}*(sqrt(($5-threshold)**2))):($5) \\\n")
                out_f.write( "     ");
                out_f.write(f" w circles fill transparent solid {transp} noborder")
                out_f.write(" ls 100")

            if no_dispersion == False:
                out_f.write( ", \\\n")
                out_f.write(f"     '{out_data}' u 1:4 w l lw 0.2 lc 'black' \n")
            else:
                out_f.write("\n")
                
            out_f.write("unset multiplot \n")

        out_f.close()

    def write_weight_data( self, outfile, unit, e_origin, atom_id, layer_id ):
        out_f = open( outfile, "w" );

        if self.nspin == 1:
            if unit == "eV":
                out_f.write("{:s}\n".format("#    dk[Bohr-1]     energy[eV]      weight"))
            elif unit == "Ry":
                out_f.write("{:s}\n".format("#    dk[Bohr-1]     energy[Ry]      weight"))
            elif unit == "Ha":
                out_f.write("{:s}\n".format("#    dk[Bohr-1]     energy[Ha]      weight"))

            for j in range( self.nband ):
                for i in range( self.nkpt ):
                    c1 = self.bands[i][j] -e_origin
                    csum = 0.0
                    if atom_id != 'none':
                        for n in range (self.natom):
                            if self.factor_atom[n] > 0:
                                csum = csum +self.weight_atom[i][j][n]
                    elif layer_id != 'none':
                        for n in range ( self.nlayer ):
                            if self.factor_layer[n] > 0:
                                csum = csum +self.weight_layer[i][j][n]
                    out_f.write("{:15.8f}{:15.8f}{:15.8f}\n".format(self.dk_list[i], 
                                                                    c1, csum ) )
                out_f.write("\n")

        else:
            out_f.write("{:s}{:s}{:s}\n".format("#                      ",
                                                "      up spin              ",
                                                "      down spin"))
            if unit == "eV":
                out_f.write("{:s}{:s}\n".format(
                    "#   dk[Bohr-1]      energy[eV]     weight",
                    "         energy[eV]      weight")  )
            elif unit == "Ry":
                out_f.write("{:s}{:s}\n".format(
                    "#   dk[Bohr-1]      energy[Ry]     weight",
                    "         energy[eV]      weight")  )
            elif unit == "Ha":
                out_f.write("{:s}{:s}\n".format(
                    "#   dk[Bohr-1]      energy[Ha]     weight",
                    "         energy[eV]      weight")  )

            nkpt2 = int( self.nkpt /2 )

            for j in range( self.nband ):
                for i in range( nkpt2 ):
                    c1 = self.bands[2*i][j] -e_origin
                    c2 = self.bands[2*i+1][j] -e_origin
                    csum1 = 0.0
                    csum2 = 0.0
                    if atom_id != 'none':
                        for n in range ( self.natom ):
                            if self.factor_atom[n] > 0:
                                csum1 = csum1 +self.weight_atom[2*i][j][n]
                                csum2 = csum2 +self.weight_atom[2*i+1][j][n]
                    elif layer_id != 'none':
                        for n in range ( self.nlayer ):
                            if self.factor_layer[n] > 0:
                                csum1 = csum1 +self.weight_layer[2*i][j][n]
                                csum2 = csum2 +self.weight_layer[2*i+1][j][n]
                    out_f.write(
                        "{:15.8f}{:15.8f}{:15.8f}{:15.8f}{:15.8f}\n".format(
                            self.dk_list[2*i], c1, csum1, c2, csum2 ) )
                out_f.write("\n")

        out_f.close()

    def set_factor( self, atom_id, layer_id ):
        if atom_id != 'none':
            count  = [ 0 for i in range(self.natom) ]
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
            for i in range(self.natom):
                if count[i] > 0:  
                    self.factor_atom[i] = 1

        elif layer_id != 'none':
            count  = [ 0 for i in range(self.nlayer) ]
            for i in range(len(layer_id)):
                itmp = layer_id[i].split('-');   inum = len(itmp)
                if inum == 2:
                    layer_s = int(itmp[0])-1;  layer_e = int(itmp[1])
                    for j in range(layer_s,layer_e):
                        if j < 0 or j > self.nlayer-1:
                            pass
                        else:
                            count[j] = count[j] +1
                else:
                    layer_s = int(itmp[0]) -1;    j = layer_s
                    if j < 0 or j > self.nlayer-1:
                        pass
                    else:
                        count[j] = count[j] +1
            for i in range(self.nlayer):
                if count[i] > 0:  self.factor_layer[i] = 1

    def conv_unit( self, unit ):
        if unit == "eV":
            self.fermi = self.fermi *27.2116
        elif unit == "Ry":
            self.fermi = self.fermi *2.0

        for i in range(self.nkpt):
            for j in range(self.nband):
                if unit == "eV":
                    self.bands[i][j] = self.bands[i][j] *27.2116
                elif unit == "Ry":
                    self.bands[i][j] = self.bands[i][j] *2

    def set_dk_list( self ):
        for i in range( self.nkpt ):
            if i == 0:
                dk = 0.0
            else:
                dx = self.kcart[i][0] -self.kcart[i-1][0]
                dy = self.kcart[i][1] -self.kcart[i-1][1]
                dz = self.kcart[i][2] -self.kcart[i-1][2]
                dk = dk +math.sqrt( dx**2 +dy**2 +dz**2 )
            self.dk_list.append( dk )

        self.dk_max = dk

    def set_kpt_cart( self ):
        for i in range( self.nkpt ):
            for j in range(3):
                self.kcart[i][j] = self.recip_vec1[j] *self.kfrac[i][0] \
                                   +self.recip_vec2[j] *self.kfrac[i][1] \
                                   +self.recip_vec3[j] *self.kfrac[i][2]

    def read_special_kpt( self, kpt_file ):
        kpt_in = open( kpt_file, "r" );
        kpt_in.readline()

        for i in range(3):
            line = kpt_in.readline()
            data = line.split();
            self.recip_vec1[i] = float(data[0]); 
            self.recip_vec2[i] = float(data[1]); 
            self.recip_vec3[i] = float(data[2]); 

        last = 0
        while True:
            line = kpt_in.readline()
            if not line:
                eof_flag = 1;  break;
                               
            data = line.split()
            if data[4] == "#":
                found = 0
                kx = float( data[0] ) /float( data[3] )
                ky = float( data[1] ) /float( data[3] )
                kz = float( data[2] ) /float( data[3] )

                data = re.split('[#\n]',line);
                for i in range( last,self.nkpt ):
                    dx = kx -self.kfrac[i][0]
                    dy = ky -self.kfrac[i][1]
                    dz = kz -self.kfrac[i][2]
                    d2 = dx**2 +dy**2 +dz**2
                    if d2 < 1.0E-6:
                        found = 1
                        self.special_kpt_id.append( i )
                        self.special_kpt_name.append( data[1] )
                        self.num_special_kpt = self.num_special_kpt +1
                        last = i;  break
        kpt_in.close()

    def read_body_lband_weight( self, lband_weight_in ):
        eof_flag = 0;
        atom_flag = 0;   layer_flag = 0

        kcount = 0;   ecount = 0
        atom_count = 0;  layer_count = 0

        while True:
            line = lband_weight_in.readline()
            if not line:
                eof_flag = 1;  break;

            if line.find('for total') > 0:  continue
            if line.find('for magnetic') > 0:  continue

            if atom_flag==1:
                data = line.split();   num_data = len(data);
                if line.find(')') > 0:
                    for j in range(1,num_data):
                        self.weight_atom[kcount][ecount][atom_count] = float(data[j])
                        atom_count = atom_count +1
                else:
                    for j in range(num_data):
                        self.weight_atom[kcount][ecount][atom_count] = float(data[j])
                        atom_count = atom_count +1

                if atom_count == self.natom:  
                    ecount = ecount +1; atom_count = 0

                if ecount == self.nband:  atom_flag = 0

            if layer_flag==1:
                data = line.split();   num_data = len(data);
                if line.find(')') > 0:
                    for j in range(1,num_data):
                        self.weight_layer[kcount][ecount][layer_count] = float(data[j])
                        layer_count = layer_count +1
                else:
                    for j in range(num_data):
                        self.weight_layer[kcount][ecount][layer_count] = float(data[j])
                        layer_count = layer_count +1

                if layer_count == self.nlayer:  
                    ecount = ecount +1; layer_count = 0

                if ecount == self.nband: layer_flag = 0

            if line.find('each atomic') > 0:
                data = line.split();   num_data = len(data);
                if self.ndim_magmom == 4:
                    jj = int( data[num_data-1] )
                    kcount = int( ( jj+1 ) /2 -1 )
                else:
                    kcount = int( data[num_data-1] ) -1
                    
                atom_flag = 1;  atom_count = 0;   ecount = 0

            if line.find('Atom decomposed') > 0:
                data = line.split();   num_data = len(data);
                if self.ndim_magmom == 4:
                    jj = int( data[num_data-1] )
                    kcount = int( ( jj+1 ) /2 -1 )
                else:
                    kcount = int( data[num_data-1] ) -1

                atom_flag = 1;  atom_count = 0;   ecount = 0

            if line.find('each layer') > 0:
                data = line.split();   num_data = len(data);
                if self.ndim_magmom == 4:
                    jj = int( data[num_data-1] )
                    kcount = int( ( jj+1 ) /2 -1 )
                else:
                    kcount = int( data[num_data-1] ) -1

                layer_flag = 1;  layer_count = 0;  ecount = 0

            if line.find('Layer decomposed') > 0:
                data = line.split();   num_data = len(data);
                if self.ndim_magmom == 4:
                    jj = int( data[num_data-1] )
                    kcount = int( ( jj+1 ) /2 -1 )
                else:
                    kcount = int( data[num_data-1] ) -1

                layer_flag = 1;  layer_count = 0;  ecount = 0

    def read_header_lband_weight( self, lband_weight_in ):
        eof_flag = 0;

        while True:
            line = lband_weight_in.readline()
            if not line:
                eof_flag = 1;  break;

            if line.find('ndim_magmom') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.ndim_magmom = int( data[ntarget] )

            if line.find('num_atoms') > 0:
                data = line.split();
                self.natom = int( data[2] )

            if line.find('num_layers') > 0:
                data = line.split();
                self.nlayer = int( data[2] )

            if line.find('Atom Info.') > 0:   break
            if line.find('Layer Info.') > 0:  break

    def read_header_band_energy( self, eband_in ):
        eof_flag = 0;

        while True:
            line = eband_in.readline()
            if not line:
                eof_flag = 1;  break;

            if line.find('num_kpoints') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.nkpt = int( data[ntarget] )

            if line.find('num_bands') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.nband = int( data[ntarget] )

            if line.find('nspin') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.nspin = int( data[ntarget] )

            if line.find('ndim_magmom') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.ndim_magmom = int( data[ntarget] )

            if line.find('band max') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.fermi = float( data[ntarget] )

            if line.find('Fermi energy') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.fermi = float( data[ntarget] )

            if line.find('nk_conv') > 0:  break;
            if line.find('kpoint list') > 0:  break;

    def read_body_band_energy( self, eband_in ):
        eof_flag = 0;    read_flag = 0;
        kcount = 0;     ecount = 0; 

        while True:
            line = eband_in.readline()

            if not line:
                eof_flag = 1;  break;

            if read_flag == 2:
                data = line.split();  num = len(data)

                for j in range(num):
                    c1 = float(data[j])
                    self.bands[kcount][ecount] = c1
                    ecount = ecount +1
                    if ecount == self.nband:
                        kcount = kcount +1; read_flag = 0

            if line.find('energy_eigen_values') > 0 or \
               line.find('energy eigenvalues') > 0:
                read_flag = 1;  ecount = 0

            if read_flag == 1 and line.find('ik =') > 0:
                data = line.split();  num = len(data)
                self.kfrac[kcount][0] = float(data[4]);
                self.kfrac[kcount][1] = float(data[5]);
                self.kfrac[kcount][2] = float(data[6]);

                read_flag = read_flag +1;


###################################### Main ########
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='plot lband')

    parser.add_argument( 'band_file', help='nfenergy.data' )
    parser.add_argument( 'kpt_file', 
                         help='k-point generation control file (bandkpt.in)' )
    parser.add_argument( 'lband_weight_file', help='wfn_local_decomp.data' )

    parser.add_argument( '--atom_id', help='list of atom ids',
                         default='none', nargs='*' )
    parser.add_argument( '--layer_id', help='list of layer ids',
                         default='none', nargs='*' )

    parser.add_argument( '--e_range', help='energy range (min and max)',
                         default='none', nargs=2 )
    parser.add_argument( '--e_inc', help='energy increment',
                         default='none', nargs=1 )

    parser.add_argument( '--unit',
                         help='unit (eV, Ry, Ha) [default:eV]',
                         default='eV' )
    parser.add_argument( '--plot_style',
                         help='plot style (1,2,3) (default:1)',
                         type=int, default=1 )
    parser.add_argument( '--circle_scale',
                         help='circle scale for plotting (default:1.0)',
                         type=float, default=1.0 )

    parser.add_argument( '--cb_range', help='color range (min and max)',
                         default='none', nargs=2 )

    parser.add_argument( '--line_width', help='window width [default:1.0]',
                         type=float, default=1.0 )
    parser.add_argument( '--window_scale', help='window scale (horizontal and vertical)',
                         default='none', nargs=2 )
    parser.add_argument( '--vertical', help='action: layout of graphs when nspin=2 vertically [default:False]',
                         action='store_true', default=False )

    parser.add_argument( '--fig_format', 
                         help='figure format (eps,png) [default:eps]',
                         default='eps' )
    parser.add_argument( '--out_file', 
                         help='output filename [default:plot_lband]',
                         default='plot_lband' )

    parser.add_argument( '--with_fermi', 
                         help='action: plot line at 0.0 eV [default:False]',
                         action='store_true', default=False )

    parser.add_argument( '--no_dispersion', 
                         help='action: dispersion is off [default:False]',
                         action='store_true', default=False )

    parser.add_argument( '--threshold',
                         help='threshold when plot_style==3 (default:0.0)',
                         type=float, default=0.0 )

    args = parser.parse_args()

    plot_lband( args.band_file, args.kpt_file, args.lband_weight_file,
                args.atom_id, args.layer_id,
                args.e_range, args.e_inc, args.unit,
                args.plot_style, args.circle_scale, args.cb_range,
                args.line_width, args.window_scale, args.vertical,
                args.fig_format, args.out_file, 
                args.with_fermi, args.no_dispersion, args.threshold )
