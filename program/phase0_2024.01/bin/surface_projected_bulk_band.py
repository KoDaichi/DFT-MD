#!/usr/bin/env python3
#
##################################################
#
# surface_projected_bulk_band.py
#   version : 1.01
#
#   contact address :  Phase System Consortium
#
##################################################

import sys, math, re;
import argparse, subprocess, shlex

class surface_projected_bulk_band:
    def __init__( self, bulk_band_file, kpt_file,
                  e_range, e_inc, unit, 
                  plot_style, cb_range, line_width, window_scale, vertical,
                  fig_format, outfile, force_normal, with_fermi,
                  surf_band_file, e_shift,
                  ndiv_erange_map, broadening_width_map, spin_sum_map ):

        self.nkpt_bulk = 0;       self.nband_bulk = 0;     self.nspin_bulk = 0;
        self.fermi_bulk = 0.0;

        self.recip_vec1 = [ 0.0 for j in range(3) ]
        self.recip_vec2 = [ 0.0 for j in range(3) ]
        self.recip_vec3 = [ 0.0 for j in range(3) ]
        self.num_special_kpt = 0
        self.special_kpt_id = [];   self.special_kpt_name = []

        bulk_band_in = open( bulk_band_file, "r" );
        self.read_header_bulk_band( bulk_band_in )

        self.bands_bulk = [ [0.0 for j in range(self.nband_bulk) ] \
                            for i in range(self.nkpt_bulk) ]
        self.kfrac_bulk = [ [0.0 for j in range(3) ] for i in range(self.nkpt_bulk) ]
        self.read_body_bulk_band( bulk_band_in )
        bulk_band_in.close()

        self.determine_nkz()
        self.nband_2d = int( self.nband_bulk *self.nkz )
        self.nkpt_2d  = int( self.nkpt_bulk /self.nkz )

        self.bands_2d = [ [0.0 for j in range(self.nband_2d) ] \
                            for i in range(self.nkpt_2d) ]
        self.kfrac_z  = [ 0.0 for i in range(self.nband_2d) ]
        self.kfrac_2d = [ [0.0 for j in range(2) ] for i in range(self.nkpt_2d) ]
        self.kcart_2d = [ [0.0 for j in range(2) ] for i in range(self.nkpt_2d) ]
        self.dk_list_2d = [];       self.dk_max_2d = 0.0

        self.conv_bulk_to_2d( self.nkz )
        
        self.read_special_kpt( kpt_file )
        self.check_kz_normal_to_surf( force_normal )

        self.set_kpt_cart_2d()
        self.set_dk_list_2d()
        self.conv_unit_bulk( unit )

        out_data_bulk = outfile +"_bulk.dat";    out_gnu  = outfile +".gnu"
        out_data_surf = outfile +"_surf.dat";    out_map = outfile +'_bulk.map'

        if surf_band_file == "none":
            with_surf_mode = 0

            e_origin_bulk = self.fermi_bulk
            self.write_band_data_2d( out_data_bulk, unit, e_origin_bulk )

            if e_range == "none":
                emin = 1.0E10;  emax = -1.0E10
                for i in range(self.nkpt_bulk):
                    for j in range(self.nband_bulk):
                        emin = min( emin, self.bands_bulk[i][j] )
                        emax = max( emax, self.bands_bulk[i][j] )
                emin = emin -e_origin_bulk
                emax = emax -e_origin_bulk
            else:
                emin = float( e_range[0] );  emax = float( e_range[1] )

        else:
            with_surf_mode = 1

            surf_band_in = open( surf_band_file, "r" );
            self.read_header_surf_band( surf_band_in )
            self.bands_surf = [ [0.0 for j in range(self.nband_surf) ] \
                                for i in range(self.nkpt_surf) ]
            self.kfrac_surf = [ [0.0 for j in range(2) ] for i in range(self.nkpt_surf) ]
            self.kcart_surf = [ [0.0 for j in range(2) ] for i in range(self.nkpt_surf) ]
            self.dk_list_surf = [];       self.dk_max_surf = 0.0

            self.read_body_surf_band( surf_band_in )
            surf_band_in.close()

            self.set_kpt_cart_surf()
            self.set_dk_list_surf()
            self.conv_unit_surf( unit )

            e_origin_surf = self.fermi_surf
            e_origin_bulk = self.fermi_surf -e_shift
            self.write_band_data_surf( out_data_surf, unit, e_origin_surf )
            self.write_band_data_2d( out_data_bulk, unit, e_origin_bulk )

            if e_range == "none":
                emin = 1.0E10;  emax = -1.0E10
                for i in range(self.nkpt_2d):
                    for j in range(self.nband_2d):
                        emin = min( emin, self.bands_2d[i][j] )
                        emax = max( emax, self.bands_2d[i][j] )
                emin = emin -e_origin_bulk;   emax = emax -e_origin_bulk
            else:
                emin = float( e_range[0] );  emax = float( e_range[1] )

            self.set_ref_map( self.nkz, unit, e_origin_bulk, emin, emax, out_map, 
                              ndiv_erange_map, broadening_width_map,
                              spin_sum_map )

        self.write_gnu_file( outfile, with_surf_mode,
                             out_gnu, out_data_bulk, out_data_surf, out_map,
                             spin_sum_map,
                             emin, emax, e_inc, unit, plot_style, 
                             cb_range, line_width, window_scale, vertical,
                             with_fermi, fig_format )

        self.run_gnuplot( out_gnu )

    def check_kz_normal_to_surf( self, force_normal ):
        c1 = math.sqrt( self.recip_vec3[0]**2 +self.recip_vec3[1]**2 )
        if c1 > 1E-8:
            print("Warning")
            print(" The 3rd reciprocal vector is not normal to surface")
            if force_normal == False:
                print(" The program abnormally exits.")
                exit(-1)
    
    def conv_bulk_to_2d( self, nkz ):
        ik = 0
        if self.nspin_bulk == 1:
            for i in range(self.nkpt_2d):
                ib = 0
                for n in range(nkz):
                    for j in range(self.nband_bulk):
                        self.bands_2d[i][ib] = self.bands_bulk[ik][j]
                        self.kfrac_z[ib]     = self.kfrac_bulk[ik][2]
                        ib = ib +1
                    self.kfrac_2d[i][0] = self.kfrac_bulk[ik][0]
                    self.kfrac_2d[i][1] = self.kfrac_bulk[ik][1]

                    ik = ik +1
        else:
            nkpt = int( self.nkpt_2d /2 )
            for i in range(nkpt):
                ib = 0
                for n in range(nkz):
                    for j in range(self.nband_bulk):
                        self.bands_2d[2*i][ib]   = self.bands_bulk[2*ik][j]
                        self.bands_2d[2*i+1][ib] = self.bands_bulk[2*ik+1][j]
                        self.kfrac_z[ib]         = self.kfrac_bulk[2*ik][2]
                        ib = ib +1
                    self.kfrac_2d[2*i][0]   = self.kfrac_bulk[2*ik][0]
                    self.kfrac_2d[2*i][1]   = self.kfrac_bulk[2*ik][1]
                    self.kfrac_2d[2*i+1][0] = self.kfrac_bulk[2*ik+1][0]
                    self.kfrac_2d[2*i+1][1] = self.kfrac_bulk[2*ik+1][1]
                    
                    ik = ik +1

    def conv_unit_bulk( self, unit ):
        if unit == "eV":
            self.fermi_bulk = self.fermi_bulk *27.2116
        elif unit == "Ry":
            self.fermi_bulk = self.fermi_bulk *2.0

        for i in range(self.nkpt_bulk):
            for j in range(self.nband_bulk):
                if unit == "eV":
                    self.bands_bulk[i][j] = self.bands_bulk[i][j] *27.2116
                elif unit == "Ry":
                    self.bands_bulk[i][j] = self.bands_bulk[i][j] *2
        for i in range(self.nkpt_2d):
            for j in range(self.nband_2d):
                if unit == "eV":
                    self.bands_2d[i][j] = self.bands_2d[i][j] *27.2116
                elif unit == "Ry":
                    self.bands_2d[i][j] = self.bands_2d[i][j] *2

    def conv_unit_surf( self, unit ):
        if unit == "eV":
            self.fermi_surf = self.fermi_surf *27.2116
        elif unit == "Ry":
            self.fermi_surf = self.fermi_surf *2.0

        for i in range(self.nkpt_surf):
            for j in range(self.nband_surf):
                if unit == "eV":
                    self.bands_surf[i][j] = self.bands_surf[i][j] *27.2116
                elif unit == "Ry":
                    self.bands_surf[i][j] = self.bands_surf[i][j] *2
            
    def determine_nkz( self ):
        kx0 = self.kfrac_bulk[0][0]
        ky0 = self.kfrac_bulk[0][1]
        for i in range( self.nkpt_bulk ):
            dx = self.kfrac_bulk[i][0] -kx0
            dy = self.kfrac_bulk[i][1] -ky0
            d2 = dx**2 +dy**2
            if d2 > 1.0E-6:
                break
        self.nkz = int( i /self.nspin_bulk )

            
    def set_ref_map( self, nkz, unit, eshift, emin, emax, out_map, 
                     ndiv_erange_map, broadening_width_map, spin_sum_map ):
        count = 0;  dk = 0.0

        if ndiv_erange_map == "none":
            de = ( emax -emin ) /1000.0;
            ne_div = int( ( emax -emin )/de ) +1
        else:
            ne_div = int(ndiv_erange_map)
            de = ( emax -emin ) /ne_div

        ovp = broadening_width_map
        
        itmp = 0
        spectr = [ [0.0 for j in range(ne_div)] for i in range(self.nkpt_2d)]

        if self.nspin_bulk == 1:
            for i in range( self.nkpt_2d ):
                emin_wk = [  1.0E10 for j in range(self.nband_bulk) ]
                emax_wk = [ -1.0E10 for j in range(self.nband_bulk) ]

                for k in range( nkz ):
                    itmp = nkz *i +k
                    for j in range( self.nband_bulk ):
                        e1 = self.bands_bulk[itmp][j] -eshift
                        c1 = min( emin_wk[j], e1 )
                        c2 = max( emax_wk[j], e1 )
                        emin_wk[j] = c1;  emax_wk[j] = c2
                        
                for n in range(ne_div):
                    count = 0
                    for j in range( self.nband_bulk ):
                        e1 = emin +de *n
                        if e1 >= emin_wk[j]-de*ovp and e1 < emax_wk[j]+de*ovp:
                            count = count +1
                    if count > 0:
                        spectr[i][n] = 1
        else:
            nkpt = int( self.nkpt_2d /2 )
            for ispin in range(2):
                for i in range( nkpt ):
                    ik = 2*i +ispin
                    emin_wk = [  1.0E10 for j in range(self.nband_bulk) ]
                    emax_wk = [ -1.0E10 for j in range(self.nband_bulk) ]

                    for k in range( nkz ):
                        itmp = ( nkz *i +k )*2 +ispin
                        for j in range( self.nband_bulk ):
                            e1 = self.bands_bulk[itmp][j] -eshift
                            c1 = min( emin_wk[j], e1 )
                            c2 = max( emax_wk[j], e1 )
                            emin_wk[j] = c1;  emax_wk[j] = c2
                        
                    for n in range(ne_div):
                        count = 0
                        for j in range( self.nband_bulk ):
                            e1 = emin +de *n
                            if e1 >= emin_wk[j]-de*ovp and e1 < emax_wk[j]+de*ovp:
                                count = count +1
                        if count > 0:
                            spectr[ik][n] = 1

        out_f = open( out_map, "w" );
        if unit == "eV":
            out_f.write("{:s}\n".format("#      dk           energy[eV]"))
        elif unit == "Ry":
            out_f.write("{:s}\n".format("#      dk           energy[Ry]"))
        elif unit == "Ha":
            out_f.write("{:s}\n".format("#      dk           energy[Ha]"))

        if self.nspin_bulk == 1:
            for i in range(self.nkpt_2d):
                dk = self.dk_list_2d[i]
                for j in range(ne_div):
                    e1 = emin +de *j
                    c1 = spectr[i][j]
                    out_f.write("{:15.8f}{:15.8f}{:15.8f}\n".format(dk,e1,c1) )
                out_f.write("\n")
        else:
            nkpt_wk = int( self.nkpt_2d /2 )
            for i in range(nkpt_wk):
                dk = self.dk_list_2d[2*i]
                if spin_sum_map == True:
                    for j in range(ne_div):
                        e1 = emin +de *j
                        c1 = ( spectr[2*i][j] +spectr[2*i+1][j] )/2.0
                        out_f.write("{:15.8f}{:15.8f}{:15.8f}\n".format(dk,e1,c1) )
                    out_f.write("\n")
                else:
                    for j in range(ne_div):
                        e1 = emin +de *j
                        c1 = spectr[2*i][j]
                        c2 = spectr[2*i+1][j] 
                        out_f.write("{:15.8f}{:15.8f}{:15.8f}{:15.8f}\n".\
                                    format(dk,e1,c1,c2) )
                    out_f.write("\n")

            
        out_f.close()

    def run_gnuplot( self, out_gnu ):
        cmd = 'gnuplot '+out_gnu
        tokens = shlex.split(cmd)
        subprocess.run(tokens)

    def write_gnu_file( self, outfile, with_surf_mode, 
                        out_gnu, out_data_bulk, out_data_surf, out_map, 
                        spin_sum_map,
                        emin, emax, e_inc, unit, plot_style,
                        cb_range, line_width, window_scale, vertical,
                        with_fermi, fig_format ):

        out_f = open( out_gnu, "w" );
        out_f.write("set nokey \n");

        if unit == "eV":
            out_f.write("set ylabel 'Energy (eV)' \n")
        elif unit == "Ry":
            out_f.write("set ylabel 'Energy (Ry)' \n")
        elif unit == "Ha":
            out_f.write("set ylabel 'Energy (Ha)' \n")

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

        if with_surf_mode == 1:
            out_f.write("set palette gray \n")

        out_f.write("{:s}{:8.5f}{:s}{:8.5f}{:s}\n".format(
            "set xr [",0.0,":",self.dk_max_2d,"]") )

        out_f.write(f"emin={emin}\n" )
        out_f.write(f"emax={emax}\n" )

        out_f.write(f"set yr [emin:emax]\n" )
        if e_inc != "none":
            e1 = float( e_inc[0] )
            out_f.write(f"set ytics {e1}\n" )

        out_f.write("set gr front \n");

        if with_fermi == True:
            out_f.write(f"set arrow front from 0.0,0.0 ")
            out_f.write(f"to {self.dk_max_2d},{0.0} nohead lt -1 lw 0.4 \n")
            
        for i in range(self.num_special_kpt):
            j = self.special_kpt_id[i]
            out_f.write(f"set arrow front from {self.dk_list_2d[j]},emin ")
            out_f.write(f"to {self.dk_list_2d[j]},emax nohead lt -1 \n")

        out_f.write("{:s}".format("set xtics (") )
        for i in range(self.num_special_kpt):
            j = self.special_kpt_id[i]
            out_f.write("{:s}{:s}{:s}{:8.5f}".format(
                '"', self.special_kpt_name[i], '"' ,self.dk_list_2d[j] ) )

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
            
        if with_surf_mode == 0: 
            if plot_style == 1:
                out_f.write(f"set style line 100 lt 7 lw {lw} lc rgb 'web-blue' \n")
                if self.nspin_bulk == 1:
                    out_f.write(f"plot '{out_data_bulk}' u 1:3 w l ls 100 \n")
                else:
#                    out_f.write(f"set multiplot layout 1,2 scale {window_width}, 1.0 \n")
                    out_f.write(f"set multiplot layout {lay_x},{lay_y} \n")
                    out_f.write(f"plot '{out_data_bulk}' u 1:3 w l ls 100 \n")
                    out_f.write(f"plot '{out_data_bulk}' u 1:4 w l ls 100 \n")
                    out_f.write("unset multiplot \n")

            elif plot_style == 2:
                if cb_range == "none":
                    cbmin = 0.0;  cbmax = 0.5
                else:
                    cbmin = float( cb_range[0] );  cbmax = float( cb_range[1] )

                out_f.write("set cblabel \"|kz|\" \n");
                out_f.write(f"set cbrange [{cbmin}:{cbmax}]\n" )
                out_f.write("set palette rgb 22,13,-31 \n");
                out_f.write(f"set style line 100 lc palette lw {lw} \n")

                if self.nspin_bulk == 1:
                    out_f.write(f"plot '{out_data_bulk}' u 1:3:(sqrt($2**2)) w l ls 100 \n")
                else:
##                    out_f.write(f"set multiplot layout 1,2 scale {window_width}, 1.0 \n")
                    out_f.write(f"set multiplot layout {lay_x},{lay_y} \n")
                    out_f.write(f"plot '{out_data_bulk}' u 1:3:(sqrt($2**2)) w l ls 100 \n")
                    out_f.write(f"plot '{out_data_bulk}' u 1:4:(sqrt($2**2)) w l ls 100 \n")
                    out_f.write("unset multiplot \n")

        else:
            out_f.write(f"set cbrange [0:1] \n" )
            out_f.write(f"unset colorbox \n" )
            out_f.write(f"set style line 100 lt 7 lw {lw} lc rgb 'web-blue' \n")
            if self.nspin_surf == 1:
                if self.nspin_bulk == 2 and spin_sum_map == False:
##                    out_f.write(f"set multiplot layout 1,2 scale {window_width}, 1.0 \n")
                    out_f.write(f"set multiplot layout {lay_x},{lay_y} \n")
                    out_f.write(f"plot '{out_map}' u 1:2:(1-0.2*$3) w image, \\\n")
                    out_f.write(f"     '{out_data_surf}' u 1:2 w l ls 100 \n")
                    out_f.write(f"plot '{out_map}' u 1:2:(1-0.2*$4) w image, \\\n")
                    out_f.write(f"     '{out_data_surf}' u 1:2 w l ls 100 \n")
                    out_f.write("unset multiplot \n")
                else:
                    out_f.write(f"plot '{out_map}' u 1:2:(1-0.2*$3) w image, \\\n")
                    out_f.write(f"     '{out_data_surf}' u 1:2 w l ls 100 \n")
            else:
#                out_f.write(f"set multiplot layout 1,2 scale {window_width}, 1.0 \n")
                out_f.write(f"set multiplot layout {lay_x},{lay_y} \n")
                if self.nspin_bulk == 2 and spin_sum_map == False:
                    out_f.write(f"plot '{out_map}' u 1:2:(1-0.2*$3) w image, \\\n")
                    out_f.write(f"     '{out_data_surf}' u 1:2 w l ls 100 \n")
                    out_f.write(f"plot '{out_map}' u 1:2:(1-0.2*$4) w image, \\\n")
                    out_f.write(f"     '{out_data_surf}' u 1:3 w l ls 100 \n")
                else:
                    out_f.write(f"plot '{out_map}' u 1:2:(1-0.2*$3) w image, \\\n")
                    out_f.write(f"     '{out_data_surf}' u 1:2 w l ls 100 \n")
                    out_f.write(f"plot '{out_map}' u 1:2:(1-0.2*$3) w image, \\\n")
                    out_f.write(f"     '{out_data_surf}' u 1:3 w l ls 100 \n")

                out_f.write("unset multiplot \n")


        out_f.close()

    def set_dk_list_2d( self ):
        for i in range( self.nkpt_2d ):
            if i == 0:
                dk = 0.0
            else:
                dx = self.kcart_2d[i][0] -self.kcart_2d[i-1][0]
                dy = self.kcart_2d[i][1] -self.kcart_2d[i-1][1]
                dk = dk +math.sqrt( dx**2 +dy**2 )
            self.dk_list_2d.append( dk )

        self.dk_max_2d = dk

    def set_dk_list_surf( self ):
        for i in range( self.nkpt_surf ):
            if i == 0:
                dk = 0.0
            else:
                dx = self.kcart_surf[i][0] -self.kcart_surf[i-1][0]
                dy = self.kcart_surf[i][1] -self.kcart_surf[i-1][1]
#                dz = self.kcart_surf[i][2] -self.kcart_surf[i-1][2]
#                dk = dk +math.sqrt( dx**2 +dy**2 +dz**2 )
                dk = dk +math.sqrt( dx**2 +dy**2 )
            self.dk_list_surf.append( dk )

        self.dk_max_surf = dk

    def write_band_data_2d( self, outfile, unit, e_origin ):
        out_f = open( outfile, "w" );

        if self.nspin_bulk == 1:
            if unit == "eV":
                out_f.write("{:s}\n".\
                            format("#    dk[Bohr-1]         kz         energy[eV]"))
            elif unit == "Ry":
                out_f.write("{:s}\n".\
                            format("#    dk[Bohr-1]         kz         energy[Ry]"))
            elif unit == "Ha":
                out_f.write("{:s}\n".\
                            format("#    dk[Bohr-1]         kz         energy[Ha]"))

            for j in range( self.nband_2d ):
                for i in range( self.nkpt_2d ):
                    kz = self.kfrac_z[j]
                    c1 = self.bands_2d[i][j] -e_origin
                    out_f.write("{:15.8f}{:15.8f}{:15.8f}\n".\
                                format(self.dk_list_2d[i], kz, c1 ) )
                out_f.write("\n")

        else:
            if unit == "eV":
                out_f.write("{:s}\n".\
                            format("#    dk[Bohr-1]         kz         energy[eV]"))
            elif unit == "Ry":
                out_f.write("{:s}\n".\
                            format("#    dk[Bohr-1]      kz      energy[Ry]"))
            elif unit == "Ha":
                out_f.write("{:s}\n".\
                            format("#    dk[Bohr-1]      kz      energy[Ha]"))

            out_f.write("{:38s}{:s}\n".format("#","up            down"))
            nkpt = int( self.nkpt_2d /2 )

            for j in range( self.nband_2d ):
                for i in range( nkpt ):
                    kz = self.kfrac_z[j]
                    c1 = self.bands_2d[2*i][j] -e_origin
                    c2 = self.bands_2d[2*i+1][j] -e_origin

                    out_f.write("{:15.8f}{:15.8f}{:15.8f}{:15.8f}\n".format( 
                        self.dk_list_2d[2*i], kz, c1, c2 ) )
                out_f.write("\n")
                        
        out_f.close()

    def write_band_data_surf( self, outfile, unit, e_origin ):
        out_f = open( outfile, "w" );

        if self.nspin_surf == 1:
            if unit == "eV":
                out_f.write("{:s}\n".format("#    dk[Bohr-1]     energy[eV]"))
            elif unit == "Ry":
                out_f.write("{:s}\n".format("#    dk[Bohr-1]     energy[Ry]"))
            elif unit == "Ha":
                out_f.write("{:s}\n".format("#    dk[Bohr-1]     energy[Ha]"))

            for j in range( self.nband_surf ):
                for i in range( self.nkpt_surf ):
                    c1 = self.bands_surf[i][j] -e_origin
                    out_f.write("{:15.8f}{:15.8f}\n".format(self.dk_list_surf[i], c1 ) )
                out_f.write("\n")

        else:
            if unit == "eV":
                out_f.write("{:s}\n".format("#    dk[Bohr-1]            energy[eV]"))
            elif unit == "Ry":
                out_f.write("{:s}\n".format("#    dk[Bohr-1]            energy[Ry]"))
            elif unit == "Ha":
                out_f.write("{:s}\n".format("#    dk[Bohr-1]            energy[Ha]"))

            out_f.write("{:s}\n".format("#                       up            down"))
            nkpt = int( self.nkpt_surf /2 )

            for j in range( self.nband_surf ):
                for i in range( nkpt ):
                    c1 = self.bands_surf[2*i][j] -e_origin
                    c2 = self.bands_surf[2*i+1][j] -e_origin

                    out_f.write("{:15.8f}{:15.8f}{:15.8f}\n".format( 
                        self.dk_list_surf[2*i], c1, c2 ) )
                out_f.write("\n")
                        
        out_f.close()

    def set_kpt_cart_surf( self ):
        for i in range( self.nkpt_surf ):
            for j in range(2):
                self.kcart_surf[i][j] = self.recip_vec1[j] *self.kfrac_surf[i][0] \
                                       +self.recip_vec2[j] *self.kfrac_surf[i][1] 

    def set_kpt_cart_2d( self ):
        for i in range( self.nkpt_2d ):
            for j in range(2):
                self.kcart_2d[i][j] = self.recip_vec1[j] *self.kfrac_2d[i][0] \
                                      +self.recip_vec2[j] *self.kfrac_2d[i][1]

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
                for i in range( last,self.nkpt_2d ):
                    dx = kx -self.kfrac_2d[i][0]
                    dy = ky -self.kfrac_2d[i][1]
                    d2 = dx**2 +dy**2
                    if d2 < 1.0E-6:
                        found = 1
                        self.special_kpt_id.append( i )
                        self.special_kpt_name.append( data[1] )
                        self.num_special_kpt = self.num_special_kpt +1
                        last = i;  break
        kpt_in.close()

    def read_header_surf_band( self, eband_in ):
        eof_flag = 0;

        while True:
            line = eband_in.readline()
            if not line:
                eof_flag = 1;  break;

            if line.find('num_kpoints') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.nkpt_surf = int( data[ntarget] )

            if line.find('num_bands') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.nband_surf = int( data[ntarget] )

            if line.find('nspin') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.nspin_surf = int( data[ntarget] )

            if line.find('band max') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.fermi_surf = float( data[ntarget] )

            if line.find('Fermi energy') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.fermi_surf = float( data[ntarget] )

            if line.find('nk_conv') > 0:  break;

    def read_body_surf_band( self, eband_in ):
        eof_flag = 0;    read_flag = 0;
        kcount = 0;     ecount = 0

        while True:
            line = eband_in.readline()
            if not line:
                eof_flag = 1;  break;

            if read_flag == 2:
                data = line.split();  num = len(data)
                for j in range(num):
                    c1 = float(data[j])
                    self.bands_surf[kcount][ecount] = c1

                    ecount = ecount +1
                    if ecount == self.nband_surf:
                        kcount = kcount +1; read_flag = 0

            if line.find('energy_eigen_values') > 0 or \
               line.find('energy eigenvalues') > 0:
                read_flag = 1;  ecount = 0

            if read_flag == 1 and line.find('ik =') > 0:
                data = line.split();  num = len(data)
                self.kfrac_surf[kcount][0] = float(data[4]);
                self.kfrac_surf[kcount][1] = float(data[5]);
#                self.kfrac_surf[kcount][2] = float(data[6]);

                read_flag = read_flag +1;

    def read_header_bulk_band( self, eband_in ):
        eof_flag = 0;

        while True:
            line = eband_in.readline()
            if not line:
                eof_flag = 1;  break;

            if line.find('num_kpoints') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.nkpt_bulk = int( data[ntarget] )

            if line.find('num_bands') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.nband_bulk = int( data[ntarget] )

            if line.find('nspin') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.nspin_bulk = int( data[ntarget] )

            if line.find('band max') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.fermi_bulk = float( data[ntarget] )

            if line.find('Fermi energy') > 0:
                data = line.split();   ntarget = len(data)-1;
                self.fermi_bulk = float( data[ntarget] )

            if line.find('nk_conv') > 0:  break;

    def read_body_bulk_band( self, eband_in ):
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
                    self.bands_bulk[kcount][ecount] = c1
                    ecount = ecount +1
                    if ecount == self.nband_bulk:
                        kcount = kcount +1; read_flag = 0

            if line.find('energy_eigen_values') > 0 or \
               line.find('energy eigenvalues') > 0:
                read_flag = 1;  ecount = 0

            if read_flag == 1 and line.find('ik =') > 0:
                data = line.split();  num = len(data)
                self.kfrac_bulk[kcount][0] = float(data[4]);
                self.kfrac_bulk[kcount][1] = float(data[5]);
                self.kfrac_bulk[kcount][2] = float(data[6]);

                read_flag = read_flag +1;


###################################### Main ########
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='surface projected bulk band structure')

    parser.add_argument( 'bulk_band_file', help='nfenergy.data' )
    parser.add_argument( 'kpt_file', 
                         help='k-point generation control file (bandkpt.in)' )

    parser.add_argument( '--e_range', help='energy range (min and max)',
                         default='none', nargs=2 )
    parser.add_argument( '--e_inc', help='energy increment',
                         default='none', nargs=1 )
    parser.add_argument( '--unit',
                         help='unit (eV, Ry, Ha) [default:eV]',
                         default='eV' )

    parser.add_argument( '--plot_style',
                         help='plot style (1,2) [default:1]',
                         type=int, default=1 )

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
                         help='output filename [default:energy_band]',
                         default='energy_band' )
    parser.add_argument( '--force_normal', 
                         help='action: proceed even when the third reciprocal vector is not normal to surface [default:False]',
                         action='store_true', default=False )
    parser.add_argument( '--with_fermi', 
                         help='action: plot line at 0.0 eV [default:False]',
                         action='store_true', default=False )

    parser.add_argument( '--surf_band_file', 
                         help='surface band file',
                         default='none' )
    parser.add_argument( '--e_shift', 
                         help='energy shift of surface bands relative to bulk bands [default:0.0]',
                         type=float, default=0.0 )
    parser.add_argument( '--ndiv_erange_map', 
                         help='divison of energy range for mapping bulk bands',
                         default='none' )
    parser.add_argument( '--broadening_width_map', 
                         help='broadening width for mapping bulk bands',
                         type=int, default=5 )
    parser.add_argument( '--spin_sum_map', 
                         help='action: generate map by summing two spin components [default:False]',
                         action='store_true', default=False )

    args = parser.parse_args()

    surface_projected_bulk_band( args.bulk_band_file, args.kpt_file,
                                 args.e_range, args.e_inc, args.unit,
                                 args.plot_style, args.cb_range,
                                 args.line_width, args.window_scale, args.vertical,
                                 args.fig_format, args.out_file, 
                                 args.force_normal, args.with_fermi,
                                 args.surf_band_file, args.e_shift, 
                                 args.ndiv_erange_map, 
                                 args.broadening_width_map, args.spin_sum_map )

 
