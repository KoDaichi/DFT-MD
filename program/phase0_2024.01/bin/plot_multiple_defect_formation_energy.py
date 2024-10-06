#!/usr/bin/env python3

import sys, math;
import argparse
import subprocess
import shlex

class multiple_plot:
    def __init__( self, fname_in, fname_out, emin, emax, vmin, vmax, image_format,
                  keypos_h, keypos_v ):
        self.num_defect_file = 0;
        self.defect_file = [];    self.defect_title = [];
        self.emin = emin;         self.emax = emax;
        self.vmin = vmin;         self.vmax = vmax;
        self.band_gap = 0.0;

        in_f = open( fname_in, "r" );
        self.read_input( in_f );
        in_f.close();

        fname_out_gnu = fname_out +'.gnu'
        out_f = open( fname_out_gnu, "w" );
        self.write_gnuplot_file( out_f, fname_out, image_format, keypos_h, keypos_v )
        out_f.close();

        cmd = 'gnuplot '+fname_out_gnu
        tokens = shlex.split(cmd)
        subprocess.run(tokens)

    def write_gnuplot_file( self, out_f, fname_out, image_format, 
                            keypos_h, keypos_v ):
        out_f.write('{:s}\n'.format("# Formation energy"));
        out_f.write('{:s}\n'.format("se gr"));
        out_f.write('{:s}{:s}{:s}{:s}\n'.format("set key opaque box ",
                                                keypos_h," ", keypos_v) );
        out_f.write('{:s}\n'.format("set st da line"));

        if self.emax < 3.0:
            out_f.write('{:s}\n'.format("set format x '%3.1f'"))

        out_f.write('{:s}\n'.format("set xlabel 'Fermi energy (eV)'"));
        out_f.write('{:s}\n'.format("set ylabel 'Formation energy (eV)'"));

        if self.emin > -100 and self.emax < 100:
            out_f.write('{:s}{:8.3f}{:s}{:8.3f}{:s}\n'.
                        format("set xr [",self.emin,":",self.emax,"]"));
        elif self.emin > -100:
            out_f.write('{:s}{:8.3f}{:s}{:s}\n'
                        .format("set xr [",self.emin,":","]"));
        elif self.emax < 100:
            out_f.write('{:s}{:s}{:8.3f}{:s}\n'.
                        format("set xr [",":",self.emax,"]"));

        if self.vmin > -100 and self.vmax < 100:
            out_f.write('{:s}{:8.3f}{:s}{:8.3f}{:s}\n'.
                        format("set yr [",self.vmin,":",self.vmax,"]"));
        elif self.vmin > -100:
            out_f.write('{:s}{:8.3f}{:s}{:s}\n'
                        .format("set yr [",self.vmin,":","]"));
        elif self.vmax < 100:
            out_f.write('{:s}{:s}{:8.3f}{:s}\n'.
                        format("set yr [",":",self.vmax,"]"));

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
        for j in range( self.num_defect_file ):
            out_f.write('  {:s}{:s}{:s}'.format("'",self.defect_file[j],"'"))
            out_f.write('{:s}'.format(" us 1:2 ") )
            out_f.write('{:s}{:s}{:s}'.format("title '", 
                                              self.defect_title[j]
                                              .replace('_',"\_"),"'"));
            if j < self.num_defect_file -1:
                out_f.write('{:s}\n'.format(", \\"));

    def read_input( self, in_f ):
        while True:
            line = in_f.readline()
            if not line: break;
                
            if line.find("&band_gap") == 0 :
                line = in_f.readline()
                self.band_gap = float(line)

            if line.find("&List") == 0 :
                while True:
                    line = in_f.readline()
                    data = line.split();   num_data = len(data)
                    if num_data == 0: break

                    self.defect_title.append( data[0] )
                    self.defect_file.append( data[1] +'.min' )
                    self.num_defect_file = self.num_defect_file +1

###################################### Main ########
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Multiple plot of formation energy of defects');

    parser.add_argument( 'input', help='input file');
    parser.add_argument( '-o', '--outfile',
                         help='output file name (default:result_all)',
                         default='result_all');
    parser.add_argument( '--emin', 
                         help='minimum of Fermi energy',
                         type=float, default = -1E4 );
    parser.add_argument( '--emax', 
                         help='maximum of Fermi energy',
                         type=float, default = 1E4 );
    parser.add_argument( '--vmin', 
                         help='minimum of Formation energy',
                         type=float, default = -1E4 );
    parser.add_argument( '--vmax', 
                         help='maximum of Formation energy',
                         type=float, default = 1E4 );
    parser.add_argument( '--keypos_h', 
                         help='horizontal key position (left/center/right)',
                         default = 'right' );
    parser.add_argument( '--keypos_v', 
                         help='vertical key position (top/center/bottom)',
                         default = 'top' );
    parser.add_argument( '--image_format', 
                         help='image format of figure (png/eps, default:png)',
                         default='png');

    args = parser.parse_args();

    multiple_plot( args.input, args.outfile, args.emin, args.emax, 
                   args.vmin, args.vmax, args.image_format,
                   args.keypos_h, args.keypos_v);


