#!/usr/bin/env python3

import math
import sys
import optparse
import glob
import os

class Cube:

    real_part = None
    imag_part = None
    ndiv      = None
    div       = None
    coord     = None
    comments  = None
    sq        = None
    molid     = None
    natm      = 0
    append    = False

    def __init__(self, cube):
        self.cube = cube
        self.parse_cube()

    def sqcube(self):
        if len(self.real_part) != len(self.imag_part):
            print('invalid molecular orbital file (len(realpart) != len(imagpart))')
            sys.exit()
        nelem = len(self.real_part)
        if self.sq is None:
            self.sq = []
        for i in range(nelem):
            self.sq.append(0.0)
        for i in range(nelem):
            self.sq[i] = self.real_part[i]*self.real_part[i] + \
                         self.imag_part[i]*self.imag_part[i]
    def set_append_mode(self, append):
        self.append = append
                    
    def export_sq_cube(self,cbfile):
        if not self.append:
            cfile = open(cbfile,'w')
        else:
            cfile = open(self.cube, 'a')

        if not self.append:
            for comm in self.comments:
                cfile.write(comm+'\n')
            cfile.write('{:>10}'.format(abs(self.natm))+' '+'{:>12.5f}'.format(0.0)+' '+'{:>12.5f}'.format(0.0)+\
                    ' '+'{:>12.5f}'.format(0.0)+'\n')

            cfile.write('{:>10}'.format(self.ndiv[0])+' '+'{:>12.5f}'.format(self.div[0][0])+' '+'{:>12.5f}'.format(self.div[0][1])+' '\
                    +'{:>12.5f}'.format(self.div[0][2])+'\n')
            cfile.write('{:>10}'.format(self.ndiv[1])+' '+'{:>12.5f}'.format(self.div[1][0])+' '+'{:>12.5f}'.format(self.div[1][1])+' '\
                    +'{:>12.5f}'.format(self.div[1][2])+'\n')
            cfile.write('{:>10}'.format(self.ndiv[2])+' '+'{:>12.5f}'.format(self.div[2][0])+' '+'{:>12.5f}'.format(self.div[2][1])+' '\
                    +'{:>12.5f}'.format(self.div[2][2])+'\n')
            for i in range(abs(self.natm)):
                cfile.write('{:>12.6f}'.format(self.coord[i][0])+' '+'{:>12.6f}'.format(self.coord[i][1])+' '+\
                            '{:>12.6f}'.format(self.coord[i][2])+' '+'{:>12.6f}'.format(self.coord[i][3])+' '+\
                            '{:>12.6f}'.format(self.coord[i][4])+'\n')

        self.sqcube()
        nelem = len(self.sq)
        for ielem in range(nelem):
            cfile.write('{:.5e}'.format(self.sq[ielem])+' ')
            if (ielem+1)%6==0:
                cfile.write('\n')
        cfile.close()

    def parse_cube(self):
        cfile = open(self.cube)
        lineno = 0
        atm_read = False
        nelem = 0
        for line in cfile:
            line = line.strip()
            lineno += 1
            if self.comments is None:
                self.comments = []
            if lineno<=2:
                self.comments.append(line)
                continue
            words = line.split()
            if lineno == 3:
                self.natm = int(words[0])
                if self.natm>=0:
                    print('error! natm must be a negative value for\
                            cube file containing molecular orbital data')
                    sys.exit()
                continue
            if lineno>3 and lineno<=6:
                if self.ndiv is None:
                    self.ndiv = []
                self.ndiv.append(int(words[0]))
                if self.div is None:
                    self.div = []
                self.div.append((float(words[1]), float(words[2]), float(words[3])))
                continue
            if lineno>6 and not atm_read:
                if self.coord is None:
                    self.coord = []
                self.coord.append((float(words[0]),float(words[1]),\
                              float(words[2]),float(words[3]),float(words[4])))
                atm_read = len(self.coord) == abs(self.natm)
                continue

            if lineno == 7+abs(self.natm):
                nelem = self.ndiv[0]*self.ndiv[1]*self.ndiv[2]
                self.molid = line
                continue

            if self.real_part is None:
                self.real_part = []

            if len(self.real_part)<nelem:
                for word in words:
                    self.real_part.append(float(word))
                continue

            if self.imag_part is None:
                self.imag_part = []

            if len(self.imag_part)<nelem:
                for word in words:
                    self.imag_part.append(float(word))
                continue
        cfile.close()

def get_inputs(inputstr):
    print('inputstr '+inputstr)
    inputs = []
    cands = inputstr.split(',')
    for cand in cands:
        cbs = glob.glob(cand)
        if len(cbs)==0:
            print('WARN : no matching file '+cand)
        for cb in cbs:
            inputs.append(cb)
    return inputs

def get_state_str_from(cbfile):
    ar = cbfile.split('.')
    if len(ar)<3:
        return None
    return ar[1]

if __name__ == '__main__':
    parser = optparse.OptionParser(usage="%prog [options]",description="output square of wave functions")
    parser.add_option('-i', '--input', dest='input',type=str,default='*.cube')
    parser.add_option('-o', '--output_prefix', dest='output_prefix', type=str, default='nfwfsq')
    parser.add_option('-a', '--append', dest='append', action='store_true', default=False)
    (options,args) = parser.parse_args()

    inputs = get_inputs(options.input)
    idefault = 0
    for inp in inputs:
        print('target wavefunction file : '+inp)
        cube = Cube(inp)
        cube.set_append_mode(options.append)
        outfile = None
        if not options.append:
            state = get_state_str_from(inp)
            if state is None:
                state = str(idefault)
                idefault += 1
            outfile = options.output_prefix+'.'+state+'.cube'
            print('square of the wavefunctions output to : '+outfile)
        else:
            print('square of the wavefunctions are appended to : '+inp)
        cube.export_sq_cube(outfile)

