#!/usr/bin/env python3

import optparse
import sys
import math

b2ang=0.529177

def log_intpl(z, z1, z2, y1, y2):
    if y1<0 or y2<0:
        return 0.0
    lyy1 = math.log10(z/z1)*math.log(y2/y1)/math.log10(z2/z1)
    return y1*10**lyy1

def extract(cubeobj, zaxis, zindex, zval, outfile):
    nx = cubeobj['meshx']['ndiv']
    ny = cubeobj['meshy']['ndiv']
    nz = cubeobj['meshz']['ndiv']
    dx = cubeobj['meshx']['div']
    dy = cubeobj['meshy']['div']
    dz = cubeobj['meshz']['div']
    charge = cubeobj['charge']
    fo = open(outfile, 'w')
    if zaxis=='1':
        for i in range(ny):
            for j in range(nz):
                z = log_intpl(zval, zindex*dx[0], (zindex+1)*dx[0], charge[zindex][i][j], charge[zindex+1][i][j])
                fo.write(str(i*dy[1])+' '+str(j*dz[2])+' '+str(z)+'\n')
            fo.write('\n')
    if zaxis=='2':
        for i in range(nx):
            for j in range(nz):
                z = log_intpl(zval, zindex*dy[1], (zindex+1)*dy[1], charge[i][zindex][j], charge[i][zindex+1][j])
                fo.write(str(i*dx[0])+' '+str(j*dz[2])+' '+str(z)+'\n')
            fo.write('\n')
    if zaxis=='3':
        for i in range(nx):
            for j in range(ny):
                z = log_intpl(zval, zindex*dz[2], (zindex+1)*dz[2], charge[i][j][zindex], charge[i][j][zindex+1])
                fo.write(str(i*dx[0])+' '+str(j*dy[1])+' '+str(z)+'\n')
            fo.write('\n')
    fo.close()
    return

def parse_cubefile(cube):
    f = open(cube)
    f.readline()
    f.readline()
    f.readline()
    cubeobj = {}
    wordsx = f.readline().split()
    wordsy = f.readline().split()
    wordsz = f.readline().split()
    ndivx = int(wordsx[0])
    divx  = (float(wordsx[1])*b2ang, float(wordsx[2])*b2ang, float(wordsx[3])*b2ang)
    ndivy = int(wordsy[0])
    divy  = (float(wordsy[1])*b2ang, float(wordsy[2])*b2ang, float(wordsy[3])*b2ang)
    ndivz = int(wordsz[0])
    divz  = (float(wordsz[1])*b2ang, float(wordsz[2])*b2ang, float(wordsz[3])*b2ang)
    cubeobj['meshx'] = {'ndiv':ndivx, 'div':divx}
    cubeobj['meshy'] = {'ndiv':ndivy, 'div':divy}
    cubeobj['meshz'] = {'ndiv':ndivz, 'div':divz}
    charge_list = [[[0.0 for i in range(ndivz)] for j in range(ndivy)] for k in range(ndivx)] # [x, y, z]
    icount = 0
    for line in f:
        words = line.split()
        try:
            int(words[0])
            continue
        except ValueError:
            pass
        for w in words:
            val = float(w)
            mod0 = icount%(ndivz*ndivy)
            ix = int(icount/(ndivz*ndivy))
            iy = int(mod0/ndivz)
            iz = mod0%ndivz
            charge_list[ix][iy][iz] = val
            icount+=1
    f.close()
    cubeobj['charge'] = charge_list
    return cubeobj

def check_index(indx, cubeobj, zaxis):
    if indx<0:
        return False
    if zaxis == '1':
        nmesh = cubeobj['meshx']['ndiv']
    if zaxis == '2':
        nmesh = cubeobj['meshy']['ndiv']
    if zaxis == '3':
        nmesh = cubeobj['meshz']['ndiv']
    return indx< nmesh-1
    
def find_pos(cubeobj, zval, zaxis):
    if zaxis == '1':
        nmesh = cubeobj['meshx']['ndiv']
        mesh  = cubeobj['meshx']['div'][0]
    if zaxis == '2':
        nmesh = cubeobj['meshy']['ndiv']
        mesh  = cubeobj['meshy']['div'][1]
    if zaxis == '3':
        nmesh = cubeobj['meshz']['ndiv']
        mesh  = cubeobj['meshz']['div'][2]
    return int(zval/mesh)

def get_zval(ind, cubeobj, zaxis):
    if zaxis == '1':
        mesh  = cubeobj['meshx']['div'][0]
    if zaxis == '2':
        mesh  = cubeobj['meshy']['div'][1]
    if zaxis == '3':
        mesh  = cubeobj['meshz']['div'][2]
    return ind*mesh

def run():
    parser = optparse.OptionParser(usage='%prog [options]',description='extract slice data from a cube file ')
    parser.add_option('-c','--cube',type=str,dest='cube',help='the target cube file. defaults to nfchr.cube', default='nfchr.cube')
    parser.add_option('-z','--zval', type=float, dest='zval', help='the height from which the slice data shall be extracted.', default=None)
    parser.add_option('-a','--zaxis',type='choice', choices=('1','2','3'), dest='zaxis',help=\
    'specify which direction is considered as the z-axis. 1 stands for the a-vector, 2 stands for the b-vector, and 3 stands for the c-vector. defaults to 3',default='3')
    parser.add_option('-i','--zindex', type=int, dest='zindex', help=\
    'the index of the height from which the slice data shall be extracted. if --zval and --zindex are both specified, the former will be preferred.', default=0)
    parser.add_option('-o','--output',type=str,dest='output',help='the file to which the extracted results are output', default='slice.dat')
    (options,args) = parser.parse_args()

    cubeobj = parse_cubefile(options.cube)
    ind = options.zindex
    if options.zval is not None:
        ind = find_pos(cubeobj, options.zval, options.zaxis)
        zval = options.zval
    else:
        zval = get_zval(ind, cubeobj, options.zaxis)

    if not check_index(ind, cubeobj, options.zaxis):
        print('index '+str(ind)+' out of range')
        sys.exit()
    extract(cubeobj, options.zaxis, ind, zval, options.output)

if __name__ == '__main__':
    run()

