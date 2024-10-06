''' provides an interface for the Open Babel program '''
import logging
logger = logging.getLogger('piou.config.obabel')

from piou import config
from piou.util import pyutil
from piou.config import cif
import os
import subprocess

class OpenBabelInterface(config.AtomConfigGenerator):
    ''' the interface coords generator for the Open Babel program. will use the CIF format as 
    an intermediate file format '''

    bfile = None

    def get_ftype(self,mode=None):
        babelh = os.popen('babel -H')
        found=False
        for b in babelh:
            if '[Write-only]' in b:
                found=True
                break
        if not found:
            babelh = os.popen('babel -L formats')
        ft = []
        for b in babelh:
            bb = b.split('--')
            if len(bb)==2:
                if mode is not None and \
                    ('Write-only' in bb[1] and mode=='r' or \
                     'Read-only'  in bb[1] and mode=='w'):
                    continue
                if '[Write-only]' in bb[1]:
                    tmp = bb[1].strip()
                    bb[1] = tmp[0:len(tmp)-13].strip()
                if '[Read-only]' in bb[1]:
                    tmp = bb[1].strip()
                    bb[1] = tmp[0:len(tmp)-12].strip()
                boo=[]
                for bbb in bb:
                    if len(bbb.strip())==0 or bbb.startswith('-'):
                        continue
                    boo.append(bbb.strip())
                if len(boo)==2:
                    ft.append(boo)
        return ft

    def get_all_subtypes(self,mode='w'):
        fts = self.get_ftype(mode=mode)
        ret=[]
        for ft in fts:
            ret.append(ft[0])
        return ret

    def select_subtype(self,mode=None):
        ft = self.get_ftype(mode=mode)
        (lines,cols) = pyutil.getTerminalSize()
        ncol=int(cols/55)
        if ncol==0:
            ncol=1
        icount=0
        ff=[]
        if mode is None:
            print('supported file types of the Open Babel program :')
        elif mode == 'r':
            print('importable file types of the Open Babel program :')
        elif mode == 'w':
            print('exportable file types of the Open Babel program :')
        line = ''
        for f in ft:
            icount += 1
            if icount%ncol==0:
                icount=0
                line += '\n'
            line += f[0].ljust(9)+' '+f[1].ljust(45)
            ff.append(f[0])
        print(line)
        while True:
            ret = pyutil.interactive(msg="the file type to be passed to the Open Babel program",\
            condition=lambda f:True,typ=str, default='')
            if len(ret)==0:
                continue
            if not ret in ff:
                print('invalid selection : ',ret,)
            else:
                break
        return ret

    def get_defaultname(self,arg):
        if self.subtype==None:
            return 'obabel.'+arg
        else:
            return 'obabel.'+str(self.subtype)

    def get_subtype_selector(self):
        return subtype_selector

    def check_support(self):
        p = subprocess.Popen('babel -H', shell=True, 
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        (sin,sout,serr) = (p.stdin, p.stdout, p.stderr)
        ii = 0
        for st in sout:
            ii += 1
        if ii==0:
            return False
        return True

    def __init__(self,coord_file):
        super().__init__()
        self.bfile = coord_file

    tmpfile='phase_io_util.cif'
    def export_atomic_configuration(self,atomic_coordinates,to_file=None,frame_no=None,all_frames=True):
        logger.info( 'subtype : '+str(self.subtype))
        cifgen = cif.CIF(coord_file=None)
        cifgen.export_atomic_configuration(atomic_coordinates,to_file=self.tmpfile,frame_no=frame_no,all_frames=all_frames)
        babel = 'babel -icif '+self.tmpfile+' -o'+self.subtype+' '+str(to_file)
        if to_file is None:
            babel = 'babel -icif '+self.tmpfile+' -o'+self.subtype
        p = subprocess.Popen(babel, shell=True, 
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        (imsg,msg,emsg) = (p.stdin, p.stdout, p.stderr)
        for m in msg:
            logger.info(m.split('\n')[0])
        for m in emsg:
            logger.info(m.split('\n')[0])
        if os.path.exists(self.tmpfile):
            os.remove(self.tmpfile)

    def _gen_atomic_configuration(self):
        babel = 'babel -i'+self.subtype+' '+self.bfile+' -ocif '+self.tmpfile
        p = subprocess.Popen(babel, shell=True, 
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        (imsg,msg,emsg) = (p.stdin, p.stdout, p.stderr)
        for m in msg:
            logger.info(m.split('\n')[0])
        for m in emsg:
            logger.info(m.split('\n')[0])
        if not os.path.exists(self.tmpfile):
            return None
        cifgen = cif.CIF(coord_file=self.tmpfile)
        ret = cifgen._gen_atomic_configuration()
        os.remove(self.tmpfile)
        return ret

    def get_name(self):
        return str(self.get_subtype())+' format via the Open Babel program'

