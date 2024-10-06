#!/usr/bin/env python3

import tkinter
import tkinter.filedialog
from matplotlib.backends.backend_tkagg import (
     FigureCanvasTkAgg, NavigationToolbar2Tk)
from tkinter import *
from tkinter.ttk import *
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import logging
import re
import optparse
import os
import sys
import math
from collections import OrderedDict
from piou.util import pyutil
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

bohr2ang=0.529177

logformat  = '%(asctime)8s %(name)15s:%(levelname)5s: %(message)s'
logformat_dbg = '%(asctime)8s %(name)20s %(lineno)6d %(levelname)10s: %(message)s'
logtimefmt = '%y/%d/%b %H:%M:%S'
logger = logging.getLogger('dos')
def init_logger(loglevel):
    if loglevel=="2":
        logging.basicConfig(level=logging.DEBUG,
                    format=logformat_dbg,
                    datefmt=logtimefmt
                    )
    elif loglevel=="1":
        logging.basicConfig(level=logging.INFO,
                    format=logformat,
                    datefmt=logtimefmt
                    )
    else:
        logging.basicConfig(level=logging.WARN,
                    format=logformat,
                    datefmt=logtimefmt
                    )

output000_p = re.compile('output\d\d\d')
def get_mostrescent_outputxxx(directory='.'):
    """ get the most rescent PHASE/0 log file (outputxxx file) under the specified directory"""
    files = os.listdir(directory)
    ret = None
    mod = -1
    for f in files:
        if output000_p.match(f) is not None:
            t = os.path.getmtime(directory+'/'+f)
            if t>mod:
                ret = f
    if ret==None:
        return None
    return directory+'/'+ret

def copy_str(a1,retval_if_None):
    if a1 is None:
        return retval_if_None
    return str(a1)

class DosAttributes:
    elements = None
    xyz = None
    layers = None
    def __init__(self,output000):
        self.output000 = output000
        self.__parse()

    def get_elements(self):
        if self.elements is None or self.xyz is None:
            return None
        if len(self.elements) != len(self.xyz):
            return None
        return self.elements

    def get_xyz(self):
        if self.elements is None or self.xyz is None:
            return None
        if len(self.elements) != len(self.xyz):
            return None
        return self.xyz

    def get_num_atoms(self):
        return len(self.elements)

    def __parse(self):
        if self.output000 is None:
            logger.warning('log file unspecified')
        elif not os.path.exists(self.output000):
            logger.warning('log file :'+str(self.output000)+' does not exist')
            return
        out = open(self.output000)
        reading_coordinates = False
        ihead = 0
        for line in out:
            line = line.strip()
            if 'Atomic coordinates' in line:
                reading_coordinates = True
                continue
            if not reading_coordinates:
                continue
            if ihead<2:
                ihead += 1
                continue
            words = line.split()
            if len(words) != 12:
                break
            try:
                atom = (float(words[5]),float(words[6]),float(words[7]))
                if self.xyz is None:    
                    self.xyz = [] 
                self.xyz.append(atom)
            except ValueError:
                logger.error('invalid atomic coordinate line : '+line)
                break
            if self.elements is None:
                self.elements = []
            self.elements.append(words[11])
        out.close()

        out = open(self.output000)
        reading_layer = False
        for line in out:
            line = line.strip()
            if '!!ldos     no,        min,           max' in line:
                reading_layer = True
                continue
            if not reading_layer:
                continue
            words = line.split()
            if len(words) != 4:
                break
            try:
                layer = (float(words[2]),float(words[3]))
                if self.layers is None:
                    self.layers = []
                self.layers.append(layer)
            except ValueError:
                logger.error('invalid layer line : '+line)
                break
        out.close()

    def __str__(self):
        if self.elements is None or self.xyz is None:
            return 'None coordinates'
        if len(self.elements) != len(self.xyz):
            return 'number of coordinates and elements do not match' 
        ret = 'number of atoms defined : '+str(len(self.elements))+'\n'
        for i in range(len(self.elements)):
            atm = self.elements[i]+' '+str(self.xyz[i][0])+' '+str(self.xyz[i][1])+' '+str(self.xyz[i][2])+'\n'
            ret += atm
        if self.layers is not None:
           ret += 'number of layers defined : '+str(len(self.layers))+'\n'
           for layer in self.layers:
               ret += '('+str(layer[0])+', '+str(layer[1])+')\n'
        return ret
            
class DosData:
    energy = []
    tdos = None
    aldos = None
    layerdos = None
    pdos = None
    spin = False
    dos_attributes = None
    def __init__(self,dosdata_file, output000):
        if not os.path.exists(dosdata_file):
            logger.error('dos file : '+str(dosdata_file)+' does not exist')
            sys.exit()
        self.dosdata_file = dosdata_file
        if output000 is None or not os.path.exists(output000):
            logger.warning('log file : '+str(output000)+' does not exist')
        else:
            self.dos_attributes = DosAttributes(output000) 
            logger.debug('dos attributes')
            logger.debug(self.dos_attributes)
        self.__parse()

    def get_num_layers(self):
        if self.dos_attributes is not None and self.dos_attributes.layers is not None:
            return len(self.dos_attributes.layers)
        nl = 0
        if self.layerdos is not None:
            keys = self.layerdos.keys()
            for key in keys:
                il=self.get_il_from_layerdoskey(key)
                if il>nl:
                    nl = il
        return nl

    def get_layer_positions(self):
        if self.dos_attributes is None or self.dos_attributes.layers is None:
            n = self.get_num_layers()
            ret = []
            for i in range(n):
                ret.append((float(i),0))
            return ret
        ret = []
        for l in self.dos_attributes.layers:
            ret.append(l)
        return ret

    def get_num_atoms(self):
        if self.dos_attributes is not None and self.dos_attributes.elements is not None:
            return self.dos_attributes.get_num_atoms()
        na=0
        if self.aldos is not None:
            keys = self.aldos.keys()
            for key in keys:
                ia=self.get_ia_from_aldoskey(key)
                if ia>na:
                    na = ia
        if self.pdos is not None:
            keys = self.pdos.keys()
            for key in keys:
                (ia,l,m,t)=self.get_almt_from_pdoskey(key)
                if ia>na:
                    na = ia
        return na

    def __parse(self):
        dosfile = open(self.dosdata_file)
        maxnat=1
        maxnlayer=1
        for line in dosfile:
            line = line.strip()
            words = line.split() 
            if (len(words)>=4 and words[1].startswith('num_atom')):
                try:
                    nat = int(words[3])
                    if nat>maxnat:
                        maxnat = nat
                except ValueError:
                    pass
            if (len(words)>=3 and words[1].startswith('ia=')):
                try:
                    nat = int(words[2])
                    if nat>maxnat:
                        maxnat = nat
                except ValueError:
                    pass
            if (len(words)>=4 and words[1].startswith('num_layer')):
                try:
                    nl = int(words[3])
                    if nl>maxnlayer:
                        maxnlayer = nl
                except ValueError:
                    pass
        dosfile.close()
        aketa = int(math.log10(maxnat))+1
        lketa = int(math.log10(maxnlayer))+1

        dosfile = open(self.dosdata_file)
        self.spin = False
        in_dos = False
        tdos_read=False
        linecount=0
        for line in dosfile:
            line = line.strip()
            words = line.split()
            linecount += 1
            if words[0].startswith('No') and not tdos_read:
                if not (len(words)==6 or len(words)==10):
                    self.__print_err_and_exit(self.dosdata_file,line,linecount)
                    sys.exit()
                in_dos = True
                if len(words)==10:
                    self.spin = True
                self.energy = []
                if self.tdos is None:
                    self.tdos = OrderedDict()
                tdoskey = 'total'
                self.tdos[tdoskey] = OrderedDict()
                self.tdos[tdoskey]['up'] = []
                self.tdos[tdoskey]['intup'] = []
                if self.spin:
                    self.tdos[tdoskey]['down'] = []
                    self.tdos[tdoskey]['intdown'] = []
                    self.tdos[tdoskey]['inttotal'] = []
                continue
            if words[0].startswith('No') and tdos_read:
                continue

            if words[0].startswith('ALDOS'):
                if self.aldos is None:
                    self.aldos = OrderedDict()
                if len(words)<4:
                    self.__print_err_and_exit(self.dosdata_file,line,linecount)
                atmkey = 'atom'+str(words[3]).zfill(aketa)
                self.aldos[atmkey] = OrderedDict()
                try:
                    self.aldos[atmkey]['iat'] = int(words[3])
                except ValueError:
                    logger.warning('invalid atom ID : '+str(words[3]))
                    self.aldos[atmkey]['iat'] = 0

                self.aldos[atmkey]['up'] = []
                self.aldos[atmkey]['intup'] = []
                if self.spin:
                    self.aldos[atmkey]['down'] = []
                    self.aldos[atmkey]['intdown'] = []
                    self.aldos[atmkey]['inttotal'] = []
                thedict = self.aldos
                thekey = atmkey
                in_dos = True
                continue
            if words[0].startswith('LAYER'):
                if self.layerdos is None:
                    self.layerdos = OrderedDict()
                if len(words)<4:
                    self.__print_err_and_exit(self.dosdata_file,line,linecount)
                layerkey = 'layer'+str(words[3]).zfill(lketa)
                self.layerdos[layerkey] = OrderedDict()
                try:
                    self.layerdos[layerkey]['ilayer'] = int(words[3])
                except ValueError:
                    logger.warning('invalid layer ID : '+str(words[3]))
                    self.layerdos[layerkey]['ilayer'] = 0
                self.layerdos[layerkey]['up'] = []
                self.layerdos[layerkey]['intup'] = []
                if self.spin:
                    self.layerdos[layerkey]['down'] = []
                    self.layerdos[layerkey]['intdown'] = []
                    self.layerdos[layerkey]['inttotal'] = []
                thedict = self.layerdos
                thekey = layerkey
                in_dos = True
                continue
            if words[0].startswith('PDOS'):
                if self.pdos is None:
                    self.pdos = OrderedDict()
                if len(words)<9:
                    self.__print_err_and_exit(self.dosdata_file,line,linecount)
                pdoskey = 'atom'+str(words[2]).zfill(aketa)+'_l'+str(words[4])+'_m'+str(words[6])+'_t'+str(words[8])
                self.pdos[pdoskey] = OrderedDict()
                try:
                    self.pdos[pdoskey]['iat'] = int(words[2])
                except ValueError:
                    logger.warning('invalid pdos atom ID : '+words[2])
                    self.pdos[pdoskey]['iat'] = 0
                try:
                    self.pdos[pdoskey]['l'] = int(words[4])
                except ValueError:
                    logger.warning('invalid pdos l ID : '+words[4])
                    self.pdos[pdoskey]['l'] = -1
                try:
                    self.pdos[pdoskey]['m'] = int(words[6])
                except ValueError:
                    logger.warning('invalid pdos m ID : '+words[6])
                    self.pdos[pdoskey]['m'] = -1
                try:
                    self.pdos[pdoskey]['t'] = int(words[8])
                except ValueError:
                    logger.warning('invalid pdos t ID : '+words[8])
                    self.pdos[pdoskey]['t'] = -1

                self.pdos[pdoskey]['up'] = []
                self.pdos[pdoskey]['intup'] = []
                if self.spin:
                    self.pdos[pdoskey]['down'] = []
                    self.pdos[pdoskey]['intdown'] = []
                    self.pdos[pdoskey]['inttotal'] = []
                thedict = self.pdos
                thekey = pdoskey
                in_dos = True
                continue
            if words[0].startswith('END'):
                in_dos = False
                if not tdos_read:
                    tdos_read = True
                continue
            if not tdos_read and in_dos:
                try:
                    if not self.spin:
                        self.energy.append(float(words[3]))
                        self.tdos[tdoskey]['up'].append(float(words[4]))
                        self.tdos[tdoskey]['intup'].append(float(words[5]))
                    else:
                        self.energy.append(float(words[4]))
                        self.tdos[tdoskey]['up'].append(float(words[5]))
                        self.tdos[tdoskey]['down'].append(float(words[6]))
                        self.tdos[tdoskey]['intup'].append(float(words[7]))
                        self.tdos[tdoskey]['intdown'].append(float(words[8]))
                        self.tdos[tdoskey]['inttotal'].append(float(words[9]))
                except ValueError:
                    self.__print_err_and_exit(self.dosdata_file,line,linecount)
            if tdos_read and in_dos:
                try:
                    if not self.spin:
                        thedict[thekey]['up'].append(float(words[4]))
                        thedict[thekey]['intup'].append(float(words[5]))
                    else:
                        thedict[thekey]['up'].append(float(words[5]))
                        thedict[thekey]['down'].append(float(words[6]))
                        thedict[thekey]['intup'].append(float(words[7]))
                        thedict[thekey]['intdown'].append(float(words[8]))
                        thedict[thekey]['inttotal'].append(float(words[9]))
                except ValueError:
                    self.__print_err_and_exit(self.dosdata_file,line,linecount)
        dosfile.close()

    def get_il_from_layerdoskey(self,layerdoskey):
        il=0
        try:
            il = int(layerdoskey[5:len(layerdoskey)])
        except ValueError:
            logger.error('invalid layerdos key '+layerdoskey)
        return il

    def get_ia_from_aldoskey(self,aldoskey):
        ia=0
        try:
            ia = int(aldoskey[4:len(aldoskey)])
        except ValueError:
            logger.error('invalid aldos key '+aldoskey)
        return ia

    def get_almt_from_pdoskey(self,pdoskey):
        ia=0
        l=0
        m=0
        t=0
        arr = pdoskey.split('_')

        try:
            ia = int(arr[0][4:len(arr[0])])
            l  = int(arr[1][1:len(arr[1])])
            m  = int(arr[2][1:len(arr[2])])
            t  = int(arr[3][1:len(arr[3])])
        except ValueError:
            logger.error('invalid pdos key '+pdoskey)

        return (ia,l,m,t)

    def get_num_data(self):
        if self.energy is None:
            return 0
        return len(self.energy)

    def __str__(self):
        if self.energy is None or self.tdos is None:
            return ''
        if len(self.energy) != len(self.tdos['total']['up']):
            logger.error('the number of data for energy and tdos is inconsistent')
            return  ''
        ret = 'TOTAL DOS\n'
        ret += self.get_dos_str(self.tdos['total'])

        if self.aldos is not None:
            aldos = sorted(self.aldos.keys())
            for ald in aldos:
                ret += 'ALDOS for '+ald+'\n'
                ret += self.get_dos_str(self.aldos[ald])

        if self.layerdos is not None:
            layerdos = sorted(self.layerdos.keys())
            for layerd in layerdos:
                ret += 'LAYERDOS for '+layerd+'\n'
                ret += self.get_dos_str(self.layerdos[layerd])

        if self.pdos is not None:
            pdos = sorted(self.pdos.keys())
            for pd in pdos:
                ret += 'PDOS for '+pd+'\n'
                ret += self.get_dos_str(self.pdos[pd])
        return ret

    def get_dos_str(self,dos):
        nume = len(self.energy)
        ret = ''
        for i in range(nume):
            line = str(self.energy[i])+ ' ' + str(dos['up'][i])
            if self.spin:
                line += ' '+str(dos['down'][i])
            ret += line+'\n'
        return ret
        
    def __print_err_and_exit(self,dosfile,line,linecount):
        logger.error('invalid line found in '+str(dosfile)+' : '+line+'; line no. '+str(linecount))
        sys.exit()

    def get_summed_dos(self,dosids,mode):
        if mode=='total':
            logger.warning('can\'t sum total dos')
            return None
        thedos = self.get_target_dos(mode)
        keys = list(thedos.keys())
        retup = []
        retdown = []
        for i in range(self.get_num_data()):
            retup.append(0.0)
        if self.spin:
            for i in range(self.get_num_data()):
                retdown.append(0.0)
        for dosid in dosids:
            dos = thedos[keys[dosid-1]]
            for i in range(self.get_num_data()):
                retup[i] += dos['up'][i]
            if self.spin:
                for i in range(self.get_num_data()):
                    retdown[i] += dos['down'][i]
        return (retup,retdown)

    def get_target_dos(self,mode):
        thedos = self.tdos
        if mode=='atom':
            thedos = self.aldos
        if mode=='layer':
            thedos = self.layerdos
        if mode=='projected':
            thedos = self.pdos
        return thedos

    def get_energy(self):
        return self.energy

    def get_total_dos(self):
        return self.tdos

    def get_aldos(self):
        return self.aldos

    def get_layer_dos(self):
        return self.layerdos

    def get_pdos(self):
        return self.pdos

    def get_doskeys(self):
        ret = []
        if self.tdos is not None:
            ret.append('total DOS')
        if self.aldos is not None:
            keys = self.aldos.keys()
            for key in keys:
                ret.append(key)
        if self.layerdos is not None:
            keys = self.layerdos.keys()
            for key in keys:
                ret.append(key)
        if self.pdos is not None:
            keys = self.pdos.keys()
            for key in keys:
                ret.append(key)
        return ret

    @staticmethod
    def get_intidarray(idstr,nelem):
        ret = []
        if idstr is None or len(idstr)==0:
            return ret
        if idstr.lower() == 'all':
            for i in range(nelem):
                ret.append(i+1)
            return ret
        idstr = idstr.replace(' ','')
        if idstr.endswith(','):
            idstr = idstr[0:len(idstr)-1]
        css = idstr.split(',')
        for cs in css:
            hss = cs.split('-')
            if len(hss)>2:
                logger.error('invalid idstr : '+idstr)
            if len(hss)==1:
                try:
                    ret.append(int(hss[0]))
                except ValueError:
                    logger.error('invalid idstr : '+idstr)
            elif len(hss)==2:
                try:
                    mini = int(hss[0])
                    maxi = int(hss[1])
                    if mini>maxi:
                        m = mini
                        mini = maxi
                        maxi = m
                    for i in range(maxi-mini+1):
                        ret.append(i+mini)
                except ValueError:
                    logger.error('invalid idstr : '+idstr)
        return ret

    def __filter_by_elemid(self,idstr,elemid,ndat):
        aid = copy_str(idstr,'all')
        ai = self.get_intidarray(aid,ndat)
        if elemid is not None and self.dos_attributes is not None:
            elemid = elemid.replace(' ','')
            elemids = elemid.split(',')
            elements = self.dos_attributes.get_elements()
            if elements is None:
                logger.error('atomic configuration unspecified')
            else:
                aid = ''
                for a in ai:
                    for e in elemids:
                        if elements[a-1] == e:
                            aid += str(a)+','
                aid = aid[0:-1]
        return self.get_intidarray(aid,ndat)

    def get_idarray(self,mode,dosidstr,atomid,layerid,elemid,lid,mid,tid,nelem):
        if dosidstr is not None:
            return self.get_intidarray(dosidstr,nelem)

        if self.dos_attributes is None:
            logger.warning('output000 file is not specified or invalid')

        if mode == 'total':
            logger.error('dosid is not specified')
            return None

        if mode=='layer':
           li = self.__filter_by_elemid(layerid,elemid,self.get_num_layers())
           return li

        ai = self.__filter_by_elemid(atomid,elemid,self.get_num_atoms())
        if mode=='atom':
            return ai

        if self.pdos is None:
            logger.error('mode==pdos, but pdos is not allocated')
            return None
        ret = []

        if lid is not None:
            lid =lid.replace('s','0')
            lid =lid.replace('p','1')
            lid =lid.replace('d','2')
            lid =lid.replace('f','3')

        icount = 1 #dosid begins from 1
        li = self.get_intidarray(copy_str(lid,'all'),4)
        mi = self.get_intidarray(copy_str(mid,'all'),7)
        ti = self.get_intidarray(copy_str(tid,'all'),2)
        pdoskeys = self.pdos.keys()
        for pdoskey in pdoskeys:
            (ia,l,m,t) = self.get_almt_from_pdoskey(pdoskey)
            if ia in ai and l in li and m in mi and t in ti:
                ret.append(icount)
            icount += 1
        return ret

    def get_aldos_heatmap_data(self):
        if self.aldos is None:
            logger.error('heatmap support is for the layer/atomic local dos only')
            return None

        keys = list(self.aldos.keys())
        xvals = []
        for key in keys:
            xvals.append(self.get_ia_from_aldoskey(key))

        yvals = self.energy
        icount = 0
        retz = []
        if self.spin:
            retzd = []
        for j in range(len(yvals)):
            rowz = []
            if self.spin:
                rowzd = []
            for i in range(len(xvals)):
                rowz.append(0.0)
                if self.spin:
                    rowzd.append(0.0)
            retz.append(rowz)
            if self.spin:
                retzd.append(rowzd)

        for j in range(len(yvals)):
            for i in range(len(xvals)):
                dos = self.aldos[keys[i]]
                retz[j][i] = dos['up'][j]
                if self.spin:
                    retzd[j][i] = dos['down'][j]

        if self.spin:
            return (xvals,yvals,retz,retzd)
        else:
            return (xvals,yvals,retz)

    def get_layerdos_heatmap_data(self):
        if self.layerdos is None:
            logger.error('heatmap support is for the layer dos only')
            return None

        lpos = self.get_layer_positions()
        xvals = []
        for x in lpos:
            xvals.append(x[0]*bohr2ang)
        yvals = self.energy
        keys = list(self.layerdos.keys())
        icount = 0
        retz = []
        if self.spin:
            retzd = []
        for j in range(len(yvals)):
            rowz = []
            if self.spin:
                rowzd = []
            for i in range(len(xvals)):
                rowz.append(0.0)
                if self.spin:
                    rowzd.append(0.0)
            retz.append(rowz)
            if self.spin:
                retzd.append(rowzd)

        for j in range(len(yvals)):
            for i in range(len(xvals)):
                dos = self.layerdos[keys[i]]
                retz[j][i] = dos['up'][j]
                if self.spin:
                    retzd[j][i] = dos['down'][j]
        

        if self.spin:
            return (xvals,yvals,retz,retzd)
        else:
            return (xvals,yvals,retz)

class DosGUI:
    dos = None
    options = None
    canvas = None
    fig = None
    ax = None
    aldos_var = None
    aldos_vars = None
    layerdos_var = None
    layerdos_vars = None
    cmap_vars = None
    heatmap_spin_vars = None
    heatmap_spin_vars_atm = None
    pdos_var = None
    pdos_vars = None
    currdosdata = None
    currdosid = None
    currdosids = None
    currsum_mode = None
    in_heatmap = False
    root = None
    nb = None
    heatmap_mode = None
    dosopts = {
    'erange':'',
    'drange':'',
    'lrange':'',
    'arange':'',
    'with_fermi':True,
    'cmap':'viridis',
    'level':100
    }

    def __init__(self,dos,options):
        self.dos = dos
        self.options = options 
        self.fig, self.ax = plt.subplots()
        self.__init_gui__()

    def __init_gui__(self):
        try:
            self.root = tkinter.Tk()
        except tkinter.TclError:
            logger.error('failed to initiate Tk; check your X environment')
            sys.exit()
        self.root.wm_title(self.options.file)

        self.gen_fig('total',-1,False,[])
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        toolbar = NavigationToolbar2Tk(self.canvas, self.root)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

        self.set_menu(self.root)
        self.build_Notebook(self.root)

    def on_tab_click(self,event):
        mode = self.nb.tab(self.nb.select(), "text")
        clicked_tab = self.nb.tk.call(self.nb._w, "identify", "tab", event.x, event.y)
        if mode == 'total':
            self.gen_fig('total',-1,None,None)
        if mode == 'atom':
            self.gen_fig('atom',self.aldos_var.get(),False,None)
        if mode == 'layer':
            self.gen_fig('layer',self.layerdos_var.get(),False,None)
        if mode == 'projected':
            self.gen_fig('projected',self.pdos_var.get(),False,None)
        if mode == 'atom (sum)':
            dosids = []
            icount=0
            for var in self.aldos_vars:
                if var.get()==1:
                    dosids.append(icount)
                icount += 1
            self.gen_fig('atom',None,True,dosids)
        if mode == 'layer (sum)':
            dosids = []
            icount=0
            for var in self.layerdos_vars:
                if var.get()==1:
                    dosids.append(icount)
                icount += 1
            self.gen_fig('layer',None,True,dosids)
        if mode == 'projected (sum)':
            dosids = []
            icount=0
            for var in self.pdos_vars:
                if var.get()==1:
                    dosids.append(icount)
                icount += 1
            self.gen_fig('projected',None,True,dosids)
        if mode == 'layer (heatmap)':
            self.heatmap_mode = 'layer'
            self.gen_heatmap()
        if mode == 'atom (heatmap)':
            self.heatmap_mode = 'atom'
            self.gen_heatmap(mode='atom')
        self.canvas.draw()

    def build_Notebook(self,root):
        self.nb = Notebook(root)
        tdosframe = Frame(self.nb)
        self.nb.add(tdosframe, text='total')
        rb = tkinter.Radiobutton(tdosframe,text='total')
        rb.pack(side=tkinter.LEFT)

        if self.dos.get_aldos() is not None:
            def update(event):
                self.update_plot()
            self.aldos_var = tkinter.IntVar()
            self.set_dospane(self.nb,'atom',self.aldos_var,self.aldos_rb,list(self.dos.get_aldos().keys()),8)
            self.aldos_vars = []
            for i in range(len(self.dos.get_aldos().keys())):
                self.aldos_vars.append(tkinter.IntVar())
            self.set_sumdospane(self.nb,'atom (sum)',self.aldos_vars,self.aldos_cb,self.clear_atom,list(self.dos.get_aldos().keys()),8)
            heatmapframe = Frame(self.nb)
            self.nb.add(heatmapframe, text='atom (heatmap)')
            rb = tkinter.Radiobutton(heatmapframe,text='heatmap')
            rb.grid(row=0,column=0)
            lbspin = Label(heatmapframe,text='spin')
            lbspin.grid(row=0,column=1)
            self.heatmap_spin_vars_atm = StringVar()
            cbspin = Combobox(heatmapframe,textvariable=self.heatmap_spin_vars_atm)
            values = ('up')
            if self.dos.spin:
                values = ('up','down')
            cbspin['values'] = values
            cbspin.grid(row=0,column=2)
            self.heatmap_spin_vars_atm.set('up')
            cbspin.bind('<<ComboboxSelected>>', update)

        if self.dos.get_layer_dos() is not None:
            def update(event):
                self.update_plot()
            self.layerdos_var = tkinter.IntVar()
            self.set_dospane(self.nb,'layer',self.layerdos_var,self.layerdos_rb,list(self.dos.get_layer_dos().keys()),8)
            self.layerdos_vars = []
            for i in range(len(self.dos.get_layer_dos().keys())):
                self.layerdos_vars.append(tkinter.IntVar())
            self.set_sumdospane(self.nb,'layer (sum)',self.layerdos_vars,self.layerdos_cb,self.clear_layer,list(self.dos.get_layer_dos().keys()),8)
            heatmapframe = Frame(self.nb)
            self.nb.add(heatmapframe, text='layer (heatmap)')
            rb = tkinter.Radiobutton(heatmapframe,text='heatmap')
            rb.grid(row=0,column=0)
            lbspin = Label(heatmapframe,text='spin')
            lbspin.grid(row=0,column=1)
            self.heatmap_spin_vars = StringVar()
            cbspin = Combobox(heatmapframe,textvariable=self.heatmap_spin_vars)
            values = ('up')
            if self.dos.spin:
                values = ('up','down')
            cbspin['values'] = values
            cbspin.grid(row=0,column=2)
            self.heatmap_spin_vars.set('up')
            cbspin.bind('<<ComboboxSelected>>', update)

        if self.dos.get_pdos() is not None:
            self.pdos_var = tkinter.IntVar()
            self.set_dospane(self.nb,'projected',self.pdos_var,self.pdos_rb,list(self.dos.get_pdos().keys()),5)
            self.pdos_vars = []
            for i in range(len(self.dos.get_pdos().keys())):
                self.pdos_vars.append(tkinter.IntVar())
            self.set_sumdospane(self.nb,'projected (sum)',self.pdos_vars,self.pdos_cb,self.clear_projected,list(self.dos.get_pdos().keys()),5)

        self.set_optpane(self.nb)

        self.nb.bind('<ButtonRelease-1>', self.on_tab_click)
        self.nb.pack(expand=1, fill="both")
        tkinter.mainloop()

    def pdos_rb(self):
        self.gen_fig('projected',self.pdos_var.get(),False,None)
        self.canvas.draw()

    def pdos_cb(self):
        dosids = []
        icount=0
        for var in self.pdos_vars:
            if var.get()==1:
                dosids.append(icount)
            icount += 1
        self.gen_fig('projected',None,True,dosids)
        self.canvas.draw()

    def clear_projected(self):
        for var in self.pdos_vars:
            var.set(0)
        self.pdos_cb()

    def aldos_rb(self):
        self.gen_fig('atom',self.aldos_var.get(),False,None)
        self.canvas.draw()

    def aldos_cb(self):
        dosids = []
        icount=0
        for var in self.aldos_vars:
            if var.get()==1:
                dosids.append(icount)
            icount += 1
        self.gen_fig('atom',None,True,dosids)
        self.canvas.draw()

    def clear_atom(self):
        for var in self.aldos_vars:
            var.set(0)
        self.aldos_cb()

    def layerdos_rb(self):
        self.gen_fig('layer',self.layerdos_var.get(),False,None)
        self.canvas.draw()

    def layerdos_cb(self):
        dosids = []
        icount=0
        for var in self.layerdos_vars:
            if var.get()==1:
                dosids.append(icount)
            icount += 1
        self.gen_fig('layer',None,True,dosids)
        self.canvas.draw()

    def clear_layer(self):
        for var in self.layerdos_vars:
            var.set(0)
        self.layerdos_cb()

    def set_dospane(self,nb,text,var,command,keys,ncol):
        keys = list(keys)
        icol = 0
        irow = 0
        nrow = int(len(keys)/ncol)+1
        layerdosframe = Frame(nb)
        nb.add(layerdosframe, text=text)
        canvas1 = Canvas(layerdosframe,width=800,height=100)
        scroll = Scrollbar(layerdosframe, command=canvas1.yview)
        canvas1.config(yscrollcommand=scroll.set, scrollregion=(0,0,0,19*nrow))
        canvas1.pack(side=LEFT, fill=BOTH, expand=True)
        scroll.pack(side=RIGHT, fill=Y)
        frameOne = Frame(canvas1, width=800, height=450)
        canvas1.create_window(350, nrow*10, window=frameOne)
        for i in range(len(keys)):
            rb = tkinter.Radiobutton(frameOne,value=i,variable=var,text=keys[i],command=command)
            rb.grid(row=irow,column=icol)
            icol += 1
            if icol==ncol:
                icol = 0
                irow += 1

    def set_sumdospane(self,nb,text,variables,command,clear_command,keys,ncol):
        keys = list(keys)
        icol = 0
        irow = 0
        nrow = int(len(keys)/ncol)+1
        layerdosframe = Frame(nb)
        nb.add(layerdosframe, text=text)
        canvas1 = Canvas(layerdosframe,width=800,height=100)
        scroll = Scrollbar(layerdosframe, command=canvas1.yview)
        canvas1.config(yscrollcommand=scroll.set, scrollregion=(0,0,0,19*(nrow+1)))
        canvas1.pack(side=LEFT, fill=BOTH, expand=True)
        scroll.pack(side=RIGHT, fill=Y)
        frameOne = Frame(canvas1, width=800, height=450)
        canvas1.create_window(350, (nrow+1)*10, window=frameOne)
        for i in range(len(keys)):
            rb = tkinter.Checkbutton(frameOne,variable=variables[i],text=keys[i],command=command)
            rb.grid(row=irow,column=icol)
            icol += 1
            if icol==ncol and not i==len(keys)-1:
                icol = 0
                irow += 1
        btn = Button(frameOne,text='clear selection',command=clear_command)
        btn.grid(row=irow+1,columns=ncol)

    emin=None
    emax=None
    dmin=None
    dmax=None
    lmin=None
    lmax=None
    with_fermi_var=None
    tflevel=None
    def set_eminmax(self):
        minv = self.emin.get()
        maxv = self.emax.get()
        if len(minv)==0 and len(maxv)==0:
            return
        try:
            if len(minv)>0:
                float(minv)
            if len(maxv)>0:
                float(maxv)
            self.dosopts['erange'] = minv+','+maxv
        except ValueError:
            return

    def set_dminmax(self):
        minv = self.dmin.get()
        maxv = self.dmax.get()
        if len(minv)==0 and len(maxv)==0:
            return
        try:
            if len(minv)>0:
                float(minv)
            if len(maxv)>0:
                float(maxv)
            self.dosopts['drange'] = minv+','+maxv
        except ValueError:
            return

    def set_lminmax(self):
        minv = self.lmin.get()
        maxv = self.lmax.get()
        if len(minv)==0 and len(maxv)==0:
            return
        try:
            if len(minv)>0:
                float(minv)
            if len(maxv)>0:
                float(maxv)
            self.dosopts['lrange'] = minv+','+maxv
        except ValueError:
            return

    def set_aminmax(self):
        minv = self.amin.get()
        maxv = self.amax.get()
        if len(minv)==0 and len(maxv)==0:
            return
        try:
            if len(minv)>0:
                float(minv)
            if len(maxv)>0:
                float(maxv)
            self.dosopts['arange'] = minv+','+maxv
        except ValueError:
            return

    def set_with_fermi(self):
        self.dosopts['with_fermi'] = self.with_fermi_var.get()
        
    def update_plot(self):
        self.set_eminmax()
        self.set_dminmax()
        if self.in_heatmap:
            self.set_lminmax()
            self.gen_heatmap(mode=self.heatmap_mode)
        else:
            self.set_with_fermi()
            self.gen_fig(None,None,False,None)
        self.canvas.draw() 

    def cmap_selected(self,event):
        self.dosopts['cmap'] = self.cmap_vars.get()
        self.update_plot()

    def set_menu(self,root):
        def open_dos():
            fTyp = [("","*")]
            iDir = os.path.abspath(os.path.dirname(__file__))
            thefile = tkinter.filedialog.askopenfilename(filetypes = [("","*")],initialdir = '.')
            path = Path(thefile)
            output000 = get_mostrescent_outputxxx(directory=str(path.parent))
            self.dos = DosData(path,output000)
            self.nb.destroy() 
            self.build_Notebook(self.root)
            self.canvas.delete()
            plt.clf()
            self.gen_fig('total',None,None,None)
            self.canvas.draw()

        def save_data():
            thefile = tkinter.filedialog.asksaveasfilename(filetypes = [("","*")],initialdir = '.')
            if len(thefile)==0:
                print('canceled')
                return

        def exit_command():
            sys.exit()

        menu_file = tkinter.Menu()
        #menu_file.add_command(label='Open',command=open_dos)
        #menu_file.add_command(label='export data to file',command=save_data)
        #menu_file.add_command(label='export graphics to file')
        menu_file.add_command(label='exit',command=exit_command)
        file_menu = tkinter.Menu()
        file_menu.add_cascade(menu=menu_file,label='File')
        root['menu'] = file_menu

    def set_optpane(self,nb):
        def keyreleased(event):
            self.update_plot()
        frame = Frame(nb)
        nb.add(frame,text='figure options')
        lblemin = Label(frame,text='emin')
        lblemax = Label(frame,text='emax')
        lblemin.grid(row=0,column=0)
        lblemax.grid(row=0,column=2)
        self.emin = tkinter.Entry(frame,width=10)
        self.emax = tkinter.Entry(frame,width=10)
        self.emax.bind('<KeyRelease>',keyreleased)
        self.emin.bind('<KeyRelease>',keyreleased)
        self.emin.grid(row=0,column=1)
        self.emax.grid(row=0,column=3)
        lbldmin = Label(frame,text='dmin')
        lbldmax = Label(frame,text='dmax')
        lbldmin.grid(row=1,column=0)
        lbldmax.grid(row=1,column=2)
        self.dmin = tkinter.Entry(frame,width=10)
        self.dmax = tkinter.Entry(frame,width=10)
        self.dmin.bind('<KeyRelease>',keyreleased)
        self.dmax.bind('<KeyRelease>',keyreleased)
        self.dmin.grid(row=1,column=1)
        self.dmax.grid(row=1,column=3)
        self.with_fermi_var = tkinter.BooleanVar()
        if self.dosopts['with_fermi']:
            self.with_fermi_var.set(True)
        else:
            self.with_fermi_var.set(False)
        with_fermi = Checkbutton(frame,var=self.with_fermi_var,text='with_fermi',command=self.update_plot)
        with_fermi.grid(row=2,column=0)
#        btn = Button(frame,text='update plot',command=self.update_plot)
#        if self.dos.get_layer_dos() is None:
#            btn.grid(row=3,column=3)
#        else:
#            btn.grid(row=3,column=10)

        if self.dos.get_layer_dos() is not None or self.dos.get_aldos() is not None:
            lbldum = Label(frame,text='       ')
            lbldum.grid(row=0,column=4)
            lblheat = Label(frame,text='heatmap options')
            lblheat.grid(row=0,column=5)
            lbllevel = Label(frame,text='level')
            lbllevel.grid(row=1,column=5)
            self.tflevel = tkinter.Entry(frame,width=10)
            self.tflevel.insert(1,str(self.dosopts['level']))
            self.tflevel.bind('<KeyRelease>',keyreleased)
            self.tflevel.grid(row=1,column=6)

            lbllmin = Label(frame,text='x min')
            lbllmax = Label(frame,text='x max')
            lbllmin.grid(row=1,column=7)
            lbllmax.grid(row=1,column=9)
            self.lmin = tkinter.Entry(frame,width=10)
            self.lmax = tkinter.Entry(frame,width=10)
            self.lmin.bind('<KeyRelease>',keyreleased)
            self.lmax.bind('<KeyRelease>',keyreleased)
            self.lmin.grid(row=1,column=8)
            self.lmax.grid(row=1,column=10)

            lbllevel = Label(frame,text='colormap')
            lbllevel.grid(row=2,column=5)
            self.cmap_vars = StringVar()
            cbcmap = Combobox(frame,textvariable=self.cmap_vars)
            cbcmap['values'] = \
            ( 
             'viridis', 'plasma', 'inferno', 'magma', 'cividis',
             'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
             'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
             'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
             'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
             'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
             'hot', 'afmhot', 'gist_heat', 'copper',
              'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
             'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
             'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'
            )
            cbcmap.grid(row=2,column=6)
            cbcmap.bind('<<ComboboxSelected>>', self.cmap_selected)
            self.cmap_vars.set('viridis')

    def gen_fig(self,mode,dosid,sum_mode,dosids):
        self.in_heatmap = False
        x  = None
        yu = None
        yd = None
        if mode is not None:
            dosdata = self.dos.get_target_dos(mode)
            self.currdosdata = dosdata
            self.currdosid = dosid
            self.currdosids = dosids
            self.currsum_mode = sum_mode
        else:
            dosdata = self.currdosdata
            dosid = self.currdosid
            dosids = self.currdosids
            sum_mode = self.currsum_mode
        if mode=='total':
            (x,yu,yd) = self.get_x_yu_yd(dosdata,'total',sum_mode,dosids)
        else:
            key = None
            if not sum_mode:
                key = list(dosdata.keys())[dosid]
            (x,yu,yd) = self.get_x_yu_yd(dosdata,key,sum_mode,dosids)
        self.ax.clear()
        self.ax.plot(x,yu,color='red',linewidth=1)
        if self.dos.spin:
            self.ax.plot(x,yd,color='green',linewidth=1) 

        plt.xlim(self.ax.get_xbound()[0],self.ax.get_xbound()[1])
        set_range('erange',self.dosopts['erange'],plt.xlim)
        plt.ylim(self.ax.get_ybound()[0],self.ax.get_ybound()[1])
        default = None
        if not self.dos.spin:
            default = '0,'
        set_range('drange',self.dosopts['drange'],plt.ylim,default_if_None=default)
        x0 = np.array([self.ax.get_xbound()[0],self.ax.get_xbound()[1]])
        y0 = np.array([0,0])
        self.ax.plot(x0,y0,linestyle='-',color='black',linewidth=0.5)
        if self.dosopts['with_fermi']:
            x0 = np.array([0,0])
            y0 = np.array([self.ax.get_ybound()[0],self.ax.get_ybound()[1]])
            self.ax.plot(x0,y0,linestyle='--',color='black',linewidth=0.5)

        plt.xlabel('Energy (eV)')
        plt.ylabel('DOS (states/eV)')

    def gen_heatmap(self,mode='layer'):
        self.in_heatmap = True
        self.ax.clear()
        ldosar = None
        if mode=='layer':
            ldosar = self.dos.get_layerdos_heatmap_data()
        elif mode=='atom':
            ldosar = self.dos.get_aldos_heatmap_data()
        x = np.array(ldosar[0])
        y = np.array(ldosar[1])
        sv = self.heatmap_spin_vars
        if mode=='atom':
            sv = self.heatmap_spin_vars_atm
        if sv.get() == 'up':
            z = np.array(ldosar[2])
        else:
            z = np.array(ldosar[3])
        minv = None
        maxv = None
        rng = self.dosopts['drange']
        if rng is not None and len(rng)>0:
            try:
                if rng.endswith(','):
                    minv = float(rng[0:-1])
                elif rng.startswith(','):
                    maxv = float(rng[1:len(rng)])
                else:
                    minmax = rng.split(',')
                    minv = float(minmax[0].strip())
                    maxv = float(minmax[1].strip())
            except (ValueError,IndexError):
                logger.warning('invalid dosmin/dosmax; will be ignored')
        self.ax.contourf(x,y,z,self.dosopts['level'],cmap=self.dosopts['cmap'],vmin=minv,vmax=maxv)
        set_range('erange',self.dosopts['erange'],plt.ylim)
        set_range('lrange',self.dosopts['lrange'],plt.xlim)
        if mode=='layer':
            plt.xlabel(r'distance along the axis perpendicular to the layers ($\mathrm{\AA}$)')
        elif mode=='atom':
            plt.xlabel(r'atom ID')
        plt.ylabel('Energy (eV)')

    def get_x_yu_yd(self,dosdata,key,sum_mode,dosids):
        x = self.dos.get_energy()
        yu = None
        yd = None
        if not sum_mode:
            yu = np.array(dosdata[key]['up'])
            yd = None
            if self.dos.spin:
                down = []
                for d in dosdata[key]['down']:
                    down.append(-1.0*d)
                yd = np.array(down)
        else:
            keys = list(dosdata.keys())
            yu = []
            for d in dosdata[keys[0]]['up']:
                yu.append(0.0)
            if self.dos.spin:
                yd = []
                for d in dosdata[keys[0]]['down']:
                    yd.append(0.0)
            for dosid in dosids:
                dosd = dosdata[keys[dosid]]['up']
                for i in range(len(dosd)):
                    yu[i] += dosd[i]
            if self.dos.spin:
                for dosid in dosids:
                    dosd = dosdata[keys[dosid]]['down']
                    for i in range(len(dosd)):
                        yd[i] -= dosd[i]
        return (x,yu,yd)
        

def boot_gui(dos,options):
    gui = DosGUI(dos,options)
                        
def get_options():
    parser = optparse.OptionParser(usage='\n%prog [OPTION] ',version='%prog 1.0')
    parser.add_option('-i','--interactive',action='store_true',dest='interactive',\
    help='specify this option in order to run this script interactively',default=False)
    parser.add_option('-g','--gui',action='store_true',dest='gui',\
    help='specify this option in order to boot the GUI',default=False)
    parser.add_option('-l','--loglevel',dest='loglevel',type='choice', choices=['1','2','3'], \
    help='specify the loglevel. 1: INFO, 2: DEBUG, 3: WARN; defaults to 1',default='1')
    parser.add_option('-f','--file',dest='file',type=str,help='the DOS file [default : dos.data]',default='dos.data')
    parser.add_option('--output000',dest='output000',type=str,help='the outputxxx file [default : the most recent outputxxx file]',default=None)

    parser.add_option('-m','--mode',dest='mode',type=str,\
    help='specify the run mode. one of total, atom, layer, projected. multiple \
    modes can be specified by a comma-separated string. defaults to total.',default='total')
    parser.add_option('-a','--action',dest='action',type=str,\
    help='specify the \'action\'. action can be one of analyze, sum, or split. multiple actions may be specified by a comma-seperated string.',\
    default='split')
    parser.add_option('--heatmap',dest='heatmap',action='store_true',\
    help='specify this option in order to create a heatmap for the layer DOS.',default=False)
    parser.add_option('-o','--output_action',dest='output_action',type='choice', choices=['genfig','storedata','both'], \
    help='specify the output action. one of genfig, storedata or both. defaults to both',default='both')

    filteropts = optparse.OptionGroup(parser,'filter options','options used to select a subset of the DOS')
    parser.add_option_group(filteropts)
    filteropts.add_option('--dosid',dest='dosid',\
    help='specify the ids of the DOS to sum; for example, 3-5,7 will be interpereted as 3,4,5,7. does not make sense when --mode=total ()',default=None)
    filteropts.add_option('--atomid',dest='atomid',\
    help='specify the atom ids. this is the same as dosid for aldos, but not for pdos.',default=None)
    filteropts.add_option('--layerid',dest='layerid',\
    help='specify the layer ids. this is the same as dosid for layerdos.',default=None)
    filteropts.add_option('--elemid',dest='elemid',help='specify the element ids.',default=None)
    filteropts.add_option('--lid',dest='lid', \
    help='specify the l ids. can be either an integer or a string; namely 0 or s, 1 or p, 2 or d and 3 or f.',default=None)
    filteropts.add_option('--mid',dest='mid',help='specify the m ids.',default=None)
    filteropts.add_option('--tid',dest='tid',help='specify the t ids.',default=None)

    figopts = optparse.OptionGroup(parser,'figure options','options for configuring the output figure')
    parser.add_option_group(figopts)
    figopts.add_option('-e','--erange',dest='erange',\
    help='specify the energy range in the form emin,emax',default=None)
    figopts.add_option('--einc',dest='einc',help='energy tics',type=float,default=None)
    figopts.add_option('-d','--drange',dest='drange',\
    help='specify the DOS range in the form dmin,dmax',default=None)
    figopts.add_option('--dinc',dest='dinc',help='DOS tics',type=float,default=None)
    figopts.add_option('--lrange',dest='lrange',\
    help='specify the layer range in the form lmin,lmax',default=None)
    figopts.add_option('--arange',dest='arange',\
    help='specify the atom ID range in the form amin,amax',default=None)
    figopts.add_option('--linc',dest='linc',help='layer tics',type=float,default=None)
    figopts.add_option('--with_fermi',dest='with_fermi',action='store_true',\
    help='specify this option in order to draw a vertical line at the Fermi level',default=False)
    figopts.add_option('--title',dest='title',action='store_true',\
    help='specify this option in order to set a title in a plot',default=False)
    figopts.add_option('--level',dest='level',type=int,help='specify the contour level for the heatmap',default=100)
    figopts.add_option('--cmap',dest='cmap', \
    help='specify the color map to be used in the layer dos heatmap plot.\
    refer to https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html for possible choices',default='viridis')
    figopts.add_option('--imgtype',dest='imgtype',type='choice', choices=['eps','ps','png','jpg','pdf','svg'],\
    help='specify the output figure file type. one of : eps, ps, png, jpg, pdf or svg. defaults to eps',default='eps')

    (options,args) = parser.parse_args()
    return options

def print_dos(dos,filter=None):
    dkeys = dos.keys()
    (lines,cols) = pyutil.getTerminalSize()
    ncol=int(cols/45)
    if ncol==0:
        ncol=1
    icount=0
    thestr = ''
    iketa = int(math.log10(len(dkeys)))+1 
    dketa = 0
    for dkey in dkeys:
        if len(dkey)>dketa:
            dketa = len(dkey)

    i = 0
    ip = 1
    retmap = None
    if filter is not None:
        retmap = {}
    for dkey in dkeys:
        if filter is not None and not ip in filter:
            ip += 1
            continue
        if filter is not None:
            retmap[i+1] = ip
        thestr += (str(i+1)+') ').rjust(iketa+2)+dkey.ljust(dketa+1)
        icount += 1
        if icount%ncol==0:
            icount=0
            thestr += '\n'
        i += 1
        ip += 1
    print(thestr)
    return retmap

def split_dos(dos,mode):
    energy = dos.get_energy()
    thedos = dos.get_target_dos(mode)
    keys = thedos.keys()
    ret = []
    for key in keys:
        ds = thedos[key]
        if 'down' in ds:
            ret.append((key,key,ds['up'],ds['down']))
        else:
            ret.append((key,key,ds['up']))
    return ret

def sum_dos(mode,dos,dosid,atomid,layerid,elemid,lid,mid,tid):
    energy = dos.get_energy()
    nelem = dos.get_num_data()
    dosids = dos.get_idarray(mode,dosid,atomid,layerid,elemid,lid,mid,tid,nelem)
    logger.debug('dosids : '+str(dosids))
    (up,down) = dos.get_summed_dos(dosids,mode)
    return ('summed_'+mode,mode+' '+str(dosids),up,down)

def analyze_dos(dos):
    if dos.dos_attributes is not None:
        logger.info('DOS attributes')
        print(str(dos.dos_attributes))
    tdos = dos.get_total_dos()
    aldos = dos.get_aldos()
    ldos = dos.get_layer_dos()
    pdos = dos.get_pdos()
    if tdos is not None:
        logger.info('Total DOS')
    if aldos is not None:
        logger.info('atomic local DOS')
        print_dos(aldos)
    if ldos is not None:
        logger.info('layer DOS')
        print_dos(ldos)
    if pdos is not None:
        logger.info('projected DOS')
        print_dos(pdos)

def output_layerdos_heatmap(ldos,spin):
    x = ldos[0]
    y = ldos[1]
    z = ldos[2]
    fname = 'layerdos_heatmap.data'
    if spin:
        fname = 'layerdos_heatmap_up.data'
    f=open(fname,'w')
    for i in range(len(x)):
        for j in range(len(y)):
            f.write(str(x[i])+' '+str(y[j])+' '+str(z[j][i])+'\n')
        f.write('\n')
    f.close()
    if spin:
        z = ldos[3]
        fname = 'layerdos_heatmap_down.data'
        f=open(fname,'w')
        for i in range(len(x)):
            for j in range(len(y)):
                f.write(str(x[i])+' '+str(y[j])+' '+str(z[j][i])+'\n')
            f.write('\n')
        f.close()

def output_aldos_heatmap(ldos,spin):
    x = ldos[0]
    y = ldos[1]
    z = ldos[2]
    fname = 'aldos_heatmap.data'
    if spin:
        fname = 'aldos_heatmap_up.data'
    f=open(fname,'w')
    for i in range(len(x)):
        for j in range(len(y)):
            f.write(str(x[i])+' '+str(y[j])+' '+str(z[j][i])+'\n')
        f.write('\n')
    f.close()
    if spin:
        z = ldos[3]
        fname = 'aldos_heatmap_down.data'
        f=open(fname,'w')
        for i in range(len(x)):
            for j in range(len(y)):
                f.write(str(x[i])+' '+str(y[j])+' '+str(z[j][i])+'\n')
            f.write('\n')
        f.close()

def output_dos(energy,dosdata,spin,fname=None):
    name = dosdata[0]
    header = dosdata[1]
    up = dosdata[2]
    if spin:
        down = dosdata[3]
    nelem = len(energy)
    if fname is None:
        f = open('dos_'+name+'.data','w')
    else:
        f = open(fname,'w')
    if header is not None:
        f.write('#'+str(header)+'\n')
    if spin:
        for i in range(nelem):
            f.write(str(energy[i])+' '+str(up[i])+' '+str(down[i])+'\n')
    else:
        for i in range(nelem):
            f.write(str(energy[i])+' '+str(up[i])+'\n')
    f.close()

def set_range(name,rng,lim,default_if_None=None):
    if rng is None or len(rng)==0:
        rng = default_if_None
    if rng is not None and len(rng)>0:
        try:
            if rng.endswith(','):
                minv = float(rng[0:-1])
                lim([minv,lim()[1]])
            elif rng.startswith(','):
                maxv = float(rng[1:len(rng)])
                lim([lim()[0],maxv])
            else:
                minmax = rng.split(',')
                minv = float(minmax[0].strip())
                maxv = float(minmax[1].strip())
                lim([minv,maxv])
        except (ValueError,IndexError):
            logger.warning('invalid '+name+'; will be ignored')

def set_ticks(inc,ticks,lim):
    if inc is not None:
        ticks(np.arange(lim()[0], lim()[1], inc))

def plot_heatmap_sub(fname,x,y,z,options,mode='layer'):
    minv = None
    maxv = None
    rng = options.drange
    if rng is not None:
        try:
            if rng.endswith(','):
                minv = float(rng[0:-1])
            elif rng.startswith(','):
                maxv = float(rng[1:len(rng)])
            else:
                minmax = rng.split(',')
                minv = float(minmax[0].strip())
                maxv = float(minmax[1].strip())
        except (ValueError,IndexError):
            logger.warning('invalid dosmin/dosmax; will be ignored')
    fig, ax = plt.subplots()
    if options.level is None:
        ax.contourf(x,y,z,cmap=options.cmap,vmin=minv,vmax=maxv)
    else:
        ax.contourf(x,y,z,options.level,cmap=options.cmap,vmin=minv,vmax=maxv)

    set_range('erange',options.erange,plt.ylim)
    if mode=='layer':
        set_range('lrange',options.lrange,plt.xlim)
    elif mode=='atom':
        set_range('arange',options.arange,plt.xlim)
    set_ticks(options.einc,plt.yticks,plt.ylim)
    set_ticks(options.linc,plt.xticks,plt.xlim)
    if options.title and mode=='layer':
        plt.title('layer DOS heatmap')
    elif options.title and mode=='atom':
        plt.title('atomic local DOS heatmap')
    if mode=='layer':
        plt.xlabel(r'distance along the axis perpendicular to the layers ($\mathrm{\AA}$)')
    elif mode=='atom':
        plt.xlabel(r'atom ID')
    plt.ylabel('Energy (eV)')
    ftype = options.imgtype
    plt.savefig(fname+'.'+ftype,format=ftype)

def plot_heatmap(ldos,spin,options,mode='layer'):
    fname = 'layerdos_heatmap'
    if mode=='atom':
        fname = 'aldos_heatmap'
    if spin:
        fname = 'layerdos_heatmap_up'
        if mode=='atom':
            fname = 'aldos_heatmap_up'
    x = np.array(ldos[0])
    y = np.array(ldos[1])
    z = np.array(ldos[2])
    plot_heatmap_sub(fname,x,y,z,options,mode)
    if spin:
        plt.clf()
        z = np.array(ldos[3])
        fname = 'layerdos_heatmap_down'
        if mode=='atom':
            fname = 'aldos_heatmap_down'
        plot_heatmap_sub(fname,x,y,z,options,mode)

def plot_dos(energy, dosdata, spin, options):
    x = np.array(energy)
    yup = np.array(dosdata[2])
    if spin:
        down = []
        for d in dosdata[3]:
            down.append(-1.0*d)
        ydown = np.array(down)
    fig, ax = plt.subplots()
    ax.plot(x,yup,color='red',linewidth=1)
    if spin:
        ax.plot(x,ydown,color='green',linewidth=1) 

    plt.xlabel('Energy (eV)')
    plt.ylabel('DOS (states/eV)')

    set_range('erange',options.erange,plt.xlim)
    plt.ylim(ax.get_ybound()[0],ax.get_ybound()[1])
    default = None
    if not spin:
        default = '0,'
    set_range('drange',options.drange,plt.ylim,default)

    set_ticks(options.einc,plt.xticks,plt.xlim)
    set_ticks(options.dinc,plt.yticks,plt.ylim)

    if options.title:
        plt.title(dosdata[1])
    if options.with_fermi:
        x0 = np.array([0,0])
        y0 = np.array([ax.get_ybound()[0],ax.get_ybound()[1]])
        ax.plot(x0,y0,linestyle='--',color='black',linewidth=0.5)
    ftype = options.imgtype
    name = dosdata[0]
    plt.savefig('dos_'+name+'.'+ftype,format=ftype)
    plt.close()

def valid_idstr(idstr,size):
    if idstr is None or len(idstr.strip())==0:
        return False
    idstr = idstr.strip()
    if idstr.lower()=='all':
        return True
    idstr = idstr.replace(' ','')
    if idstr.endswith(','):
        idstr = idstr[0:len(idstr)-1]
    css = idstr.split(',')
    for cs in css:
        hss = cs.split('-')
        if len(hss)>2:
            return False
        if len(hss)==1:
            try:
                it = int(hss[0])
                if it<1 or it>size:
                    return False
            except ValueError:
                return False
        elif len(hss)==2:
            try:
                mini = int(hss[0])
                maxi = int(hss[1])
                if mini<1:
                    return False 
                if maxi>size:
                    return False
            except ValueError:
                return False
    return True

def get_real_idstr(mapid,cand_ids,idstr):
    ret = ''
    if idstr == 'all':
        for i in cand_ids:
            ret += str(i)+','
        return ret[0:len(ret)-1]
    intids = DosData.get_intidarray(idstr,-1)
    for i in intids:
        ret += str(mapid[i])+','
    return ret[0:len(ret)-1]

def range_check(x):
    if len(x.strip())==0:
        return True
    x = x.strip()
    if x.startswith(','):
        try:
            float(x[1:len(x)])
            return True
        except ValueError:
            return False
    elif x.endswith(','):
        try:
            float(x[0:len(x)-1])
            return True
        except ValueError:
            return False
    else:
        try:
            minmax = x.replace(' ','').split(',')
            minv = float(minmax[0])
            maxv = float(minmax[1])
            if maxv>=minv:
                return True
        except (ValueError, IndexError):
            return False 
        return False
        
def interactive(dos,options):
    mchoices = ['total']
    pd = dos.get_pdos()
    if dos.get_aldos() is not None:
        mchoices.append('atom')
    if dos.get_layer_dos() is not None:
        mchoices.append('layer')
    if pd is not None:
        mchoices.append('projected')

    mode = pyutil.interactive('select mode',default=0,choice=mchoices)
    action = 0
    idstr = None
    dosdata = []
    if mode != 0:
        achoices = ['split','sum']
        if mchoices[mode] == 'layer' or mchoices[mode] == 'atom':
            achoices = ['split','sum','heatmap']
        action = pyutil.interactive('select action',default=0,choice=achoices)

    if action == 0:
        splitted = split_dos(dos,mchoices[mode])
        for spl in splitted:
            dosdata.append(spl)
    if action == 1:
        cand_ids = dos.get_idarray(mchoices[mode],None,options.atomid,options.layerid,\
                       options.elemid,options.lid,options.mid,options.tid,-1)
        mapid = None
        if mchoices[mode]=='atom':
            mapid = print_dos(dos.get_aldos(),filter=cand_ids)
        if mchoices[mode]=='layer':
            mapid = print_dos(dos.get_layer_dos(),filter=cand_ids)
        if mchoices[mode]=='projected':
            mapid = print_dos(dos.get_pdos(),filter=cand_ids)
        idstr = pyutil.interactive('DOS to sum ',default='',condition=lambda x: valid_idstr(x,len(mapid.keys())))
        idstr = get_real_idstr(mapid,cand_ids,idstr)
        dosdata.append(sum_dos(mchoices[mode],dos,idstr,options.atomid,options.layerid,options.elemid,\
                                       options.lid,options.mid,options.tid))

    if options.output_action == 'both' or options.output_action == 'genfig':
        erange = pyutil.interactive('the energy range emin,emax',default='',condition=lambda x: range_check(x))
        options.erange = erange 
        drange = pyutil.interactive('the dos range dmin,dmax',default='',condition=lambda x: range_check(x))
        options.drange = drange 
        if action != 2:
            fermi = pyutil.interactive('should a verticle line be drawn at the Fermi lelvel?',choice=['yes','no'],default=0)
            if fermi==0: 
                options.with_fermi = True
            else:
                options.with_fermi = False
        if action == 2 and mchoices[mode] == 'atom':
            arange = pyutil.interactive('the atom ID range amin,amax',default='',condition=lambda x: range_check(x))
            options.arange = arange
        if action == 2 and mchoices[mode] == 'layer':
            lrange = pyutil.interactive('the layer range lmin,lmax',default='',condition=lambda x: range_check(x))
            options.lrange = lrange
 
        imchoice = ['eps','ps','png','jpg','pdf','svg']
        imgtype = pyutil.interactive('the image type',default=0,choice=imchoice)
        options.imgtype = imchoice[imgtype]

    if options.output_action == 'both' or options.output_action == 'storedata':
        if action == 2:
            if mchoices[mode] == 'layer':
                ldosar = dos.get_layerdos_heatmap_data()
                output_layerdos_heatmap(ldosar,dos.spin)
            elif mchoices[mode] == 'atom':
                ldosar = dos.get_aldos_heatmap_data()
                output_aldos_heatmap(ldosar,dos.spin)
        else:
            for dosd in dosdata:
                output_dos(dos.get_energy(),dosd,dos.spin)

    if options.output_action == 'both' or options.output_action == 'genfig':
        if action == 2:
            if mchoices[mode] == 'layer':
                ldosar = dos.get_layerdos_heatmap_data()
                plot_heatmap(ldosar,dos.spin,options,mode=mchoices[mode])
            elif mchoices[mode] == 'atom':
                ldosar = dos.get_aldos_heatmap_data()
                plot_heatmap(ldosar,dos.spin,options,mode=mchoices[mode])
        else:
            for dosd in dosdata:
                plot_dos(dos.get_energy(),dosd,dos.spin,options)

def run():
    options = get_options()
    init_logger(options.loglevel)
    output000 = options.output000
    if output000 is None:
        path = Path(options.file)
        output000 = get_mostrescent_outputxxx(directory=str(path.parent))
    dos = DosData(options.file,output000)
    logger.debug(str(dos))
    if options.interactive:
        interactive(dos,options)
        sys.exit()

    if options.gui:
        boot_gui(dos,options)
        sys.exit()
 
    actions = options.action.split(',')
    modes = options.mode.split(',')
    dosdata = []
    vmodes = ['total','atom','layer','projected']
    for mode in modes:
        if not mode in vmodes:
            logger.error('invalid mode : '+mode)
            continue
        for action in actions:
            if action == 'analyze':
                logger.info('analyzing '+options.file)
                analyze_dos(dos)
            elif action == 'split':
                logger.info('split DOS data and output the results to a file')
                logger.info('mode : '+mode)
                splitted = split_dos(dos,mode)
                for spl in splitted:
                    dosdata.append(spl)
            elif action == 'sum':
                if mode == 'total':
                    logger.error('the \'sum\' action can\'t be used with mode=total')
                else:
                    logger.info('sum DOS data and output the results to a file')
                    logger.info('mode  : '+mode)
                    dosdata.append(sum_dos(mode,dos,options.dosid,options.atomid,options.layerid,options.elemid,\
                                       options.lid,options.mid,options.tid))
            else:
                logger.error('unknown action : '+action)

        if options.output_action == 'both' or options.output_action == 'storedata':
            if options.mode=='layer' and options.heatmap:
                ldosar = dos.get_layerdos_heatmap_data()
                output_layerdos_heatmap(ldosar,dos.spin)
            elif options.mode=='atom' and options.heatmap:
                aldosar = dos.get_aldos_heatmap_data()
                output_aldos_heatmap(aldosar,dos.spin)
            else:
                for dosd in dosdata:
                    output_dos(dos.get_energy(),dosd,dos.spin)

        if options.output_action == 'both' or options.output_action == 'genfig':
            if options.mode=='layer' and options.heatmap:
                ldosar = dos.get_layerdos_heatmap_data()
                plot_heatmap(ldosar,dos.spin,options)
            elif options.mode=='atom' and options.heatmap:
                aldosar = dos.get_aldos_heatmap_data()
                plot_heatmap(aldosar,dos.spin,options,mode='atom')
            else:
                for dosd in dosdata:
                    plot_dos(dos.get_energy(),dosd,dos.spin,options)

if __name__ == '__main__':
    run()

