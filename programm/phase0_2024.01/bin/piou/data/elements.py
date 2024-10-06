''' module for manipulating element info. '''
import logging
import re
from . import constants
from piou.util import pyutil
from os.path import join

logger = logging.getLogger("piou.data.elements")

class ElementInfo:
    __metaclass__ = pyutil.Singleton
    element_info = []
    to_au_mass = constants.atomic_mass/constants.au_mass

    def __init__(self):
        if len(self.element_info)==0:
            self.load_element_info()

    def load_element_info(self):
        dir = join(pyutil.map_package_to_dir('piou'),'data')
        eleminfo = dir+"/elementinfo"
        f = open(eleminfo,'r')
        ii = 1
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            words = line.split()
            if len(words)<7:
                continue
            elemdict={}
            elemdict['name'] = words[0]
            elemdict[ii]=words[0]
            elemdict['rname'] = re.compile(words[0]+"\d*$")
            elemdict['atomic_number'] = ii
            ii += 1
            try:
                elemdict['covalent_radius'] = float(words[2])*constants.Angstrom_2_Bohr
                elemdict['mass'] = float(words[6]) * self.to_au_mass
                elemdict['color'] = (float(words[3]),float(words[4]),float(words[5]))
                elemdict['radius'] = float(words[1])
            except ValueError:
                pass
            self.element_info.append(elemdict)
        f.close()
        logger.debug("loaded element info from "+eleminfo)
        for elem in self.element_info:
            logger.debug("element : "+elem['name']+" mass : "+str(elem['mass']))
    
    def get_element_attribute(self,element_name, attribute_name):
        ''' get the specified attribute of the specified element. Will return None 
            if the specified element or the corresponding attribute is undefined. '''
        elem = self.get_element_info(element_name)
        if elem is None:
            return None

        if not attribute_name in elem.keys():
            return None
        else:
            return elem[attribute_name]

    def set_element_attribute(self,element_name, attribute_name,value):
        ''' set the value of the specified attribute of the specified element. '''
        elem = self.get_element_info(element_name)
        if elem is not None:
            elem[attribute_name] = value

    def get_element_info(self,element_name):
        ''' get the specified element info. Will return None if no corresponding element info exist.'''
        if element_name is None:
            return None
        for elem in self.element_info:
            if elem['rname'].match(element_name):
                return elem

    def element_exists(self,element_name):
        ''' test whether the specified element exists. '''
        for elem in self.element_info:
            if elem['rname'].match(element_name):
                return True
        return False

    def get_element_name_from_atomic_number(self,atomic_number):
        for elem in self.element_info:
            if atomic_number in elem:
                return elem[atomic_number]

