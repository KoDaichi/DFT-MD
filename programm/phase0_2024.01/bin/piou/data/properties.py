import logging
import os

from os.path import join
from os.path import expanduser
from piou.util import pyutil

logger = logging.getLogger("piou.data.properties")

class Property:
    __metaclass__ = pyutil.Singleton
    props = {}
    propfile_from_arg=''
    def __init__(self):
        pass

    def set_propfile_from_args(self,p):
        self.propfile_from_arg=p

    def set_property(self,property_name,val):
        ''' set the specified property '''
        if len(self.props.keys())==0:
            self.load_props()
            logger.debug("properties defined:")
            for k,v in self.props.items():
                logger.debug(" "+k+" = "+v)
        self.props[property_name] = val

    def get_property(self,property_name,valtype="string"):
        ''' obtain the specified property. the properties can be defined in :
           ~./chase/chase.ini, ~/.chase/acvconfig or ~/.piou

            the method will try to convert the obtained string according to 'valtype'.
            'valtype' can be one of None or 'string' (no conversion), 'float', 'int', or 'bool'
            if conversion fails, will return None for 'float' and 'int', False for 'bool'.  '''
        if len(self.props.keys())==0:
            self.load_props()
            logger.debug("properties defined:")
            for k,v in self.props.items():
                logger.debug(" "+k+" = "+v)
        if property_name in self.props:
            if valtype is None or valtype=='string':
                return self.props[property_name]
            elif valtype == 'float' or valtype == float:
                try:
                    return float(self.props[property_name])
                except ValueError:
                    logger.error("property : ["+property_name+" = "+\
                                 self.props[property_name]+"] is not a float")
                    return None
            elif valtype == 'int' or valtype == int:
                try:
                    return int(self.props[property_name])
                except ValueError:
                    logger.error("property : ["+property_name+" = "+\
                                 self.props[property_name]+"] is not an int")
                    return None
            elif valtype == 'bool' or valtype==bool:
                return pyutil.parse_bool(self.props[property_name])
            
        return None

    def load_props(self):
        dir = join(pyutil.map_package_to_dir('piou'),'data')
        propfile = dir+"/props.properties"
        self.load_props_sub(propfile)
    
        homedir = expanduser("~")
        self.load_props_sub(homedir+"/.chase/chase.ini")
        self.load_props_sub(homedir+"/.chase/acvconfig")
        self.load_props_sub(homedir+"/.piou")
        if len(self.propfile_from_arg)!=0:
            self.load_props_sub(self.propfile_from_arg)

    def load_props_sub(self,filename):
        if not os.path.exists(filename):
            logger.debug("property file : "+filename+" does not exist.")
            return
        f = open(filename,'r')
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            words=line.split('=')
            if len(words)>=2:
                self.props[words[0].strip()] = words[1].strip()
        f.close()

