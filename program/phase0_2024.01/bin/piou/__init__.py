
import logging
logger=logging.getLogger('piou')
import sys

if sys.hexversion >= 0x02030000: #Python 2.3 and over
    import logging
    import optparse
else:
    from piou.ext import logging
    from piou.ext import optparse
import os

import piou.data.properties

def initialize_logger_and_opts(pioutil):
    parser = optparse.OptionParser(\
    usage="%prog [options]",\
    description=\
    "The PHASE I/O utility script."+ "\n"+\
    "This script will enable you to :\n"+\
    " 1. validate your PHASE input, \n"+\
    " 2. generate default input file from atomic coordinates and \n"+\
    " 3. convert PHASE input/output data for use with other programs, and vice-versa.\n"+\
    "refer to the user's manual for details.",\
    version="%prog 0.5")

    opts = get_opts(parser,pioutil)
    init_logger(opts.loglevel)
    if opts.propfile is not None:
        piou.data.properties.Property().set_propfile_from_args(opts.propfile)
    if opts.ppdir is not None:
        piou.data.properties.Property().set_property("pp.default_pp_dir",opts.ppdir)
    logger.debug("pp directory : "+piou.data.properties.Property().get_property("pp.default_pp_dir"))

    return opts

def get_opts(parser,pioutil):
    parser.add_option("-l", "--loglevel", type="choice",choices=("0","1","2"),\
                      dest="loglevel", help="specify the log level (one of 0,1,2)",default="1")
    parser.add_option("", "--prop",dest="propfile",help=\
   "specify the property file. this will overload the default settings",default=None)
    parser.add_option("--ppdir",dest="ppdir",help="specify the directory where the pp files reside."\
    +" if unspecified, the string specified in the property file will be used.",default=None)
    parser.add_option("-b","--batch",action="store_true",\
                      help="specify this option in order to perform this script in batch mode (for geninp and conv). ",\
                      dest="batch",default=False)

    pioutil.fill_opts(parser)

    (options, args) = parser.parse_args()
    return options

logformat  = '%(levelname)10s: %(message)s'
logformat_dbg = '%(asctime)8s %(name)20s %(lineno)6d %(levelname)10s: %(message)s'
logtimefmt = '%a %d %b %H:%M:%S'
logger  = logging.getLogger("main")

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

def run_main(pioutil):
    inst = pioutil()
    print(inst.get_description())
    print('Copyright (C) PHASE System Consortium')
    opts=initialize_logger_and_opts(inst)
    inst.initialize(opts)
    inst.run()

class PHASE_IO_utility:
    ''' the base class for the main script of the various phase io utility '''

    def __init__(self):
        pass

    def get_name(self):
        ''' return the name of this utility '''
        raise NotImplementedError()

    def get_description(self):
        ''' return the description for this utility'''
        raise NotImplementedError()

    def fill_opts(self,parser):
        ''' fill the optsions specific for this phase io utility '''
        raise NotImplementedError()

    def initialize(self,opts):
        ''' initialize using the specified options '''
        raise NotImplementedError()

    def run(self):
        ''' run the main script'''
        raise NotImplementedError()

