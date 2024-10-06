#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import piou

from piou.config import phase
import logging
import os

logger = logging.getLogger("inpcheck")

class InputValidator(piou.PHASE_IO_utility):

    dirs_checked = []
    specfile_specified=False

    def fill_opts(self,parser):
        parser.add_option("-s","--specfile",dest="specfile",\
        help="specify the name of the specification file.",default=None)
        parser.add_option("-r","--recursive",action="store_true",dest="recursive",\
        help="set this option in order to apply inpcheck recusively.",default=False)

    def initialize(self,opts):
        self.recursive=opts.recursive

    def __init__(self):
        self.program = 'phase'

    def do_inpcheck(self):
        logger.info("-- running the input validator --")
        logger.info("specfile : "+self.program+".spec")
        self._inpcheck(".",self.recursive)
        if len(self.dirs_checked)>=2:
            logger.debug("checked a total of : "+str(len(self.dirs_checked))+" directories")
            for dir in self.dirs_checked:
                logger.debug(" "+dir)

    def _inpcheck(self, dir, recur):
        os.chdir(dir)
        logger.info("checking directory : "+os.getcwd())
        files=os.listdir(".")
        sfile=None
        if not self.specfile_specified:
            sfile=self.program+'.spec'
            specfile_specified=True
        #if 'file_names.data' in files:
        self.dirs_checked.append(os.getcwd())
        inp = phase.Input(sfile,coord_file=None)
        if inp.inpfile_exists():
            try:
                inp.validate()
            except phase.InputValidationError:
                logger.error("input validation FAILED!")
                #pass
                #logger.error("input validation FAILED!")
            nwarn = inp.get_nwarn()
            nerror = inp.get_nerror()
            if nwarn == 0 and nerror == 0:
                logger.info("no errors/warnings were found in "+os.getcwd())
            elif nerror==0:
                logger.warn("found "+str(nwarn)+" warnings in "+os.getcwd())
                logger.warn("check the log for details.")
            else:
                logger.error("found "+str(nerror)+" errors and "+str(nwarn)+" warnings in "+os.getcwd())
                logger.error("check the log for details.")
        elif not recur:
            logger.error(" input file ["+inp.get_file_name()+"] does not exist under "+os.getcwd())

        if recur:
            for file in files:
                if os.path.isdir(file):
                    self._inpcheck(file,True)
            os.chdir("..")
            logger.info("")

    def get_description(self):
        return \
        "input data validator utility for PHASE"

    def get_name(self):
        return 'the input file validator'

    def run(self):
        self.do_inpcheck()

if __name__ == "__main__":
    piou.run_main(InputValidator)

