''' various utility functions and classes '''

import logging
logger = logging.getLogger("piou.util.pyutil")

import os
import sys

try:
    import readline
except ImportError:
    pass
else:
    readline.parse_and_bind("tab: complete")

try:
    input = raw_input
except NameError:
    pass

class Callable:
    ''' use this class to generate so-called class methods.
        here's an example:

        Class Something
            ...
            ...
            def do_something():
                ...
                ...

            from psf.util import pyutil
            do_something = pyutil.Callable(do_something)

        Something.do_something() will then behave like a class method in Java. '''
    def __init__(self, anycallable):
        self.__call__ = anycallable

class Singleton(type):
    ''' python-way of realizing the Singleton pattern. By doing
        class Boo():
            __metaclass__ = pyutil.Singleton
            def ...
        class Boo will behave as a Singleton. '''
    def __init__(self, *args):

        type.__init__(self, *args)

        self._instance = None

    def __call__(self, *args):
        if self._instance is None :
            self._instance = type.__call__(self, *args)
        return self._instance

def my_import(name,ignore_error=False):
    ''' get module under the specified package. '''
    mod = None
    try:
        mod = __import__(name)
    except Exception:
        if not ignore_error:
            logger.error('import error!')
        return None

    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

def map_package_to_dir(package_name):
    '''
    obtain the path to a specific package.
    will return None if no such package exists.
    '''
    boo = my_import(package_name,True)
    if boo != None:
        return boo.__path__[0]
    logger.error('  could not find package : '+package_name)
    return None

def get_all_modules(packageName):
    ''' obtain a list of all modules under the specified package
        parameter : packageName the name of the desired package
        will return None if packageName does not exist. '''
    import glob
    import os,re
    dir = map_package_to_dir(packageName)
    if dir==None:
        return None
    files = os.listdir(dir)
    reg = re.compile(".*\.py$",re.IGNORECASE)
    files = filter(reg.search, files)
    mods = []
    for file in files:
        mods.append(file[:-3])

    ret = []
    for mod in mods:
        if mod!='__init__':
            ret.append(my_import(packageName+'.'+mod,True))
    return ret

def get_all_methods(mod,prefix=''):
    ''' get all methods defined in the specified module.
        parameter : mod the desired module object '''
    ret =[]
    for method in dir(mod):
        attr = getattr(mod,method)
        if callable(attr):
            if type(attr)==type(get_all_methods):
                if method.startswith(prefix):
                    ret.append(attr)
    return ret

def execute(commands,currdir=None,log_offset=''):
    ''' execute commands according to the instructions given by the arg.
        the instruction object is a simple list. here's an example:

        from psf.util import pyutil
        commands=[['environ', PATH boo],['command', bar],['command', hoge]]
        pyutil.execute(commands)

        this example will first set the environment variable PATH to boo,
        then execute the command bar, and will finally execute the command hoge.
        The standard output will be redirected
        to the logger's info method, while the standard error output will be
        redirected to the logger's error method. '''
    import os
    import time
    import subprocess
    for command in commands:
        if command[0].strip()== 'environ':
            line = command[1].split()
            if len(line)==2:
                os.environ[line[0]] = line[1]
                logger.info(log_offset+"environment variable "+line[0]+" was set to "+line[1])

    for command in commands:
        if command[0].strip() == 'command':
            if currdir==None:
                logger.info(log_offset+"executing command : ["+command[1]+']')
            else:
                logger.info(log_offset+'executing command : ['+command[1]+'] under '+str(currdir))
            try:
                popen=subprocess.Popen(command[1],cwd=currdir,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
                while 1:
                    line = popen.stderr.readline()
                    if not line: break
                    logger.error(log_offset+line)
                while 1:
                    line=popen.stdout.readline()
                    if not line: break
                    logger.info(log_offset+line)

            except OSError:
                logger.error(log_offset+'failed execution of : '+command[1])

def get_all_classes(mod):
    ''' get all classes defined under the specified module.
        parameter : mod the desired module object '''
    ret=[]
    for clas in dir(mod):
        attr = getattr(mod,clas)
        if callable(attr):
            if type(getattr(mod,clas))==type(OSError):
                ret.append(attr)
    return ret

def get_homedir():
    ''' obtain the home directory of the system.
        will return $HOME for UNIX, %userprofile% for WinNT.
        won't work for Win9x (but who cares!) '''
    import os
    if os.name=='nt':
        return os.environ["userprofile"]
    return os.environ["HOME"]

def parse_bool(str,default=None):
    ''' convert string value to a Boolean. The rools are:
        None or len(str)==0 -> False if 'default' is None, or 'default'
        'true','yes','on','1' -> True
        'false','no','off','0' -> False
        defaults to False. '''
    if isinstance(str,bool):
        return str
    if str==None or len(str)==0:
        if default is None or not isinstance(default,bool):
            return False
        else:
            return default
    if str.lower()== 'yes':
        return True
    if str.lower()=='no':
        return False
    if str.lower()=='on':
        return True
    if str.lower()=='off':
        return False
    if str=='0':
        return False
    if str=='1':
        return True
    if str.lower()=='true':
        return True
    if str.lower()=='false':
        return False
    return False

def is_bool(str):
    ''' determine whether the string given in arg
        is of type boolean.
        parameter : str the string to be tested.'''
    str = str.lower()
    return str=='yes' or str=='no' or str=='0' or str=='1' or str=='on' or str=='off' or str=='true' or str=='false'

def get_ws(length):
    ''' return white-space characters '''
    ret =''
    for i in range(length):
        if i==0:
            continue
        ret += '  '
    return ret
    
def pathsplit(p, rest=[]):
    (h,t) = os.path.split(p)
    if len(h) < 1: return [t]+rest
    if len(t) < 1: return [h]+rest
    return pathsplit(h,[t]+rest)

def commonpath(l1, l2, common=[]):
    if len(l1) < 1: return (common, l1, l2)
    if len(l2) < 1: return (common, l1, l2)
    if l1[0] != l2[0]: return (common, l1, l2)
    return commonpath(l1[1:], l2[1:], common+[l1[0]])

def relpath(p1, p2):
    if p1.endswith(os.path.sep):
        pl = pl[:len(pl)-1]
    if p2.endswith(os.path.sep):
        p2 = p2[:len(p2)-1]
    (common,l1,l2) = commonpath(pathsplit(p1), pathsplit(p2))
    p = []
    if len(l1) > 0:
        p = [ '../' * len(l1) ]
    p = p + l2
    return os.path.join( *p )

def interactive(msg,condition=None,typ=None,default=None,choice=None,return_index=True,skippable=False,\
    msg_invalid_input=None,ask_mode=False,ask_mode_msg=None):
    default_str='['+str(default)+']: '
    if len(str(default))==0:
        default_str = ': '
    while True:
        print
        line = None
        if choice is not None:
            print(msg)
            prestr="("
            for ind in range(len(choice)):
                choi = choice[ind]
                print(str(str(ind)+". ").ljust(4)+choi)
                if ind <len(choice)-1:
                    prestr+=str(ind)+"/"
                else:
                    if not skippable:
                        prestr+=str(ind)+"/x) "+default_str
                    else:
                        prestr+=str(ind)+"/a/x) "+default_str
            if skippable:
                print('a.'.ljust(4)+'apply the default settings here after')
            print('x.'.ljust(4)+'Exit')
            #line = raw_input('Please enter a selection '+prestr)
            line = input('Please enter a selection '+prestr)
            if str(line).strip()=='x':
                sys.exit()
            if skippable and str(line).strip()=='a':
                return None
            if len(str(line).strip())==0:
                if not return_index:
                    return choice[default]
                else:
                    return default
            ichoi=0
            try:
                ichoi=int(str(line).strip())
            except ValueError:
                print('Enter an integer')
                continue
            if ichoi<0 or ichoi>=len(choice):
                print('Invalid choice : '+str(ichoi))
                continue
            if not return_index:
                return choice[ichoi]
            else:
                return ichoi
        else:
            msgstr = "Please enter "+msg+", or type x to exit. "+default_str 
            if skippable:
                msgstr = "Please enter "+msg+\
                ", or type a to apply the default settings here after, or type x to exit. "+default_str 
            #line = raw_input(msgstr)
            line = input(msgstr)
            print(line)
            if str(line).strip()=='x':
                sys.exit()
            elif skippable and str(line).strip()=='a':
                return None
            if (len(line)==0 and default is not None) or condition is not None and condition(line):
                if len(line)==0 and default is not None:
                    line = default
                if not condition(line):
                    if msg_invalid_input is None:
                        print('invalid input : '+line)
                    else:
                        print(msg_invalid_input+line)
                        if ask_mode:
                            ll=''
                            if ask_mode_msg==None:
                                #ll=raw_input('overwrite? (y/n) [n] ')
                                ll=input('overwrite? (y/n) [n] ')
                            else:
                                #ll=raw_input(ask_mode_msg)
                                ll=input(ask_mode_msg)
                            if ll.lower().startswith('y'):
                                return str(line).strip()
                    continue
                if typ==int:
                    try:
                        val=int(str(line).strip())
                    except ValueError:
                        print('enter an int')
                        continue
                elif typ==float:
                    try:
                        val=float(str(line).strip())
                    except ValueError:
                        print('enter a float')
                        continue
                elif typ==bool:
                    val=parse_bool(line)
                    if val is None:
                        print('enter a bool')
                        continue
                else:
                    val=str(line).strip()
                return val
            if msg_invalid_input is None:
                print('invalid input : '+line)
            else:
                print(msg_invalid_input+line)
                if ask_mode:
                    ll=''
                    if ask_mode_msg==None:
                        #ll=raw_input('overwrite? (y/n) [n] ')
                        ll=input('overwrite? (y/n) [n] ')
                    else:
                        #ll=raw_input(ask_mode_msg)
                        ll=input(ask_mode_msg)
                    if ll.lower().startswith('y'):
                        return str(line).strip()
def getTerminalSize():
    def ioctl_GWINSZ(fd):
        try:
            import fcntl, termios, struct, os
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return None
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.envron['COLUMNS'])
        except:
            cr = (25, 160)
    return int(cr[0]), int(cr[1])

from threading import Thread
class TailF(Thread):
    ''' class which has the same functionality as 'tail -f'.
        !!!UNFINISHED!!! '''
    def __init__(self,filename):
        Thread.__init__(self)
        self.filename = filename
        self.bufsize = 5
        self.lines = []

    def run(self):
        import time, os
        import os.path
        file = open(self.filename,'rb')
        logger.debug("starting tail -f "+self.filename)
        while 1:
            where = file.tell()
            line = file.readline()
            if not line:
                time.sleep(1)
                file.seek(where)
            else:
                print(line,) # already has newline
                #self.lines.append[line]
                #if(len(self.lines)==bufsize):
                #    for l in self.lines:
                #        print l,
                #    self.lines = []
