#!/usr/bin/env python3

# Make dependencies for use in Makefile for all Fortran source files within list of directories

import sys
import glob
import os


class FortranFile:
    """Simple class to store names of dependencies"""
    def __init__(self, filename):
        self.filename = filename
        self.program_name = None
        self.units = None
        self.depends = None
        self.file_depends = set()


def warp_filename(srcfile, ext=''):
    base = os.path.basename(srcfile)
    return base.rpartition('.')[0] + ext

if len(sys.argv) == 1:
    dir_list = [os.getcwd()]
elif len(sys.argv) == 2:
    vpath = sys.argv[1]
    dir_list = vpath.split(":")
else:
    raise ValueError('Incorrect number of command line arguments!')

# Loop through directories identifying fortran files

filelist = []
for d in dir_list:
    filelist.extend(glob.glob(os.path.join(d, '*.F90')))
    filelist.extend(glob.glob(os.path.join(d, '*.f90')))
ff_list = []

for filename in filelist:
    """Scan through file, searching for modules or programs"""
    cur_module = None
    program_name = None
    local_modules = set()
    local_depends = set()
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line_start = line.lstrip().rstrip('\n').lower()
        if line_start.startswith('program ') or line_start=='program':
            # We only deal with one program unit
            if program_name is not None:
                print('File {}'.format(filename))
                raise ValueError('File has more than one program unit!')
            if cur_module is not None:
                print('File {}'.format(filename))
                raise ValueError('Cannot have programs inside of modules!')
            if '&' in line:
                print('File {}'.format(filename))
                raise ValueError('Continuation lines on program lines not supported...')
            cur_program = line.lstrip().split()[1]
            program_name = cur_program
        elif (line_start.startswith('end program ') or
              line_start.startswith('endprogram ') or
              line_start=='end program' or
              line_start=='endprogram'):
            if cur_program is None:
                print('File {}'.format(filename))
                raise ValueError('File uses "end program" without starting "program"!')
            if '&' in line:
                print('File {}'.format(filename))
                raise ValueError('Continuation lines on end program lines not supported...')
        elif line_start.startswith('module '):
            # We are entering a module
            if '&' in line:
                print('File {}'.format(filename))
                raise ValueError('Continuation lines on module lines not supported...')
            this_mod = line.lstrip().split()[1]
            cur_module = this_mod
            local_modules.add(this_mod)
        elif line_start.startswith('end module ') or line_start.startswith('endmodule '):
            if not cur_module:
                print('File {}'.format(filename))
                raise ValueError('File uses "end module" without starting "module"!')
            if '&' in line:
                print('File {}'.format(filename))
                raise ValueError('Continuation lines on end module lines not supported...')
            cur_module = None
        elif line_start.startswith('use ') or line_start.startswith('use,'):
            # We have a potential dependency; first check for intrinsic module
            next_part = line_start[3:].lstrip(' ,')
            if next_part.startswith('intrinsic'):
                # Ignore intrinsic modules
                continue
            if next_part.startswith('non_intrinsic'):
                # trim this bit
                next_part = line_start.partition('non_intrinsic')[2].lstrip(' :')
            else:
                next_part = line_start[3:].lstrip(' :')
            if ',' in next_part:
                next_part = next_part.partition(',')[0]
            # Add dependency
            mod_name = next_part.strip()
            local_depends.add(mod_name)
            
    # Scan through dependencies removing 'local' modules
    local_depends -= local_modules
    
    # FortranFile object
    ff = FortranFile(filename)
    ff.units = local_modules
    ff.depends = local_depends
    if program_name is not None:
        ff.program_name = program_name
    ff_list.append(ff)
    
# Now scan through list of FortranFile objects without programs, matching each file's dependencies to other files
for mod in ff_list:
    for dep in mod.depends:
        for mod2 in ff_list:
            if mod is mod2:
                continue
            if dep in mod2.units:
                mod.file_depends.add(mod2.filename)

#for mod in ff_list:
    #print('ff object: ', mod.filename)
    #print('program name: ', mod.program_name)
    #print('units: ', mod.units)
    #print('dependencies: ', mod.depends)
    #print('file dependencies: ', mod.file_depends)
    #print('----------------')

with open('fort.dep','w') as f:
    f.write('\n')
    for mod in ff_list:
        if mod.program_name is None:
            item = warp_filename(mod.filename, '.o')
        else:
            item = warp_filename(mod.filename)
        if not mod.file_depends:
            continue
        depends_str = ' '.join((warp_filename(fname, '.o') for fname in mod.file_depends))
        f.write('{} : {}\n\n'.format(item, depends_str))
