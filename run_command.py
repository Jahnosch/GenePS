#!/usr/bin/env python3
import subprocess
import shutil
import sys
import contextlib
import os
import tempfile as tmp

'''
COMMENTS
- run_command.py does more than just run commands, give it a more general name
'''

def run_cmd(**kwargs):
    cmd = kwargs['command']
    if type(cmd) != list:
        cmd = kwargs['command'].split() # sanitation
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True,
                         bufsize=-1) # buffersize of system
    wait = kwargs.get('wait', False)
    if wait:
        p.wait()
        if p.returncode == 0:
            pass
    else:
        return iter(p.stdout.readline, '')


def which(program):
    if shutil.which(program):
        return program
    else:
        print(program + " not Found\n")
        sys.exit(1)


def check_programs(*arg):
    for program in arg:
        which(program)


@contextlib.contextmanager
def cd(newdir, cleanup=lambda: True):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
        cleanup()


@contextlib.contextmanager
def tempdir():
    dirpath = tmp.mkdtemp()

    def cleanup():
        shutil.rmtree(dirpath)
    with cd(dirpath, cleanup):
        yield dirpath
