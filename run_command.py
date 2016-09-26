#!/usr/bin/env python3
import subprocess
import shutil
import sys
import contextlib
import os
import tempfile as tmp


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


