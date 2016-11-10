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
        return False


def check_programs(*arg):
    error_list = []
    for program in arg:
        if which(program) is False:
            error_list.append("\t[!] {} not found! Please install and add to it PATH variable".format(program))
    if error_list:
        print("\n".join(error_list))


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
