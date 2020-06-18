#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 07:27:17 2020

@author: christophe
"""

import distutils.spawn

class Collector:
    def __init__(self, sras = [], outdir = "sequences"):
        self._sras = sras
        self._outdir = outdir


    def _check_for_tools(self):
        name = 'fastq-dump'
        if distutils.spawn.find_executable(f"{name}.exe") != None:
            self._logger.warning(f'{name} found in system directory')
            if name == 'fastq-dump':
                self._fastq_dump = 'fastq-dump.exe'
            elif name == 'blastn':
                self._blastn = 'blastn.exe'
            else:
                self._makeblastdb = 'makeblastdb.exe'
        elif distutils.spawn.find_executable(name) != None:
            self._logger.warning(f'{name} found in system directory')
            if name == 'fastq-dump':
                self._fastq_dump = 'fastq-dump'
            elif name == 'blastn':
                self._blastn = 'blastn'
            else:
                self._makeblastdb = 'makeblastdb'
        else:
            if sys.platform == 'win32':
                if name == 'fastq-dump':
                    self._fastq_dump = 'bin\windows\fastq-dump.exe'
                elif name == 'blastn':
                    self._blastn = 'bin\windows\blastn.exe'
                else:
                    self._makeblastdb = 'bin\windows\makeblastdb.exe'
            elif sys.platform == 'linux':
                if name == 'fastq-dump':
                    self._fastq_dump = 'bin/linux/./fastq-dump'
                elif name == 'blastn':
                    self._blastn = 'bin/linux/./blastn'
                else:
                    self._makeblastdb = 'bin/linux/./makeblastdb'
            elif sys.platform in ['darwin']:
                if name == 'fastq-dump':
                    self._fastq_dump = 'bin/mac/./fastq-dump'
                elif name == 'blastn':
                    self._blastn = 'bin/mac/./blastn'
                else:
                    self._makeblastdb = 'bin/mac/./makeblastdb'


if __name__ == '__main__':
    collect = Collector(sras = ['ERR037527'])
    collect.get()