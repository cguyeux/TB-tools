#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:21:13 2020

@authors: Christophe Guyeux, Guislaine Refrégier, Christophe Sola
"""

import argparse
import distutils.spawn
import logging
import logging.config
import os.path
import pathlib
import shutil
import subprocess as sp
import sys


class CRISPRbuilderTB:
    def __init__(self):
        # Parse arguments
        self._parse_args()
        self._linux_or_mac = sys.platform.startswith('linux') or sys.platform.startswith('darwin')
        if not self._linux_or_mac:
            self._logger.warning(f'{sys.platform} OS not recommended (sharp reduction in speed')
        self._run()


    def _parse_args(self):
        """
        Parses the arguments provided to the CRISPRbuilder-TB
        """
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-sra",
                            help="accession number to deal with",
                            default='',
                            type=str)
        parser.add_argument('-d',
                            '--debug',
                            help="debugging mode of logs",
                            action="store_const",
                            dest="loglevel",
                            const=logging.DEBUG,
                            default=logging.WARNING
                            )
        parser.add_argument('-v',
                            '--verbose',
                            help="verbose logs",
                            action="store_const",
                            dest="loglevel",
                            const=logging.INFO
                            )
        parser.add_argument("-out",
                            "--output_directory",
                            default='sequences',
                            help="directory of outputs",
                            type=str)
        parser.add_argument("--num_threads",
                            default='1',
                            help="number of threads",
                            type=str)
        parser.add_argument("-evalue",
                            default='1e-7',
                            help="evalue when blasting spacers, DRs, etc.",
                            type=str)
        args = parser.parse_args()
        logging.basicConfig(level=args.loglevel)
        self._logger = logging.getLogger()
        self._sra = args.sra
        self._output_dir = pathlib.Path(args.output_directory)
        self._num_threads = args.num_threads
        self._evalue = args.evalue


    def _run(self):
        self._check_for_tools()
        self._prepare_directory()
        self._collect_and_prepare_sra()
        self._make_blast_db()
        self._sequences_of_interest()



    def _prepare_directory(self):
        """
        Prepare the structure of sequence directory        
        """
        self._logger.debug('Preparing directory structure')
        self._dir = self._output_dir / self._sra
        self._dir.mkdir(exist_ok=True, parents=True)


    def _check_for_tools(self):
        for name in ['fastq-dump', 'blastn', 'makeblastdb']:
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


    def _collect_and_prepare_sra(self):
        self._sra_shuffled = self._dir / (self._sra + '_shuffled')
        self._sra_shuffled = self._sra_shuffled.with_suffix('.fasta')
        if not os.path.isfile(self._sra_shuffled):
            sp.run([self._fastq_dump,
                    '--split-files',
                    '--fasta',
                    '-O', self._dir,
                    self._sra
                    ])
        if self._linux_or_mac:
            self._logger.debug('Renaming reads')  
            files = [k for k in os.listdir(self._dir) 
                     if k.startswith(self._sra) and k.endswith('.fasta')]
            for fic in [k.lstrip(self._sra).rstrip('.fasta') for k in files]:
                sp.run(["sed",
                        "-i", f"s/{self._sra}./{self._sra}{fic}./g",
                        (self._dir / (self._sra + fic)).with_suffix('.fasta')
                        ])
            self._logger.warning('Concatenating SRA files')   
            files = [self._dir / k for k in files]
            with open(self._sra_shuffled, "w+") as f:
                sp.call(["cat", *files], stdout=f)
        else:
            self._logger.info(f'We will only use one fasta file (OS={sys.platform})')
            files = [k for k in os.listdir(self._dir) 
                     if k.startswith(self._sra) and k.endswith('.fasta')]
            files = [self._dir / k for k in files]
            shutil.copyfile(files[0], self._sra_shuffled)
            

    def _make_blast_db(self):
        if not any([p.suffix == '.nin' for p in pathlib.Path(self._dir).iterdir()]):            
            completed = sp.run([self._makeblastdb,
                                '-in', self._sra_shuffled,
                                '-dbtype', 'nucl',
                                '-title', self._sra,
                                '-out', self._dir / self._sra])
            assert completed.returncode == 0


    def _sequences_of_interest(self):
        completed = sp.run([self._blastn,
                            '-num_threads', self._num_threads,
                            '-query', pathlib.Path('data') / 'fastas' / 'crispr_patterns.fasta',
                            '-evalue', self._evalue,
                            '-task', "blastn",
                            '-db', self._dir / self._sra,
                            '-outfmt',
                            '10 sseqid sstart send',
                            '-out', self._dir / (self._sra + '_crispr_patterns.blast')])
        assert completed.returncode == 0



if __name__ == "__main__":
    CRISPRbuilderTB()
