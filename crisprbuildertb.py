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
        args = parser.parse_args()
        logging.basicConfig(level=args.loglevel)
        self._logger = logging.getLogger()
        self._sra = args.sra
        self._output_dir = pathlib.Path(args.output_directory)


    def _run(self):
        self._check_for_tools()
        self._prepare_directory()
        self._collect_sra()
        self._prepare_sra()
        self._make_blast_db()
        self._sequences_of_interest()
        self._clean_directory()


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


    def _collect_sra(self):
        """
        Download the SRA files if needed.
        """
        self._sra_file_1 = self._dir / (self._sra + '_1')
        self._sra_file_1 = self._sra_file_1.with_suffix('.fasta')
        if not os.path.isfile(self._sra_file_1):
            sp.run([self._fastq_dump,
                    '--split-files',
                    '--fasta',
                    '-O', self._dir,
                    self._sra
                    ])
        self._sra_file_2 = self._dir / (self._sra + '_2')
        self._sra_file_2 = self._sra_file_2.with_suffix('.fasta')
        if not os.path.isfile(self._sra_file_1) or not os.path.isfile(self._sra_file_2):
            raise NotImplementedError("Fasta stem suffices are not _1 and _2")
        self._logger.warning('SRA files found in sequence directory')


    def _prepare_sra(self):
        """
        Merge the two downloaded fasta files of the paired Illumina WGS data.
        Valid only for linux or mac. For windows, we only take the first file.
        """
        if self._linux_or_mac:
            self._sra_shuffled = self._dir / (self._sra + '_shuffled')
            self._sra_shuffled = self._sra_shuffled.with_suffix('.fasta')
            if not os.path.isfile(self._sra_shuffled):
                self._logger.debug('Renaming reads')        
                for fic in ['_1', '_2']:
                    sp.run(["sed",
                            "-i", f"s/{self._sra}./{self._sra}{fic}./g",
                            (self._dir / (self._sra + fic)).with_suffix('.fasta')
                            ])
                self._logger.warning('Shuffling SRA files')      
                with open(self._sra_shuffled, "w+") as f:
                    sp.call(["cat",
                             self._sra_file_1,
                             self._sra_file_2],
                            stdout=f
                            )
        else:
            self._logger.info(f'We will only use one fasta file (OS={sys.platform})')
            self._sra_shuffled = self._sra_file_1
            

    def _make_blast_db(self):
        if not any([p.suffix == '.nin' for p in pathlib.Path(self._dir).iterdir()]):            
            completed = sp.run([self._makeblastdb,
                                '-in', self._sra_shuffled,
                                '-dbtype', 'nucl',
                                '-title', self._sra,
                                '-out', self._dir / self._sra])
            assert completed.returncode == 0


    def _sequences_of_interest(self):
        txt = ''
        for k in dicofind:
            txt += '>'+k
            txt += os.linesep
            txt += dicofind[k]
            txt += os.linesep
        with open('sequences.fasta','w') as f:
            f.write(txt)
        completed = sp.run([self._blastn,
                            '-num_threads', self._num_threads,
                            '-query', self._dir_data / 'spoligo_vitro.fasta',
                            '-evalue', self._spol_evalue,
                            '-task', "blastn",
                            '-db', self._dir / self._sra,
                            '-outfmt',
                            '10 qseqid sseqid sstart send qlen length score evalue',
                            '-out', self._dir / (self._sra + '_vitro_old-spol.blast')])


completed = subprocess.run("blastn -query /tmp/interet.fasta -task blastn -evalue 1e-7 -db "+item+" -num_threads 8 -outfmt '10 sseqid sstart send' -out /tmp/"+item, shell = True)
assert completed.returncode == 0


    def _clean_directory(self):
        os.remove('sequences.fasta')


    


if __name__ == "__main__":
    CRISPRbuilderTB()
