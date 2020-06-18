#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:21:13 2020

@authors: Christophe Guyeux, Guislaine Refrégier, Christophe Sola
"""

import argparse
import Bio.pairwise2
import Bio.SeqIO
import copy
import distutils.spawn
import logging
import logging.config
import os.path
import pathlib
import shutil
import subprocess as sp
import sys

from collections import deque
from networkx import DiGraph
import pydot
from networkx.drawing.nx_pydot import write_dot


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
            self._suffices = [k.lstrip(self._sra).rstrip('.fasta') for k in files]
            for fic in self._suffices:
                sp.run(["sed",
                        "-i", f"s/{self._sra}./{self._sra}{fic}./g",
                        (self._dir / (self._sra + fic)).with_suffix('.fasta')
                        ])
            self._logger.warning('Concatenating SRA files')   
            self._sra_files = [self._dir / k for k in files]
            with open(self._sra_shuffled, "w+") as f:
                sp.call(["cat", *self._sra_files], stdout=f)
        else:
            self._logger.info(f'We will only use one fasta file (OS={sys.platform})')
            files = [k for k in os.listdir(self._dir) 
                     if k.startswith(self._sra) and k.endswith('.fasta')]
            files = [self._dir / k for k in files]
            shutil.copyfile(files[0], self._sra_shuffled)
            self._sra_files = [files[0]]
            
            

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
        seqs = {}
        with open(self._dir / (self._sra + '_crispr_patterns.blast'), 'r') as f:
            for u in f.read().split('\n')[:-1]:
                nom, deb,fin = u.split(',')
                deb,fin = eval(deb),eval(fin)
                seqs[nom] = deb<fin
        fasta_sequences = {}
        self._SEQS = []
        for u in self._sra_files:
            fasta_sequences[u.stem] = Bio.SeqIO.parse(open(u),'fasta')
        for sens in fasta_sequences:
            for fasta in fasta_sequences[sens]:
                if fasta.id in seqs:
                    if seqs[fasta.id]:
                        self._SEQS.append(str(fasta.seq))
                    else:
                        self._SEQS.append(str(fasta.seq.reverse_complement()))


    def __get_chaine(self, G, node):
        chaine = node.split('*')[0]
        successor = list(G.neighbors(node))
        suite = ''
        score = 0
        while len(successor) == 1:
            score += G[node][successor[0]]['weight']
            chaine += '*'+successor[0].split('*')[0]
            suite = '*'.join(successor[0].split('*')[1:])
            node = successor[0]
            successor = list(G.neighbors(successor[0]))
        if suite != '':
            chaine += '*'+suite   
        return (chaine, successor, score)


    def __formate_contigs(self, chaine):
        chaine = chaine.replace('AAAACCCCGAGAGGGGACGGAAAC', '*DRb2')
        chaine = chaine.replace('motif_fin1*motif_fin2*motif_fin3', 'motif_fin')
        chaine = chaine.replace('**', '*')
        cpt, txt = 0, ''
        for k in chaine.split('*'):
            txt += k+'*'
            cpt += len(k)+1
            if cpt > 75:
                cpt = 0
                txt += '\n'
        return txt


    def __get_contigs(self, triplets, sensibility = 200):
        arcs = {}
        for k in triplets:
            for l in triplets:
                if k[1:] == l[:taille_tuple-1]:
                    cle = '*'.join(k)+'+'+'*'.join(l)
                    if cle not in arcs:
                        arcs[cle] = 1
                    else:
                        arcs[cle] += 1
        G = DiGraph()
        for arc in arcs:
            k,l = arc.split('+')
            if arcs[arc] >=sensibility:
                G.add_edge(k,l, weight = arcs[arc])
        write_dot(G, item+'.dot')
        
        investigated = []
        roots = deque([k for k in G.nodes() if list(G.predecessors(k)) == []])
        contigs = []
        while len(roots) != 0:
            node = roots.popleft()
            investigated.append(node)
            chaine, suite, score = self.__get_chaine(G, node)
            if suite != []:
                roots.extend([u for u in suite if u not in investigated])
                for k in suite:
                    if chaine+k not in contigs:
                        contigs.append((chaine+'*'+k, score))
            else:
                if chaine not in [u[0] for u in contigs]:
                    contigs.append((chaine, score))
    
        if len(contigs) == 2:
            if contigs[0][0].endswith('IS6110') and contigs[1][0].startswith('finIS6110'):
                txt = contigs[0][0]+contigs[1][0]
                contigs = [(txt.replace('IS6110finIS6110', 'IS6110'), contigs[0][1]+contigs[1][1])]
            elif contigs[1][0].endswith('IS6110') and contigs[0][0].startswith('finIS6110'):
                txt = contigs[1][0]+contigs[0][0]
                contigs = [(txt.replace('IS6110finIS6110', 'IS6110'), contigs[0][1]+contigs[1][1])]
                
        return contigs

    def _parse_matches(self):
        dicofind, dico_cas = {}, {}
        with open(pathlib.Path('data') / 'fastas' / 'crispr_patterns.fasta') as f:
            txt = f.read()
        for k in txt.split('>')[1:]:
            dicofind[k.split('\n')[0]] = k.split('\n')[1]
        with open(pathlib.Path('data') / 'fastas' / 'cas_patterns.fasta') as f:
            txt = f.read()
        for k in txt.split('>')[1:]:
            dico_cas[k.split('\n')[0]] = k.split('\n')[1]
        S=[]
        for k in self._SEQS:
            s = copy.deepcopy(k)
            for l in dicofind:
                s = s.replace(dicofind[l],'*'+l+'*')
            for l in dico_cas:
                s = s.replace(dico_cas[l],'*'+l+'*')
            s = s.replace('**','*')
            s = s.replace('*GTCGTCAGACCCAAAACCC*','*rDRa1*')
            s = s.replace('*CCCCGAGAGGGGACGGAAAC*', '*DRb1*')
            s = s.replace('*AAAACCCCGAGAGGGGACGGAAAC*', '*DRb2*')
            if '*' in s:
                S.append(s.split('*')[1:-1])
        
        taille_tuple = 3
        liste = [u for u in [[k[l:l+taille_tuple] for l in range(len(k)-(taille_tuple-1))] for k in S if len(k) >= taille_tuple-1] if u != []]
        triplets = [item for sublist in liste for item in sublist]
        
        occurrences_minimales = 3
        sommets = [item for sublist in triplets for item in sublist]
        sommets = sorted(set([k for k in sommets if sommets.count(k) >= occurrences_minimales]))

        triplets = [k for k in triplets if all([u in sommets for u in k])]




if __name__ == "__main__":
    CRISPRbuilderTB()
