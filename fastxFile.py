#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 08:24:01 2017

@author: Koen Herten

This is pbunmap v1. A tool to compare the output of multiple read mappers

Copyright 2017, Koen Herten, All rights reserved

This file is part of aftermerge.

pbunmap is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pbunmap is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pbunmap.  If not, see <http://www.gnu.org/licenses/>.

 
"""

import fastx

class FastxFile:
    
    def __init__(self, file_name, isfasta = False):
        '''
        initiate all variables
        '''
        self._isfasta = isfasta
        if isfasta:
            self._file = FastaFile(file_name)
        else:
            self._file = FastqFile(file_name)
        
    @property
    def isfasta(self):
        return self._isfasta
        
    def nextFastx(self):
        return self._file.nextFasta()
    
    def close(self):
        self._file.close()

class FastqFile:
    
    def __init__(self, file_name):
        '''
        initiate all variables
        '''
        self._file = open(file_name)
        
        
    def nextFastq(self):
        line1 = self._file.readline()
        if (line1 is None or line1 is ""):
            return None
        line1 = line1.rstrip()
        line2 = self._file.readline()
        line2 = line2.rstrip()
        line3 = self._file.readline()
        line3 = line3.rstrip()
        line4 = self._file.readline()
        line4 = line4.rstrip()
        return fastx.FastqRead(line1, line2, line3, line4)
        
    def nextFasta(self):
        fastq = self.nextFastq()
        if fastq is None:
            return None
        return fastq.toFastaRead()
        
    def close(self):
        self._file.close()
        
        
class FastaFile:
    
    def __init__(self, file_name):
        '''
        initiate all variables
        '''
        self._file = open(file_name)
        self._prevline = None
        
        
    def nextFasta(self):
        if self._prevline is None:
            self._prevline = self._file.readline()
            if (self._prevline is None or line is ""):
                return None
            self._prevline = self._prevline.rstrip()
        line = self._file.readline()
        seq = ""
        while not str(line).startwith(">") or line is None or line is "":
            line = line.rstrip()
            seq = "{}{}".format(seq, line)
            line = self._file.readline()
        fasta = fastx.FastaRead(self._prevline, seq)
        self._prevline = line
        return fasta
        
        
    def close(self):
        self._file.close()