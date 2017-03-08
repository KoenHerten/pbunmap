#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 08:17:07 2017

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


import sys

class FastqRead:
    
    def __init__(self, queryname, sequence, plus, quality):
        '''
        initiate all variables
        '''
        self._queryname = queryname
        self._sequence = sequence
        self._plus = plus
        self._quality = quality
        if (self._queryname is None or self._sequence is None or self._plus is None or self._quality is None) or (self._queryname is "" or self._sequence is "" or self._plus is "" or self._quality is ""):
            print("NO valid fastq block found", file=sys.stderr)
    
    @property
    def queryname(self):
        return str(self._queryname)
        
    @property
    def name(self):
        return str(self._queryname).lstrip("@")
        
    @property
    def sequence(self):
        return str(self._sequence)
        
    @property
    def quality(self):
        return str(self._sequence)
        
    @property
    def sequenceLength(self):
        return len(str(self._sequence))
        
    def toFastaRead(self):
        return FastaRead(">{}".format(self.name), self.sequence)
        
        
class FastaRead:
    
    def __init__(self, queryname, sequence):
        '''
        initiate all variables
        '''
        self._queryname = queryname
        self._sequence = sequence
        if (self._queryname is None or self._sequence is None) or (self._queryname is "" or self._sequence is ""):
            print("NO valid fasta block found", file=sys.stderr)
        
    
    @property
    def queryname(self):
        return str(self._queryname)
        
    @property
    def name(self):
        return str(self._queryname).lstrip(">")
        
    @property
    def sequence(self):
        return str(self._sequence)
        
    @property
    def sequenceLength(self):
        return len(str(self._sequence))