#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 07:59:15 2017

This is pbunmap v1. A tool to unmap PacBio reads who were to short

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

from argparse import ArgumentParser
from argparse import FileType
import sys
import samRead
import samFile
import fastx
import fastxFile 
import os.path



def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if __name__ == '__main__':
    #add all arguments
    parser = ArgumentParser(description='pbunmap v1.0 unmapping wrongly (short) mapped PacBio reads')
    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('-fasta', '--isfasta', action="store_true", dest="isfasta", help="The read file is a fasta file (else fastq file)")
    parser.add_argument('readfile', help="The fastq or fasta file with the raw data")
    parser.add_argument('-ccs', '--ccs', action="store_true", help="These are PacBio CCS reads (fastq description should end with ccs (some mappers add extra numbers, these are removed to get the original query name))")
    parser.add_argument('-l', '--length', help="The minimum length of a mapped fragment", default=0, type=int)
    parser.add_argument('-p', '--percent', help="The percentage of the raw read that has to be mapped to the reference (between 0-1)", default=0, type=float)
    parser.add_argument('samfile', help="The mapped data of the readfile")
    
    
    
    args = parser.parse_args()
    #if verbose, print the parameters
    if args.verbose:
        eprint("Verbose: {}".format(args.verbose))
        eprint("Read File is fasta: {}".format(args.isfasta))
        eprint("Reads are CCS reads: {}".format(args.ccs))
        eprint("Minimum length to be mapped: {}".format(args.length))
        eprint("Minimum percentage to be mapped: {}".format(args.percent))
        eprint("Read File: {}".format(args.readfile))
        eprint("Sam file: {}".format(args.samfile))
        
    #check values
    max_length = 999999
    if args.length < 0 or args.length > max_length:
        eprint("The minimum length must be between 0 and {}".format(max_length))
        exit(1)
    if (args.percent < 0 or args.percent > 1):
        eprint("The percentage must be between 0 and 1")
        exit(1)
    
    #check files
    if not os.path.exists(args.readfile):
        eprint("Read File {} does not exist".format(args.readfile))
        exit(1)
    if not os.path.exists(args.samfile):
        eprint("Sam file: {} does not exist".format(args.samfile))
        exit(1)
        
        
    #parsing the fastq/fasta file to a length directory
    readfile = fastxFile.FastxFile(args.readfile)

    read_length_dict = {}
    
    fasta = readfile.nextFasta()
    fasta_count = 0
    while (fasta is not None):
        fasta_count = fasta_count + 1
        if (args.verbose and fasta_count%100==0):
            eprint("reads processed: {}".format(fasta_count))
        read_length_dict[fasta.name] = fasta.sequenceLength
        fasta = readfile.nextFasta()
        
    readfile.close()
    
    if (args.verbose):
        eprint("Processed {} reads in the read file".format(fasta_count))
        
    #go over sam file, and correct
    samfile = samFile.SamFile(args.samfile, args.ccs)
    
    samread = samfile.nextSamRead()
    
    sam_count = 0
    previous_samread = None
    second_count = 0
    suplement_count = 0
    toshort = 0
    while (samread is not None):
        sam_count = sam_count + 1
        if (args.verbose and sam_count%100==0):
            eprint("mapped processed: {}".format(sam_count))
        #print(samread.line)
        raw_length = read_length_dict[samread.qname]
        if (raw_length is not None or raw_length is not ""):
            #fragment found
            mapped_length = samread.getMappedLength()
            #check length
            if (samread.isSecondaryAlignment()):
                #is a secondary alignment
                #print("SECONDARY")
                second_count = second_count +1
                samread.changeToUnmapped()
            if (samread.issuplementary()):
                #is mapped to multiple places
                suplement_count = suplement_count + 1
            if ((mapped_length < args.length) or mapped_length < (raw_length * args.percent)):
                #print("TO SHORT")
                toshort = toshort + 1
                samread.changeToUnmapped()
        #print(samread.line)
        previous_samread = samread
        samread = samfile.nextSamRead()
    
    samfile.close()
    
    if (args.verbose):
        eprint("Processed {} reads in the mapped file".format(sam_count))
        eprint("Found {} secondary mapped reads".format(second_count))
        eprint("Found {} suplementary mapped reads".format(suplement_count))
        eprint("Found {} mapped reads that were to short".format(toshort))
    