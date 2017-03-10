# pbunmap 

This tool unmap reads that are mapped, but where the length of the mapped subsequence is lower than a certain length, or is lower than a certain percentage of the complete sequence. Ideal if the read is mapping uniquely to a location, but for instance only for 200bp, while the actual sequence is 2000bp. In this case, the read is actually mapped incorrectly, and should be unmapped. However it is still possible that the mapping quality is high, and will not be filtered out in any way.

##Licence

All parts of this tool is licenced under GPLv3.  
A copy of this licence is included under LICENSE.

##Dependencies

pbunmap was developed and tested usint Python 3.5.2

##Help
samcompare is really simple to use:
The only thing needed is a list of name sorted sam files. The output of the new sam file is on the standard output, metrics are printed to the standard error output.

Example of usage:
```bash

python pbunmap.py sample.fastq sample.sam

```

##Parameters

For the latest information use:
```bash
python pbunmap.py -h
```

positional arguments:
| readfile | The fastq or fasta file with the raw data |
| samfile | The mapped data of the readfile |

optional arguments:
| -h, --help | show this help message and exit |
| -v, --verbose | |
| -fasta, --isfasta | The read file is a fasta file (else fastq file) |
| -ccs, --ccs | These are PacBio CCS reads (fastq description should end with ccs (some mappers add extra numbers, these are removed to get the original query name)) |
| -l LENGTH, --length LENGTH | The minimum length of a mapped fragment (default: 200) |
| -p PERCENT, --percent PERCENT | The percentage of the raw read that has to be mapped to the reference (between 0-1) (default: 0.5) |
| -ns, --not-correct-supplementary | If there are supplementary reads passing the filter, do not correct the primary |

