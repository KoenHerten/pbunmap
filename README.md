# pbunmap 

This tool unmap reads that are mapped, but where the length of the mapped subsequence is lower than a certain length, or is lower than a certain percentage of the complete sequence. Ideal if the read is mapping uniquely to a location, but for instance only for 200bp, while the actual sequence is 2000bp. In this case, the read is actually mapped incorrectly, and should be unmapped. However it is still possible that the mapping quality is high, and will not be filtered out in any way.

##Licence

All parts of this tool is licenced under GPLv3.  
A copy of this licence is included under LICENSE.

##Help
samcompare is really simple to use:
The only thing needed is a list of name sorted sam files.

Example of usage:
```bash

python pbunmap.py sample.fastq sample.sam

```

