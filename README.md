## EDSM

Elastic Degenerate String Matching.

If you make use of this software, please cite the following:

> Robert Grossi, Costas S. Iliopoulos, Chang Liu, Nadia Pisanti, Solon P. Pissis, Ahmad Retha, Giovanna Rosone, Fatima Vayani and Luca Versari; On-line pattern matching on similar texts; (in preparation).

### How to Install

`~$ ./pre-install.sh && make`

### How to Execute

`~$ ./edsm seq.txt pattern`

or

`~$ ./edsm reference.fasta variants.vcf pattern`

The `pattern` may be a string or the path to a file containing a single pattern.

If you want to use a compressed vcf file (*.vcf.gz), please make sure its accompanying tbi file is also present in the same directory. You can also use `Tabix` to generate a tbi file.

### License

GNU GPLv3 License; Copyright (C) 2017 Chang Liu, Solon P. Pissis, Ahmad Retha and Fatima Vayani.

