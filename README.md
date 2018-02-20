## EDSM

Elastic Degenerate String Matching.

If you make use of this software, please cite the following:

> R. Grossi, C. S. Iliopoulos, C. Liu, N. Pisanti, S. P. Pissis, A. Retha, G. Rosone, F. Vayani, L. Versari, "On-Line Pattern Matching on Similar Texts", in 28th Annual Symposium on Combinatorial Pattern Matching (CPM 2017), J. R. Juha Kärkkäinen, W. Rytter, Eds., Dagstuhl, Germany: Schloss Dagstuhl--Leibniz-Zentrum fuer Informatik, pp. 9:1-9:14.

### How to Install

You may need to install one or more libraries before compiling EDSM. Please read INSTALL.md.

### How to Execute

`~$ ./edsm seq.txt pattern`

or

`~$ ./edsm reference.fasta variants.vcf pattern`

The `pattern` may be a string or the path to a file containing a single pattern.

If you want to use a compressed vcf file (*.vcf.gz), please make sure its accompanying tbi file is also present in the same directory. You can also use `Tabix` to generate a tbi file.

### License

GNU GPLv3 License; Copyright (C) 2017 Chang Liu, Solon P. Pissis, Ahmad Retha and Fatima Vayani.

