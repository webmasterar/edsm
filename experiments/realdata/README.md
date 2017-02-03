## Obtaining realdata experiment files and instructions

The realdata experiment files are from the 1000 Genomes Project. You will need to download three sets of files to run the experiment.

The reference sequences were downloaded from: ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/

For example, we download the sequence file for chromosome Y:

`~$ wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.Y.fa.gz`

The reference sequence files must be extracted to work with EDSM:

`~$ gzip -d Homo_sapiens.GRCh37.75.dna.chromosome.Y.fa.gz`

The variant call files (vcf.gz) and their associated tbi files were downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

Download the vcf and tbi files for chromosome Y:

`~$ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz`

`~$ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi`

Then you will be able to run EDSM for the random pattern of length 8:

`~$ ./../../edsm Homo_sapiens.GRCh37.75.dna.chromosome.Y.fa ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz ./randomPatterns/8.txt`

### Problem solving

Sometimes you will get an error like:

`[tabix++] the index file is older than the vcf file. Please use '-f' to overwrite or reindex.`

Be sure to have downloaded the tbi file after the vcf file in future. Try using `touch` to update the tbi file's timestamp:

`~$ touch ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi`

Now try running EDSM again. If it still does not work then you need to recreate the tbi file using Tabix. You can usually install this application through your regular package manager, e.g.:

`~$ sudo apt-get install tabix`

Then you can delete the old tbi file and create the new one by entering:

`~$ rm ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi`

`~$ tabix -p vcf ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz`
