# circos
Collection of random scritps, for circos plot of a genome.
A good tutorial to follow:
https://dbsloan.github.io/TS2019/exercises/circos.html

1. After the denovo assembly. sort the contigs by length using BBMap package:
```
#!/bin/bash

#SBATCH -p test # Partition to submit to (comma separated)
#SBATCH -J bbmap
#SBATCH -c 6 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=10gb
#SBATCH --open-mode=append
#SBATCH -t 0:20:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module load GCC/7.3.0-2.30 OpenMPI/3.1.1 BBMap/38.35
sortbyname.sh in=TdomV1.asm.bp.hap1.p_ctg.fasta out=sorted.length.TdomV1.asm.bp.hap1.p_ctg.fa length descending
```
2. Install seqkit tool through bioconda to prepare tab delimited file for the scaffold/chromosome coordingates in the genome asembly
```
conda install -c bioconda seqkit
```
then
```
seqkit fx2tab --length --name --header-line --gc --gc-skew --base-count AT --base-count GC TdomV1.asm.bp.hap1.p_ctg.fasta >seqkit.fx2tab-length-name-header-line-gc-AT-GC.out
```
3. Prepare karyotype file for the assembly, column should be like:`chr - ID LABEL START END COLOR`

```
awk '{print "chr", "-", $1, $1, "0", $2}' sorted.seqkit.fx2tab-length-name-header-line-gc-AT-GC.out > karyotype.thermobia.txt
```
add color column at the seventh column
```
awk '{print $0, (NR%2==0)?"blue":"red"}' karyotype.thermobia.txt > karyotype.thermobia.color.txt
```
by default circos only allows 200 rows at max for karyotype files so extracting just 200 scaffolds to try at first.
```
head -n 36 karyotype.thermobia.color.txt > karyotype.thermobia.color.first36.txt
```
may be replace the levels with 1-2-3 continuous numbers
```
awk '{$4=NR}1' karyotype.thermobia.color.first36.txt >karyotype.thermobia.color.first36.levelsedited.txt
```




After preparing circos.conf file start an interactive session in canon cluster and load circos module
```
salloc -p test -t 60 --mem 4000
```
```
module load GCCcore/7.3.0 Circos/0.69-6-Perl-5.28.0
```


