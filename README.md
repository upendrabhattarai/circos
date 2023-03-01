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
by default circos only allows 200 rows at max for karyotype files, but I am extracting only longest 36 scaffolds to try at first.
```
head -n 36 karyotype.thermobia.color.txt > karyotype.thermobia.color.first36.txt
```
may be replace the levels with 1-2-3 continuous numbers
```
awk '{$4=NR}1' karyotype.thermobia.color.first36.txt >karyotype.thermobia.color.first36.levelsedited.txt
```
extract the longest 36 sequences to calculate gc skewness
```
awk '/^>/ && ++c>36 {exit} {print}'sorted.length.TdomV1.asm.bp.hap1.p_ctg.fa > sorted.first36.sequences.fa
```
Now we will calculate the gc skewness and level with color using a script
`gc_skew_for_circos.pl`

```
#!/usr/bin/perl 

#A program to analyze GC skew in a (circular) genome with a given window and step size.
#I have modified this to produce output for circos

use warnings;
use strict;
use POSIX;

my $usage = "\nUSAGE: $0 fastaFile windowSize stepSize chromosomeName[forCircos] posColor[forCircos] negColor[forCircos]\n\n";

my $file = shift or die ($usage);
my $windowSize = shift or die ($usage);
my $stepSize = shift or die ($usage);
my $chromosomeName = shift or die ($usage);
my $posColor = shift or die ($usage);
my $negColor = shift or die ($usage);
my $leftSize = ceil(($windowSize - 1)/2);
my $rightSize = floor (($windowSize - 1)/2);


my ($headerRef, $seqRef) = get_fasta_names_and_seqs($file);
my @seqArray = @$seqRef;

my $seq = $seqArray[0];
my $length = length ($seq);

for (my $i = 0; $i < $length; $i += $stepSize){
        my $gcount = 0;
        my $ccount = 0;
        
        my $left = $i - $leftSize;
        my $right = $i + $rightSize;
        
        
        if ($left < 0){
                $left = $length + $left;
        } 
        if ($right >= $length){
                $right = $right - $length;
        }
        

        if ($left < $right){
                my $substr = uc(substr ($seq, $left, $right - $left + 1));
                $ccount += $substr =~ tr/C/C/;
                $gcount += $substr =~ tr/G/G/;          
        }else{
                my $substr1 = uc(substr ($seq, 0, $right + 1));
                my $substr2 = uc(substr ($seq, $left));
                $ccount += $substr1 =~ tr/C/C/;
                $gcount += $substr1 =~ tr/G/G/;         
                $ccount += $substr2 =~ tr/C/C/;
                $gcount += $substr2 =~ tr/G/G/;         
        }       
        print "$chromosomeName\t", $i+1, "\t", $i+1, "\t", ($gcount - $ccount)/($gcount + $ccount), "\tfill_color=";
        if (($gcount - $ccount)/($gcount + $ccount) > 0){
                print $posColor, "\n";
        }else{
                print $negColor, "\n";
        }

}

exit;



sub get_fasta_names_and_seqs {
        use strict;
        use warnings;

        my ($inputfilename) = @_;
        my @fasta_names = ();
        my @fasta_seqs= ();

                   
        unless ( open(FILEDATA, $inputfilename) ) {
                print STDERR "Cannot open file \"$inputfilename\"\n\n"; #print error message
                exit; #exit the program
        }       

        my @filedata = <FILEDATA>; #Read the lines of the file into an array
        close FILEDATA;
        
        my $seq_count = 0; #this will be used to keep track of the number of sequences
        foreach my $line (@filedata){
                chomp $line;
                if ($line =~ /^\s*$/) {next;} #ignore line if it is blank
                
                elsif ($line =~ /^>/) { #if the line is a header line (begins with ">")...
                        if ($line =~ /^>.*[\w]+/){
                                my $partialLine = substr ($&, 1);
                                push (@fasta_names, $partialLine); #add that line to an array of fasta names
                                push (@fasta_seqs, ""); #and add a new blank element to an array of sequences
                                ++$seq_count; #also increment our counter which keeps track of sequence number
                        }
                }       
                
                else { #if the line's not blank or a header, add it to the current sequence 
                        $fasta_seqs[$seq_count-1] .= $line;
                }
                
                $fasta_seqs[$seq_count-1] =~s/\s//g; #remove all whitespace from the current  sequence
        }
        
        return (\@fasta_names, \@fasta_seqs);

}
```



After preparing circos.conf file start an interactive session in canon cluster and load circos module
```
salloc -p test -t 60 --mem 4000
```
```
module load GCCcore/7.3.0 Circos/0.69-6-Perl-5.28.0
```


