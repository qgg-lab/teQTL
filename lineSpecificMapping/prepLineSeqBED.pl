#! /usr/bin/perl -w
# process VCF file to obtain a bed file
# for a given DGRP line for subsequent lift
# input is VCF and a line number
# output a bed file with sites to modify
# ============================================================

use Getopt::Long;
use strict;
use warnings;
use List::Util qw(first sum);

my $vcf;
my $dgrp;
my @samples;
my @line;

# option
# ============================================================

GetOptions("vcf=s" => \$vcf,
           "dgrp=s" => \$dgrp);

# ============================================================

open VCF, "<$vcf";

while ( <VCF> ) {
  
  chomp $_;
  if ($_ =~ m/^\#\#/) {
    next;
  } elsif ($_ =~ m/^\#CHROM/) {
    @line = split /\t/, $_;
    @samples = @line[9..$#line];
    last;
  }
  
}

my $idx = first { $samples[$_] eq $dgrp } 0..$#samples;
if (!defined($idx)) {
	print STDERR "Cannot find $dgrp in the VCF file.\n";
	exit;
}

# initiate variables
# ============================================================

my $pre_chr = 0;
my $pre_start = 0;
my $pre_end = 0;
my $pre_id = "";
my $pre_ref = "";
my $pre_alt = "";
my @thisLineGeno = ();
my @alleles = ();
my @thisLineAlleles = ();

# chr, start, end, ref, alt, id or any other info
# if the line carries homo ref, nothing to be done
# if the line carries homo alt, output the alt allele
# if the line carries het, output the first allele
#    in this case, if the one of the alleles is ref
#    it will the first allele, otherwise, output one of the two
# ============================================================

# find the first site that contains both non-ref allele

while ( <VCF> ) {
  
  chomp $_;
  @line = split /\t/, $_;
  @thisLineGeno = split /:/, $line[$idx + 9];
  @thisLineAlleles = split /\//, $thisLineGeno[0];
  @alleles = split /,/, $line[3].",".$line[4];
  
  if ($thisLineAlleles[0] eq "0" || $thisLineAlleles[1] eq "0" || 
  	  $alleles[$thisLineAlleles[0]] eq "*" || $alleles[$thisLineAlleles[1]] eq "*") {
	next;
  } else {
  	$pre_start = $line[1] - 1;
  	$pre_end = $line[1] - 1 + length($line[3]);
  	$pre_chr = $line[0];
		$pre_id = $line[2];
		$pre_ref = $line[3];
		$pre_alt = $alleles[$thisLineAlleles[0]];
		last;
  }

}

# start output

while ( <VCF> ) {
  
  chomp $_;
  @line = split /\t/, $_;
  @thisLineGeno = split /:/, $line[$idx + 9];
  @thisLineAlleles = split /\//, $thisLineGeno[0];
  @alleles = split /,/, $line[3].",".$line[4];
  
  if (($line[0] ne $pre_chr || $line[1] - 1 >= $pre_end)) {
		
	  if ($thisLineAlleles[0] ne "0" && $thisLineAlleles[1] ne "0" && 
	  	  $alleles[$thisLineAlleles[0]] ne "*" && $alleles[$thisLineAlleles[1]] ne "*") {
		
			# output previous line
			print $pre_chr, "\t", $pre_start, "\t", $pre_end, "\t", $pre_ref, "\t", 
						$pre_alt, "\t", $pre_id, "\n";
		
    	$pre_start = $line[1] - 1;
    	$pre_end = $line[1] - 1 + length($line[3]);
    	$pre_chr = $line[0];
  		$pre_id = $line[2];
  		$pre_ref = $line[3];
  		$pre_alt = $alleles[$thisLineAlleles[0]];
		
		}
  
	}
	
}
  
# last item in the buffer
# ============================================================

print $pre_chr, "\t", $pre_start, "\t", $pre_end, "\t", $pre_ref, "\t", 
			$pre_alt, "\t", $pre_id, "\n";

close VCF;
