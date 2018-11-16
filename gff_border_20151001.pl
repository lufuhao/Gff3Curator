#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage:  in.gff.gz out.gff.gz

Extract GFF border from GFF to let IGV display

v20150710

EOH
die USAGE if (scalar(@ARGV) !=2 or [0] eq '-h' or [0] eq '--help');


my $inputfile=$ARGV[0];
die "Error: invalid gff file\n" unless (defined $inputfile and -s $inputfile);
my $outputfile=$ARGV[1];
unlink $outputfile if (defined $outputfile and -s $outputfile);
if ($inputfile=~/\.gff\d*$/i) {
	open (GFFIN, $inputfile) || die "Error: can not open input GFF file: $inputfile\n";
}
elsif ($inputfile=~/\.gff\d*\.gz$/i) {
	open (GFFIN, "zcat $inputfile | ") || die "Error: can not open input GFF file: $inputfile\n";
}
else {
	die "Error: unknown to open GFF file\n";
}
open (GFFOUT, " | bgzip > $outputfile") || die "Error: can not output GFF file: $outputfile\n";
while (my $line=<GFFIN>) {
	chomp $line;
	next if ($line=~/^#/);
	my @arr=split(/\t/, $line);
	next unless (scalar(@arr)>=8);
	print GFFOUT $arr[0], "\t", $arr[1], "\tBORDER\t", $arr[3], "\t", $arr[3], "\t", $arr[5], "\t", $arr[6], "\t", $arr[7], "\t", "ID=$arr[3]\n";
	print GFFOUT $arr[0], "\t", $arr[1], "\tBORDER\t", $arr[4], "\t", $arr[4], "\t", $arr[5], "\t", $arr[6], "\t", $arr[7], "\t", "ID=$arr[4]\n";
}
close GFFIN;
close GFFOUT;
