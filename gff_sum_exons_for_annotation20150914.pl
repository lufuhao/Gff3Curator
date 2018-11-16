#!/usr/bin/env perl
use strict;
use warnings;
use Storable qw/dclone/;
use Data::Dumper qw /Dumper/;
use FuhaoPerl5Lib::FastaKit qw/IndexFasta Codon2AA SeqRevComp/;
use SVG;
use Bio::DB::Fasta;
use Term::ANSIColor;
use constant USAGE =><<EOH;

usage: $0 fasta region mincount *.gff[.gz]

v20160126

EOH
die USAGE if (scalar(@ARGV) !=1 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');
my $error_color_code="Red on_bright_yellow";


### Input
my @fastas=glob "1.*.fasta";
die "Error: cann not dete4ct fasta file\n" if (scalar(@fastas)!=1);
my $fastafile=shift @fastas;
my $region=shift @ARGV;
my @gfffiles=glob "*.gz";
die "Error: can not find GFF.gz file\n" if (scalar(@gfffiles)==0);
#my $mincount=shift @ARGV;
#unless ($mincount=~/^\d+$/) {
#	die "Error: invalid minimum count: $mincount\n";
#}
my $mincount=3;



### Define global variables
my %genenotes=('0' => 'pseudogene', '1' => 'pseudogene(internal_stop)', '2' => 'pseudogene(no_stop)', '3' => 'pseudogene(no_start)', '4' => 'pseudogene(frameshift)', '5' => 'Part5missing', '6' => 'Part3missing', '7' => 'PartMmissing', '8' => 'Part53missing');
my $seqid='';
my $start;
my $end;
my $strand;
my $startcodon;
my $stopcodon;
my %featureLR=();
my %featureRL=();
##Format: %border=($seqid => pos => count)
my %border=();
my $minpos='';
my $maxpos='';
my $standinput;
my $WAprintnotes='';
##Format: %privisionalmrna=('mrna1' => ('leftutr' => ($start, $end), 
##                        			    'rightUTR' => ($start, $end),
##                                      'utr5' => ($start, $end),
##                                      'utr3' => ($start, $end),
##                                      'exons' => (($start, $end), 
##                                                  ($start, $end), ...)
##                                      'startcodon' => ($start, $end)
##                                      'stopcodon' => ($start, $end)
##                         'mrna2' => 
my %privisionalmrna=();
my $mrna_startcode='01';
my $gff2col='wheat3DL_FuhaoLu';
my %fastalength=();
my $test=0;
my %abinitiogff=();
my %evmgff=();


### Index fasta
unless (-s $fastafile) {
	die "Error: invalid fasta file\n";
}
unless (IndexFasta($fastafile)) {
	die "Error: index fasta failed: $fastafile\n";
}
unless (open(FASTAINDEX, "< $fastafile.fai")) {
	die "Error: failed to read fasta index: $fastafile.fai\n";
}
while (my $line=<FASTAINDEX>) {
	chomp $line;
	if ($line=~/^(\S+)\s+(\d+)/) {
		$fastalength{$1}=$2;
	}
}
close FASTAINDEX;
my $fastaDB=Bio::DB::Fasta->new("$fastafile");
#my $seq = $fastaDB->seq('lyrata', 1 => 108);



### Get region
if ($region=~/:/) {
	my @arr1=split(/:/, $region);
	$seqid=shift @arr1;
	if ($arr1[0]=~/^\S+$/) {
		if ($arr1[0]=~/-/) {
			my @arr2=split(/-/, $arr1[0]);
			($start, $end)=@arr2;
		}
		else {
			$start=$arr1[0];
		}
	}
}
if ($seqid=~/^\S+$/) {
	print "\n\n\n### SUMMART ###\n";
	print "\tSeqID: $seqid\n";
}
else {
	die "Error: invalid SeqID\n";
}
if ($start=~/^\d+$/) {
	print "\tStart: $start\n";
}
else {
	print "\tStart: undefined\n";
}
if ($end=~/^\d+$/) {
	print "\tEnd  : $end\n";
}
else {
	print "\tEND  : undefined\n";
}
print "### SUMMARY ###\n\n\n";



### Retrieve exon borders
foreach my $indgff (@gfffiles) {
	unless (-s $indgff) {
		print "Warnings: can not find GFF file: $indgff\n";
		next;
	}
	print "Info: GFF: $indgff\n";
	close GFFIN if (defined fileno(GFFIN));
	if ($indgff=~/\.gz$/i) {
		unless (open(GFFIN, "zcat $indgff |")) {
			die "Error: failed to open GZip $indgff\n";
		}
	}
	elsif ($indgff=~/\.gff\d{0,1}$/i) {
		unless (open(GFFIN, "cat $indgff |")) {
			die "Error: failed to open GFF $indgff\n";
		}
	}
	else {
		die "Error: failed to guess GFF format $indgff\n";
	}
	my $linenum=0;
	while (my $line=<GFFIN>) {
		$linenum++;
		chomp $line;
		next if ($line=~/^#/);
		my @arr3=split(/\t/, $line);
		if (scalar(@arr3)!=9) {
			die "Error: col !=9 at line $linenum in GFF $indgff\n";
		}
		next if ($arr3[0] ne $seqid);
		next if ($arr3[2] =~ /^(gene)|(mRNA)$/i);
		if ($start=~/^\d+$/) {
			next if ($arr3[3]<$start);
		}
		if ($end=~/^\d+$/) {
			next if ($arr3[4]>$end);
		}
		if ($minpos=~/^\d+$/) {
			if ($arr3[3]<$minpos) {
				$minpos=$arr3[3];
			}
		}
		else {
			$minpos=$arr3[3];
		}
		if ($maxpos=~/^\d+$/) {
			if ($arr3[4]>$maxpos) {
				$maxpos=$arr3[4];
			}
		}
		else {
			$maxpos=$arr3[4];
		}
		if ($arr3[2] =~ /^exon$/i) {
			$featureLR{$arr3[0]}{$arr3[3]}{$arr3[4]}{'exon'}++;
			$featureRL{$arr3[0]}{$arr3[4]}{$arr3[3]}{'exon'}++;
			$border{$arr3[0]}{$arr3[3]}++;
			$border{$arr3[0]}{$arr3[4]}++;
			if ($indgff=~/\.abinitio\./i) {
				$abinitiogff{$arr3[0]}{$arr3[3]}{$arr3[4]}++;
			}
			if ($indgff=~/\.evm\./i) {
				$evmgff{$arr3[0]}{$arr3[3]}{$arr3[4]}++;
			}
		}
		elsif ($arr3[2] =~ /^CDS$/i) {
			$featureLR{$arr3[0]}{$arr3[3]}{$arr3[4]}{'CDS'}++;
			$featureRL{$arr3[0]}{$arr3[4]}{$arr3[3]}{'CDS'}++;
		}
	}
	close GFFIN;
}

#print "Test \%featureLR\n";print Dumper \%featureLR;print "\n"; ### For test ###
#print "Test \%featureRL\n";print Dumper \%featureRL;print "\n"; ### For test ###
#print "Test \%border\n";print Dumper \%border;print "\n"; ### For test ###



###Edit Gene start and end
print "\n\n\nMin: $minpos\tMAX: $maxpos\n";
print "Need to edit GeneStart (MIN: $minpos) ? [NEW number] or [return] :";
$test=0;
while ($test==0) {
	$standinput=&GetStandin;
	if ($standinput=~/^\d+$/) {
		$minpos=$standinput;
		$test=1;
		print "Input: Gene Start accepted : $standinput\n\n";
	}
	elsif ($standinput eq '' or $standinput eq 'q') {
		$test=1;
		print "Input: Gene Start accepted : $minpos\n\n";
	}
}

print "\nNeed to edit GeneEnd (MAX: $maxpos) ? [NEW number] or [return] :";
$test=0;
while ($test==0) {
	$standinput=&GetStandin;
	if ($standinput=~/^\d+$/) {
		$maxpos=$standinput;
		$test=1;
		print "Input: Gene End accepted : $standinput\n";
	}
	elsif ($standinput eq '' or $standinput eq 'q') {
		$test=1;
		print "Input: Gene End accepted : $maxpos\n";
	}
}
$privisionalmrna{$mrna_startcode}{'min'}=$minpos;
$privisionalmrna{$mrna_startcode}{'max'}=$maxpos;



### Define strandness
$strand='+';
$test=0;
while ($test==0) {
	print "\n\n\nStrand ([+]/-) : ";
	$strand=&GetStandin;
	if ($strand eq '-') {
		$strand='-';
	}
	elsif ($strand eq '+') {
		$strand='+';
	}
	if ($strand=~/^(-)|(\+)$/) {
		$test=1;
	}
}
print "Strand accepted: $strand\n";
$privisionalmrna{$mrna_startcode}{'exons'}=[];



### Filter border
my @tempseqids=keys %border;
foreach my $seqname (@tempseqids) {
	my @arr5=keys %{$border{$seqname}};
	foreach (@arr5) {
		if ($_<$minpos or $_>$maxpos) {
			delete $border{$seqname}{$_} if (exists $border{$seqname}{$_});
		}
		if (exists $border{$seqname}{$_} and $border{$seqname}{$_} < $mincount) {
			delete $border{$seqname}{$_};
		}
	}
	if (scalar(keys %{$border{$seqname}})<1) {
		delete $border{$seqname};
	}
}
#print "Test \%border\n";print Dumper \%border;print "\n"; ### For test ###



### Cleaning %featureLR and %featureRL hash
my @arr5=keys %featureLR;
foreach my $seqname (@arr5) {
	my @arr6=keys %{$featureLR{$seqname}};
	foreach (@arr6) {
		if (exists $border{$seqname} and exists $border{$seqname}{$_}) {
			my @arr6=keys %{$featureLR{$seqname}{$_}};
			foreach my $temp (@arr6) {
				delete $featureLR{$seqname}{$_}{$temp} unless (exists $border{$seqname}{$temp});
			}
			if (scalar(keys %{$featureLR{$seqname}{$_}})<1) {
				delete $featureLR{$seqname}{$_};
			}
		}
		else {
			delete $featureLR{$seqname}{$_};
		}
	}
	if (scalar(keys %{$featureLR{$seqname}})<1) {
		delete $featureLR{$seqname};
	}
	
	@arr6=();
	@arr6=keys %{$featureRL{$seqname}};
	foreach (@arr6) {
		if (exists $border{$seqname} and exists $border{$seqname}{$_}) {
			my @arr6=keys %{$featureRL{$seqname}{$_}};
			foreach my $temp (@arr6) {
				delete $featureRL{$seqname}{$_}{$temp} unless (exists $border{$seqname}{$temp});
			}
			if (scalar(keys %{$featureRL{$seqname}{$_}})<1) {
				delete $featureRL{$seqname}{$_};
			}
		}
		else {
			delete $featureRL{$seqname}{$_};
		}
	}
	if (scalar(keys %{$featureRL{$seqname}})<1) {
		delete $featureRL{$seqname};
	}
}
#print "Test \%featureLR\n";print Dumper \%featureLR;print "\n"; ### For test ###
#print "Test \%featureRL\n";print Dumper \%featureRL;print "\n"; ### For test ###



##Format: %exons2keep=($seqID => ($start => ($end => ++))
my %exons2keep=();
my %predictionselect=();
print "\n\nABINITIO GFF: \n";
if (exists $abinitiogff{$seqid} and scalar(keys %{$abinitiogff{$seqid}})>0) {
	$predictionselect{1}='ABINITIO';
	my $abinitiono=0;
	foreach my $abinitiostart (sort {$a<=> $b} keys %{$abinitiogff{$seqid}}) {
		foreach  my $abinitioend (sort {$a<=> $b} keys %{$abinitiogff{$seqid}{$abinitiostart}}) {
			print "\t$abinitiono.\t$abinitiostart-$abinitioend\tCount: ";
			if (exists $featureLR{$seqid} and exists $featureLR{$seqid}{$abinitiostart} and exists $featureLR{$seqid}{$abinitiostart}{$abinitioend} and exists $featureLR{$seqid}{$abinitiostart}{$abinitioend}{'exon'}) {
				print $featureLR{$seqid}{$abinitiostart}{$abinitioend}{'exon'}."\n";
			}
			else {
				print "?\n";
			}
			$abinitiono++;
		}
	}
}
print "\n\nEVM GFF: \n";
if (exists $evmgff{$seqid} and scalar(keys %{$evmgff{$seqid}})>0) {
	$predictionselect{2}='EVM';
	my $abinitiono=0;
	foreach my $abinitiostart (sort {$a<=> $b} keys %{$evmgff{$seqid}}) {
		foreach  my $abinitioend (sort {$a<=> $b} keys %{$evmgff{$seqid}{$abinitiostart}}) {
			print "\t$abinitiono.\t$abinitiostart-$abinitioend\tCount: ";
			if (exists $featureLR{$seqid} and exists $featureLR{$seqid}{$abinitiostart} and exists $featureLR{$seqid}{$abinitiostart}{$abinitioend} and exists $featureLR{$seqid}{$abinitiostart}{$abinitioend}{'exon'}) {
				print $featureLR{$seqid}{$abinitiostart}{$abinitioend}{'exon'}."\n";
			}
			else {
				print "?\n";
			}
			$abinitiono++;
		}
	}
}



my $orfmin=0;
my $orfcount=0;
my $orfmax=0;
if (scalar(keys %predictionselect)>0) {
	$test=0;
	while ($test==0) {
		$orfmin=0;
		$orfcount=0;
		$orfmax=0;
		print "\n\n\nABINITIO or EVM annotations: \n";
		foreach (sort {$a <=>$b} keys %predictionselect) {
			print "\t", $_, ".\t", $predictionselect{$_}, "\n";
		}
		print "Keep ABINITIO or EVM annotations? [0=NO ...] : ";
		my $predictor=&GetStandin;
		$predictor=2 if ($predictor eq '');
		if ($predictor eq '0') {
			print "Info: Ignore predicted GFF, user would define exons\n";
			$test=1;
		}
		elsif (exists $predictionselect{$predictor}) {
			if ($predictionselect{$predictor} eq 'ABINITIO') {
				foreach my $abinitiostart (sort {$a<=> $b} keys %{$abinitiogff{$seqid}}) {
					if ($orfcount==0) {
						$orfmin=$abinitiostart;
						$orfcount++;
					}
					foreach my $abinitioend (sort {$a<=> $b} keys %{$abinitiogff{$seqid}{$abinitiostart}}) {
						$exons2keep{$seqid}{$abinitiostart}{$abinitioend}++;
						if ($abinitioend>$orfmax) {
							$orfmax=$abinitioend;
						}
					}
				}
				$orfmax-=2;
				$test=1;
			}
			elsif ($predictionselect{$predictor} eq 'EVM') {
				foreach my $abinitiostart (sort {$a<=> $b} keys %{$evmgff{$seqid}}) {
					if ($orfcount==0) {
						$orfmin=$abinitiostart;
						$orfcount++;
					}
					foreach  my $abinitioend (sort {$a<=> $b} keys %{$evmgff{$seqid}{$abinitiostart}}) {
						$exons2keep{$seqid}{$abinitiostart}{$abinitioend}++;
						if ($abinitioend>$orfmax) {
							$orfmax=$abinitioend;
						}
					}
				}
				$orfmax-=2;
				$test=1;
			}
		}
		else {
			print "Info: KEY not defined, input KEY NUMBER again ?\n";
		}
	}
}
%abinitiogff=();
%evmgff=();



### Define start codon
my $test4=0;
while ($test4==0) {
	$test=0;
	while ($test==0) {
		print "\n\n\nStart codon ";
		if ($strand eq '+') {
			print "(DEFAULT: $orfmin) " if ($orfmin>0);
		}
		elsif ($strand eq '-') {
			print "(DEFAULT: $orfmax) " if ($orfmax>0);
		}
		print " [y]=use deaflt | [INT]= custom | 0/? = unknown : ";
		$startcodon=&GetStandin;
		if ($startcodon eq '' or $startcodon eq 'y') {
			$startcodon=($strand eq '+') ? $orfmin : $orfmax;
			$test=1;
		}
		elsif ($startcodon=~/^\d+$/ and $startcodon>0) {
			$test=1;
		}
		else {
			print "Start codon undefined ? [y]|n : ";
			my $teststart=&GetStandin;
			if ($teststart eq '' or $teststart eq 'y') {
				$startcodon='?';
				$test=1;
			}
		}
		if ($startcodon=~/^\d+$/) {
			my $startseq=$fastaDB->seq($seqid, $startcodon => $startcodon+2);
			$startseq=SeqRevComp($startseq) if ($strand eq '-');
			if (defined $startseq and $startseq =~/^atg$/i) {
				print "Start codon accepted: $startcodon\n";
				$test=1;
			}
			else {
				print STDERR colored("Error: invalid StartCodon: $startseq", "$error_color_code"), "\n";
				$test=0;
			}
		}
	}

	### Define stop codon
	$test=0;
	while ($test==0) {
		print "\n\n\nStop codon : ";
		if ($strand eq '+') {
			print "(DEFAULT: $orfmax) " if ($orfmax>0);
		}
		elsif ($strand eq '-') {
			print "(DEFAULT: $orfmin) " if ($orfmin>0);
		}
		print " [y]=use deaflt | [INT]= custom | 0/? = unknown : ";
		$stopcodon=&GetStandin;
		if ($stopcodon eq '' or $stopcodon eq 'y') {
			$stopcodon=($strand eq '+') ? $orfmax : $orfmin;
			$test=1;
		}
		elsif ($stopcodon=~/^\d+$/ and $stopcodon>0) {
			$test=1;
		}
		else {
			print "Stop codon undefined ? [y]|n : ";
			my $teststop=&GetStandin;
			if ($teststop eq '' or $teststop eq 'y') {
				$stopcodon='?';
				$test=1;
			}
		}
		if ($stopcodon=~/^\d+$/) {
			my $stopseq=$fastaDB->seq($seqid, $stopcodon => $stopcodon+2);
			$stopseq=SeqRevComp($stopseq) if ($strand eq '-');
			if (defined $stopseq and $stopseq =~/^(tga)|(tag)|(taa)$/i) {
				print "Stop codon accepted: $stopcodon\n\n\n";
				$test=1;
			}
			else {
				print STDERR colored("Error: invalid StopCodon : $stopseq", "$error_color_code"), "\n";
				$test=0;
			}
		}
	}
	
	
	if ($startcodon=~/^\d+$/ and $stopcodon=~/^\d+$/) {
		if ($strand eq '-') {
			if ($startcodon > $stopcodon) {
				$test4=1;
			}
			elsif ($startcodon < $stopcodon) {
				print "\n***Warnings: startcodon $startcodon > stopcodon $stopcodon on strand $strand, Switch it ? [y]|n :\n";
				$standinput=&GetStandin;
				if ($standinput eq '' or $standinput eq 'y') {
					my $tempstart=$startcodon;
					$startcodon=$stopcodon;
					$stopcodon=$tempstart;
					print "Info: NEW startcodon: $startcodon\nInfo: NEW stopcodon $stopcodon\n";
					$test4=1;
				}
			}
			else {
				print STDERR "\n***", colored("Error: invalid StartCodon ($startcodon) or Stopcodon ($stopcodon) on strand $strand", "$error_color_code"), "\n";
			}
		}
		elsif ($strand eq '+') {
			if ($startcodon > $stopcodon) {
				print STDERR "\n", colored("***Warnings: startcodon $startcodon < stopcodon $stopcodon on strand $strand, Switch it ? [y]|n :", "$error_color_code"), "\n";
				$standinput=&GetStandin;
				if ($standinput eq '' or $standinput eq 'y') {
					my $tempstart=$startcodon;
					$startcodon=$stopcodon;
					$stopcodon=$tempstart;
					print "Info: NEW startcodon: $startcodon\nInfo: NEW stopcodon $stopcodon\n";
					$test4=1;
				}
			}
			elsif ($startcodon < $stopcodon) {
				$test4=1;
			}
			else {
				print STDERR "\n***", colored("Error: invalid StartCodon ($startcodon) or Stopcodon ($stopcodon) on strand $strand", "$error_color_code"), "\n";
			}
		}
	}
	else {
		$test4=1;
	}
}



### Define UTR
print "###### UTR summary #####\n";
if ($strand eq '+') {
	if ($startcodon=~/^\d+$/) {
		print "Start codon: $startcodon - ".($startcodon+2)."\n";
		$privisionalmrna{$mrna_startcode}{'startcodon'}=[$startcodon, $startcodon+2];
		if ($minpos<($startcodon-1)) {
			print "5UTR: $minpos - ".($startcodon-1)."\n";
			my @temparr=($minpos, $startcodon-1);
			push (@{$privisionalmrna{$mrna_startcode}{'utr5'}}, \@temparr);
			push (@{$privisionalmrna{$mrna_startcode}{'leftutr'}}, \@temparr);
		}
	}
	if ($stopcodon=~/^\d+$/) {
		print "Stop codon: $stopcodon - ".($stopcodon+2)."\n";
		$privisionalmrna{$mrna_startcode}{'stopcodon'}=[$stopcodon, $stopcodon+2];
		if (($stopcodon+3)<$maxpos) {
			print "3UTR: ".($stopcodon+3)." - $maxpos\n";
			my @temparr=($stopcodon+3, $maxpos);
			push (@{$privisionalmrna{$mrna_startcode}{'utr3'}}, \@temparr);
			push (@{$privisionalmrna{$mrna_startcode}{'rightutr'}}, \@temparr);
		}
	}
}
elsif ($strand eq '-') {
	if ($startcodon=~/^\d+$/) {
		print "Start codon: $startcodon - ".($startcodon+2)."\n";
		$privisionalmrna{$mrna_startcode}{'startcodon'}=[$startcodon, $startcodon+2];
		if (($startcodon+3)<$maxpos) {
			print "5UTR: ".($startcodon+3)." - $maxpos\n";
			my @temparr=($startcodon+3, $maxpos);
			push (@{$privisionalmrna{$mrna_startcode}{'utr5'}}, \@temparr);
			push (@{$privisionalmrna{$mrna_startcode}{'rightutr'}}, \@temparr);
		}
	}
	if ($stopcodon=~/^\d+$/) {
		print "Stop codon: $stopcodon - ".($stopcodon+2)."\n";
		$privisionalmrna{$mrna_startcode}{'stopcodon'}=[$stopcodon, $stopcodon+2];
		if ($minpos<($stopcodon-1)) {
			print "3UTR: $minpos - ".($stopcodon-1)."\n";
			my @temparr=($minpos, $stopcodon-1);
			push (@{$privisionalmrna{$mrna_startcode}{'utr3'}}, \@temparr);
			push (@{$privisionalmrna{$mrna_startcode}{'leftutr'}}, \@temparr);
		}
	}
}

my $exonno=0;
foreach my $start2 (sort {$a <=>$b} keys %{$exons2keep{$seqid}}) {
	foreach my $end2 (sort {$a<=>$b} keys %{$exons2keep{$seqid}{$start2}}) {
		print "Exon: \t$exonno\t$start2-$end2\n";
		$exonno++;
	}
}
print "\n\nAdd UTR:\n";
if (scalar(keys %{$exons2keep{$seqid}})==1) {
	my @arrleft=keys %{$exons2keep{$seqid}};
	my $geneleft=$arrleft[0];
	if ($geneleft=~/^\d+$/ and $geneleft != $minpos) {
		print "Test: Geneleft\tMIN: $minpos\n"; ### for test ###
		my @arrtight=keys %{$exons2keep{$seqid}{$geneleft}};
		if (scalar(@arrtight)==1) {
			my $generight=$arrtight[0];
			if ($generight=~/^\d+$/) {
				print "\n\n\nMIN: $minpos\tMAX: $maxpos\n";
				print "Change Exon $geneleft-$generight to $minpos-$maxpos [y/1] | n/0 : ";
				$standinput=&GetStandin;
				unless ($standinput=~/(^n$)|(^no$)|(^0$)/i) {
					delete $exons2keep{$seqid}{$geneleft};
					$exons2keep{$seqid}{$minpos}{$maxpos}++;
				}
			}
		}
	}
}
elsif (scalar(keys %{$exons2keep{$seqid}})>1) {
	my @arrleft1=sort {$a <=>$b} keys %{$exons2keep{$seqid}};
	my $geneleft=shift @arrleft1;
	my @arrleft2=sort {$a <=>$b} keys %{$exons2keep{$seqid}{$geneleft}};
	if ($geneleft=~/^\d+$/) {
		my @arrleft2=sort {$a <=>$b} keys %{$exons2keep{$seqid}{$geneleft}};
		foreach (@arrleft2) {
			print "\n\n\nMIN: $minpos\n";
			print "Change LEFT Exon $geneleft-$_ to $minpos-$_ [y/1] | n/0 : ";
			$standinput=&GetStandin;
			unless ($standinput=~/(^n$)|(^no$)|(^0$)/i) {
				delete $exons2keep{$seqid}{$geneleft}{$_};
				$exons2keep{$seqid}{$minpos}{$_}++;
			}
		}
		if (scalar(keys %{$exons2keep{$seqid}{$geneleft}}<1)) {
			delete $exons2keep{$seqid}{$geneleft};
		}
	}
	
	
	
	my $generight=0;
	foreach my $start2 (sort {$a <=>$b} keys %{$exons2keep{$seqid}}) {
		foreach my $end2 (sort {$a<=>$b} keys %{$exons2keep{$seqid}{$start2}}) {
			$generight=$end2 if ($end2>$generight);
		}
	}
	my @arrright2=();
	foreach my $start2 (sort {$a <=>$b} keys %{$exons2keep{$seqid}}) {
		foreach my $end2 (sort {$a<=>$b} keys %{$exons2keep{$seqid}{$start2}}) {
			push (@arrright2, $start2) if ($end2==$generight);
		}
	}
	foreach (@arrright2) {
		print "\n\n\nMAX: $maxpos\n";
		print "Change RIGHT Exon $_-$generight to $_-$maxpos [y/1] | n/0 : ";
		$standinput=&GetStandin;
		unless ($standinput=~/(^n$)|(^no$)|(^0$)/i) {
			delete $exons2keep{$seqid}{$_}{$generight};
			$exons2keep{$seqid}{$_}{$maxpos}++;
		}
	}
}



#print Dumper \%featureLR; ### For test ###



### Select exons
$test=0;
while ($test==0) {
	print "\n\n\nExons seltected: \n";
	$exonno=0;
	my %exonhash=();
	foreach my $start2 (sort {$a <=>$b} keys %{$exons2keep{$seqid}}) {
		foreach my $end2 (sort {$a<=>$b} keys %{$exons2keep{$seqid}{$start2}}) {
			print "Exon: \t$exonno\t$start2-$end2\t\t";
			print "Alt: ";
			if (exists $featureLR{$seqid} and exists $featureLR{$seqid}{$start2}) {
				foreach (sort {$a<=>$b} keys %{$featureLR{$seqid}{$start2}}) {
					unless (exists $exons2keep{$seqid}{$start2}{$_}) {
						print "$start2-$_";
						if (exists $featureLR{$seqid}{$start2}{$_}{'exon'} and $featureLR{$seqid}{$start2}{$_}{'exon'} =~/^\d+$/) {
							print '[', $featureLR{$seqid}{$start2}{$_}{'exon'}, ']';
						}
						print "\t";
					}
				}
			}
			if (exists $featureRL{$seqid} and exists $featureRL{$seqid}{$end2}) {
				foreach (sort {$a<=>$b} keys %{$featureRL{$seqid}{$end2}}) {
					unless (exists $exons2keep{$seqid}{$_}{$end2}) {
						print "$_-$end2";
						if (exists $featureRL{$seqid}{$end2}{$_}{'exon'} and $featureRL{$seqid}{$end2}{$_}{'exon'}=~/^\d+$/) {
							print '[', $featureRL{$seqid}{$end2}{$_}{'exon'}, ']';
						}
						print "\t";
					}
				}
			}
			print "\n";
			
			$exonhash{$exonno++}="$start2-$end2";
		}
	}
	print "Info: Min: $minpos\tMAX: $maxpos\n";
	print "Input: exon border \n\t[start-end]\n\t[start]\n\t[End]\n\t[d]index\n\t[j]index1+index2\n\t[r]index1+index2\n\t[q]uit: \n\n\nYour choice: ";
	my $position=&GetStandin;
	if ($position=~/^(\d+)-(\d+)$/) {
		if ($1<$2) {
			print "Accepted exon: $position\n";
			$exons2keep{$seqid}{$1}{$2}++;
		}
		else {
			print STDERR colored("Warnings: invalid numbers: $position", "$error_color_code"), "\n";
		}
	}
	elsif ($position=~/^\d+$/) {
		my $num2choose=0;
		my %hash2choose=();
		if (exists $featureLR{$seqid} and exists $featureLR{$seqid}{$position}) {
			foreach (sort {$a<=>$b} keys %{$featureLR{$seqid}{$position}}) {
				print "Start: \t$num2choose\t$position-$_\tCount: $featureLR{$seqid}{$position}{$_}{'exon'}\n";
				$hash2choose{$num2choose++}="$position-$_";
			}
		}
		if (exists $featureRL{$seqid} and exists $featureRL{$seqid}{$position}) {
			foreach (sort {$a<=>$b} keys %{$featureRL{$seqid}{$position}}) {
				print "End: \t$num2choose\t$_-$position\tCount: $featureRL{$seqid}{$position}{$_}{'exon'}\n";
				$hash2choose{$num2choose++}="$_-$position";
			}
		}
		if ($num2choose==0) {
			print STDERR colored("NO suggestion for this exon border, TRY [start-end]", "$error_color_code"), "\n";
			next;
		}
		else {
			print "Select which exons to keep (dot delimited) : ";
			$standinput=&GetStandin;
			my @choice=split(/\./, $standinput);
			foreach my $indchoice (@choice) {
				if (exists $hash2choose{$indchoice}) {
					print "Accepted exon: $hash2choose{$indchoice}\n";
					my @arrexon=split(/-/, $hash2choose{$indchoice});
					$exons2keep{$seqid}{$arrexon[0]}{$arrexon[1]}++;
				}
			}
		}
	}
	elsif ($position=~/^d(\d+)$/) {
		if (exists $exonhash{$1}){
			print "Delete exon $1: $exonhash{$1}\n";
			my @arrexon=split(/-/, $exonhash{$1});
			if (exists $exons2keep{$seqid}{$arrexon[0]} and $exons2keep{$seqid}{$arrexon[0]}{$arrexon[1]}) {
				delete $exons2keep{$seqid}{$arrexon[0]}{$arrexon[1]};
				if (scalar(keys %{$exons2keep{$seqid}{$arrexon[0]}})==0) {
					delete $exons2keep{$seqid}{$arrexon[0]};
				}
			}
		}
		
	}
	elsif ($position=~/^j(\d+)\+(\d+)$/) {
		my @exonarr1=split(/-/, $exonhash{$1});
		my $exonadd1=shift @exonarr1;
		my @exonarr2=split(/-/, $exonhash{$2});
		my $exonadd2=pop @exonarr2;
		$exons2keep{$seqid}{$exonadd1}{$exonadd2}++;
	}
	elsif ($position=~/^r(\d+)\+(\d+)$/) {
		my ($exonadd1, $exonlft1)=split(/-/, $exonhash{$1});
		my ($exonadd2, $exonlft2)=split(/-/, $exonhash{$2});
		$exons2keep{$seqid}{$exonadd1}{$exonlft2}++;
		delete $exons2keep{$seqid}{$exonadd1}{$exonlft1};
		delete $exons2keep{$seqid}{$exonadd1} if (scalar(keys %{$exons2keep{$seqid}{$exonadd1}})<1);
		delete $exons2keep{$seqid}{$exonadd2}{$exonlft2};
		delete $exons2keep{$seqid}{$exonadd2} if (scalar(keys %{$exons2keep{$seqid}{$exonadd2}})<1);
	}
	else {
		print "Quit exons editing? [y]|n : ";
		$standinput=&GetStandin;
		unless ($standinput eq 'n') {
			if (scalar(keys %{$exons2keep{$seqid}})>0) {
				$test=1;
			}
			else {
				print STDERR colored("Error: no exons selected", "$error_color_code"), "\n";
			}
		}
	}
}




#print "Test: \%exons2keep\n";print Dumper \%exons2keep; print "\n"; ### For test ###
%featureLR=();
%featureRL=();


### join Exons to provisional mRNA
my $randomctrl=1;
foreach my $seqname (keys %exons2keep) {
	foreach my $start4 (sort {$a<=>$b} keys %{$exons2keep{$seqname}}) {
		my @endnum=sort {$a <=> $b} keys %{$exons2keep{$seqname}{$start4}};
#		print "Test: SeqID:Start $seqname:$start4 ends: @endnum\n"; ### For test ###
#		if ($start4==1387) {print "Test: \%privisionalmrna\n"; print Dumper \%privisionalmrna; print "\n";} ### For test ###
		my @mrnagroups=keys %privisionalmrna;
		foreach my $group (@mrnagroups) {
			if (exists $privisionalmrna{$group}{'exons'} and scalar(@{$privisionalmrna{$group}{'exons'}})>0 and defined $privisionalmrna{$group}{'exons'}[-1] and defined $privisionalmrna{$group}{'exons'}[-1][-1] and $privisionalmrna{$group}{'exons'}[-1][-1]=~/^\d+$/ and $start4<=$privisionalmrna{$group}{'exons'}[-1][-1]) {
				next;
			}
			my @endnum2=@endnum;
			my @expandgroups=();
			if (scalar(@endnum)>1) {
				push (@expandgroups, $group);
				for (my $i=1; $i<scalar(@endnum); $i++) {
					$privisionalmrna{"$group.$randomctrl"}=dclone \%{$privisionalmrna{$group}};
					push (@expandgroups, "$group.$randomctrl");
					$randomctrl++;
				}
			}
			elsif (scalar(@endnum)==1) {
				push (@expandgroups, $group);
			}
			else {
				print STDERR colored("Warnings: SeqID:Start $seqname:$start4 has no ends", "$error_color_code"), "\n";
			}
#			print "Test: expand groups: @expandgroups\n"; ### For test ###
			foreach my $indgrp (@expandgroups) {
				my $thisend=shift @endnum2;
#				print "Test: start - end : $start4 - $thisend\n"; ### For test ###
#				print "Test: Group: $indgrp, Add array $start4 - $thisend\n"; ### For test ###
				push (@{$privisionalmrna{$indgrp}{'exons'}}, [$start4, $thisend]);
			}
		}
#		if ($start4==1387) {print "Test: \%privisionalmrna\n"; print Dumper \%privisionalmrna; print "\n";} ### For test ###
	}
}
#print "Test: \%privisionalmrna\n"; print Dumper \%privisionalmrna; print "\n"; ### For test ###


unless (&PlotGeneStructure(\%privisionalmrna)) {
	print STDERR colored("Error: Draw PlotGeneStructure failed", "$error_color_code"), "\n";
	exit 1;
}

#print "Test: \%privisionalmrna\n"; print Dumper \%privisionalmrna; print "\n"; ### For test ###

$test=0;
while ($test==0) {
	print "Ploting  ... \n";
	print "\nNeed to edit? [y]|n : ";
	$standinput=&GetStandin;
	if ($standinput eq 'n') {
		$test=1;
		last;
	}
	else {
		my $test2=0;
		while ($test2==0) {
			print "########## EDIT MODE ##########\nAll mRNAs:\n";
			print "Ploting  ... \n";
			&PrintAnnotation(\%privisionalmrna);
			my $test3=0;
			my $test4_edit=0;
			while ($test3==0) {
				print "Input the mRNA ID ("; 
				map {print $_.' ' } (sort keys %privisionalmrna);
				print "or [q]uit if ID not existed) (DEFAULT: q): ";
				my $mrnaid2edit=&GetStandin;
				if ($mrnaid2edit eq '' or $mrnaid2edit=~/(^q$)|(^quit$)/i) {
					$test3=1;
					last;
				}
				unless (exists $privisionalmrna{$mrnaid2edit}) {
					print STDERR "\n\n\n", colored("Error: mRNA ID not existed : $mrnaid2edit", "$error_color_code"), "\n\n\n";
					last;
				}
				print "Actions: [c]opy, [e]dit, [d]eletion, [q]uit (DEFAULT: e) ? ";
				$standinput=&GetStandin;
				if ($standinput eq 'c' or $standinput eq 'copy') {
					$privisionalmrna{"$mrnaid2edit.$randomctrl"}=dclone $privisionalmrna{$mrnaid2edit};
					print "mRNA ID created: $mrnaid2edit.$randomctrl\n";
					$randomctrl++;
					last;
				}
				elsif ($standinput eq 'e' or $standinput eq 'edit' or $standinput eq '') {
					$privisionalmrna{$mrnaid2edit}=&EditAnnotation(\%{$privisionalmrna{$mrnaid2edit}});
					last;
				}
				elsif ($standinput eq 'd' or $standinput eq 'delete') {
					delete $privisionalmrna{$mrnaid2edit};
					last;
				}
				elsif ($standinput eq 'q' or $standinput eq 'quit') {
					$test3=1;
					$test4_edit=1;
					last;
				}
			}
			print "New Plotting ...";
			unless (&PlotGeneStructure(\%privisionalmrna)) {
				print STDERR colored("Error: Draw PlotGeneStructure failed", "$error_color_code"), "\n";
				exit 1;
			}
			print "finished\n";
			if ($test4_edit==0) {
				print "\nStay in EDIT MODE, y|[n] ?  : ";
				$standinput=&GetStandin;
				if ($standinput eq 'y') {
					next;
				}
				elsif ($standinput eq '' or $standinput eq 'n') {
					print "########## DISPLAY MODE ##########\n";
					$test2=1;
				}
			}
		}
	}
}



&WriteAnnotation(\%privisionalmrna);







my $exonnum=0;
my %startcodonlist=();
my %stopcodonlist=();
foreach my $sepid (keys %privisionalmrna) {
#	print "mRNA: $sepid";
	if (exists $privisionalmrna{$sepid}{'exons'} and scalar(@{$privisionalmrna{$sepid}{'exons'}})>$exonnum) {
		$exonnum=scalar(@{$privisionalmrna{$sepid}{'exons'}});
	}
	if (exists $privisionalmrna{$sepid}{'startcodon'} and defined $privisionalmrna{$sepid}{'startcodon'}[0] and $privisionalmrna{$sepid}{'startcodon'}[0]=~/^\d+$/) {
#		print "\tStartcodon\t$privisionalmrna{$sepid}{'startcodon'}[0]";
		$startcodonlist{$privisionalmrna{$sepid}{'startcodon'}[0]}++;
	}
	if (exists $privisionalmrna{$sepid}{'stopcodon'} and defined $privisionalmrna{$sepid}{'stopcodon'}[0] and $privisionalmrna{$sepid}{'stopcodon'}[0]=~/^\d+$/) {
#		print "\tStopcodon\t$privisionalmrna{$sepid}{'stopcodon'}[0]";
		$stopcodonlist{$privisionalmrna{$sepid}{'stopcodon'}[0]}++;
	}
	print "\n";
}
print "\n\n\nmin\tmax\tstartcodon\tstopcodon\tMaxExons\tNote\n";
print "\n$minpos\t$maxpos\t";
if (scalar(keys %startcodonlist)>0){
	print join('/', sort {$a<=>$b} keys %startcodonlist), "\t";
}
else {
	print "NaN\t";
}
if (scalar(keys %stopcodonlist)>0){
	print join('/', sort {$a<=>$b} keys %stopcodonlist), "\t";
}
else {
	print "NaN\t";
}
print $exonnum, "\t", $WAprintnotes, "\n";

exit 0;



#####################################################################
########################### SUB functions ###########################
#####################################################################
sub ChooseHash {
	my $CHhashin=shift @_;
	
	my $CHsubinfo='SUB(ChooseHash)';
	my $CHstartcode=0;
	my $CHtestfinish=0;
	my %CHhashselector=();
	my %tempquick=();
	my %temphash=();
	my %finalhash=();
	
	while ($CHtestfinish==0) {
		$CHstartcode=0;
		%CHhashselector=();
		%tempquick=();
		foreach my $CHstart (sort {$a <=> $b} keys %{$CHhashin}) {
			foreach my $CHend (sort {$a <=>$b} keys %{${$CHhashin}{$CHstart}}) {
				$CHhashselector{$CHstartcode}="$CHstart  ---  $CHend";
				$tempquick{$CHstartcode}{'start'}=$CHstart;
				$tempquick{$CHstartcode}{'end'}=$CHend;
				$CHstartcode++;
			}
		}
		my $CHtestchooser=0;
		while ($CHtestchooser==0) {
			%temphash=%CHhashselector;
			print "\n\n\n", $CHsubinfo, ": Choose Exons to REMOVE (dot delimited list): \n";
			foreach (sort {$a <=>$b} keys %temphash) {
				print "\t", $_, ".\t". $temphash{$_}. "\n";
			}
			print "\n", $CHsubinfo, ": Input selection: ";
			my $CHoptions=&GetStandin;
			my @CHarr=split(/\./, $CHoptions);
			foreach (@CHarr) {
				if (exists $temphash{$_}) {
					delete $temphash{$_};
				}
				else {
					print "\n", $CHsubinfo, "Warnings: Non exists option: $_\n";
				}
			}
			print "\n", $CHsubinfo, ": Exons after clean: \n";
			foreach (sort {$a <=>$b} keys %temphash) {
				print "\t", $_, ".\t". $temphash{$_}. "\n";
			}
			my $CHteminitor=0;
			print "\n", $CHsubinfo, ": Is this final selection? [y]|n :";
			$CHteminitor=&GetStandin;
			if ($CHteminitor eq '' or $CHteminitor eq 'y') {
				print $CHsubinfo, ": Final selection: accepted\n";
				$CHtestchooser=1;
			}
		}
		my $CH2ndsure=0;
		print $CHsubinfo, ": Sure? [y]|n :";
		$CH2ndsure=&GetStandin;
		if ($CH2ndsure eq '' or $CH2ndsure eq 'y') {
			print "\n", $CHsubinfo, ": Exit selection\n";
			$CHtestfinish=1;
		}
	}
	foreach (sort {$a<=>$b} keys %temphash) {
		if (exists $tempquick{$_}{'start'} and exists $tempquick{$_}{'end'}) {
			$finalhash{$tempquick{$_}{'start'}}{$tempquick{$_}{'end'}}++;
		}
	}
	return \%finalhash;
}


### Global: $minpos, $maxpos, $seqid, $strand
sub PlotGeneStructure {
	my $PGShash=shift;
	
	my $PGSsubinfo='SUB(PlotGeneStructure)';
	my $PGSsvg;
	my $PGSwidthmargin=100;### width margin=200
	my $PGSheightmargin=100; ### height margin = 200
	my $geneheight=50;
	my $exonheight=30;
	my $nummrna=scalar(keys %{$PGShash});
	my $textmargin=10;
	my $textfontsize=14;
	
	if ($nummrna<1) {
		print STDERR $PGSsubinfo, colored("Error: empty hash to plot", "$error_color_code"), "\n";
		return 0;
	}
	
	my $PGSx0 = $PGSwidthmargin - $minpos;
	my $PGSy0 = $PGSheightmargin;
	my $PGSwidth=$maxpos - $minpos + 2 * $PGSwidthmargin;
	my $PGSheight=$nummrna * $geneheight + 2 * $PGSheightmargin;
#	print $PGSsubinfo, "Test: Plot size: width:height $PGSwidth:$PGSheight\n"; ### For test ###
	$PGSsvg = SVG->new(width=>$PGSwidth, height=>$PGSheight);
	my $PGSwhichgene=0;
	foreach my $PGSindgrp (sort keys %{$PGShash}) {
		my $PGStextheight=$geneheight * $PGSwhichgene + $PGSy0 + $exonheight / 2;
		$PGSsvg->text(x => $PGSwidthmargin, y => $PGStextheight, width => 10, height => 10, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=>$textfontsize, "-cdata" => "$PGSindgrp");
		
		if (exists ${$PGShash}{$PGSindgrp}{'exons'} and scalar(@{${$PGShash}{$PGSindgrp}{'exons'}})>0) {
			my $PGSlastexonRight=0;;
			foreach (my $PGSi=0; $PGSi<scalar(@{${$PGShash}{$PGSindgrp}{'exons'}}); $PGSi++) {
				next unless (scalar(@{${$PGShash}{$PGSindgrp}{'exons'}[$PGSi]})==2);
				my ($PGSexonLeft, $PGSexonRight)=@{${$PGShash}{$PGSindgrp}{'exons'}[$PGSi]};
				unless (defined $PGSexonLeft and $PGSexonLeft=~/^\d+$/ and defined $PGSexonRight and $PGSexonRight=~/^\d+$/) {
					print STDERR $PGSsubinfo, colored("Warnings: invalid exon coordinates: $PGSexonLeft - $PGSexonRight for mRNA $PGSindgrp", "$error_color_code"), "\n";
					next;
				}
				my $PGSrectX=$PGSexonLeft+$PGSx0;
				my $PGSrectY=$geneheight * $PGSwhichgene + $PGSy0;
				my $PGSexonwidth=$PGSexonRight-$PGSexonLeft;
				my $PGSexonid=$PGSindgrp."_$PGSexonLeft-$PGSexonRight";
				$PGSsvg->rectangle(x 		=> $PGSrectX, 
								   y 		=> $PGSrectY, 
								   width  	=> $PGSexonwidth, 
								   height => $exonheight, 
								   id=> "$PGSexonid",
								   style => {'fill' => 'rgb(0, 255, 0)',
                                             'stroke'         => 'black',
                                             'stroke-width'   =>  0,
                                             'stroke-opacity' =>  0.50,
                                             'fill-opacity'   =>  1,
                                    		},
                                );
				$PGStextheight=$PGSrectY+ $exonheight + $textmargin;
				$PGSsvg->text(x => $PGSrectX, y => $PGStextheight, width => 10, height => 10, "font-family"=>"Arial", "text-anchor"=>"start","font-size"=>$textfontsize, "-cdata" => "$PGSexonLeft");
				my $PGStestX2=$PGSexonRight+$PGSx0;
				$PGSsvg->text(x => $PGStestX2, y => $PGStextheight, width => 10, height => 10, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=>$textfontsize, "-cdata" => "$PGSexonRight");
				
				if ($PGSi>0) {
					my $PGSthisexonLeft=$PGSrectX-1;
					my $PGSthisexonhight=$geneheight * $PGSwhichgene + $PGSy0 + $exonheight / 2;
#					print $PGSsubinfo, "Test: intron: $PGSlastexonRight - $PGSthisexonLeft\n"; ### For test ###
					$PGSsvg->line(x1 => $PGSlastexonRight, y1 => $PGSthisexonhight, x2 => $PGSthisexonLeft, y2 => $PGSthisexonhight, stroke=>'black',"stroke-width"=>3); 
					###Intron GT-AG
					my $PGSintronleftstart=$PGSlastexonRight-$PGSx0;
					my $PGSintronleftend=$PGSintronleftstart+1;
					my $PGSleftseq = $fastaDB->seq($seqid, $PGSintronleftstart => $PGSintronleftend);
					my $PGSintronrightend=$PGSthisexonLeft-$PGSx0;
					my $PGSintronrightstart=$PGSintronrightend-1;
					my $PGSrightseq= $fastaDB->seq($seqid, $PGSintronrightstart => $PGSintronrightend);
					if ($strand eq '-') {
						$PGSleftseq=reverse $PGSleftseq;
						$PGSleftseq=~tr/actgACTG/TGACTGAC/;
						$PGSrightseq=reverse $PGSrightseq;
						$PGSrightseq=~tr/actgACTG/TGACTGAC/;
					}
					$PGSsvg->text(x => $PGSlastexonRight, y => $PGSthisexonhight-5, width => 10, height => 10, "font-family"=>"Arial", "text-anchor"=>"start","font-size"=>$textfontsize, "-cdata" => "$PGSleftseq");
					$PGSsvg->text(x => $PGSthisexonLeft, y => $PGSthisexonhight-5, width => 10, height => 10, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=>$textfontsize, "-cdata" => "$PGSrightseq");
				}
				$PGSlastexonRight=$PGSexonRight + $PGSx0 + 1;
#				print $PGSsubinfo, "Test: intron: $PGSlastexonRight\n"; ### For test ###
			}
		}
		
		LEFTUTR: {if (exists ${$PGShash}{$PGSindgrp}{'leftutr'} and scalar(@{${$PGShash}{$PGSindgrp}{'leftutr'}[0]})==2) {
			my ($PGSleftutrX1, $PGSleftutrX2, $PGSleftutrY1, $PGSleftutrY2)=(0, 0, 0, 0);
			($PGSleftutrX1, $PGSleftutrX2)=@{${$PGShash}{$PGSindgrp}{'leftutr'}[0]};
			unless (defined $PGSleftutrX1 and $PGSleftutrX1=~/^\d+$/ and defined $PGSleftutrX2 and $PGSleftutrX2=~/^\d+$/) {
				print STDERR $PGSsubinfo, colored("Warnings: invalid leftUTR number: $PGSleftutrX1 - $PGSleftutrX2", "$error_color_code"), "\n";
				last LEFTUTR;
			}
			$PGSleftutrX1+=$PGSx0;
			$PGSleftutrX2+=$PGSx0;
			$PGSleftutrY1=$geneheight * $PGSwhichgene + $PGSy0 + $exonheight / 2;
			$PGSleftutrY2=$PGSleftutrY1;
#			print $PGSsubinfo, "Test: LeftUTR: X1:Y1:X2:Y2 $PGSleftutrX1:$PGSleftutrY1:$PGSleftutrX2:$PGSleftutrY2\n"; ### For test ###
			$PGSsvg->line(x1 => $PGSleftutrX1, y1 => $PGSleftutrY1, x2 => $PGSleftutrX2, y2 => $PGSleftutrY2, stroke=>'red',"stroke-width"=>5);
			$PGStextheight=$PGSleftutrY1+ $exonheight / 2 + $textmargin;
			$PGSsvg->text(x => $PGSleftutrX1, y => $PGStextheight, width => 10, height => 10, "font-family"=>"Arial", "text-anchor"=>"start","font-size"=>$textfontsize, "-cdata" => "${$PGShash}{$PGSindgrp}{'leftutr'}[0][0]");
		}
		else {
#			print $PGSsubinfo, "Test: mRNA $PGSindgrp has no LeftUTR\n"; ### For test ###
		}}###LEFTUTR
		
		RIGHTUTR: {if (exists ${$PGShash}{$PGSindgrp}{'rightutr'} and scalar(@{${$PGShash}{$PGSindgrp}{'rightutr'}[0]})==2) {
			my ($PGSrightutrX1, $PGSrightutrX2, $PGSrightutrY1, $PGSrightutrY2)=(0, 0, 0, 0);
			($PGSrightutrX1, $PGSrightutrX2)=@{${$PGShash}{$PGSindgrp}{'rightutr'}[0]};
			unless (defined $PGSrightutrX1 and $PGSrightutrX1=~/^\d+$/ and defined $PGSrightutrX2 and $PGSrightutrX2=~/^\d+$/) {
				print STDERR $PGSsubinfo, colored("Warnings: invalid rightUTR number: $PGSrightutrX1 - $PGSrightutrX2", "$error_color_code"), "\n";
				last RIGHTUTR;
			}
			$PGSrightutrX1+=$PGSx0;
			$PGSrightutrX2+=$PGSx0;
			$PGSrightutrY1=$geneheight * $PGSwhichgene + $PGSy0 + $exonheight / 2;
			$PGSrightutrY2=$PGSrightutrY1;
#			print $PGSsubinfo, "Test: RightUTR: X1:Y1:X2:Y2 $PGSrightutrX1:$PGSrightutrY1:$PGSrightutrX2:$PGSrightutrY2\n"; ### For test ###
			$PGSsvg->line(x1 => $PGSrightutrX1, y1 => $PGSrightutrY1, x2 => $PGSrightutrX2, y2 => $PGSrightutrY2, stroke=>'red',"stroke-width"=>5);
			$PGStextheight=$PGSrightutrY1+ $exonheight / 2 + $textmargin;
			$PGSsvg->text(x => $PGSrightutrX2, y => $PGStextheight, width => 10, height => 10, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=>$textfontsize, "-cdata" => "${$PGShash}{$PGSindgrp}{'rightutr'}[0][1]");
		}
		else {
#			print $PGSsubinfo, "Test: mRNA $PGSindgrp has no RightUTR\n"; ### For test ###
		}}###RIGHTUTR
		
		
		$PGSwhichgene++;
	}
	my $PGSout = $PGSsvg->xmlify;
	open SVGFILE, ">temp.svg";
	print SVGFILE $PGSout;
	close SVGFILE;
	return 1;
}



sub PrintAnnotation {
	my $PShash=shift;
	
	my $PSsubinfo='SUB(PrintAnnotation)';
	
	foreach my $PSmrnaid (sort keys %{$PShash}) {
		print "\n\n\n", '######  ', $PSmrnaid, ":\n";
		print "\tStrand: ", $strand, "\n";
		my $PStest_start=0;
		my $PStest_end=0;
		my $PStest_exon=0;
		if (exists ${$PShash}{$PSmrnaid}{'leftutr'}) {
			print "\tLeftUTR:  ${$PShash}{$PSmrnaid}{'leftutr'}[0][0] - ${$PShash}{$PSmrnaid}{'leftutr'}[0][1]\n";
		}
		if (exists ${$PShash}{$PSmrnaid}{'rightutr'}) {
			print "\tRightUTR: ${$PShash}{$PSmrnaid}{'rightutr'}[0][0] - ${$PShash}{$PSmrnaid}{'rightutr'}[0][1]\n";
		}
		if (exists ${$PShash}{$PSmrnaid}{'startcodon'} and defined ${$PShash}{$PSmrnaid}{'startcodon'}[0]) {
			print "\t", colored("StartCodon: ${$PShash}{$PSmrnaid}{'startcodon'}[0] - ${$PShash}{$PSmrnaid}{'startcodon'}[1]", "$error_color_code"), "\n";
			$PStest_start=1;
		}
		if (exists ${$PShash}{$PSmrnaid}{'stopcodon'} and defined ${$PShash}{$PSmrnaid}{'stopcodon'}[0]) {
			print "\t", colored("StopCodon:  ${$PShash}{$PSmrnaid}{'stopcodon'}[0] - ${$PShash}{$PSmrnaid}{'stopcodon'}[1]", "$error_color_code"), "\n";
			$PStest_end=1;
		}
		if (exists ${$PShash}{$PSmrnaid}{'exons'}) {
			$PStest_exon=1 if (scalar(@{${$PShash}{$PSmrnaid}{'exons'}})>0);
			print "\tExons (", scalar(@{${$PShash}{$PSmrnaid}{'exons'}}), ") : ";
			foreach (@{${$PShash}{$PSmrnaid}{'exons'}}) {
				print "$_->[0] - $_->[1], ";
			}
			print "\n";
		}
		
		my $PStestphase=0;
		my $PScdsarr=[];
		my $PSphase=[];
		if ($PStest_exon==1) {
			if ($PStest_start==1 and $PStest_end==1) {
				if ($strand eq '+') {
					print "\tGetting CDS: Start: ${$PShash}{$PSmrnaid}{'startcodon'}[0], End: ${$PShash}{$PSmrnaid}{'stopcodon'}[1], Strand: $strand\n";
					$PScdsarr=&GetCDS(\@{${$PShash}{$PSmrnaid}{'exons'}}, ${$PShash}{$PSmrnaid}{'startcodon'}[0], ${$PShash}{$PSmrnaid}{'stopcodon'}[1]);
				}
				elsif ($strand eq '-') {
					print "Getting CDS: Start: ${$PShash}{$PSmrnaid}{'stopcodon'}[0], End: ${$PShash}{$PSmrnaid}{'startcodon'}[1], Strand: $strand\n";
					$PScdsarr=&GetCDS(\@{${$PShash}{$PSmrnaid}{'exons'}}, ${$PShash}{$PSmrnaid}{'stopcodon'}[0], ${$PShash}{$PSmrnaid}{'startcodon'}[1]);
				}
				if (scalar(@{$PScdsarr})>0) {
					@{${$PShash}{$PSmrnaid}{'cds'}}=@{$PScdsarr};
					($PStestphase, $PSphase)=&GetCdsPhase(\@{$PScdsarr}, $PSmrnaid, 5);
					if ($PStestphase>0 and scalar(@{$PSphase}) == scalar(@{$PScdsarr})) {
						@{${$PShash}{$PSmrnaid}{'cdsphase'}}=@{$PSphase};
						$PStestphase=1;
					}
				}
			}
			elsif ($PStest_start==1) {
				if ($strand eq '+') {
					print "\tGetting CDS: Start: ${$PShash}{$PSmrnaid}{'startcodon'}[0], End: NONE, Strand: $strand\n";
					$PScdsarr=&GetCDS(\@{${$PShash}{$PSmrnaid}{'exons'}}, ${$PShash}{$PSmrnaid}{'startcodon'}[0], ' ');
				}
				elsif ($strand eq '-') {
					print "\tGetting CDS: Start: NONE, End: ${$PShash}{$PSmrnaid}{'startcodon'}[1], Strand: $strand\n";
					$PScdsarr=&GetCDS(\@{${$PShash}{$PSmrnaid}{'exons'}}, ' ', ${$PShash}{$PSmrnaid}{'startcodon'}[1]);
				}
				if (scalar(@{$PScdsarr})>0) {
					@{${$PShash}{$PSmrnaid}{'cds'}}=@{$PScdsarr};
					($PStestphase, $PSphase)=&GetCdsPhase(\@{$PScdsarr}, $PSmrnaid, 5);
					if ($PStestphase>0 and scalar(@{$PSphase}) == scalar(@{$PScdsarr})) {
						@{${$PShash}{$PSmrnaid}{'cdsphase'}}=@{$PSphase};
						$PStestphase=1;
					}
				}
			}
			elsif ($PStest_end==1) {
				if ($strand eq '+') {
					print "\tGetting CDS: Start: NONE, End: ${$PShash}{$PSmrnaid}{'stopcodon'}[1], Strand: $strand\n";
					$PScdsarr=&GetCDS(\@{${$PShash}{$PSmrnaid}{'exons'}}, ' ', ${$PShash}{$PSmrnaid}{'stopcodon'}[1]);
				}
				elsif ($strand eq '-') {
					print "\tGetting CDS: Start: ${$PShash}{$PSmrnaid}{'stopcodon'}[0], End: NONE, Strand: $strand\n";
					$PScdsarr=&GetCDS(\@{${$PShash}{$PSmrnaid}{'exons'}}, ${$PShash}{$PSmrnaid}{'stopcodon'}[0], ' ');
				}
				if (scalar(@{$PScdsarr})>0) {
					@{${$PShash}{$PSmrnaid}{'cds'}}=@{$PScdsarr};
					($PStestphase, $PSphase)=&GetCdsPhase(\@{$PScdsarr}, $PSmrnaid, 3);
					if ($PStestphase>0 and scalar(@{$PSphase}) == scalar(@{$PScdsarr})) {
						@{${$PShash}{$PSmrnaid}{'cdsphase'}}=@{$PSphase};
						$PStestphase=1;
					}
				}
			}
			else {
				my $PSestseq=&GetEst(\@{${$PShash}{$PSmrnaid}{'exons'}});
				print "mRNA $PSmrnaid EST: \n$PSestseq\n";
			}
		}
#		print $PSsubinfo, "Test: CDS\n"; print Dumper $PScdsarr; print "\n"; ### For test ###
#		print $PSsubinfo, "Test: Phase\n"; print Dumper $PSphase; print "\n"; ### For test ###
#		print $PSsubinfo, "Test: CDS2\n"; print Dumper ${$PShash}{$PSmrnaid}{'cds'}; print "\n"; ### For test ###
#		print $PSsubinfo, "Test: Phase2\n"; print Dumper ${$PShash}{$PSmrnaid}{'cdsphase'}; print "\n"; ### For test ###
		if ($PStestphase) {
			print "\tCDS: ";
			for (my $PSi=0; $PSi<scalar(@{$PScdsarr}); $PSi++) {
				print "(${$PShash}{$PSmrnaid}{'cdsphase'}[$PSi]) ${$PShash}{$PSmrnaid}{'cds'}[$PSi][0] - ${$PShash}{$PSmrnaid}{'cds'}[$PSi][1], ";
			}
			print "\n";
			print "Translation: ";
			my $PSaa='';
			if ($PStest_start==1) {
				$PSaa=&CdsTranslation(\@{${$PShash}{$PSmrnaid}{'cds'}}, $PSmrnaid, 5);
			}
			elsif ($PStest_end==1) {
				$PSaa=&CdsTranslation(\@{${$PShash}{$PSmrnaid}{'cds'}}, $PSmrnaid, 3);
			}
			print $PSaa."\n";
		}
	}
}



### GetCDS from exons array
### GetCDS
### Global: 
### Dependency: 
### Note: 
sub GetCDS {
	my ($GCexonarr, $GCstart, $GCend)=@_;
	
	my $GCsubinfo='SUB(GetCDS)';
	my @GCCDSarr=();
	
	foreach my $GCexon (@{$GCexonarr}) {
		if (defined $GCstart and $GCstart=~/^\d+$/ and defined $GCend and $GCend=~/^\d+$/) {
			if ($GCexon->[0]<=$GCstart and $GCstart <= $GCexon->[1]) {
				if ($GCend<=$GCexon->[1]) {
					push (@GCCDSarr, [$GCstart, $GCend]);
					last;
				}
				elsif ($GCexon->[1] < $GCend) {
					push (@GCCDSarr, [$GCstart, $GCexon->[1]]);
				}
			}
			elsif ($GCexon->[0]>=$GCstart) {
				if ($GCend<=$GCexon->[1]) {
					push (@GCCDSarr, [$GCexon->[0], $GCend]);
					last;
				}
				else {
					push (@GCCDSarr, [$GCexon->[0], $GCexon->[1]]);
				}
			}
		}
		elsif (defined $GCend and $GCend=~/^\d+$/) {
			if ($GCexon->[1]<=$GCend) {
				push (@GCCDSarr, [$GCexon->[0], $GCexon->[1]]);
			}
			elsif ($GCexon->[0]<$GCend and $GCexon->[1]>=$GCend) {
				push (@GCCDSarr, [$GCexon->[0], $GCend]);
			}
		}
		elsif (defined $GCstart and $GCstart=~/^\d+$/) {
			if ($GCexon->[0]<=$GCstart and $GCexon->[1]>=$GCstart) {
				push (@GCCDSarr, [$GCstart, $GCexon->[1]]);
			}
			elsif ($GCexon->[0]>=$GCstart) {
				push (@GCCDSarr, [$GCexon->[0], $GCexon->[1]]);
			}
		}
	}
	return \@GCCDSarr;
}



### GetCDS phase
### GetCdsPhase($CDSarr, 5/3)
### Global: $strand
### Dependency:
### Note: 
sub GetCdsPhase {
	my ($GCPcdsarr, $GCPmrnaid, $GCPwhichend)=@_;
	
	my $GCPsubinfo='SUB(GetCdsPhase)';
	my @GCPphase=();
	$GCPwhichend=5 unless (defined $GCPwhichend and $GCPwhichend==3);
	$GCPmrnaid='mRNA' unless (defined $GCPmrnaid);
	
#	print Dumper $GCPcdsarr; ### For test ###
	my $phase=0;
	if ($strand eq '+') {
		if ($GCPwhichend==5) {
			foreach my $GCPindcds (@{$GCPcdsarr}) {
				push (@GCPphase,  $phase);
				$phase=3-($GCPindcds->[1]-$GCPindcds->[0]+1-$phase)%3;
				$phase=0 if ($phase==3);
			}
		}
		elsif ($GCPwhichend==3) {
			for (my $GCPi=(scalar(@{$GCPcdsarr})-1); $GCPi>=0; $GCPi--) {
	#			print "Start: ${$GCPcdsarr}[$GCPi][0]\n"; ### For test ###
	#			print "End: ${$GCPcdsarr}[$GCPi][1]\n"; ### For test ###
				$phase=(${$GCPcdsarr}[$GCPi][1]-${$GCPcdsarr}[$GCPi][0] + 1 - (3-$phase) )%3;
				$phase=0 if ($phase==3);
				unshift (@GCPphase,  $phase);
			}
			
		}
		$phase=$GCPphase[0];
	}
	elsif ($strand eq '-') {
		if ($GCPwhichend==5) {
			for (my $GCPi=(scalar(@{$GCPcdsarr})-1); $GCPi>=0; $GCPi--) {
				unshift (@GCPphase,  $phase);
	#			print "Start: ${$GCPcdsarr}[$GCPi][0]\n"; ### For test ###
	#			print "End: ${$GCPcdsarr}[$GCPi][1]\n"; ### For test ###
				$phase=3- ((${$GCPcdsarr}[$GCPi][1]-${$GCPcdsarr}[$GCPi][0] + 1 - $phase)%3);
				$phase=0 if ($phase==3);
			}
		}
		elsif ($GCPwhichend==3) {
			foreach my $GCPindcds (@{$GCPcdsarr}) {
				$phase=3-($GCPindcds->[1]-$GCPindcds->[0]+1-$phase)%3;
				$phase=0 if ($phase==3);
				push (@GCPphase,  $phase);
			}
		}
		$phase=$GCPphase[-1];
	}
	
	if ($phase != 0) {
		print STDERR "\n", $GCPsubinfo, colored("Warnings: mRNA ( ID: $GCPmrnaid ) might be pseudogene (NOT 3*)", "$error_color_code"), "\n";
#		print "($GCPsubinfo)Info: strand: $strand; phase: @GCPphase\n"; ### For test ###
		return ($phase, \@GCPphase);
	}
	elsif (scalar(@GCPphase)>0) {
#		print "($GCPsubinfo)Info: strand: $strand; phase: @GCPphase\n"; ### For test ###
		return (3, \@GCPphase);
	}
	else {
		return 0;
	}
}



### CDS translation
### &CdsTranslation($CDSarr, 5/3)
### Global: 
### Dependency: 
### Note: 
sub CdsTranslation {
	my ($CTcdsarr, $CTmrnaid, $CTwhichend)=@_;
	
	my $CTsubinfo='SUB(CdsTranslation)';
	my $CTaaseq='';
	$CTwhichend=5 unless (defined $CTwhichend and $CTwhichend==3);
	my $CTsuffix='';
	my $CTlastcodonpart='';
	my $CTnumstart=0;
	my $CTnumstop=0;
	
	my ($CTtestphase, $CTphasearr)=&GetCdsPhase($CTcdsarr, $CTmrnaid, $CTwhichend);
#	print "CDS phase\n"; print Dumper $CTphasearr; print "\n"; ### For test ###
	if ($CTtestphase==3) {
		$CTsuffix='-----Intact-----';
	}
	elsif ($CTtestphase==0) {
		$CTsuffix='-----Error-----';
	}
	else {
		$CTsuffix="-----Pseudogene$CTtestphase-----";
	}
	
	if ($strand eq '+') {
		if ($CTwhichend==5) {
			for (my $CTi=0; $CTi<scalar(@{$CTcdsarr}); $CTi++) {
				$CTaaseq.="---Exon".($CTi+1)."Phase";
				$CTaaseq.= ($CTtestphase>0) ? ${$CTphasearr}[$CTi] : '?';
				$CTaaseq.='---';
				my $CTj=${$CTcdsarr}[$CTi][0];
				while ($CTj<=${$CTcdsarr}[$CTi][1]) {
					$CTlastcodonpart='' if (length($CTlastcodonpart)==3);
					my $CTcodonend=$CTj+(3-length($CTlastcodonpart)-1);
					if ($CTcodonend>${$CTcdsarr}[$CTi][1]) {
						$CTcodonend=${$CTcdsarr}[$CTi][1];
					}
					$CTlastcodonpart.=$fastaDB->seq($seqid, $CTj => $CTcodonend);
					if (length($CTlastcodonpart)==3) {
						my $CTthisaa=&Codon2AA($CTlastcodonpart);
	###					print "Codon: $CTlastcodonpart\tAA: $CTthisaa\n"; ### For test ###
						$CTaaseq.=$CTthisaa;
						if ($CTthisaa eq '*') {
							$CTnumstop++;
							$CTaaseq.="[$CTj-$CTcodonend]";
						}
						elsif ($CTthisaa =~/^m$/i) {
							$CTaaseq.="[$CTj-$CTcodonend]";
							$CTnumstart++;
						}
						$CTlastcodonpart='';
					}
					$CTj=$CTcodonend+1;
				}
			}
		}
		elsif ($CTwhichend==3) {
			my $CTcodonend=0;
			my $CTsodonstarted=0;
			for (my $CTi=(scalar(@{$CTcdsarr})-1); $CTi>=0; $CTi--) {
#				print "Exon".($CTi+1).": ${$CTcdsarr}[$CTi][0] - ${$CTcdsarr}[$CTi][1]\n"; ### For test ###
				my $CTj=${$CTcdsarr}[$CTi][1];
				while ($CTj>=${$CTcdsarr}[$CTi][0]) {
					$CTlastcodonpart='' if (length($CTlastcodonpart)==3);
					my $CTcodonstart=$CTj-(3-length($CTlastcodonpart)-1);
					if ($CTcodonstart<${$CTcdsarr}[$CTi][0]) {
						$CTcodonstart=${$CTcdsarr}[$CTi][0];
					}
					$CTlastcodonpart=$fastaDB->seq($seqid, $CTcodonstart => $CTj).$CTlastcodonpart;
					if (length($CTlastcodonpart)==3) {
						my $CTthisaa=&Codon2AA($CTlastcodonpart);
	#					print "Codon: $CTlastcodonpart\tAA: $CTthisaa\n"; ### For test ###
						if ($CTsodonstarted==0) {
							$CTcodonend=$CTj;
						}
						if ($CTthisaa eq '*') {
							$CTnumstop++;
							$CTaaseq="[$CTcodonstart-$CTcodonend]".$CTaaseq;
						}
						elsif ($CTthisaa =~/^m$/i) {
							$CTnumstart++;
							$CTaaseq="[$CTcodonstart-$CTcodonend]".$CTaaseq;
						}
						$CTaaseq=$CTthisaa.$CTaaseq;
						$CTlastcodonpart='';
						$CTsodonstarted=0;
					}
					else {
						if ($CTsodonstarted==0) {
							$CTcodonend=$CTj;
							$CTsodonstarted=1;
						}
					}
					$CTj=$CTcodonstart-1;
				}
				$CTaaseq='---'.$CTaaseq;
				$CTaaseq= (($CTtestphase>0) ? ${$CTphasearr}[$CTi] : '?').$CTaaseq;
				$CTaaseq="---Exon".($CTi+1)."Phase".$CTaaseq;
			}
		}
	}
	elsif ($strand eq '-') {
		if ($CTwhichend==5) {
			my $CTcodonend=0;
			my $CTsodonstarted=0;
			for (my $CTi=(scalar(@{$CTcdsarr})-1); $CTi>=0; $CTi--) {
#				print "Exon".($CTi+1).": ${$CTcdsarr}[$CTi][0] - ${$CTcdsarr}[$CTi][1]\n"; ### For test ###
				my $CTj=${$CTcdsarr}[$CTi][1];
				$CTaaseq.="---Exon".($CTi+1)."Phase";
				$CTaaseq.= ($CTtestphase>0) ? ${$CTphasearr}[$CTi] : '?';
				$CTaaseq.='---';
				
				while ($CTj>=${$CTcdsarr}[$CTi][0]) {
					$CTlastcodonpart='' if (length($CTlastcodonpart)==3);
					my $CTcodonstart=$CTj-(3-length($CTlastcodonpart)-1);
					if ($CTcodonstart<${$CTcdsarr}[$CTi][0]) {
						$CTcodonstart=${$CTcdsarr}[$CTi][0];
					}
					$CTlastcodonpart=$fastaDB->seq($seqid, $CTcodonstart => $CTj).$CTlastcodonpart;
					if (length($CTlastcodonpart)==3) {
						$CTlastcodonpart=SeqRevComp($CTlastcodonpart);
						my $CTthisaa=&Codon2AA($CTlastcodonpart);
	#					print "Codon: $CTlastcodonpart\tAA: $CTthisaa\n"; ### For test ###
						if ($CTsodonstarted==0) {
							$CTcodonend=$CTj;
						}
						$CTaaseq.=$CTthisaa;
						if ($CTthisaa eq '*') {
							$CTnumstop++;
							$CTaaseq.="[$CTcodonstart-$CTcodonend]";
						}
						elsif ($CTthisaa =~/^m$/i) {
							$CTnumstart++;
							$CTaaseq.="[$CTcodonstart-$CTcodonend]";
						}
						
						$CTlastcodonpart='';
						$CTsodonstarted=0;
					}
					else {
						if ($CTsodonstarted==0) {
							$CTcodonend=$CTj;
							$CTsodonstarted=1;
						}
					}
					$CTj=$CTcodonstart-1;
				}
			}
		}
		elsif ($CTwhichend==3) {
			for (my $CTi=0; $CTi<scalar(@{$CTcdsarr}); $CTi++) {
				my $CTj=${$CTcdsarr}[$CTi][0];
				while ($CTj<=${$CTcdsarr}[$CTi][1]) {
					$CTlastcodonpart='' if (length($CTlastcodonpart)==3);
					my $CTcodonend=$CTj+(3-length($CTlastcodonpart)-1);
					if ($CTcodonend>${$CTcdsarr}[$CTi][1]) {
						$CTcodonend=${$CTcdsarr}[$CTi][1];
					}
					$CTlastcodonpart.=$fastaDB->seq($seqid, $CTj => $CTcodonend);
					if (length($CTlastcodonpart)==3) {
						$CTlastcodonpart=SeqRevComp($CTlastcodonpart);
						my $CTthisaa=&Codon2AA($CTlastcodonpart);
	###					print "Codon: $CTlastcodonpart\tAA: $CTthisaa\n"; ### For test ###
						
						if ($CTthisaa eq '*') {
							$CTnumstop++;
							$CTaaseq="[$CTj-$CTcodonend]".$CTaaseq;
						}
						elsif ($CTthisaa =~/^m$/i) {
							$CTaaseq="[$CTj-$CTcodonend]".$CTaaseq;
							$CTnumstart++;
						}
						
						$CTaaseq=$CTthisaa.$CTaaseq;
						$CTlastcodonpart='';
					}
					$CTj=$CTcodonend+1;
				}
				$CTaaseq='---'.$CTaaseq;
				$CTaaseq= (($CTtestphase>0) ? ${$CTphasearr}[$CTi] : '?').$CTaaseq;
				$CTaaseq="---Exon".($CTi+1)."Phase".$CTaaseq;
			}
		}
		
	
	}
	$CTaaseq.=$CTsuffix;
	$CTaaseq.="start $CTnumstart stop $CTnumstop";
	return $CTaaseq;
}



### Edit each feature for a specified mRNA
### EditAnnotation(\%annot)
### Global: $standinput
### Dependency:
### Note:
sub EditAnnotation {
	my $EAannot=shift;
	my $EAexit=0;
	
	my $EAsubinfo='SUB(EditAnnotation)';
	my $EAedit=dclone $EAannot;
	my %EAexcludekeys=('cds' => 1, 'cdsphase' => 1, 'leftutr' => 1, 'rightutr' => 1, 'utr3' => 1, 'utr5' => 1);
	
	while ($EAexit==0) {
		&OutAnnot($EAedit);
		my @EAarrkey=();
		my $EAitems2edit;
		my %EAitemshash=();
		$EAitemshash{'q'}='QUIT';
		$EAitemshash{'s'}='strand';
		my $EAi=0;
		print "\n\n\n", $EAsubinfo, "Info: select items to edit (CDS feature auto-created): \n";
		foreach (sort keys %{$EAedit}) {
			unless (exists $EAexcludekeys{$_}) {
				print "$EAi:\t$_\n";
				$EAitemshash{$EAi++}=$_;
			}
		}
		print "s:\tNEW\n";
		$EAitemshash{'s'}='NEW';
		
		print $EAsubinfo, "Input: select feature to edit or [q]uit : ";
		$standinput=&GetStandin;
		if (exists $EAitemshash{$standinput}) {
			print $EAsubinfo, "Info: Accepted to edit: $standinput: $EAitemshash{$standinput}\n";
		}
		elsif ($standinput eq '') {
			print STDERR $EAsubinfo, colored("Error: NO feature selected: return ...", "$error_color_code"), "\n";
			$EAexit=1;
			next;
		}
		else {
			print STDERR $EAsubinfo, colored("Error: index not existed: $standinput", "$error_color_code"), "\n";
			next;
		}
		if ($EAitemshash{$standinput} eq 'leftutr') {
			print $EAsubinfo, "Info: LeftUTR: ${$EAedit}{'leftutr'}[0][0] - ${$EAedit}{'leftutr'}[0][1]\n";
			print $EAsubinfo, "Select: [e]dit, [d]elete, or [q]uit? : ";
			$standinput=&GetStandin;
			if ($standinput =~/(^e$)|(^edit$)/i) {
				my $EAinput=&StandInput;
				${$EAedit}{'leftutr'}[0][0]=${$EAinput}[0];
				${$EAedit}{'leftutr'}[0][1]=${$EAinput}[1];
			}
			elsif ($standinput =~/(^d$)|(^delete$)/i) {
				delete ${$EAedit}{'leftutr'};
			}
			elsif ($standinput =~/(^q$)|(^quit$)/i) {
				next;
			}
		}
		elsif ($EAitemshash{$standinput} eq 'rightutr') {
			print $EAsubinfo, "Info: RightUTR: ${$EAedit}{'rightutr'}[0][0] - ${$EAedit}{'rightutr'}[0][1]\n";
			print $EAsubinfo, "Select: [e]dit, [d]elete, or [q]uit? : ";
			$standinput=&GetStandin;
			if ($standinput =~/(^e$)|(^edit$)/i) {
				my $EAinput=&StandInput;
				${$EAedit}{'rightutr'}[0][0]=${$EAinput}[0];
				${$EAedit}{'rightutr'}[0][1]=${$EAinput}[1];
			}
			elsif ($standinput =~/(^d$)|(^delete$)/i) {
				delete ${$EAedit}{'rightutr'};
			}
			elsif ($standinput =~/(^q$)|(^quit$)/i) {
				next;
			}
		}
		elsif ($EAitemshash{$standinput} eq 'startcodon') {
			print $EAsubinfo, "Info: StartCodon: ${$EAedit}{'startcodon'}[0] - ${$EAedit}{'startcodon'}[1]\n";
			print $EAsubinfo, "Input: [e]dit, [d]elete, or [q]uit? : ";
			$standinput=&GetStandin;
			if ($standinput =~/(^e$)|(^edit$)/i) {
				my $EAteststart=0;
				while ($EAteststart==0) {
					print $EAsubinfo, "Input: NEW startcodon (NOW ${$EAedit}{'startcodon'}[0] - ${$EAedit}{'startcodon'}[1]) [INT] | [q]uit : ";
					my $EAinput=&GetStandin;
					if ($EAinput=~/^\d+$/) {
						${$EAedit}{'startcodon'}[0]=$EAinput;
						${$EAedit}{'startcodon'}[1]=$EAinput+2;
						$EAteststart=1;
						print "Info: NEW startcodon defined: ${$EAedit}{'startcodon'}[0] - ${$EAedit}{'startcodon'}[1]\n";
					}
					elsif ($EAinput =~ /^(q)|(quit)$/i) {
						$EAteststart=1;
					}
					else {
						print "Try again ...\n";
					}
				}
			}
			elsif ($standinput =~/(^d$)|(^delete$)/i) {
				delete ${$EAedit}{'startcodon'};
			}
			elsif ($standinput =~/(^q$)|(^quit$)/i) {
				next;
			}
		}
		elsif ($EAitemshash{$standinput} eq 'stopcodon') {
			print $EAsubinfo, "Info: StopCodon: ${$EAedit}{'stopcodon'}[0] - ${$EAedit}{'stopcodon'}[1]\n";
			print $EAsubinfo, "Input: [e]dit, [d]elete, or [q]uit? : ";
			$standinput=&GetStandin;
			if ($standinput =~/(^e$)|(^edit$)/i) {
				my $EAtestend=0;
				while ($EAtestend==0) {
					print $EAsubinfo, "Input: NEW stopcodon (NOW ${$EAedit}{'stopcodon'}[0] - ${$EAedit}{'stopcodon'}[1]) [INT] : ";
					my $EAinput=&GetStandin;
					if ($EAinput=~/^\d+$/) {
						${$EAedit}{'stopcodon'}[0]=$EAinput;
						${$EAedit}{'stopcodon'}[1]=$EAinput+2;
						$EAtestend=1;
						print "Info: NEW stopcodon defined: ${$EAedit}{'stopcodon'}[0] - ${$EAedit}{'stopcodon'}[1]\n";
					}
					elsif ($EAinput =~ /^(q)|(quit)$/i) {
						$EAtestend=1;
					}
					else {
						print "Try again ...\n";
					}
				}
			}
			elsif ($standinput =~/(^d$)|(^delete$)/i) {
				delete ${$EAedit}{'stopcodon'};
			}
			elsif ($standinput =~/(^q$)|(^quit$)/i) {
				next;
			}
		}
		elsif ($EAitemshash{$standinput} eq 'exons') {
			print $EAsubinfo, "Info: Exons: \n\n";
			my $EAtestexon=0;
			my $EAexons=dclone \@{${$EAedit}{'exons'}};
			while ($EAtestexon==0) {
				my %EAexonhash=();
				my $EAj=0;
				foreach (@{$EAexons}) {
					print "Exon: \t$EAj.\t$_->[0]-$_->[1]\n";
					$EAexonhash{$EAj++}="$_->[0]-$_->[1]";
				}
				print "\n\n", "$EAsubinfo: MIN: $minpos MAX: $maxpos\n";
				print $EAsubinfo, "Select exon or [i]nsert r [q]uit : ";
				my $EAexon2edit=&GetStandin;
				if ($EAexon2edit eq '' or $EAexon2edit =~/(^q$)|(^quit$)/i) {
					$EAtestexon=1;
					next;
				}
				
				elsif ($EAexon2edit =~/(^i$)|(^insert$)/i) {
					print "\nInfo: NEW exon, Now input new [INT-INT] : ";
					my $EAinput=&StandInput;
					$EAexonhash{$EAj++}="${$EAinput}[0]-${$EAinput}[1]";
				}
				elsif ($EAexon2edit =~ /^\d+$/ and exists $EAexonhash{$EAexon2edit}) {
					print $EAsubinfo, "Select: [e]dit, [i]nsert, [d]elete, or [q]uit? : ";
					$standinput=&GetStandin;
					$standinput='q' if ($standinput eq '');
					if ($standinput =~/(^e$)|(^edit$)/i) {
						print "\nInfo: OLD exon ($EAexonhash{$EAexon2edit}), Now input new [INT-INT] : ";
						my $EAinput=&StandInput;
						$EAexonhash{$EAexon2edit}="${$EAinput}[0]-${$EAinput}[1]";
					}
					elsif ($standinput =~/(^d$)|(^delete$)/i) {
						delete $EAexonhash{$EAexon2edit};
					}
					elsif ($standinput =~/(^i$)|(^insert$)/i) {
						my $EAinput=&StandInput;
						$EAexonhash{$EAj++}="${$EAinput}[0]-${$EAinput}[1]";
					}
					elsif ($standinput =~/(^q$)|(^quit$)/i) {
						$EAtestexon=1;
						next;
					}
				}
				else {
					print STDERR $EAsubinfo, colored("Warnings: invalid selection: $EAexon2edit", "$error_color_code"), "\n";
					next;
				}
				$EAexons=SortExons(\%EAexonhash);
			}
			${$EAedit}{'exons'}=$EAexons;
		}
		elsif ($EAitemshash{$standinput} eq 'min') {
			print $EAsubinfo, "Info: MIN: ${$EAedit}{'min'}, Input NEW number of [q]uit: ";
			$standinput=&GetStandin;
			if ($standinput =~/^\d+$/i) {
				${$EAedit}{'min'}=$standinput;
				$minpos=$standinput if ($standinput< $minpos);
			}
			else {
				next;
			}
		}
		elsif ($EAitemshash{$standinput} eq 'max') {
			print $EAsubinfo, "Info: MAX: ${$EAedit}{'max'}, Input NEW number of [q]uit: ";
			$standinput=&GetStandin;
			if ($standinput =~/^\d+$/i) {
				${$EAedit}{'max'}=$standinput;
				$maxpos=$standinput if ($standinput>$maxpos);
			}
			else {
				next;
			}
		}
		elsif ($EAitemshash{$standinput} eq 'strand') {
			print $EAsubinfo, "Input: strand (current: $strand) or [q]uit: ";
			$standinput=&GetStandin;
			if ($standinput =~/^(\+)|(-)$/) {
				$strand=$standinput;
				print $EAsubinfo, "Info: strand accepted: $strand\n";
			}
			else {
				next;
			}
		}
		elsif ($EAitemshash{$standinput} eq 'NEW') {
			my %EAaddfeature=(1=>'startcodon', 2=>'stopcodon');
			print $EAsubinfo, "Info: available feature: 1. startcodon, 2. stopcodon\n";
			print $EAsubinfo, "Input: feature key: ";
			my $EAnewfeaturenum=&GetStandin;
			unless (exists $EAaddfeature{$EAnewfeaturenum}){
				print STDERR $EAsubinfo, colored("Error: feature key NOT existed: $EAnewfeaturenum", "$error_color_code"), "\n";
				print STDERR $EAsubinfo, colored("Error: start again ...", "$error_color_code"), "\n";
				next;
			}
			my $EAnewfeature=$EAaddfeature{$EAnewfeaturenum};
			unless (exists ${$EAedit}{$EAnewfeature}) {
				if ($EAnewfeature eq 'startcodon' or $EAnewfeature eq 'stopcodon') {
					my $EAtestfeature=0;
					while ($EAtestfeature==0) {
						print $EAsubinfo, "Input: feature value [INT]|[q]uit: ";
						my $EAfeaturevalue=&GetStandin;
						if ($EAfeaturevalue=~/^(\d+)$/) {
							${$EAedit}{$EAnewfeature}[0]=$1;
							${$EAedit}{$EAnewfeature}[1]=${$EAedit}{$EAnewfeature}[0] + 2;
							print $EAsubinfo, "Info: NEW feature $EAnewfeature: $1 - ".(${$EAedit}{$EAnewfeature}[0] + 2)." accepted\n";
							$EAtestfeature=1;
						}
						elsif ($EAfeaturevalue eq 'q' or $EAfeaturevalue eq 'quit') {
							$EAtestfeature=1;
						}
					}
				}
			}
			else {
				print STDERR $EAsubinfo, colored("Error: feature key existed: $EAnewfeature", "$error_color_code"), "\n";
			}
		}
		elsif ($EAitemshash{$standinput} eq 'QUIT') {
			$EAexit=1;
		}
		else {
			$EAexit=1;
		}
	}
	return $EAedit;
}



### Input start-end from STDIN
### StandInput()
### Global: 
### Dependency:
### Note: 
sub StandInput {
	
	my $SIsubinfo='SUB(StandInput)';
	my $SIexit=0;
	my @SIarr=();
	
	while ($SIexit==0) {
		print $SIsubinfo, "Input: NEW(start-end): ";
		my $SIstandin=&GetStandin;
		@SIarr=split(/-/, $SIstandin);
		unless ($SIarr[0]=~/^\d+$/ and $SIarr[1]=~/^\d+$/ and $SIarr[0] < $SIarr[1]) {
			print STDERR "\n\n\n$SIsubinfo: ", colored("invalid number: $SIarr[0] - $SIarr[1]", "$error_color_code"), "\n\n\n";
			next;
		}
		print $SIsubinfo, "Accepted: $SIarr[0] - $SIarr[1]\n";
		print $SIsubinfo, "Exit_this_feature_edit: [y]|n: ";
		$SIstandin=&GetStandin;
		if ($SIstandin eq '' or $SIstandin eq 'y') {
			$SIexit=1;
		}
	}
	return \@SIarr;
}



###

sub OutAnnot {
	my $OAannot=shift;
	
	print "\n\n\n";
	if (exists ${$OAannot}{'leftutr'}) {
		print "LeftUTR: ${$OAannot}{'leftutr'}[0][0] - ${$OAannot}{'leftutr'}[0][1]\n";
	}
	else {
		print "LeftUTR: [undefined]\n";
	}
	if (exists ${$OAannot}{'rightutr'}) {
		print "RightUTR: ${$OAannot}{'rightutr'}[0][0] - ${$OAannot}{'rightutr'}[0][1]\n";
	}
	else {
		print "RightUTR: [undefined]\n";
	}
	if (exists ${$OAannot}{'utr5'}) {
		print "5'UTR: ${$OAannot}{'utr5'}[0][0] - ${$OAannot}{'utr5'}[0][1]\n";
	}
	else {
		print "5'UTR: [undefined]\n";
	}
	if (exists ${$OAannot}{'utr3'}) {
		print "LeftUTR: ${$OAannot}{'utr3'}[0][0] - ${$OAannot}{'utr3'}[0][1]\n";
	}
	else {
		print "3'UTR: [undefined]\n";
	}
	if (exists ${$OAannot}{'startcodon'}) {
		print "StartCodon: ${$OAannot}{'startcodon'}[0] - ${$OAannot}{'startcodon'}[1]\n";
	}
	else {
		print "StartCodon: [undefined]\n";
	}
	if (exists ${$OAannot}{'stopcodon'}) {
		print "StopCodon: ${$OAannot}{'stopcodon'}[0] - ${$OAannot}{'stopcodon'}[1]\n";
	}
	else {
		print "StopCodon: [undefined]\n";
	}
	if (exists ${$OAannot}{'exons'}) {
		my $OAnum=0;
		foreach (@{${$OAannot}{'exons'}}) {
			$OAnum++;
			print "Exon($OAnum): $_->[0] - $_->[1]\n";
		}
	}
	else {
		print "Exons: [undefined]\n";
	}
}



### SortExons
### &SortExons()
### Global: 
### Dependency: 
### Note:
sub SortExons {
	my $SEexon=shift;
	
	my $SEsubinfo='SUB(SortExons)';
	my %SEtempexon=();
	my $SEretarr=[];
	my $SElastend=0;
	
	foreach (keys %{$SEexon}) {
		if (${$SEexon}{$_}=~/^(\d+)-(\d+)$/) {
			if ($2>$1) {
				$SEtempexon{$1}{$2}++;
			}
			else {
				print STDERR $SEsubinfo, colored("Error: invalid Exon ($_) End>=Start: ${$SEexon}{$_}", "$error_color_code"), "\n";
			}
		}
	}
	
	foreach my $SEstart (sort {$a <=>$b} keys %SEtempexon) {
		unless ($SEstart>$SElastend) {
			print STDERR $SEsubinfo, colored ("Error: Exon start < last end $SEstart <= $SElastend", "$error_color_code"), "\n";
			next;
		}
		my @SEendarr=(sort {$a <=>$b} keys %{$SEtempexon{$SEstart}});
		if (scalar(@SEendarr)==1) {
			push (@{$SEretarr}, [$SEstart, $SEendarr[0]]);
			$SElastend=$SEendarr[0];
		}
		elsif (scalar(@SEendarr)>1) {
			foreach (@SEendarr) {
				push (@{$SEretarr}, [$SEstart, $_]);
				$SElastend=$_;
			}
		}
		else {
			print STDERR $SEsubinfo, colored("Error: Exon start has no ends: $SEstart", "$error_color_code"), "\n";
		}
	}
	return $SEretarr;
}



### Write Annotation to file
### WriteAnnotation(\%provisionalmrna)
### Global: $seqid, $start, $end, $gff2col, %fastalength, $strand
### Dependency:
### Note: 
sub WriteAnnotation {
	my $WAmrna=shift;
	
	local *WAOUTGFF;
	my $WAsubinfo='SUB(WriteAnnotation)';
	my $WAgff3out=$seqid.'_'.$start.'-'.$end.'.gff3';
	my $WAgeneid=$seqid.'_'.$start.'-'.$end.'.gene';
	my $WAmrnaid='trans0001';
	my $WAexonid='exon0001';
	my $WAtest=0;
	my $WAgenescore;
	
	while ($WAtest==0) {
		print "\n\n\n", $WAsubinfo, "Input: Gene Score [INT]: ";
		$WAgenescore=&GetStandin;
		if ($WAgenescore=~/^\d+$/) {
			$WAtest=1;
		}
		else {
			print STDERR $WAsubinfo, colored("Error: invalid Gene Score number, try again...", "$error_color_code"), "\n";
		}
	}
	
	close WAOUTGFF if (defined fileno(WAOUTGFF));
	
	print "\n\n\n", $WAsubinfo, "Selection: notes: \n";
	foreach (sort {$a <=> $b} keys %genenotes) {
		print "\t$_\t$genenotes{$_}\n";
	}
	
	
	$WAtest=0;
	while ($WAtest==0) {
		print "\n\n\n", $WAsubinfo, "Input: [N/ENTER]: \n";
		my $WAinputnote=&GetStandin;
		if (exists $genenotes{$WAinputnote}) {
			$WAprintnotes=$genenotes{$WAinputnote};
			$WAtest=1;
		}
		elsif ($WAinputnote eq '') {
			$WAtest=1;
		}
	}
	
	unless (open(WAOUTGFF, "> $WAgff3out")) {
		 print STDERR $WAsubinfo, colored("Error: can not write GFF3 file: $WAgff3out", "$error_color_code"), "\n";
		 return 0;
	}
	print WAOUTGFF "##gff-version 3\n";
	print WAOUTGFF "##sequence-region   $seqid $minpos $maxpos\n";
	print WAOUTGFF $seqid,"\t", $gff2col, "\tgene\t", $minpos, "\t", $maxpos, "\t", $WAgenescore, "\t", $strand, "\t", '.', "\tID=$WAgeneid";
	if ($WAprintnotes ne '') {
		print WAOUTGFF ";Note=$WAprintnotes\n";
	}
	else {
		print WAOUTGFF "\n";
	}
	foreach my $WAindmrnaid (sort keys %{$WAmrna}) {
		$WAexonid='exon0001';
		my @WArrnanote=();
		if (exists ${$WAmrna}{$WAindmrnaid}{'startcodon'}) {
			push (@WArrnanote, 'StartCodon');
		}
		else {
			push (@WArrnanote, 'missingStartCodon');
		}
		if (exists ${$WAmrna}{$WAindmrnaid}{'utr5'}) {
			push (@WArrnanote, 'UTR5');
		}
		else {
			push (@WArrnanote, 'missingUTR5');
		}
		if (exists ${$WAmrna}{$WAindmrnaid}{'stopsodon'}) {
			push (@WArrnanote, 'StopCodon');
		}
		else {
			push (@WArrnanote, 'missingStopCodon');
		}
		if (exists ${$WAmrna}{$WAindmrnaid}{'utr3'}) {
			push (@WArrnanote, 'UTR3');
		}
		else {
			push (@WArrnanote, 'missingUTR3');
		}
		print WAOUTGFF $seqid,"\t", $gff2col, "\tmRNA\t", ${$WAmrna}{$WAindmrnaid}{'min'}, "\t", ${$WAmrna}{$WAindmrnaid}{'max'}, "\t", $WAgenescore, "\t", $strand, "\t", '.', "\tID=$WAgeneid.$WAmrnaid;Parent=$WAgeneid;Note=". join('.', @WArrnanote) ."\n";
		if (exists ${$WAmrna}{$WAindmrnaid}{'startcodon'}) {
			print WAOUTGFF $seqid,"\t", $gff2col, "\tstartcodon\t", ${$WAmrna}{$WAindmrnaid}{'startcodon'}[0], "\t", ${$WAmrna}{$WAindmrnaid}{'startcodon'}[1], "\t", $WAgenescore, "\t", $strand, "\t", '.', "\tID=$WAgeneid.$WAmrnaid.startcodon;Parent=$WAgeneid.$WAmrnaid\n";
		}
		if (exists ${$WAmrna}{$WAindmrnaid}{'stopcodon'}) {
			print WAOUTGFF $seqid,"\t", $gff2col, "\tstopcodon\t", ${$WAmrna}{$WAindmrnaid}{'stopcodon'}[0], "\t", ${$WAmrna}{$WAindmrnaid}{'stopcodon'}[1], "\t", $WAgenescore, "\t", $strand, "\t", '.', "\tID=$WAgeneid.$WAmrnaid.stopcodon;Parent=$WAgeneid.$WAmrnaid\n";
		}
		foreach (@{${$WAmrna}{$WAindmrnaid}{'exons'}}) {
			print WAOUTGFF $seqid,"\t", $gff2col, "\texon\t", $_->[0], "\t", $_->[1], "\t", $WAgenescore, "\t", $strand, "\t", '.', "\tID=$WAgeneid.$WAmrnaid.$WAexonid;Parent=$WAgeneid.$WAmrnaid\n";
			$WAexonid++;
		}
		if (exists ${$WAmrna}{$WAindmrnaid}{'cds'}) {
			for (my $WAi=0; $WAi<scalar(@{${$WAmrna}{$WAindmrnaid}{'cds'}}); $WAi++) {
				print WAOUTGFF $seqid,"\t", $gff2col, "\tCDS\t", ${$WAmrna}{$WAindmrnaid}{'cds'}[$WAi][0], "\t", ${$WAmrna}{$WAindmrnaid}{'cds'}[$WAi][1], "\t", $WAgenescore, "\t", $strand, "\t", ${$WAmrna}{$WAindmrnaid}{'cdsphase'}[$WAi], "\tID=$WAgeneid.$WAmrnaid.cds;Parent=$WAgeneid.$WAmrnaid\n";
				$WAexonid++;
			}
		}
		$WAmrnaid++;
	}
	close WAOUTGFF;
	return 1;
}



sub GetEst {
	my $GEexons=shift;
	
	my $GEsubinfo='SUB(GetEst)';
	my $GEestseq='';
	
	foreach (@{$GEexons}) {
#		print $GEsubinfo, "Test: EST from exon $_->[0] - $_->[1]\n"; ### For test ###
		my $GEstart=$_->[0];
		my $GEend=$_->[1];
#		print $GEsubinfo, "Test: EST from exon $GEstart - $GEend\n"; ### For test ###
		$GEestseq.= $fastaDB->seq($seqid, $GEstart => $GEend);
	}
	
	$GEestseq=SeqRevComp($GEestseq) if ($strand eq '-');
#	print $GEsubinfo, "Test: EST: $GEestseq\n";### For test ###
	return $GEestseq;
}



### input from STDIN
sub GetStandin {
	my $GSinput=<STDIN>;
	chomp $GSinput;
	return $GSinput;
}
