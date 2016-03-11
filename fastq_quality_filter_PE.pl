#!/usr/bin/perl

use strict;

# filters out paired-end reads and corresponding barcodes if AVERAGE Q score in any of the three is below the specified threshold

my $usage = "perl fastq_quality_filter_PE.pl FASTQ1 FASTQ2 BARCODE MIN-QUALITY MIN-PCT-BASES OUT-FASTQ1 OUT-FASTQ2 OUT-BARCODE\n";
my $fastq1_file = shift or die $usage;
my $fastq2_file = shift or die $usage;
my $barcode_file = shift or die $usage;
my $min_quality = shift or die $usage;
my $min_pct_bases = shift or die $usage;
my $fastq1_outfile = shift or die $usage;
my $fastq2_outfile = shift or die $usage;
my $barcode_outfile = shift or die $usage;

#@NS500273:4:H0NHFAGXX:1:11201:6165:1036 1:N:0:GTTGGTCAATAT
#NAAAAAGAGCGAGAAAAAAGTCAAAAAAAATAAAACGAAAAAGTCTTATAAATATCGTGTCGCATCTGCTGAGAAAATCCAGGCAATAAATTTGGGATAGCAGTGGGTTAAGGCATGCTAGAGTGCAGTAGGGGAGATTGGAATTCCTGGA
#+
##))<A<F.<FF.<<FFF<FFA7FF<AA77FFF7F.A<FFFFFF...F)A<A)F))FFFFFFF.))7)FFF.AF)F<FF<7FFF)F<7)F.FFFFFAAF<7FFFF<FFF7..F<<.FF.77FFFFFF)FAFF<AAFFF7A...F7<.AAFF.

#==> Haiti_Run1_Read2.fastq <==
#@NS500273:4:H0NHFAGXX:1:11201:6165:1036 2:N:0:GTTGGTCAATAT
#AAAAAATATATCAAGCAACAAAAAGAATTAACATAATACGCACTGAAACAAGATTACAAAAATATCAACGAAACACTTGAAAAAAATGAAAAACAAAAAAAAAACAGGAAAAAAAGGCTAACCTACTGCAAAAAAACATGCACAAACCCCC
#+
#<A7.7FFFFFF.)FFFA)F.<F.FFF7F<F7A))77)<..FFFFF<..F<F.).)A7.)<F)FFF7.<FF7FFA)7F.<FFFFAFFF.FFFAF)F<AAFFFAFAFFF.FFF)FFF7.)77)7)..AA)7F<FF<<F<7<<))..F..)<)F

open(FASTQ1, "<$fastq1_file") || die "Unable to open: $!";
open(FASTQ2, "<$fastq2_file") || die "Unable to open: $!";
open(BARCODE, "<$barcode_file") || die "Unable to open: $!";
open(FASTQ1_OUT, ">$fastq1_outfile") || die "Unable to open: $!";
open(FASTQ2_OUT, ">$fastq2_outfile") || die "Unable to open: $!";
open(BARCODE_OUT, ">$barcode_outfile") || die "Unable to open: $!";

my $head_f1;
my $head_f2;
my $head_b;
my $seq_f1;
my $seq_f2;
my $seq_b;
my $qual_f1;
my $qual_f2;
my $qual_b;

my $num_failed_reads;
my $num_total_reads;
my $num_failed_bases;
my $num_total_bases;

while (my $line = <FASTQ1>) {
	# get header lines
	chomp $line; $head_f1 = $line;
	$line = <FASTQ2>; chomp $line; $head_f2 = $line;
	$line = <BARCODE>; chomp $line; $head_b = $line;
	# get seq lines
	$line = <FASTQ1>; chomp $line; $seq_f1 = $line;
	$line = <FASTQ2>; chomp $line; $seq_f2 = $line;
	$line = <BARCODE>; chomp $line; $seq_b = $line;
	# skip a line
	<FASTQ1>; <FASTQ2>; <BARCODE>;
	# get quality lines
	$line = <FASTQ1>; chomp $line; $qual_f1 = $line;
	$line = <FASTQ2>; chomp $line; $qual_f2 = $line;
	$line = <BARCODE>; chomp $line; $qual_b = $line;
	
	# quality filtering
	my $pass = 1;
	my $q_sum = 0;
	my $num_above = 0;
	my $num_total = 0;
	foreach (split //, $qual_f1) {
		my $q = ord($_)-33;
		$q_sum += $q;
		$num_total++;
		if ($q < $min_quality) {
			$num_failed_bases++;
		}
		else {
			$num_above++;
		}
		$num_total_bases++;
	}
##	my $q_avg = ($q_sum / $num_total);
##	if ($q_avg < $min_quality) {
##		$pass = 0;
###		print STDERR "failed F1 with q_sum=$q_sum num_total=$num_total q_avg=$q_avg min_qual=$min_quality\n";
##	}
##	my $q_sum = 0;
##	my $num_total = 0;
	foreach (split //, $qual_f2) {
		my $q = ord($_)-33;
		$q_sum += $q;
		$num_total++;
		if ($q < $min_quality) {
			$num_failed_bases++;
		}
		else {
			$num_above++;
		}
		$num_total_bases++;
	}
##	my $q_avg = ($q_sum / $num_total);
##	if ($q_avg < $min_quality) {
##		$pass = 0;
###		print STDERR "failed F2 with q_sum=$q_sum num_total=$num_total q_avg=$q_avg min_qual=$min_quality\n";
##	}
##	my $q_sum = 0;
##	my $num_total = 0;
	foreach (split //, $qual_b) {
		my $q = ord($_)-33;
		$q_sum += $q;
		$num_total++;
		if ($q < $min_quality) {
			$num_failed_bases++;
		}
		else {
			$num_above++;
		}
		$num_total_bases++;
	}
	my $q_avg = ($q_sum / $num_total);
##	if ($q_avg < $min_quality) {
##		$pass = 0;
###		print STDERR "failed B with q_sum=$q_sum num_total=$num_total q_avg=$q_avg min_qual=$min_quality\n";
##	}
	if ((100*($num_above / $num_total)) < $min_pct_bases) {
		$pass = 0;
#		print STDERR "failed with num_above=$num_above num_total=$num_total min_pct_bases=$min_pct_bases\n";
	}
	if ($pass == 1) {
		# print R1, R2, and barcode
		print FASTQ1_OUT "$head_f1\n$seq_f1\n+\n$qual_f1\n";
		print FASTQ2_OUT "$head_f2\n$seq_f2\n+\n$qual_f2\n";
		print BARCODE_OUT "$head_b\n$seq_b\n+\n$qual_b\n";
#		print STDERR "passed with num_above=$num_above num_total=$num_total min_pct_bases=$min_pct_bases\n";
	}
	else {
		$num_failed_reads++;
	}
	$num_total_reads++;
	
	if (($num_total_reads % 100000) == 0) {
		my $num_passed_reads = $num_total_reads - $num_failed_reads;
		my $num_passed_bases = $num_total_bases - $num_failed_bases;
		print STDERR "$num_passed_reads\/$num_total_reads ($num_passed_bases\/$num_total_bases bases) passed so far...\n";
	}
	
}

print STDERR "Q$min_quality filtering removed $num_failed_reads\/$num_total_reads reads\n";

close(FASTQ1);
close(FASTQ2);
close(BARCODE);
close(FASTQ1_OUT);
close(FASTQ2_OUT);
close(BARCODE_OUT);


