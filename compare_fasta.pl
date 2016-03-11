#!/usr/bin/perl

use strict;

# convert GTF file to FASTA file

my $f1 = shift;
my $f2 = shift;

sub rc {
my $dna = shift;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

open(F1, "<$f1") || die "Unable to open: $!";
open(F2, "<$f2") || die "Unable to open: $!";


my %seqs1;
my %seqs2;
my %strands;
my $id = "";
my $seq = "";
my $strand;
# store entries from tophat
while (my $line = <F1>) {
	chomp $line;
	if ($line =~ /^>\d+ (\S+) \S+([+-])/) {
		if ($id ne "") {
			$seqs1{$id} = $seq;
			$strand = $2;
		}
		$id = $1;
		$seq = "";
		$strands{$id} = $strand;
#		print "setting strand of id=$id = $strand\n";
	}
	else {
		$seq .= uc($line);
	}
}
$seqs1{$id} = $seq;
# reverse complement minus strand seqs
foreach my $key ( keys %seqs1 ) {
  if ($strands{$key} eq "-") {
  	$seqs1{$key} = rc($seqs1{$key});
  }
}

# store entries from gtf2fasta
$id = "";
$seq = "";
$strand = "";
# store entries from tophat
while (my $line = <F2>) {
	chomp $line;
	if ($line =~ /^>(\S+)/) {
		if ($id ne "") {
			$seqs2{$id} = $seq;
		}
		$id = $1;
		$seq = "";
	}
	else {
		$seq .= $line;
	}
}
$seqs2{$id} = $seq;

close(F1);
close(F2);

# compare sequences
foreach my $key ( keys %seqs1 ) {
	if ($seqs1{$key} eq $seqs2{$key}) {
		#print "$key ... OKAY\n";
	}
	else {
		die "$key does not match!\n1: $seqs1{$key}\n2: $seqs2{$key}\n";
	}
}



