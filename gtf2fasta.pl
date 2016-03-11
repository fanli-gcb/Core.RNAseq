#!/usr/bin/perl

use strict;
use Getopt::Long qw( :config posix_default no_ignore_case );
use Pod::Usage;

# convert GTF file to FASTA file

my $selecttype = "exon";
my $fieldname = "transcript_id";
my $gtf_input_file;
my $genome_fasta_file;
my $output_fasta_file;
my $man = 0;
my $help = 0;
# Parse options and print usage if there is a syntax error,
# or if usage was explicitly requested.
GetOptions('type=s' => \$selecttype, 'field=s' => \$fieldname, 'help|?' => \$help, 'man' => \$man) or pod2usage(-verbose => 1);
pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

my $remaining = @ARGV;
if ($remaining < 3) {
        printf STDERR "\nERROR: Missing required arguments(s).\n\n";
        pod2usage(-verbose => 1);
}
$gtf_input_file = shift;
$genome_fasta_file = shift;
$output_fasta_file = shift;
if (@ARGV) {
        printf STDERR "\nWarning: Extra arguments/options detected. Exiting!\n\n";
        pod2usage(-verbose => 1);
}

#chr1    unknown exon    11874   12227   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS16107";
#chr1    unknown exon    12613   12721   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS16107";
#chr1    unknown exon    13221   14409   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS16107";
#chr1    unknown exon    14362   14829   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS8151";
#chr1    unknown exon    14970   15038   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS8151";


sub rc {
my $dna = shift;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}


open(INPUT, "<$gtf_input_file") || die "Unable to open: $!";
open(GENOME, "<$genome_fasta_file") || die "Unable to open: $!";
open(OUTPUT, ">$output_fasta_file") || die "Unable to write: $!";

# store genome in memory
my %genome;
my $id = "";
my $seq = "";
while (my $line = <GENOME>) {
	chomp $line;
	if ($line =~ /^>(\S+)/) {
		if ($id ne "") {
			$genome{$id} = $seq;
			print STDERR "storing $id\n";
		}
		$id = $1;
		$seq = "";
	}
	else {
		$seq .= $line;
	}
}
$genome{$id} = $seq;
print STDERR "storing $id\n";

# process GTF file
my %seqs;
my %strands;
my $key;
my @seqordering;

while (my $line = <INPUT>) {
	chomp $line;
	my ($chrom, $source, $type, $start, $end, $score, $strand, $frame, $attr_string) = split(/\t/, $line);
	next unless ($type eq $selecttype);
	if ($attr_string =~ /$fieldname \"(\S+)\";/) {
		$key = $1;
	}
	else {
		die "ERROR: Failed to parse $attr_string for field $fieldname!\n";
	}
	
	# store sequence
	if (!(defined($seqs{$key}))) {
		$seqs{$key} = "";
		$strands{$key} = $strand;
		push @seqordering, $key;
	}
	my $seq = substr $genome{$chrom}, $start-1, $end-$start+1;
	$seqs{$key} .= uc($seq);
}


# reverse complement minus strand seqs
foreach my $key ( @seqordering ) {
  if ($strands{$key} eq "-") {
  	$seqs{$key} = rc($seqs{$key});
  }
  
  print OUTPUT ">$key\n$seqs{$key}\n";
}

close(GENOME);
close(INPUT);
close(OUTPUT);


__END__

=head1 NAME
 
gtf2fasta - Convert GTF file to FASTA format
 
=head1 SYNOPSIS
 
gtf2fasta [options] GTF_input_file genome_fasta_file output_fasta_file

=head1 OPTIONS

=over 8
 
=item B<-t/--type TYPE>

Extract features matching this TYPE (3rd column of GTF file). [Default: exon]

=item B<-f/--field FIELD>

Field identifier used to separate sequences (from the GTF attribute string). [Default: transcript_id]

=back

=head1 ARGUMENTS

=over 8

=item B<GTF_input_file>

Input GTF file

=item B<genome_fasta_file>

Genome sequence file in FASTA format. Please check that the chromosome names match...

=item B<output_fasta_file>

Output file in FASTA format

=back

=head1 DESCRIPTION

Convert GTF file to FASTA format, extracting sequence using the provided genome FASTA file

=cut

