#!/usr/bin/perl

use Bio::SeqIO;
use Bio::Seq;
use strict;
use Getopt::Std;
$| = 1;

my %opt;
# r - reverse complement
# s - starting coord
# f - fasta file
# i - seq identifier (required if multifasta)
# d - debug
&getopts('rs:f:i:d', \%opt);
our $DEBUG = $opt{'d'};

# read in fasta file
if (!defined $opt{'f'}) { die "Need to specify sequence file with -f\n";}
our %SEQ;
my $sfo = Bio::SeqIO->new(-file => $opt{'f'});
while (my $sro = $sfo->next_seq) {
    warn "Loading sequence " . $sro->display_id . "\n" if ($DEBUG);
    $SEQ{$sro->display_id} = $sro;
}
$sfo->close;

# get the sequence to be rotated
my $so;
if (scalar(keys %SEQ) == 1) {
    my @k = keys %SEQ;
    $so = $SEQ{$k[0]};
} elsif (! defined $opt{'i'}) {
    die "Multifasta input detected. Need to specify which molecule is to be rotated with -i\n";
} else {
    $so = $SEQ{$opt{'i'}}
}

# rotate sequence
if (! defined $opt{'s'}) { die "Need to provide starting coord with -s. (Can also specify reverse complement with -r.)\n";}
my $roto;
if ($opt{'r'}) {
    my $first_part = $so->subseq(1, $opt{'s'});
    my $second_part = $so->subseq($opt{'s'}+1, $so->length);
    my $rot_seq = $second_part . $first_part;
    $roto = Bio::Seq->new(-display_id => $so->display_id,
			   -seq => $rot_seq);
    $roto = $roto->revcom;
} else {
    my $first_part = $so->subseq($opt{'s'}, $so->length);
    my $second_part = $so->subseq(1, $opt{'s'}-1);
    my $rot_seq = $first_part . $second_part;
    $roto = Bio::Seq->new(-display_id => $so->display_id,
			   -seq => $rot_seq);
}

if ($roto->length != $so->length) { die "Pasta fa zool! The length's don't match!\n"; }

# print it out
my $outfile = $opt{'f'} . ".rotated.fa";
my $outo = Bio::SeqIO->new(-file => ">$outfile",
			   -format => 'fasta');
$outo->write_seq($roto);

exit();
