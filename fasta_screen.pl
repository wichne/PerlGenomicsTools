#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use Getopt::Std;

my %args;
# i is for input file path
# f is for format (not implemented yet)
&getopts('i:l:L:g:G:h',\%args);

if ($args{h}) {
    print "
USAGE:

fasta_screen.pl -i fastafile -l minimum_length -L maximum_length -g min_GC -G max_GC

Output is to STDOUT

";
}
my $in;
if ($args{i}) {
    $in = new Bio::SeqIO(-format => 'fasta',
			    -file => $args{i});
} else {
    die "Need to provide fastafile with -i";
}
my $out = new Bio::SeqIO(-format => 'fasta');

while (my $seq = $in->next_seq()) {
    if (defined $args{l} && $seq->length < $args{l}) { next }
    if (defined $args{L} && $seq->length > $args{L}) { next }
    if (defined $args{g} || defined $args{G}) {
	my %CHAR;
	while ($seq->seq =~ /(.)/g) { $CHAR{$1}++ }
	my $gc = sprintf "%.2f",(($CHAR{G} + $CHAR{C} + $CHAR{g} + $CHAR{c})/$seq->length);
	if (defined $args{g} && $gc < $args{g}) { next }
	if (defined $args{G} && $gc > $args{G}) { next }
    }

    $out->write_seq($seq);
}
exit();
