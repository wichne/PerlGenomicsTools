#!/usr/bin/perl

use Bio::SeqIO;
use Bio::SeqIO::fastq;
use Getopt::Std;

my %arg;
&getopts('f:t:q:o:',\%arg);
my $file = $arg{'f'};
my $end5_trim = $arg{'t'};

my $in = Bio::SeqIO->new(-file => $file);
my $outfile = $file . ".trimmed";
$outfile = $arg{'o'} if ($arg{'o'});
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -file => ">$outfile");

while (my $seqo = $in->next_seq) {
    my $seq = $seqo->seq;
    my $qual = $seqo->qual;
    my $id = $seqo->id;
#    print ">$id $desc\n$seq\n" . join(":",@$qual) . "\n\n";

    my $i;
    for ($i=0; $i<@$qual;$i++) {
	if ($qual->[$i] == 2) { last }
    }

    my $qseq = substr($seq, $end5_trim, $i);
    if (!$qseq) { next }
    my $qqual = [@{$qual}[0..$i-1]];

    my $qseqo = Bio::Seq::Quality->new
	( -qual => $qqual,
	  -seq =>  $qseq,
	  -id  => $id,
	);
    $out->write_seq($qseqo);
}
exit;
