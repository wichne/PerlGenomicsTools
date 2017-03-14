#!/usr/bin/perl

use Bio::SeqIO;

my $infile = $ARGV[0];

my $ino = new Bio::SeqIO(-file => $infile);

my $outfile1 = $infile . ".1";
my $outo1 = new Bio::SeqIO(-file => ">$outfile1", -format => $ino->format);
my $outfile2 = $infile . ".2";
my $outo2 = new Bio::SeqIO(-file => ">$outfile2", -format => $ino->format);

my @out = ($outo1, $outo2);

my $lasto;
while (my $seqo = $ino->next_seq()) {
    if ($lasto && $seqo->display_id ne $lasto->display_id) {
	warn $lasto->display_id . " apparently has no match (" . $seqo->display_id . "). Not writing.\n";
	$lasto = $seqo;
	next;
    } elsif ($lasto && $seqo->display_id eq $lasto->display_id) {
	$outo1->write_seq($lasto);
	$outo2->write_seq($seqo);
	$lasto = "";
    } else {
	$lasto = $seqo;
    }
}
