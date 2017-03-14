#!/usr/bin/perl

use Bio::AlignIO;

my $infile = $ARGV[0];
my $outfile = $infile . ".sto";

my $in = Bio::AlignIO->new(-file => $infile);
my $out = Bio::AlignIO->new(-file => ">$outfile",
			    -format => 'stockholm');
while ( my $aln = $in->next_aln() ) {
    $out->write_aln($aln);
}

exit();
