#!/usr/bin/perl

use Bio::AlignIO;
use Getopt::Std;

my %arg;
&getopts('i:f:t:', \%arg);
my $input_file = $arg{i} or die "Need to provide input file with -i";
my $from_format = $arg{f};
my $to_format = $arg{t} or die "Need to provide output format with -t";

$output_file = $input_file;
$output_file =~ s/\.\w+$/\.$to_format/;

$in  = Bio::AlignIO->new(-file   => $input_file,
			 -format => $from_format);
$out = Bio::AlignIO->new(-file   => ">$output_file" ,
			 -format => $to_format);

while ( my $aln = $in->next_aln() ) {
    $out->write_aln($aln);
}
