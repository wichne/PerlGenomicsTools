#!/usr/bin/perl
#Script used to parse out abundance counts for MaxBin bins
#usage:
#process_maxbin.pl -s prefix.abund prefix.index.fasta


use Getopt::Std;

my %arg;
&getopt('s:', \%arg);
my $stats_file = $arg{'s'};
if (! $stats_file) { die "Provide scaffold/contig fastaStats output with -s\n";}

use strict;
foreach my $file (@ARGV) {
    $file =~ /(\w+)\.(\d+)\.fasta/;
    my $prefix = $1;
    my $index = $2;
    my $out = "$prefix.$index.csv";
    open my $OUT, ">$out";
    my $abund_file = $prefix . ".abund";
    open my $in, $file;
    while (my $line = <$in>) {
	if ($line =~ /^>(\S+)/) {
	    my $scaf = $1;

	    my $stats = `grep -w $scaf $stats_file`;
	    $stats =~ /length\=(\d+).*G\+C\=([\d\.]+)/;
	    my $len = $1; my $gc = $2;
	    
	    my $abundl = `grep -w $scaf $abund_file`;
	    chomp $abundl;
	    my ($id, $abund) = split/\t/,$abundl;
	    print $OUT "$scaf\t$len\t$gc\t$abund\n";
	}
    }
}
