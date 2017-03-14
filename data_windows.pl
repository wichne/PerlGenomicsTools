#!/usr/bin/perl
use Getopt::Std;
use Statistics::Basic qw(:all);

&getopts('w:s:');
my $step = $opt_s;
my $window = $opt_w;

if (! $step) { die "Provied step with -s" }
if (! $window) { die "Provide window with -w" }

# read in data
while (my $l = <STDIN>) {
    chomp $l;
    my ($seqacc, $pos, $val) = split/\s+/, $l;
    $DATA{$seqacc}->[$pos] = $val;
}

foreach my $acc (keys %DATA) {
    my $ref = $DATA{$acc};
    for (my $i=0; $i+$window<@{$ref}; $i+=$step) {
	my $m = mean(@{$ref}[$i..$i+$window]);
	my $pos = int($i+$window/2);
	print "$acc\t$pos\t$m\n";
    }
}
