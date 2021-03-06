#!/usr/bin/perl
#Edited 11/11/2015 by JMo to incorporate shorter scaffolds

use strict;
use Statistics::Basic qw (:all);
use Number::Format;

my $statsfile = shift @ARGV; # output from fastaStats
my @file = @ARGV; # per-base-coverage file generated by samtools depth command


my (%DATA, %L);
open my $stats, $statsfile or die "bam_to_coverage.pl fastaStats_file by_base_cov_file\n";
# set the DATA structure length to the length of the molecule. Critical for correct coverage determination.
while (my $line=<$stats>) {
    $line =~ />(\S+)\s.*length\=(\d+)/;
#    if ($2 < 999) { next }
    $DATA{$1}->[$2] = 0; # set array to length of sequence
    $L{$1}->{length} = $2; # store length of sequence
}
close $stats;

# This just converts the per-base coverage file to a data structure.
foreach my $file (@file) {
    open my $in, $file;
    while (my $line = <$in>) {
	chomp $line;
	my ($contig, $pos, $cov) = split/\s+/,$line;

	# mystery contig?
	if (! defined($DATA{$contig})) { print STDERR "Where did $contig come from? It's not in the stats file.\n"; }

	# data problem?
	if ($DATA{$contig}->[$pos] != 0) { die "Why does contig $contig position $pos already have a value?\n"; }

	# set coverage value for position
	$DATA{$contig}->[$pos] = $cov;

	# keep track of minimum and maximum coverages
	$L{$1}->{min_cov} = $cov if ($cov < $L{$1}->{min_cov} || ! $L{$1}->{min_cov});
	$L{$1}->{max_cov} = $cov if ($cov > $L{$1}->{max_cov} || ! $L{$1}->{max_cov});
    }
}

# there is a problem that where we have 16S, we get strong cov values that can skew results
# screen for value >3*stddev of set
foreach my $c (sort { $#{$DATA{$b}} <=> $#{$DATA{$a}} } keys %DATA) {
    my @current = @{$DATA{$c}};
    my $mean = &mean(@{$DATA{$c}});
    my $stdev = &stddev(@{$DATA{$c}});
    my @new = grep { $_ > ($mean - 3 * $stdev) && $_ < ($mean + 3 * $stdev) } @current;
    printf "# Molecule %s: length=%d; min_cov=%d; max_cov=%d; good_range=%d-%d; good_bases=%d\n",
    ($c, $L{$c}->{length}, $L{$c}->{min_cov}, $L{$c}->{max_cov}, ($mean - (3 * $stdev)), ($mean + (3 * $stdev)), scalar(@new));
    $DATA{$c} = \@new;
}

print "# scafId scaflen covAvg covStdev basesAnalyzed first3rdCovAvg first3rdCovStdev first3rdLength last3rdCovAvg last3rdCovStdev last3rdLength\n";
foreach my $c (sort { $L{$b}->{'length'} <=> $L{$a}->{'length'} } keys %DATA) {
    # whole thing
#    printf "$c\t%.2f\n", (&mean(@{$DATA{$c}}));
#    printf "$c\t%.2f\t%.2f\t%i\n", (&mean(@{$DATA{$c}}), &stddev(@{$DATA{$c}}), $#{$DATA{$c}});

    # the 0 position is currently empty because coords go from one,
    # so shift the array
#    shift @{$DATA{$c}};
    # first third
    my $len = scalar(@{$DATA{$c}});
    my @first = @{$DATA{$c}}[0..int($len/3)];
#    printf "$c\t%.2f\t%.2f\t%i\n", (&mean(@first), &stddev(@first), $#first);
    # last third
    my @last = @{$DATA{$c}}[($len - int($len/3))..$len];
    printf "$c\t%i\t%.2f\t%.2f\t%i\t%.2f\t%.2f\t%i\t%.2f\t%.2f\t%i\n",
    ($L{$c}->{'length'},
     &mean(@{$DATA{$c}}), &stddev(@{$DATA{$c}}), $#{$DATA{$c}},
     &mean(@first), &stddev(@first), $#first,
     &mean(@last), &stddev(@last), $#last);
}
