#!/usr/bin/perl
use strict;
use Statistics::Basic qw (:all);
use Number::Format;

my $statsfile = shift @ARGV;
my @covfiles = @ARGV;

my %REG;
open (my $stats, $statsfile) or die "bam_to_coverage.pl coords_file(id low high) by_base_cov_file \n";
# set the DATA structure length to the length of the molecule.
while (my $line=<$stats>) {
    chomp $line;# =~ 
    my ($id, $lo, $hi, @therest);
    if ($line =~ />(\S+)\s.*length\=(\d+)/) {
	($id, $lo, $hi) = ($1, 1, $2);
    } else {
	($id, $lo, $hi, @therest) = split /\s+/, $line;
    }
#    if (($hi-$lo) < 2000) { next }
    push @{$REG{$id}}, [$lo, $hi];
}
close $stats;

#my $last;
my %DATA;
foreach my $file (@covfiles) {
    open (my $in, $file) or die "Couldn't open $file: $!\n";
    while (my $line = <$in>) {
	chomp $line;
	my ($contig, $pos, $cov) = split/\s+/,$line;
#    if (! defined($DATA{$contig})) { next }
#    print STDERR "$contig\n" if ($contig ne $last);
	$DATA{$contig}->[$pos] = $cov;
#    push @{$DATA{$contig}}, $cov;
#	$last = $contig;
    }
}

print "# scafId covAvg covStdev Length first3rdCovAvg first3rdCovStdev first3rdLength last3rdCovAvg last3rdCovStdev last3rdLength\n";
foreach my $c (keys %REG) {
#foreach my $c (sort { $#{$DATA{$b}} <=> $#{$DATA{$a}} } keys %DATA) {
    # whole thing
#    printf "$c\t%.2f\n", (&mean(@{$DATA{$c}}));
#    printf "$c\t%.2f\t%.2f\t%i\n", (&mean(@{$DATA{$c}}), &stddev(@{$DATA{$c}}), $#{$DATA{$c}});
    # first third
    foreach my $coord_ref (@{$REG{$c}}) {
	my ($lo, $hi) = @$coord_ref;
	my $len = $hi - $lo + 1;
#	print STDERR "$c: $lo, $hi, $len\n";
	#my $len = scalar(@{$DATA{$c}});
	my @reg = @{$DATA{$c}}[$lo..$hi];
	my @first = @{$DATA{$c}}[$lo..($lo+int($len/3)-1)];
#    printf "$c\t%.2f\t%.2f\t%i\n", (&mean(@first), &stddev(@first), $#first);
    # last third
	my @last = @{$DATA{$c}}[($hi - (int($len/3) - 1))..$hi];
	printf "$c\t$lo\t$hi\t%.2f\t%.2f\t%i\t%.2f\t%.2f\t%i\t%.2f\t%.2f\t%i\n",
	(&mean(@reg), &stddev(@reg), $len,
	 &mean(@first), &stddev(@first), $#first,
	 &mean(@last), &stddev(@last), $#last);
    }
}
