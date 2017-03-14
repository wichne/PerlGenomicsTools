#!/usr/bin/perl
use strict;
use Statistics::Basic qw (:all);
use Number::Format;

my %DATA;
while (my $file = shift @ARGV) {
#my $coordsfile = shift @ARGV;

    my $last;
    open my $in, $file;
    while (my $line = <$in>) {
	chomp $line;
	my ($contig, $pos, $cov) = split/\s+/,$line;
#    if (! defined($DATA{$contig})) { next }
#    print STDERR "$contig\n" if ($contig ne $last);
	$DATA{$contig}->[$pos] += $cov;
#    push @{$DATA{$contig}}, $cov;
	$last = $contig;
    }
}

print "# scafId covAvg covStdev Length first3rdCovAvg first3rdCovStdev first3rdLength last3rdCovAvg last3rdCovStdev last3rdLength\n";
#open my $coords, $coordsfile or die "bam_to_coverage.pl by_base_cov_file fastaStats_file\n";
while (my $line=<STDIN>) {
    chomp $line;
    my ($contig, $lo, $hi, $len) = split/\t/, $line;
    if ($len < 20) { next }
    my @reg = @{$DATA{$contig}}[$lo..$hi];
		# whole thing
#    printf "$c\t%.2f\n", (&mean(@{$DATA{$c}}));
#    printf "$c\t%.2f\t%.2f\t%i\n", (&mean(@{$DATA{$c}}), &stddev(@{$DATA{$c}}), $#{$DATA{$c}});
		# first third
#    my @first = @{$DATA{$c}}[0..int($len/3)];
#    printf "$c\t%.2f\t%.2f\t%i\n", (&mean(@first), &stddev(@first), $#first);
		# last third
#    my @last = @{$DATA{$c}}[($len - int($len/3))..$len];
		printf "$contig\t$lo\t$hi\t%.2f\t%.2f\t%i\n", #%.2f\t%.2f\t%i\t%.2f\t%.2f\t%i\n",
		(&mean(@reg), &stddev(@reg), $len);
#     &mean(@first), &stddev(@first), $#first,
#     &mean(@last), &stddev(@last), $#last);
}
