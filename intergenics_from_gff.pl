#!/usr/bin/perl

my $length_cutoff = 90;

my $file = $ARGV[0];
open (my $in, $file);
while (my $l = <$in>) {
    chomp $l;
    @f = split/\t/, $l;
    $seq_acc = $f[0];
    $lo = $f[3];
    $hi = $f[4];

    if (! $prev_acc && $lo >= 90) {
	printf "1\t%i\t$seq_acc\t%i\n", ($lo - 1, $lo - 1);
    } elsif ($prev_acc && $prev_acc ne $seq_acc) {
	printf "%s\tEND\t$prev_acc\t\n", ($prev_hi + 1);
    } if ($prev_acc eq $seq_acc && ($lo - $prev_hi + 1) >= $length_cutoff) {
	printf "%i\t%i\t$seq_acc\t%i\n", ($prev_hi + 1, $lo -1, ($lo - $prev_hi));
    }

    $prev_hi = $hi;
    $prev_acc = $seq_acc;
}
