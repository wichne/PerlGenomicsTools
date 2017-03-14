#!/usr/bin/perl



sub trinucleotide_
$seq;
$seq .= rc($seq); # append the reverse complement to the sequence to calculate the symmetrical odds ratio
$seq = uc($seq); # we don't need to deal with upper vs lower case
$seqlen = length($seq);

# collect the character in the sequence
my @s = split("", $seq);
foreach my $p (@s) {
    $F{$p}++; # this is the sum of the characters;
}

#foreach my $c (keys %S) {
#    $F{$c} = $S{$c}/$seqlen; # collect the frequencies of the individual characters
#}

foreach my $x (keys %S) {
    foreach my $y (keys %S) {
	foreach my $z (keys %S) {
	    while ($seq =~ /$x$y$z/g) {
		$F{$x.$y.$z}++;
	    }
	    $P{$x.$y.$z} = $F{$x.$y.$z}/$F{$x}$F{$y}$F{$z};
	}
    }
}
