#!/usr/bin/perl

use Bio::SearchIO;
use strict;
use Getopt::Std;

my $usage = "
-S percent similarity cutoff
-I percent identity cutoff
-B bit score cutoff
-E E-value cutoff
-L percent length of query cutoff
-l percent length of match cutoff
-a length of alignment region cutoff
-n emulate ncbi tabular output
-h print this message

OUTPUT
[0] query identifier
[1] query length
[2] query start
[3] query end
[4] hit identifier
[5] hit description
[6] length
[7] hit start
[8] hit end
[9] strand (1/-1) or frame (1/2/3/-1/-2/-3)
[10] bit score
[11] %similarity
[12] %identity
[13] E value\n";

$| = 1;

my %opt;
&getopts('S:I:B:E:L:l:a:nh', \%opt);

if ($opt{h}) { die $usage }

my $input_file = $ARGV[0];

my $in;
if ($input_file) {
    $in = new Bio::SearchIO(-format => 'blast',
			    -file => $input_file);
} else {
    $in = new Bio::SearchIO(-format => 'blast',
			    -fh => \*STDIN);
}
while (my $result = $in->next_result ) {
unless ( $opt{'n'} ) {
        print "\# " . $result->algorithm . " " . $result->algorithm_version . "\n";
    print "\# Parameters ";
    foreach ($result->available_parameters) {
        print "$_ = " . $result->get_parameter($_) . ";";
    }
    print "\n\# Query: " . $result->query_name . "\n"; # . " length=" . $result->query_length . "\n";
    print "\# Database: " . $result->database_name . "\n"; # . " entries=" . $result->database_entries . " letters=" . $result->database_letters . "\n";
#    print "\# Statistics " . $result->available_statistics . "\n";
    print "# " . $result->num_hits . " hits found\n";
}
while (my $hit = $result->next_hit) {
    while (my $hsp = $hit->next_hsp) {
	if (&filter($result, $hit, $hsp)) {
	    my $dir;
	    if ($hsp->algorithm eq "BLASTN" ||
		$hsp->algorithm eq "MEGABLAST") {
		$dir = $hsp->strand('query') == $hsp->strand('hit') ? 1 : -1;
	    } else {
		$dir = ($hsp->frame('query') + 1) * $hsp->strand("query");
	    }
	    if ($opt{'n'}) {
		    my $gap_op = scalar(split(/[\-\.]+/,$hsp->query_string)) - 1 + scalar(split(/[\-\.]+/,$hsp->hit_string)) - 1;
		    my ($end5, $end3) = $hsp->strand > 0 ? ($hsp->start('query'), $hsp->end('query')) :
			($hsp->end('query'), $hsp->start('query'));
		    print join ("\t", ($result->query_name,
				       $hit->name,
				       sprintf("%.2i", $hsp->percent_identity),
				       $hsp->length,
				       $hsp->length - $hsp->num_conserved,
				       $gap_op,
				       $end5,
				       $end3,
				       $hsp->start('subject'),
				       $hsp->end('subject'),
				       $hsp->evalue,
				       $hsp->bits));
		    print "\n";
		} else {
		    print join ("\t", ($result->query_name,
				       $result->query_length,
				       $hsp->start('query'),
				       $hsp->end('query'),
				       $hit->name,
				       $hit->description,
				       $hit->length,
				       $hsp->start('subject'),
				       $hsp->end('subject'),
				       $dir,
				       $hsp->bits,
				       sprintf("%.2i",($hsp->frac_conserved * 100)),
				       sprintf("%.2i", $hsp->percent_identity),
				       $hsp->evalue));
		    print "\n";
		}
	    }
	}
    }
}

exit();


sub filter {
    my ($result, $hit, $hsp) = @_;
    my $pass = 1;

    if (defined $opt{B} && $hsp->bits < $opt{B}) { $pass = 0 }
    if (defined $opt{S} && $hsp->frac_conserved * 100 < $opt{S}) { $pass = 0 }
    if (defined $opt{I} && $hsp->percent_identity < $opt{I}) { $pass = 0 }
    if (defined $opt{E} && $hsp->evalue > $opt{E}) { $pass = 0 }
    if (defined $opt{L} && $hsp->length('query')/$result->query_length * 100 < $opt{L}) { $pass = 0 }
    if (defined $opt{l} && $hsp->length('hit')/$hit->length * 100 < $opt{l}) { $pass = 0 }
    if (defined $opt{a} && $hsp->length('query') < $opt{a}) { $pass = 0 }
    return $pass;
}
