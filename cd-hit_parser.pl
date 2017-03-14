#!/usr/bin/perl
use Bio::SeqIO;

my $nr_file = $ARGV[0];
my $clstr_file = $nr_file . ".clstr";

open (my $clstr, $clstr_file) or die "Can't open clstr file $clstr_file: $!\n";

while (my $l = <$clstr>) {
    chomp $l;
    if ($l =~ /^\>\S+\s+(\d+)/) { # should look like '>Cluster 0'
	$c =  $1;
    } else {
	my ($i, $len, $a, $desc, $at, $perid) = split/\s+/, $l;
	$a =~ s/>//;
	$LU{$a} = $c;
	push @{$C{$c}}, $a;
    }
}

my $nr = Bio::SeqIO->new(-file => $nr_file);
my $nr_new = Bio::SeqIO->new(-file => ">$nr_file.tmp", -format => 'fasta');
while (my $seqo = $nr->next_seq) {
    my $d = $seqo->desc;
    my $a = $seqo->display_id;
    $seqo->display_id("Cluster_$LU{$a}_" . scalar(@{$C{$LU{$a}}}));
    $seqo->desc("$a $d");
    $nr_new->write_seq($seqo);
}
