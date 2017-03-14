#!/usr/bin/perl
use strict;

# This program looks for bi-directional best hits between two organisms
# Input is blastp output from an all_v_all search between two genomes.
# The top non-self hit for each gene is stored in a hash, then we go 
# through each protein from one of the genomes and check for reciprocal
# best hits. Cutoffs are >30% id over > 70% of the query length
# for the first iteration. (The length cutoff only has to be unidirectional
# to accommodate fusion proteins.) The avg and stdev for per id are
# calculated from the initial results. In the second iteration,
# the cutoff is set to avg-stdev so the top 84% of the pairs are grabbed.
#  Clear?


#require "$ENV{SGC_SCRIPTS}/sgc_library.dbi";
#use DBI;


my $file = shift @ARGV;
my @infos = @ARGV;

my $DEBUG = 1;
my $VERBOSE = 1;

our $length_cutoff = 0.7;
our $id_cutoff = 30;

our %D;

# Read the info files
our %PROT;
foreach my $info_file (@infos) {
    $info_file =~ /(\w+)\.info/;
    my $db = $1;
    open my $in, $info_file or die "Can't open $info_file: $!\n";
    while (my $line = <$in>) {
	chomp $line;
	my @f = split /\t/, $line;
#	$f[0] =~ s/^S\d+\_//;
	$PROT{$db}->{$f[0]}->{location} = $f[4];
	$PROT{$db}->{$f[0]}->{direction} = $f[3];
	$D{$f[0]}->{'genome'} = $db;
    }
}
our %DB;
foreach my $db (keys %PROT) {
    foreach my $prot (sort {$PROT{$db}->{$a}->{location} <=> $PROT{$db}->{$b}->{location}} keys %{$PROT{$db}}) {
	push @{$DB{$db}}, $prot;
	$PROT{$db}->{$prot}->{position} = scalar(@{$DB{$db}}) - 1;
    }
}

&read_file($file); # populates %D

&synteny;

my $cutoff = 0;
our %BBH;
my $per_id_ref = &look_for_bbh($cutoff);

my $index = @$per_id_ref;
my $sum = 0;
foreach (@$per_id_ref) { $sum += $_ }
my $avg = $sum/$index;
my $x = 0;
foreach (@$per_id_ref) { $x += ($_ - $avg)**2 }
my $stdev = sqrt($x/$index);
my $keep_threshold = sprintf "%.2f\n",($avg-$stdev);

if ($VERBOSE) {
    printf "index =    $index\n";
    printf "Avg per_id = %.2f\n",$avg;
    printf "St dev     = %.2f\n",$stdev;
    print "Keep threshold  = $keep_threshold";
}

open OUT, ">G1G2.blastp";

open my $IN, $file or die "Can't open file $file : $!\n";
while (my $line = <$IN>) {
    if ($line =~ /^\#/) { print OUT $line; next; }
    chomp $line;
    my @f = split/\t/, $line;
    if ($f[12] >= $keep_threshold ||
	($D{$f[0]}->{'synteny'} &&
	 $f[10] > 50)) {
	print OUT join("\t", ($f[0], $f[4], $f[12], $f[6], 0, 0, $f[2], $f[3], $f[7], $f[8], $f[13], $f[10]));
	print OUT "\n";
    }
}
#foreach my $p (sort {$a cmp $b} keys %BBH) {
#    if ($D{$p}->{'per_id'} >= $keep_threshold ||
#	$D{$p}->{'synteny'}) {
#	print OUT $D{$p}->{'header'} . $D{$p}->{'blast_data'} . "\n";
#    }
#}

exit();

sub look_for_bbh {
    my @per_id_ref;
    foreach my $p (keys %D) {
	if (! defined $D{$p}->{'tophit'}) { next }
	my $q = $D{$p}->{'tophit'}; # ie top hit from other genome
	if (! defined $D{$q}) {   # perhaps because length cutoff wasnt' met for a fusion protein...
	    print OUT $D{$p}->{$q} . "\n"; #"$p\t$q\t$D{$p}->{per_id}\n";
	    push @per_id_ref, $D{$p}->{'per_id'};
	    $BBH{$p} = {$q};
	} elsif ($D{$q}->{'tophit'} eq $p) {
	    print "$p\t$q\t$D{$p}->{per_id}\n";
	    push @per_id_ref, ($D{$p}->{'per_id'}, $D{$q}->{'per_id'});
	    $BBH{$p} = {$q};
	}
    }
    return (\@per_id_ref);
}

sub read_file {
    my $file = shift;
    open my $IN, $file or die "Can't open file $file : $!\n";

    my $header;
    my $last;
    while (my $line = <$IN>) {
	# Need the header stuff for MSOAR
	if ($line =~ /^\#/) {
	    $header .= $line;
	    next;
	}
	chomp $line;
	my @f = split/\t/, $line;
	if ($line && @f == 14) { 
	    # tab_blast data
	    # Save header and reset variable.
	    if ($f[0] eq $f[4]) { next; } # self hit
	    
	    # determine if these are from the same genome to ignore paralog hits
	    if ($D{$f[0]}->{'genome'} eq $D{$f[4]}->{'genome'}) { next }
#	my ($qg, $sg);
#	if ($f[0] =~ /(\w+)\_/) { $qg = $1 }
#	elsif ($f[0] =~ /([a-zA-Z]+)/) { $qg = $1 }
#	else {$qg = substr($f[0], 3) }
#	if ($f[4] =~ /(\w+)\_/) { $sg = $1 }
#	elsif ($f[4] =~ /([a-zA-Z]+)/){ $sg = $1 }
#	else { $sg = substr($f[4],3) }
#	if ($qg eq $sg) {next}  # This ensures paralog hits are not top hits
	    
	    # This is the meat right here. Comparison against length and id cutoffs
	    # populates the Data table.
	    
	    if ($f[0] ne $last &&  # this is the top hit
		$f[12] >= $id_cutoff &&
		$f[6]/$f[1] >= $length_cutoff) {
		$D{$f[0]} = {'header' => $header,
			     'tophit' => $f[4],
			     'per_id' => $f[12],
			     'blast_data' => $line};
	    }
	    $last = $f[0];
	    $header = "";
	} elsif ($line && @f == 12) {
	    # mtab data
	    if ($f[0] eq $f[1]) { 
		$D{$f[0]}->{self_hit} = $f[11] * 1;
		next;
	    } # self hit
	    if ($D{$f[0]}->{'genome'} eq $D{$f[1]}->{'genome'}) { next }

	    # for this data, we don't have the query and reference seq lengths,
	    # so we'll use a bit score cutoff instead
	    if ($f[0] ne $last &&  # this is the top hit
		$f[2] >= $id_cutoff &&
		$f[11] > $D{$f[0]}->{self_hit}/10) {
		$D{$f[0]} = {'header' => $header,
			     'tophit' => $f[1],
			     'per_id' => $f[2],
			     'blast_data' => $line};
	    }
	    $last = $f[0];
	    $header = "";
	} else {
	    die "Don't recognize format of tabulated blast output\n";
	}
    }
    close $IN;
}


sub synteny {
    foreach my $p (keys %D) {
	my $db = $D{$p}->{genome}; 
	if (! defined $D{$p}->{'tophit'}) { next }
	my $q = $D{$p}->{'tophit'}; # ie top hit from other genome
	my $qdb = $D{$q}->{genome};

	my $ppos = $PROT{$db}->{$p}->{'position'};
	my $p_usg = $PROT{$db}->{$p}->{'direction'} eq "+" ? $DB{$db}->[$ppos-1] : $DB{$db}->[$ppos+1];
	my $p_dsg = $PROT{$db}->{$p}->{'direction'} eq "+" ? $DB{$db}->[$ppos+1] : $DB{$db}->[$ppos-1];
	my $qpos = $PROT{$qdb}->{$q}->{'position'};
	my $q_usg = $PROT{$qdb}->{$q}->{'direction'} eq "+" ? $DB{$qdb}->[$qpos-1] : $DB{$qdb}->[$qpos+1];
	my $q_dsg = $PROT{$qdb}->{$q}->{'direction'} eq "+" ? $DB{$qdb}->[$qpos+1] : $DB{$qdb}->[$qpos-1];

	if ($D{$p_usg}->{'tophit'} eq $q_usg ||
	    $D{$p_dsg}->{'tophit'} eq $q_dsg) { $D{$p}->{'synteny'} = 1 }
    }
}
