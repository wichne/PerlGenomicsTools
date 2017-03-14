#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = "
# This program looks for bi-directional best hits between organism pairs.
# Input is blastp output from an all_v_all search between multiple genomes.
# The top non-self hit for each gene to each genome is stored in a hash, then we go 
# through each protein from one of the genomes and check for reciprocal
# best hits. Cutoffs are >30% id over > 70% of the query length
# for the first iteration. (The length cutoff only has to be unidirectional
# to accommodate fusion proteins.) The avg and stdev for %id (AAI) are
# calculated from the initial results. In the second iteration, the cutoff 
# is set to avg-2*stdev, constraining the 'best' hits to reasonable scores.
# Clear?

USAGE: multi_bbh.pl -i inputfile -n number_of_genomes \@info_files

-i - blast_tab.pl output from an all vs all blastp search.

";


#require "$ENV{SGC_SCRIPTS}/sgc_library.dbi";
#use DBI;

my %arg;
&getopts('i:n:h', \%arg);

if ($arg{'h'}) { die $usage }
my $file = $arg{'i'} or die "Need to provide input btab file with -i\n";
my @infos = @ARGV;
if (!@infos) {
    print STDERR "No info files provided. Synteny will not be evaluated.\n";
}

my $VERBOSE = 1;

our $length_cutoff = 0.7;
our $id_cutoff = 30;

our %HITS; # Lookup table: prot acc -> genome; also holds blast hit info

# Read the info files for positional information
our %POS;
foreach my $info_file (@infos) {
    # genome name is defined by info file prefix
    $info_file =~ /([^\/]*)\.info/;
    my $db = $1;
    print STDERR "Found db $db from $info_file\n";
    open my $in, $info_file or die "Can't open $info_file: $!\n";
    while (my $line = <$in>) {
	chomp $line;
	my @f = split /\t/, $line;
	if (@f > 4) { die "File $info_file is in wrong format (too many columns)\n"; }
#	$f[0] =~ s/^S\d+\_//;
	$POS{$db}->{$f[0]}->{location} = $f[1];
	$POS{$db}->{$f[0]}->{direction} = $f[2];
	$POS{$db}->{$f[0]}->{contig} = $f[3];
	$HITS{$f[0]}->{'genome'} = $db;  # relate prot acc to a genome
    }
}

# Foreach db, sort proteins by position
our %DB; # Key: genome; value: ordered array of prot accs
foreach my $db (keys %POS) {
    foreach my $prot (sort {$POS{$db}->{$a}->{location} <=> $POS{$db}->{$b}->{location}} keys %{$POS{$db}}) {
	push @{$DB{$db}}, $prot;
	$POS{$db}->{$prot}->{position} = $#{$DB{$db}}; # store position index of protein
    }
}

&read_file($file); # populates %HITS

#&synteny if (%POS);

my $cutoff = 0;
our %BBH;
my $per_id_ref = &look_for_bbh($cutoff); # populates %BBH

my @DB = sort keys %$per_id_ref;

# Roll through the pair-wise comparisons
for (my $i=0; $i<@DB; $i++) {
    my $db1 = $DB[$i];
    print STDERR "Doing $db1...\n";
    for (my $j=$i; $j<@DB; $j++) {
	my $db2 = $DB[$j];
	print STDERR "...with $db2\n";
	if (! defined $per_id_ref->{$db1}->{$db2}) { next } # If we don't have scores for these two dbs

	# Calculate the AAI for the two genomes so we can set a threshhold.
	my $index = scalar(@{$per_id_ref->{$db1}->{$db2}}); # how many bbhs are there?
	my $sum = 0;
	foreach (@{$per_id_ref->{$db1}->{$db2}}) { $sum += $_ } # sum the AAIDs
	my $avg = $sum/$index; # calculate the average

	# Calculate the stdev for bbh AAI
	my $x = 0;
	foreach (@{$per_id_ref->{$db1}->{$db2}}) { $x += ($_ - $avg)**2 }
	my $stdev = sqrt($x/$index);

	# The threshhold is 2 stdev from the average.
	my $keep_threshold = sprintf "%.2f\n",($avg-2*$stdev);
	if ($keep_threshold < $id_cutoff) { $keep_threshold = $id_cutoff }
	
	if ($VERBOSE) {
	    print "$db1 - $db2\n";
	    print "index =    $index\n";
	    printf "Avg per_id = %.2f\n",$avg;
	    printf "St dev     = %.2f\n",$stdev;
	    print "Keep threshold  = $keep_threshold\n";
	}
	
	open OUT, ">${db1}.${db2}.ortho";	

	# This old code prints out the MSOAR input, I think
	# open my $IN, $file or die "Can't open file $file : $!\n";
	# while (my $line = <$IN>) {
	#     if ($line =~ /^\#/) { print OUT $line; next; }
	#     chomp $line;
	#     my @f = split/\t/, $line;
	#     if ($f[12] >= $keep_threshold ||
	# 	($HITS{$f[0]}->{'synteny'} &&
	# 	 $f[10] > 50)) {
	# 	print OUT join("\t", ($f[0], $f[4], $f[12], $f[6], 0, 0, $f[2], $f[3], $f[7], $f[8], $f[13], $f[10]));
	# 	print OUT "\n";
	#     }
	# }

	# we want to print out because this is the multiparanoid input
	# OID     bits    'species'               'IPscore'  gene        boostrap
	# 1       1302    dros_sh_dfl_12.pep      1.000   S0_Dshi_0002    100%
        # 1       1302    bin05.pep       1.000   bin05_02438     100%

	# some of it we'll make up (like IPscore and bootstrap)

	my $index = 0;
	foreach my $p (sort {$a cmp $b} keys %BBH) {
	    if ($HITS{$p}->{'genome'} ne $db1) { next }
	    if (defined $BBH{$p}->{$db2} &&
		$HITS{$p}->{$db2}->{'per_id'} >= $keep_threshold ||
 		$HITS{$p}->{$db2}->{'synteny'}) { # synteny can overcome a low score!!!
		$index++;

		# multiParanoid input
		print OUT "$index\t$HITS{$p}->{$db2}->{bits}\t$db1\t$HITS{$p}->{$db2}->{score}\t$p\t100%\n";
		print OUT "$index\t$HITS{$BBH{$p}->{$db2}}->{$db1}->{bits}\t$db2\t$HITS{$BBH{$p}->{$db2}}->{$db1}->{score}\t$BBH{$p}->{$db2}\t100%\n";

		# not sure what this is
#		print OUT "$p\t$BBH{$p}->{$db2}\t$HITS{$p}->{$db2}->{bits}\n";
# 		print OUT $HITS{$p}->{'header'} . $HITS{$p}->{'blast_data'} . "\n";

		# Basic output
#		print OUT "$p\t$BBH{$p}->{$db2}\n";
 	    }
	}
	close OUT;
    }
}

exit();

sub look_for_bbh {
    my %per_id_ref;

    # Go through each protein
    foreach my $p (keys %HITS) {
	my $db = $HITS{$p}->{'genome'};

	# go through each db
	foreach my $qdb (keys %DB) {
	    # skip if there's no hits
	    if (! defined $HITS{$p}->{$qdb}->{'tophit'}) { next }

	    my $q = $HITS{$p}->{$qdb}->{'tophit'}; # ie top hit from other genome
	    # Check for unidirectional hits. 
	    if (! defined $HITS{$q}->{$db}) {   # perhaps because length cutoff wasnt' met for a fusion protein...
		# Not sure why this code is here. Do we really want to include 
		# unidirectional hits?
#		print OUT "$p\t$q\t$HITS{$p}->{$qdb}->{per_id}\t!!! no reciprocal hit\n";
#		push @per_id_ref, $HITS{$p}->{'per_id'};
#		$BBH{$p} = {$q};
	    }
	    # This is the magic - when two proteins are reciprocal best matches
	    elsif ($HITS{$q}->{$db}->{'tophit'} eq $p) {
#		print "$p\t$q\t$HITS{$p}->{$qdb}->{bits}\n";
		# Add per_ids to array for AAAI calculation
		push @{$per_id_ref{$db}->{$qdb}}, ($HITS{$p}->{$qdb}->{'per_id'}, $HITS{$q}->{$db}->{'per_id'});
		push @{$per_id_ref{$qdb}->{$db}}, ($HITS{$p}->{$qdb}->{'per_id'}, $HITS{$q}->{$db}->{'per_id'});
		$BBH{$p}->{$qdb} = $q;
		$BBH{$q}->{$db} = $p;
	    }
	}
    }
    return (\%per_id_ref);
}

sub read_file {
    my $file = shift;
    print STDERR "Reading file $file\n";
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

	    if ($f[0] eq $f[4]) {
		$HITS{$f[0]}->{self_hit} = $f[10];
		next;
	    } # self hit
	    # determine if these are from the same genome to ignore paralog hits
	    if (! (defined $HITS{$f[0]}->{'genome'} &&
		   defined $HITS{$f[4]}->{'genome'})) {
		die "Hey, fool. The accessions in your fasta file don't match the .info file. $f[0] $f[4]\n";
	    }
		   
	    if ($HITS{$f[0]}->{'genome'} eq $HITS{$f[4]}->{'genome'}) { next }

	    # This is the meat right here. Comparison against length and id cutoffs
	    # populates the Data table.
	    if (!defined $HITS{$f[0]}->{$HITS{$f[4]}->{'genome'}} &&  # this is the top hit to this genome
		$f[12] >= $id_cutoff &&
		$f[6]/$f[1] >= $length_cutoff) {
		my $score = $HITS{$f[0]}->{self_hit} ? $f[10]/$HITS{$f[0]}->{self_hit} : 1;
		$HITS{$f[0]}->{$HITS{$f[4]}->{'genome'}} = {'header' => $header,
							    'tophit' => $f[4],
							    'per_id' => $f[12],
							    'bits'   => $f[10],
							    'score'  => 1,
							    'blast_data' => $line};
	    }
	    $header = "";
	} elsif ($line && @f == 12) {
	    # mtab data
	    if ($f[0] eq $f[1]) { 
		$HITS{$f[0]}->{self_hit} = $f[11] * 1;
		next;
	    } # self hit
	    if ($HITS{$f[0]}->{'genome'} eq $HITS{$f[1]}->{'genome'}) { next }

	    # for this data, we don't have the query and reference seq lengths,
	    # so we'll use a bit score cutoff instead
	    if (!defined $HITS{$f[0]}->{$HITS{$f[1]}->{'genome'}} &&  # this is the top hit to this genome
		$f[2] >= $id_cutoff &&
		$f[11] > $HITS{$f[0]}->{self_hit}/10) {
		my $score = $HITS{$f[0]}->{self_hit} ? $f[11]/$HITS{$f[0]}->{self_hit} : 1;
		$HITS{$f[0]}->{$HITS{$f[1]}->{'genome'}} = {'header' => $header,
							    'tophit' => $f[1],
							    'per_id' => $f[2],
							    'bits'   => $f[11],
							    'score'  => 1,
							    'blast_data' => $line};
	    }
	    $header = "";
	} else {
	    die "Don't recognize format of tabulated blast output\n";
	}
    }
    close $IN;
}


sub synteny {
    # Here we are looking for synteny. Synteny is operationally defined as
    # the genes upstream and/or downstream of the subject gene having BBHs to
    # the genes upstream and/or downstream of the object genes. Synteny ends
    # up with a value of 0 (no synteny), 1 (synteny with either upstream or
    # downstream gene, but not both), or 2 (synteny with upstream and downstream
    # genes).
    foreach my $p (keys %HITS) {
	my $db = $HITS{$p}->{genome};
	foreach my $qdb (keys %DB) {
	    if (! defined $HITS{$p}->{$qdb}->{'tophit'}) { next }
	    my $q = $HITS{$p}->{$qdb}->{'tophit'}; # ie top hit from other genome

	    # Get position index of subject gene
	    my $ppos = $POS{$db}->{$p}->{'position'};
	    # Get upstream gene
	    my $p_usg = $POS{$db}->{$p}->{'direction'} eq "+" ? $DB{$db}->[$ppos-1] : $DB{$db}->[$ppos+1];
	    # Get downstream gene
	    my $p_dsg = $POS{$db}->{$p}->{'direction'} eq "+" ? $DB{$db}->[$ppos+1] : $DB{$db}->[$ppos-1];

	    # Get position of object gene
	    my $qpos = $POS{$qdb}->{$q}->{'position'};
	    # Get upstream gene
	    my $q_usg = $POS{$qdb}->{$q}->{'direction'} eq "+" ? $DB{$qdb}->[$qpos-1] : $DB{$qdb}->[$qpos+1];
	    # Get downstream gene
	    my $q_dsg = $POS{$qdb}->{$q}->{'direction'} eq "+" ? $DB{$qdb}->[$qpos+1] : $DB{$qdb}->[$qpos-1];
	    
	    if ($HITS{$p_usg}->{'tophit'} eq $q_usg) { $HITS{$p}->{$qdb}->{'synteny'}++ }
	    if ($HITS{$p_dsg}->{'tophit'} eq $q_dsg) { $HITS{$p}->{$qdb}->{'synteny'}++ }
	}
    }
}
