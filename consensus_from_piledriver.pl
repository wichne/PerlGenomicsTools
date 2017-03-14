#!/usr/bin/perl

use Bio::SeqIO;
use Bio::Seq;
use Getopt::Std;
use strict;
$| = 1;

my %args;
&getopts('p:f:c:D',\%args);

my $piledriver_file = $args{'p'};
my $seq_file = $args{'f'};
my $depth_cutoff = $args{'c'} ? $args{'c'} : 20;
my $DEBUG = $args{'D'};

# read reference sequences into memory
our %SEQ;
my $sfo = Bio::SeqIO->new(-file => $seq_file);
while (my $sro = $sfo->next_seq) {
    warn "Loading sequence " . $sro->display_id . "\n" if ($DEBUG);
    my @seq = split(/\s*/, $sro->seq);
    $SEQ{$sro->display_id} = \@seq;
}
$sfo->close;

# read in piledriver output
# for ease of replacement, we will want to go through these results backwards
my @pd;
warn "Reading in $piledriver_file" if ($DEBUG);
open(my $pfh, $piledriver_file) or die "Can't open piledriver file $piledriver_file: $!\n";
while (my $line = <$pfh>) {
    if ($line =~ /^chrom/) { next }  # this is the header line
    chomp $line;
    my @line = split(/\t/, $line);
    if (! defined $SEQ{$line[0]}) { die "Piledriver seq_id $line[0] not found in sequence file $seq_file.\n"; }
    if ($line[4] < $depth_cutoff) {
	warn "Depth ($line[4]) less than cutoff ($depth_cutoff): $line[0] $line[2]\n" if ($DEBUG);
    } elsif (($line[6] > $line[5]) || ($line[12]/$line[4] > 0.5)) {
	# only evaluate the ones that might represent changes (a_depth > r_depth or num_I > depth/2)
	unshift @pd, \@line;
    }
}
close $pfh;

foreach my $pdref(@pd) {
    my ($chrom,
	$start, # 0-based coord
	$end,   # 1-based coord
	$ref,
	$depth,
	$r_depth,
	$a_depth,
	$num_A,
	$num_C,
	$num_G,
	$num_T,
	$num_D,
	$num_I,
	$totQ_A,
	$totQ_C,
	$totQ_G,
	$totQ_T,
	$all_ins,
	$num_F_A,
	$num_F_C,
	$num_F_G,
	$num_F_T,
	$num_F_D,
	$num_F_I,
	$totQ_F_A,
	$totQ_F_C,
	$totQ_F_G,
	$totQ_F_T,
	$all_F_ins,
	$num_R_A,
	$num_R_C,
	$num_R_G,
	$num_R_T,
	$num_R_D,
	$num_R_I,
	$totQ_R_A,
	$totQ_R_C,
	$totQ_R_G,
	$totQ_R_T,
	$all_R_ins,
	$sample_1) = @$pdref;
    if ($num_A/$depth > 0.5 && $ref ne "A") {
	warn "Replacing $ref at $chrom $end with A ($num_A/$depth)\n" if ($DEBUG);
	&replace($chrom, $start, "A");
    } elsif ($num_C/$depth > 0.5 && $ref ne "C") {
	warn "Replacing $ref at $chrom $end with C ($num_C/$depth)\n" if ($DEBUG);
	&replace($chrom, $start, "G");
    } elsif ($num_G/$depth > 0.5 && $ref ne "G") {
	warn "Replacing $ref at $chrom $end with G ($num_G/$depth)\n" if ($DEBUG);
	&replace($chrom, $start, "G");
    } elsif ($num_T/$depth > 0.5 && $ref ne "T") {
	warn "Replacing $ref at $chrom $end with T ($num_T/$depth)\n" if ($DEBUG);
	&replace($chrom, $start, "T");
    } elsif ($num_D/$depth > 0.5) {
	warn "Deleting $ref at $chrom $end ($num_D/$depth)\n" if ($DEBUG);
	&delete($chrom, $start);
    }

    if ($num_I/$depth > 0.5) {
	# count frequency of each insertion string
	my @ins = split(/\,/, $all_ins);
	my %INS;
	foreach my $x(@ins) { $INS{$x}++; }
	# find the most common insertion
	@ins = sort {$INS{$b}<=>$INS{$a}} keys %INS;
	my ($MCI, $num_MCI) = ($ins[0], $INS{$ins[0]});
	if ($num_MCI/$num_I > 0.5) { 
	    warn "Inserting $MCI after $chrom $end ($num_MCI/$num_I/$depth)" if ($DEBUG);
	    &insert($chrom, $start, $MCI);
	} else {
	    warn "At $chrom $end, no identified insertion was dominant. No action.";
	}
    }
}

my $consensus_output = $seq_file . ".consensus.fa";
my $ofh = Bio::SeqIO->new(-file => ">$consensus_output",
			  -format => "fasta");
foreach my $seqid (keys %SEQ) {
    my $seq = join("", @{$SEQ{$seqid}});
    my $seqo = Bio::Seq->new(-display_id => $seqid,
			     -seq => $seq);
    $ofh->write_seq($seqo);
}
$ofh->close;
exit();

sub replace {
    my ($chrom, $coord, $val) = @_;
    $SEQ{$chrom}->[$coord] = $val;
}

sub delete {
    my ($chrom, $coord) = @_;
    my @seq = @{$SEQ{$chrom}};
    my $hi_idx = $#seq;
    my @new = (@seq[0..$coord-1], @seq[$coord+1..$hi_idx]);
    if (@{$SEQ{$chrom}} - @new != 1) { die "Math wrong in delete sub: " . @{$SEQ{$chrom}} . " - " . scalar(@new) . "\n"; }
    else { $SEQ{$chrom} = \@new }
}

sub insert {
    my ($chrom, $coord, $ins) = @_;
    my @ins = split(/\s*/, $ins);
    my @seq = @{$SEQ{$chrom}};
    my $hi_idx = $#seq;
    my @new = (@seq[0..$coord], @ins, @seq[$coord+1..$hi_idx]);
    if (@new - @{$SEQ{$chrom}} != scalar(@ins)) { die "Math wrong in insert sub\n"; }
    else { $SEQ{$chrom} = \@new }
}
