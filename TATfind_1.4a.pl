#!/usr/bin/perl

# This program replicates the functionality of TATFIND v 1.4 (http://http://signalfind.org/tatfind.html
# originally published as Rose, R.W., T. Brüser,. J. C. Kissinger, and M. Pohlschröder. 2002. 
# Adaptation of protein secretion to extremely high salt concentrations by extensive use of the twin 
# arginine translocation pathway. Mol. Microbiol. 5: 943-950)
# WITH ONE DIFFERENCE (thus the 1.4a, instead of 1.4):
# This program accepts KR as well as RR as the core of the TAT signal

use strict;
use Bio::SeqIO;
use Getopt::Std;

my $debug = 0;

if (@ARGV < 1) {
    print "usage: TATFIND.pl infile\n";
    die;
}

my $infile = $ARGV[0];
my $outfile = $infile;
$outfile =~ s/.*\///;
$outfile .= ".tatfind.out";
open OUTFILE, ">$outfile";
print OUTFILE "#TATFIND - version 1.4a January 2016; Nelson, adapted from  Brueser, Rose, Kissinger & Pohlschroder\nSEARCH RESULTS for $infile\:\n" if (! $infile);  
print OUTFILE "#acc\tTATrange\tTATseq\tHyPhBrange\thydrophobicity\tcleavage site\tTMH\tHydrophobic13\tFirst44\n";

#################################################################################
## Set up a hash with the hydrophobicity values                                ##
## Values are taken from: Cid et al., (1992) Hydrophobicity and structural     ##
## classes in proteins. Protein Engineering 5:373-375                          ##
#################################################################################

my %hydrophobicity = (
    E => '-1.14',
    Q => '-1.10', 
    D => '-1.04',
    S => '-0.97', 
    G => '-0.80', 
    N => '-0.77',
    T => '-0.77',
    R => '-0.42',
    K => '-0.41',
    P => '-0.09', 
    A => '0.02',
    H => '0.26',
    C => '0.77',
    M => '1.00', 
    Y => '1.11',
    V => '1.13',
    L => '1.14',
    F => '1.35',
    W => '1.71',
    I => '1.81',
    );

################################################################################
## Drink in the input fasta sequences
################################################################################
my $sfo = Bio::SeqIO->new(-file => $infile);

while (my $rec = $sfo->next_seq) {
    my ($TAT_signal, $RRplus1, $RRplus4, $hyd_sum, $preceeding_res);
    my $seq = $rec->seq;
    my @aa = split("", $rec->seq);
    my $header = $rec->display_id . " " . $rec->desc;
    print OUTFILE $rec->display_id;
    print OUTFILE "\nResults for $header\:\n" if ($debug);

################################################################################
## Search for the following pattern between residues 2 and 35 of the 
## predicted protein:
##  (X−1)R=0R+1(X+2)(X+3)(X+4) S/TRRXFLK
## where the hydrophobicity score (Cid et al. 1992) of the amino acid at position
##   X−1 ≤ 0.26; (so ARNDQEGHKPST)
##   X+2 ≤ 0.02; (so ARNDQEGKPST)
##   X+3 ≥ −0.77 (so ANCILMFPTWYV) (positively charged residues are excluded from this position);
##   and X+4 was one of ILVMF
###############################################################################
    print OUTFILE "First 60 aa: " . join("", @aa[0..59]) . "\n" if ($debug);
    my $decision = "N";
    if ($seq =~ /([ARNDQEGHKPST][KR]R([ARNDQEGKPST])[ANCHILMFPTWYV][ILVMF])(.)/ig) {
	my $TAT_sig = $1;
	$RRplus1 = $2;
	$RRplus4 = $3;

	my $RR_sig_start = (pos $seq) - 6;
	my $RR_sig_end = (pos $seq) - 1;
	print OUTFILE "\t$RR_sig_start..$RR_sig_end\t$TAT_sig";
	print OUTFILE "Found twin-arginine $1 at $RR_sig_start - $RR_sig_end\n" if ($debug);
	if (pos $seq < 100) {
	    # reset search pointer
	    pos $seq = $RR_sig_start;
###############################################################################
##  (i) search for an uncharged stretch of at least 13 residues in the 22 
##      residues following the RR
##  (ii) search for:
##       a charged residue after the twin arginine OR
##       a charged residue four residues after the twin arginine OR
##       a basic residue immediately preceeding the hydrophobic region
##  (iii) the hydrophobicity sum of the first 13 residues of the uncharged 
##        region must be < 8.0
###############################################################################	
	    if ($seq =~ /(.)([^RKDE]{12,})/ig) {
		$preceeding_res = $1;
		my $unch_region = $2;
		my $unch_start = length($`) + 2;
		my $unch_min = $unch_start + 12; # set to the minimum possible
		my $unch_end = $unch_start + length($&) - 2;
		print OUTFILE "\t$unch_start..$unch_end";
		print OUTFILE "Found an uncharged stretch of 13+ aa at ($unch_start..$unch_end): $unch_region\n" if ($debug);
		if ($unch_start > $RR_sig_start && $unch_start <= $RR_sig_start + 23) {
		    my @hyd = split("", $unch_region);
		    for (my $i=0; $i<13; $i++) {
			$hyd_sum += $hydrophobicity{$hyd[$i]};
		    }
		    print OUTFILE "\t$hyd_sum";
		    
		    if (($RRplus1 =~ /[DERK]/i || $RRplus4 =~ /[DERK]/i)
			|| ($preceeding_res =~ /[RK]/i && $unch_start > $RR_sig_end + 1)
                        ) {
			if ($hyd_sum < 8.0) {
			    $decision = "Y";

			    # test to see if there's a signal peptide cleavage site.
			    my $cleavage = "-";

			    # write a temp fasta file
			    my $tmpfile = ".tmp.$$";
			    my $seqo = Bio::Seq->new(-display_id => $rec->display_id,
						     -seq => join("", @aa[0..120]));
			    my $filo = Bio::SeqIO->new(-file => ">$tmpfile",
						       -format => 'fasta');
			    $filo->write_seq($seqo);

			    # first, look for SpI using signalp
			    my $sigp_cmd = "signalp -t gram- -c 120 $tmpfile";
			    my $sigp_result = `$sigp_cmd` || die "Can't execute $sigp_cmd\n";
			    my @sigp_results = split(/\n/, $sigp_result);
			    # first two lines are comment, so examine third line
			    my @f = split(/\s+/, $sigp_results[2]);
			    if ($f[2] >= $unch_min && $f[9] eq "Y") { $cleavage = $f[2] }

			    if ($cleavage eq "-") {
				# then look for SpII using LipoP
				my $lipop_cmd = "LipoP $tmpfile";
				my $lipop_result = `$lipop_cmd` || die "Can't execute $lipop_cmd: $!\n";
				if ($lipop_result =~ /SpII.*cleavage\=(\d+)/ && $1 >= $unch_min) {
				    $cleavage=$1;
				}
			    }
			    print OUTFILE "\t$cleavage";

			    # are there too many TMHs? How many is too many?
			    my $tmh;
			    my $tmhmm_cmd = "tmhmm $tmpfile";
			    my $tmhmm_results = `$tmhmm_cmd`;
			    if ($tmhmm_results =~ /PredHel\=(\d+)/) {
				$tmh = $1;
			    }
			    print OUTFILE "\t$tmh";

#			    print OUTFILE $rec->display_id . "\t$RR_sig_start..$RR_sig_end\t$TAT_sig\t$unch_start..$unch_end\t$hyd_sum\t$cleavage\t$tmh\t" . join("", @hyd[0..12]) . "\t" . join("", @aa[0..56]) . "\n";
			    print OUTFILE "\t" . join("", @hyd[0..12]) . "\t" . join("", @aa[0..56]) . "\n";
			} else { print OUTFILE "\t\t" . join("", @hyd[0..12]) . "\t" . join("", @aa[0..56]) . "\tHydrophobicity sum is too high\n"; }
		    } else { print OUTFILE "\t\t\t" . join("", @hyd[0..12]) . "\t" . join("", @aa[0..56]) . "\tDidn't meet charged residue criteria: $preceeding_res$RRplus1..$RRplus4\n"; }
		} else { print OUTFILE "\t\t\t\t\t" .  join("", @aa[0..56]) . "\tUncharged region is out of range.\n"; }
	    } else { print OUTFILE "\t\t\t\t\t\t" .  join("", @aa[0..56]) . "\tNo uncharged region >=13 aa\n";}
	} else { print OUTFILE "\t\t\t\t\t\t" .  join("", @aa[0..56]) . "\tTAT signal is out of range (2..35)\n"; }
    } else { print OUTFILE "\t\t\t\t\t\t\t\t". join("", @aa[0..56]) . "\tNo TAT signal found\n"; }
}
