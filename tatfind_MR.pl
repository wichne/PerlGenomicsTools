#!/usr/bin/perl
#program = TATFIND.pl
#Usage:  TATFIND.pl inputfile 
#Results are located in a file called inputfile.TATout
#This program reads a single text file (containing one or more protein sequences
#in fasta format and searches them using a series a rules to determine if any of
#them are likely to be secreted via the TAT pathway.  Putative candidates are 
#printed out along with a hydrophobicity score for the first 13 hydrophobic 
#residues of the hydrophobic stretch following the twin arginine residues.
# 
##NOTE: For this program to work correctly, the fasta formatted files MUST have at
#      least 60 characters per line.  If not, please reformat your data.
#
#COPYRIGHT Jessica C. Kissinger June 11, 2001
#

use strict;
use Bio::SeqIO;
use Getopt::Std;

if (@ARGV < 1) {
    print "usage: TATFIND.pl infile\n";
    die;
}

#################################################################################
## If you want to debug the script, change the debug value to "1"              ##
## If you want to supress the text message at the top of the results, change   ##
## the value of $info to "0".                                                  ##
#################################################################################

my $debug = 0;
my $info = 1;
 
#################################################################################
## Set up a hash with the hydrophobicity values                                ##
## Values are taken from: Cid et al., (1992) Hydrophobicity and structural     ##
## classes in proteins. Protein Engineering 5:373-375                          ##
#################################################################################

my %hydrophobicity_table = (A => '0.02',
			    R => '-0.42',
			    N => '-0.77',
			    D => '-1.04',
			    C => '0.77',
			    Q => '-1.10', 
			    E => '-1.14',
			    G => '-0.80', 
			    H => '0.26',
			    I => '1.81',
			    L => '1.14',
			    K => '-0.41',
			    M => '1.00', 
			    F => '1.35',
			    P => '-0.09', 
			    S => '-0.97', 
			    T => '-0.77',
			    W => '1.71',
			    Y => '1.11',
			    V => '1.13');

my $infile = $ARGV[0];
open OUTFILE, ">$infile.TATMRout";
print OUTFILE "TATFIND - version 1.2 January, 2002; Brueser, Rose, Kissinger & Pohlschroder\n\nSEARCH RESULTS for $infile\:\n\n" if ($info);  

#################################################################################
## Redefine the record separator so that we can read in complete sequences     ##
## Note: the sequences must be in fasta format                                 ##
## Grab the first 44 amino acids of each sequence                              ##
#################################################################################

my $sfo = Bio::SeqIO->new(-file => $infile);

while (my $rec = $sfo->next_seq) {
    my @aa = split("", $rec->seq);
    my $header = $rec->display_id . " " . $rec->desc;
    print OUTFILE "\nResults for $header\:\n" if ($debug);
    
#################################################################################
## Begin a series of tests "rules" for TAT discovery                           ##
## Rule #1  Can we find the Twin arginine motif in the first 44aa?             ##
##          If we can, locate and grab the 22aa upstream following             ##
##          the first arginine residue.                                        ##
#################################################################################

    my $first44 = join("", @aa[0..43]);
    print OUTFILE "The first 44 amino acids are $first44\n" if ($debug);
    if ($first44 =~/([HAKRNTGSDQECY][KR]R[HAKRNTGSDQEV][IFLVMAT][IFLVMAT])/ig){ 
	my $endRRxxx = pos ($first44);
	print OUTFILE "The twinarg pattern has been found in the above, it is $1\n" if ($debug);
	print OUTFILE "The Pattern RR*** ends at position $endRRxxx\n" if ($debug);
	my $beginRR = $endRRxxx - 4; #actually this omits the first "R"
	my $next22 = join("", @aa[$beginRR..($beginRR+21)]);
	print OUTFILE "The 22 AA following the first R are $next22\n" if ($debug);
	
#################################################################################
## Rule #2  Is there a hydrophobic stretch of at least 13aa within 21aa        ##
##          immediately following the "RR"?                                    ##
##          If so, gather up information for rules #3a, #3b and #4             ##
#################################################################################
	
        if ($next22 =~/([^RKDE]{13,})/ig) { 
	    my $endhydro = pos ($next22);
	    my $length = length $1;
	    my $hydro = $1;
	    print OUTFILE "The hydrophobic region is $hydro\n" if ($debug);
	    print OUTFILE "The length of the hydrophobic region is $length aa\n" if ($debug);

	    $hydro =~ /(.{13})/;
	    my $hydro13 = $1;
	    my $sum = 0;
	    @13 = split ("", $1);
	    foreach (@13) {
		if (defined ($hydrophobicity_table{$_})){
		    $sum+=$hydrophobicity_table{$_};
		}
	    }
	    
	    my $beginhydro = $endhydro - $length +1;
	    print OUTFILE "The hydrophobic region begins at $beginhydro\n" if ($debug);
	
#################################################################################
##  Rule #3a Special case:  TwinArg:single charged residue:hydrophobic region  ## 
##           check for: RR[DERK][hydrophobic stretch begins immediately]       ##
##           Tally up the hydrophobicity of the first 13 hydrophobic AA's      ##
#################################################################################
	    
	    if ($next22 =~/^R[DERK]$hydro/ig){
		print OUTFILE $rec->display_id . "\tRule3a\t$sum\t$hydro13\n$first44\n" if ($sum < 8.0);
		print OUTFILE "The 13 hydrophobic residues to be summed are: $hydro13\n" if ($debug);
		print OUTFILE "Rule 3a applies, here is the complete coding sequence:\n\n" if ($debug);
	       
#		print OUTFILE ">$header has a hydrophobicity score of $sum \n$first44\n\/\/\n" if ($sum < 8.0);
		
#################################################################################
## Rule #3b Special case: RR...[DERK][hydrophobic region is here]              ##
##          Tally up the hydrophobicity of the first 13 hydrophobic residues   ##       
#################################################################################
		
	    } elsif ($next22 =~/R...[DERK]$hydro/i){
		print OUTFILE $rec->display_id . "\tRule3b\t$sum\t$hydro13\n$first44\n" if ($sum < 8.0);
		print OUTFILE "The 13 hydrophobic residues to be summed are: $1\n" if ($debug);
		print OUTFILE "Rule 3b applies, here is the complete coding sequence:\n\n" if ($debug);
#		print OUTFILE ">$header has a hydrophobicity score of $sum \n$first44\n\/\/\n" if ($sum < 8.0);
		
		
#################################################################################
## Rule #4  If rule #3 doesn't apply, check for the presence of a basic        ##
##          residue immediately preceeding the hydrophbic region.  If not      ##
##          present, toss it out, otherwise we have a candidate.               ##
##          Tally up the hydrophobicity of the first 13 hydrophobic AA's       ##
#################################################################################
		
	    } elsif ($next22 !~/[KR][^DERK]{13,}/i){
		print OUTFILE "Rules 3and 4 didn't apply, search aborted.\n\/\/\n" if ($debug);
		next;
	    } else {
		print OUTFILE $rec->display_id . "\tRule4\t$sum\t$hydro13\n$first44\n" if ($sum < 8.0);
		print OUTFILE "The 13 hydrophobic residues to be summed are: $hydro13\n" if ($debug);
		print OUTFILE "Rule 4 applies, here is the complete coding sequence:\n\n" if ($debug);
#		print OUTFILE ">$header has a hydrophobicity score of $sum \n$line1\n$rest\n\/\/\n" if ($sum < 8.0);
	    }
	} else {
	    print OUTFILE "No hydrophobic stretch following the TAT site, search aborted.\n\/\/\n" if ($debug);
	}
    } else {
	print OUTFILE "No TAT recognition site in the first 44 aa, search aborted.\n\/\/\n" if ($debug);
    }
}
close OUTFILE;
print "All Done!  The results are located in $infile.out\n" if ($debug);











