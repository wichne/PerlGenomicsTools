#!/usr/bin/perl

# input - fasta pep file
# output - tbl format hmmscan result and csv file with analysis

# run hmmscan on pep file
# parse tbl output

use Getopt::Std;
use strict;

my $arg = {};
&getopts('f:L:o:p:g:nh', $arg);

if ($$arg{h}) {
    print STDERR "Usage est_genome_completeness.pl [-f aa_fasta_file] [-g genome_fasta_file] [-L hmm_library] [-o hmmer_output_file] [-p processors] [-n (use hmmer_output_file instead of running hmmer search)]\n";
}

my $gen_file = $$arg{g};
my $pep_file = $$arg{f};
if ($gen_file && ! $pep_file) {
    $pep_file = $gen_file . ".pep";
    $pep_file =~ s/.*\///; # get rid of path
}
my $outtbl = $$arg{o} ? $$arg{o} : $pep_file . ".tblout";
$outtbl =~ s/.*\///; # get rid of path
my $hmm_lib = $$arg{L} ? $$arg{L} : "/common/db/hmms/CSCG/single_copy_genes.TIGRFAM.hmm";
my $ncpu = $$arg{p} ? $$arg{p} : 1;

# run prodigal to get predicted proteins
if ($gen_file && ! $pep_file) {
    my $tt = $$arg{t} ? $$arg{t} : 11;
    my $prodout = $gen_file . ".prodigal.gff";
    $prodout =~ s/.*\///; # get rid of path
    my $prod_cmd = "prodigal -a $pep_file -g $tt -i $gen_file -m -o $prodout -f gff";
    system($prod_cmd);
}

# run the hmm search
unless ($$arg{o} && $$arg{n}) {
    my $hmmcmd = "hmmscan --cut_nc --tblout $outtbl --cpu $ncpu $hmm_lib $pep_file";
    print STDERR $hmmcmd, "\n";
    
    my $rawout = `$hmmcmd`;
    $rawout =~ /Target model\(s\)\:\s+(\d+)/;
}

# parse the output
my @models;
my @descs;
my $models = `grep NAME $hmm_lib`;
while ($models =~ /NAME\s+(\S+)/g) {
    my $mod = $1;
    push @models, $mod;
}
my $descs = `grep DESC $hmm_lib`;
while ($descs =~ /DESC\s+(.+)/g) {
    my $desc = $1;
    if ($desc =~ /\:\s*(.*)/) { push @descs, $1 }
    else { push @descs, $desc }
}

my $n_mod = scalar(@models);

open my $in, $outtbl;
my %DATA;
while (my $l = <$in>) {
    chomp $l;
    if ($l =~ /^#/) { next }
    else {
	my ($target,
	    $hmmacc,
	    $query,
	    $qacc,
	    $tot_evalue,
	    $tot_score,
	    $tot_bias,
	    $dom_evalue,
	    $dom_score,
	    $dom_bias,
	    $exp,
	    $reg,
	    $clu,
	    $ov,
	    $env,
	    $dom,
	    $rep,
	    $inc,
	    $hmm_desc) = split(/\s+/, $l, 19);
	if ($tot_evalue <= 1e-10) {
	    $DATA{$target}++;
	}
    }
}

my $scg = scalar(keys %DATA);

printf "Input file: $gen_file $pep_file $outtbl\nNumber of SCG models searched: %d\nNumber of SCG models hit: %d\nEstimated completeness: %.4f\n", ($n_mod, $scg, $scg/$n_mod);

my $csv = $outtbl . ".scg";
open (my $out, ">$csv") or die;
select $out;

print "file\t", join("\t", @models), "\n";
print "\t", join("\t", @descs), "\n";
print "$pep_file";
foreach my $scg (@models) {
    print "\t" . ($DATA{$scg} * 1);
}
print "\n";
close $out;
select STDOUT;
exit();
