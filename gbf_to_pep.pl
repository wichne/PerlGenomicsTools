#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use Bio::Seq;
#use Getopt::Std;
use Data::Dumper;
use strict;
use lib "/share/scripts/";
use ENV;
use vars qw ($dbh %FEAT_NAME);

my $args = {};
#&getopts('i:F:', $args);
my $usage = "ls *.gbf | gbf_to_pep.pl\n";
if ($ARGV[0] eq "-h") {
    die $usage;
}
print "@ARGV\n";

#my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
#my $format = $args->{'F'} ? $args->{'F'} : undef;
while (my $filename = <STDIN>) {
    chomp $filename;
    my $outfile = $filename;
    $outfile =~ s/\.\w+$/\.pep/;
    my $pepout = Bio::SeqIO->new(-file => ">$outfile", -format => 'fasta');

    my $in = Bio::SeqIO->new(-file => $filename);
    
    while (my $seqo = $in->next_seq) {
	my @features = $seqo->get_SeqFeatures();
	# the first feature in a genbank file is the sequence itself
	my $seqObj = shift @features;

	my $seq_id = $seqObj->seq_id;
	if (! defined $seq_id) { warn "No seq_id found in object"; next;}
	
	foreach my $featObj (@features) {
	    my ($prot, $accession, $description);
	    if ($featObj->primary_tag eq "CDS") {
		my @all_tags = $featObj->get_all_tags;
		if (grep /\btranslation\b/, @all_tags) {
		    ($prot) = $featObj->get_tag_values("translation");
		} else {
		    $prot = $featObj->seq->translate->seq unless (grep /\bpseudo\b/, @all_tags);
		}

		if ($featObj->has_tag('locus_tag')) {
		    $accession = [$featObj->get_tag_values('locus_tag')]->[0];
		}
		if ($featObj->has_tag('product')) {
		    my ($p) = $featObj->get_tag_values("product");
		    $description .= "product=$p;";
		} 
		if ($featObj->has_tag('gene')) {
		    my ($g) = $featObj->get_tag_values("gene");
		    $description .= "gene=$g;";
		    if (!$accession) { $accession = $g }
		} 
		if ($featObj->has_tag('EC_number')) {
		    my ($p) = $featObj->get_tag_values("EC_number");
		    $description .= "EC_number=$p;";
		} 

		if ($prot && $accession) {
		    my $prot_obj = Bio::Seq->new(-display_id => $accession,
						 -description => $description,
						 -seq => $prot);
		    $pepout->write_seq($prot_obj);
		} else { warn "CDS at " . $featObj->location->to_FTstring . " is missing an accession\n"; }
	    }
	}
    }
}

exit();

