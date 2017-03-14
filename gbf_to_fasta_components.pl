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

#my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
#my $format = $args->{'F'} ? $args->{'F'} : undef;
while (my $filename = <STDIN>) {
    print STDERR "Doing $filename...\n";
    chomp $filename;
    my $outfile = $filename;
    my $seqfile = $filename;
    my $nafile = $filename;
    $outfile =~ s/\.\w+$/\.pep/;
    $seqfile =~ s/\.\w+$/\.seq/;
    $nafile =~ s/\.\w+$/\.fasta/;
    my $out = Bio::SeqIO->new(-file => ">$outfile", -format => 'fasta');
    my $seqout = Bio::SeqIO->new(-file => ">$seqfile", -format => 'fasta');
    my $naout = Bio::SeqIO->new(-file => ">$nafile", -format => 'fasta');

    my $in = Bio::SeqIO->new(-file => $filename);

    while (my $seqo = $in->next_seq) {
	my @features = $seqo->get_SeqFeatures();
	$naout->write_seq($seqo);
	# the first feature in a genbank file is the sequence itself
	my $seqObj = shift @features;
	
	my $seq_id = $seqObj->seq_id;
	if (! defined $seq_id) { warn "No seq_id found in object"; next;}
	
	foreach my $featObj (@features) {
	    my ($prot, $accession, $description);
	    $seqout->write_seq($featObj->seq);
	    if ($featObj->primary_tag eq "CDS") {
		my @all_tags = $featObj->get_all_tags;
		if (grep /\btranslation\b/, @all_tags) {
		    ($prot) = $featObj->get_tag_values("translation");
		} else {
		    $prot = $featObj->seq->translate->seq;
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
		    $out->write_seq($prot_obj);
		} else { warn "CDS at " . $featObj->location->to_FTstring . " is missing an accession\n"; }
	    }
	}
    }
}

exit();

