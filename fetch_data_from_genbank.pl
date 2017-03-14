#!/usr/bin/perl
use LWP::Simple;
use Getopt::Long;
use HTML::Entities;
use strict;

my $db_type_mode = {'all' => { 'docsum' => ['xml'],
			       'uilist' => ['xml', 'text'] },
		    'bioproject' => { 'xml' => ['xml'] },
		    'biosample' => { 'full' => ['xml', 'text'] },
		    'biosystems' => { 'xml' => ['xml'] },
		    'nucest' => { '' => ['text', 'asn.1'],
				  'native' => ['xml'],
				  'acc' => ['text'],
				  'fasta' => ['text', 'xml'],
				  'seqid' => ['text'],
				  'gb' => ['text', 'xml'],
				  'gbc' => ['xml'],
				  'est' => ['text'] },
		    'gene' => { '' => ['asn.1', 'xml'],
				'gene_table' => ['text'] },
		    'gds' => { 'summary' => ['text'] },
		    'nucgss' => { '' => ['text', 'asn.1'],
				  'native' => ['xml'],
				  'acc' => ['text'],
				  'fasta' => ['text', 'xml'],
				  'seqid' => ['text'],
				  'gb' => ['text', 'xml'],
				  'gbc' => ['xml'],
				  'gss' => ['text'] },
		    'homologene' => { '' => ['asn.1', 'xml'],
				      'alignmentscores' => ['text'],
				      'fasta' => ['text'],
				      'homologene' => ['text'] },
		    'mesh' => { 'full' => ['text'] },
		    'nlmcatalog' => { '' => ['text', 'xml'] },
		    'nuccore' => { '' => ['text', 'asn.1'],
				   'native' => ['xml'],
				   'acc' => ['text'],
				   'fasta' => ['text', 'xml'],
				   'seqid' => ['text'],
				   'gb' => ['text', 'xml'],
				   'gbc' => ['xml'],
				   'ft' => ['text'],
				   'gbwithparts' => ['text'],
				   'fasta_cds_na' => ['text'],
				   'fasta_cds_aa' => ['text'] },
		    'popset' => { '' => ['text', 'asn.1'],
				  'native' => ['xml'],
				  'acc' => ['text'],
				  'fasta' => ['text', 'xml'],
				  'seqid' => ['text'],
				  'gb' => ['text', 'xml'],
				  'gbc' => ['xml'] },
		    'protein' => { '' => ['text', 'asn.1'],
				   'native' => ['xml'],
				   'acc' => ['text'],
				   'fasta' => ['text', 'xml'],
				   'seqid' => ['text'],
				   'ft' => ['text'],
				   'gp' => ['text', 'xml'],
				   'gpc' => ['xml'],
				   'ipg' => ['xml'] },
		    'pubmed' => { '' => ['asn.1', 'xml'],
				  'medline' => ['text'],
				  'abstract' => ['text'] },
		    'pmc' => { '' => ['xml'],
			       'medline' => ['text'] },
		    'snp' => { '' => ['asn.1', 'xml'],
			       'flt' => ['text'],
			       'fasta' => ['text'],
			       'rsr' => ['text'],
			       'ssexemplar' => ['text'],
			       'chr' => ['text'],
			       'docset' => ['text'],
			       'sra' => 'full' => ['xml'],
			       'taxonomy' => '' => ['xml'] },
};

my (@id, @acc, @term, $from_db, $to_db, $rettype, $retmode);
&GetOptions("id=i" => \@id,
	   "acc=s" => \@acc,
	   "term=s" => \@term,
	   "dbsearch=s" => \$from_db,
	   "dbdata=s" => \$to_db,
	   "rettype=s" => \$rettype,
	   "retmode=s" => \$retmode);

until (defined $db_type_mode->{'all'}->{$rettype} ||
       defined $db_type_mode->{$to_db}->{$rettype}) {
    print "Valid rettypes for $to_db are:\n";
    print join("\n", (keys(%{$db_type_mode->{'all'}}), keys(%{$db_type_mode->{$to_db}}))) . "\n";
    print "Enter a type: ";
    $rettype = <>;
    chomp $rettype;
}

until ((defined $db_type_mode->{'all'}->{$rettype}
	&& grep $retmode,@{$db_type_mode->{'all'}->{$rettype}}) ||
       (defined $db_type_mode->{$to_db}->{$rettype}
	&& grep $retmode,@{$db_type_mode->{$to_db}->{$rettype}})) {
    print "Valid retmodes for $to_db:$rettype are:\n";
    print join("\n", (@{$db_type_mode->{'all'}->{$rettype}}, @{$db_type_mode->{$to_db}->{$rettype}})) . "\n";
    print "Enter a mode: ";
    $retmode = <>;
    chomp $retmode;
}

my $base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

if (@term) {
    my $term = encode_entities(join("+",@term));
    my $search_url = $base . "esearch.fcgi?db=$from_db&term=$term&usehistory=y";
    my $result;
    $result = get($search_url);
    if (! $result) { die }

    # parse WebEnv and QueryKey
    my $web = $1 if ($result =~ /<WebEnv>(\S+)<\/WebEnv>/);
    my $key = $1 if ($result =~ /<QueryKey>(\d+)<\/QueryKey>/);

    my $link_url = $base . "elink.fcgi?dbfrom=$from_db&db=$to_db&WebEnv=$web&query_key=$key&usehistory=y";
    my $links;
    $links = get($link_url);
    if (!$links) { die }

    # parse ids
    while ($links =~ /<Id>(\d+)<\/Id>/g) {
	push @id, $1;
    }
    my $idstr = join(",", @id);
    my $fetch_url = $base . "efetch.fcgi?db=$to_db&id=$idstr";
    my $out3 = get($fetch_url);

    my $outfile = $from_db . "." . $term . "." . $to_db . "." . $rettype . "." . $retmode;
    open(OUT, ">$outfile");
    print OUT $out3;
    close OUT;
}

#$post_url = "epost.fcgi?db=databaseB&id=id_list2";
#$search_url = "esearch.fcgi?db=databaseB&term=query+AND+%23key&WebEnv=Webenv&usehistory=y";
#$fetch_url = "efetch.fcgi?db=databaseB&WebEnv=Webenv2&query_key=key2";
exit();
