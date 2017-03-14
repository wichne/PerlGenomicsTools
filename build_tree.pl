#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use strict;
use Getopt::Std;
use Bio::DB::GenBank;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Tools::Run::Phylo::Phylip::ProtPars;
use Bio::TreeIO;
use Bio::AlignIO;

my %arg;
# p is for the tmp file if you've already run this
# D is the database, P is the db password
# f is the input aa fasta file
# b is for the blast output

&getopts('f:b:D:P:p:h:m:t:', \%arg);

if ($arg{h}) {
    print "build_tree.pl -f input_faa -D db -P dbpass -b btab [ -p tmp_file -m alignmentfile ]

Run a blast of your-favorite-protein(s) against some database. Use tab_blast on the blast output. This is the btab file.
If you've run this program previously, you can input the .temp file and bypass the sequence retrieval.
";
    exit();
}

my $outfilebase;
if ($arg{'m'}) {
    $arg{'m'} =~ /(.*)\.\w+$/;
    $outfilebase = $1 ? $1 : $arg{'m'};
} elsif ($arg{'f'}) {
    $arg{'f'} =~ /(.*)\.\w+$/;
    $outfilebase = $1 ? $1 : $arg{'f'};
} elsif ($arg{'b'}) {
    $arg{'b'} =~ /(.*)\.\w+$/;
    $outfilebase = $1 ? $1 : $arg{'b'};
}
if (! $outfilebase) { $outfilebase = $$ }
$outfilebase;

my $alnfile = $arg{'m'} ? $arg{'m'} : "";
my ($aln, $anonaln, $anonref);
my $treefile = $arg{'t'} ? $arg{'t'} : "";
if (! $treefile) {
    if (! $alnfile) {
	# make the alignment
	my $tmp_file;
	if ($arg{p}) {
	    $tmp_file = $arg{p};
	} else {
	    $tmp_file = "$outfilebase.fa";
	    my $outfo = Bio::SeqIO->new(-file => ">$tmp_file",
					-format => 'fasta');
	    &get_seqs_from_fasta($arg{'f'}, $outfo) if ($arg{'f'});
	    &get_seqs_from_nr($arg{'b'}, $outfo) if ($arg{'b'});
	    $outfo->close;
	}

	# Now build the multiple sequence alignment
	# and the NJ tree
	my @params = ();
	my $factory = Bio::Tools::Run::Alignment::Muscle->new(@params);
	my $aln = $factory->run($tmp_file);
	$alnfile = $outfilebase . ".afa";
	my $alnfo = Bio::AlignIO->new(-file => ">$alnfile",
				      -format => "fasta");
	$alnfo->write_aln($aln);
    } else {
	my $alnfo = Bio::AlignIO->new(-file => $alnfile);
	$aln = $alnfo->next_aln
    }
    # fol-de-rol to make phylip work
    my($anonaln, $anonref) = $aln->set_displayname_safe();
    my $tmpalnfile = $alnfile;
    my $tmpalnfo = Bio::AlignIO->new(-file => ">$tmpalnfile",
    my $alnfo = Bio::AlignIO->new(-file => ">alnfile",
                      -format => 'phylip');
    $alnfo->write_aln($aln);

    my $sub_model = "PROTGAMMAJTT";
    my $bootstrap = 100;
    my $rapidBootstrapRandomSeed = 12345;
    my $parsimonyRandomSeed = 12345;
    my $phylip_algorithm = "a";
    my $cpus = 15;
    system("raxml -f $phylip_algorithm -x $rapidBootstrapRandomSeed -p $parsimonyRandomSeed -# $bootstrap -m $sub_model -T $cpus -n $outfilebase -s $alnfile");
    $treefile = "RAxML_bestTree.$outfilebase";
    my $tree_factory = Bio::Tools::Run::Phylo::Phylip::ProtPars->new();
    my $tree = $tree_factory->run($anonaln);
} else {
    my $alnfo = Bio::AlignIO->new(-file => $alnfile);
    $aln = $alnfo->next_aln;
    ($anonaln, $anonref) = $aln->set_displayname_safe();
}

my $treefo = Bio::TreeIO->new(-file => $treefile);
my $tree = $treefo->next_tree;
foreach my $node_id (keys %$anonref) {
    my @node = $tree->find_node($node_id);
    if (@node == 1) {
	$node[0]->id($anonref->{$node_id});
    } else {
	warn "$node_id show up in two places in tree\n";
    }
}
    
my $treefo = Bio::TreeIO->new(-file => ">$outfilebase.nwk",
			      -format => 'newick');
$treefo->write_tree($tree);

#unlink $tmp_file;

exit();

sub get_seqs_from_fasta {	# this is the input fastafile of proteins to align
    my $file = shift;
    my $outfo = shift;
    
    my $seqfo = Bio::SeqIO->new(-file => $file);
    
    my $host = $ENV{DBSERVER};
    my $dbh = DBI->connect("dbi:mysql:host=$host;db=$arg{D}", $ENV{USER}, $arg{P});
	
    # Get the bin names for the sequences in the fasta file
    while (my $seqo = $seqfo->next_seq) {
	my $id = $seqo->display_name;
	my $desc = $seqo->description;
	my $bin;
	if ($desc =~ /bin\=(\S+)/) {
	    $bin = $1;
	} elsif ($id =~ /(mat\-\d+\_\d+)\_scaffold_(\d+)/) {
	    my $seqidq = "SELECT seq_id FROM sequence_accessions a WHERE a.seq_accession=\"$1\"";
	    my ($seq_id) = $dbh->selectrow_array($seqidq);
	    my $binq = "SELECT ss.name FROM sequence_sets ss, seq_set_link l, sequence_accessions a WHERE a.seq_accession=\"$1\" AND l.seq_id=a.seq_id AND ss.set_id=l.set_id";
	    ($bin) = $dbh->selectrow_array($binq);
	    $bin = "none" if (! $bin);
	    $seqo->display_id("${bin}_${seq_id}_$2");
	}
	$outfo->write_seq($seqo);
    }
}


sub get_seqs_from_nr {
    my $file = shift;
    my $outfo = shift;
    # connect to Genbank to retrieve blast hits?
    my $gb = Bio::DB::GenBank->new(-retrievaltype => 'tempfile' ,
				   -format => 'Fasta',
				   -db => 'protein');

    # Add the nr hits to the fasta file
    # Make sure they're non-redundant first
    open my $nr, $file;
    my %ACC;
    while (my $line = <$nr>) {
	next if ($line =~ /^\#/);
	my @f = split/\t/, $line;
	if ($f[4] =~ /\|(\w{2}\_\w+)/) {
	    $ACC{$1}++;
	} else {
#	$ACC{$f[4]}++;
	}
    }
    my @accs = keys %ACC;
    
    use LWP::Simple;
    my $rc = getstore("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=" . join(",", @accs),"$$.seqs.fasta");
    if (is_error($rc)) { die "$arg{f} $rc\n"; }
    my $ncbi_in = Bio::SeqIO->new(-file => "$$.seqs.fasta");
    while (my $inseq = $ncbi_in->next_seq) {
	my $id = $inseq->display_id;
	$id =~ /gi\|(\d+)/;
	my $new_id = $1;
	my $desc = $inseq->desc;
	$desc =~ /\[(\S+)/;
	$new_id .= "_$1";
	$inseq->primary_id($new_id);
	$inseq->display_id($new_id);
	$outfo->write_seq($inseq);
    }
}
