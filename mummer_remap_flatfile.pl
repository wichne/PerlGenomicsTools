#!/usr/bin/perl
=head1 NAME
    
    mummer_remap.dbi
    
=head1 USAGE
    
    mummer_remap_flatfile.pl -refdb <fastafile> -O <output file> -qrydb <fastafile> -W <window size> -mumout <delta file> -reffeat <gffFile> -qryfeat<gffFile>
    
=head1 OPTIONS

    REQUIRED
    -refdb      fasta of ref sequence (map from)
    -reffeat    gff of ref features
    -qrydb      fasta of query sequence (map to)
    -qryfeat    gff of qry features
    -O          file name where output will be written

    OPTIONAL
    -W          specify window size nucmer/promer should use (defaults to 20 for nucmer, 6 for promer)
    -mumout     specify a delta file to use as basis for comparision (suppresses default run of nucmer/promer)

=head1 DESCRIPTION

    Program is designed to use nucmer/promer output to compare two or more genomes (using one as a reference), and output feature-to-feature or feature-to-sequence mappings. By default a nucmer search will be performed on the specified assemblies.
    The resulting alignments are analyzed in order of descending length (thus longer regions of synteny are mapped first). Coordinates of reference features that overlap the alignment are transformed to coordinates on the query molecule, and then a query feature that matches one or both of those coordinates is sought and if found written to the output file. If the end5 or end3 (but not both) are matched, the mapped coordinate is presented in the output comment field. No greater than a 5% difference in length is allowed for two features to map. If there is more than 5% different in length, the mapping is commented out and reported as a 'partial'. Mappings are checked for length differences, and divisibility by 3. If a query feature overlaps the transformed coordinates, but does not match either coordinate (and is in the same direction), it is listed as a 'potential' mapping and commented out.
    Each query feature can only be mapped once; reference features can map to more than one query feature. Reference features that do not map to a query feature are not listed anywhere in the output.

=head1 INPUT

    The user can specify databases, assemblies and feat_types to be used in the comparison.
    If the user specifies a mummer output file, it must be a delta file, not a coords or cluster file.
    Lockfile formats are not compatible across operating systems, and may need to be regenerated.

=head1 OUTPUT

All query features appear in the output file.

Output file columns:
    0  query database
    1  query asmbl_id
    2  feat_type
    3  feat_name (or placeholder for feature-to-sequence mappings)
    4  query end5
    5  query end3
    6  reference database
    7  reference asmbl_id
    8  reference feat_name
    9  reference end5
   10  reference end3
   11  comment

Some example output:
A perfect feature-to-feature mapping:
ntbp04  38      ORF     ORF05718        5505838 5503964 gbp1710b        451     ORFE00982       2436438 2434564 

A mapping where the length differs, but is within tolerance (note the reported end5 in the comment):
ntbp04  38      ORF     ORF06755        6473200 6473078 gbp1710b        450     ORFF02674       2540130 2539948 mapped qry end5=6473260;length_diff;

A mapping where the length difference is outside of tolerance (partial):
#ntbp04 38      ORF     ORF02876        2744440 2744619 gbp1710b        450     ORFF01344       1123593 1126454 mapped qry end5=2744443;length_diff;partial mapping;

A potential mapping:
#ntbp04 38      ORF     ORF06645        6377728 6376895 gbp1710b        450     ORFF02570       2444033 2443698 potential;mapped qry end5=6377203;length_diff;

Feature-to-sequence mapping. Note the query feature is not divisible by 3 (q%3):
ntbp04  38      ORF     NEW00255        2957967 2957853 gbp1710b        450     ORFF00104       4009251 4009364 length_diff;q%3;

=head1  CONTACT

    Bill Nelson
    wnelson@tigr.org or
    sgc@tigr.org

=begin comment
    ## legal values for status are active, inactive, hidden, unstable
  status: active
  keywords: mummer nucmer promer remapping transfer
=end comment

=cut

use strict;
use Cwd;
use Cwd 'chdir';
use Getopt::Long;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Pod::Usage;
use lib $ENV{SCRIPTS};
use Mummer_remap;
use Mummer_delta_parser;
use Pod::Usage; 

# set up options
my ($refseq,
    $reffeat,
    $qryseq,
    $qryfeat,
    $min_word_length,
    $out_file,
    $deltafile,
    $qry_lock,
    $ref_lock,
    $DEBUG,
    $VERBOSE,
    $help,
    $qryfile);

&GetOptions(
    'debug'=>\$DEBUG,
    'refdb=s'=>\$refseq,
    'reffeat=s'=>\$reffeat,
    'qrydb=s'=>\$qryseq,
    'qryfeat=s'=>\$qryfeat,
    'qryfile=s'=>\$qryfile,
    'W=i'=>\$min_word_length,
    'O=s'=>\$out_file,
    'mumout=s'=>\$deltafile,
    'h' => \$help
    );

# get rid of path info
my $refdb = $refseq;
$refdb =~ s/.*\///;
$refdb =~ s/\.*//;

my $qrydb = $qryseq;
$qrydb =~ s/.*\///;
$qrydb =~ s/\.*//;

if (! $refdb || ! $out_file) { &pod2usage( {-exitval => 1, -verbose => 2, -output => \*STDOUT} ); exit;}

if (! $min_word_length) {
    $min_word_length = $qrydb eq $refdb ? 100 : 20 }

# set up executables
my $program = "nucmer";
my $MUMMER = "/usr/local/bin/nucmer";
if (! -x $MUMMER) { die "'$MUMMER' is not executable.\n"; }

my $percent_length_difference_allowed = 0.2;

# set up log file
my $logfile = "mummer_remap_$$.log";
open LOG, ">$logfile" or warn "Log file '$logfile' cannot be written: $!\n";
select LOG; $| = 1; select STDOUT;
my $time = localtime(time);
print LOG "START: $time\n";

# set up output file
if ($out_file) {
    open OUT, ">$out_file" or die "Can't open '$out_file' for writing: $!\n";
    select OUT;
    $| = 1;
}

# Make structure that holds reference sequences
my $REFSEQ = {};
&get_sequences($refseq, $refdb, $REFSEQ);
# get features for reference sequences
my $REFFEAT = {};
&get_features($reffeat, $refdb, $REFFEAT);

# Make structure that holds query sequences
my $QRYSEQ = {};
my $QRYFEAT = {};
&get_sequences($qryseq, $qrydb, $QRYSEQ);
&get_features($qryfeat, $qrydb, $QRYFEAT);

# run nucmer(/promer?)
if (! $deltafile) {
    ($deltafile) = &run_alignment($refseq, $qryseq, $program);
}

# get alignments
my $remapper = new Mummer_remap($deltafile);
my $delta_parser = $remapper->get_delta_parser;
my @alignments = $delta_parser->get_alignments;

my $NEW = 0; # for holding new gene mappings
my %RESULTS;
my $MAPREF = {};

# go through alignments from longest to shortest
foreach my $alno (sort {$b->get_query_coord_span_length <=> $a->get_query_coord_span_length} @alignments) {
    my $refid = $alno->get_reference_accession;
    my ($ref_coord_lo, $ref_coord_hi) = $alno->get_reference_coords;

    my $qryid = $alno->get_query_accession;
    my ($qry_coord_lo, $qry_coord_hi) = $alno->get_query_coords;

    print LOG "Alignment: Reference $refid"
	. " ($ref_coord_lo-$ref_coord_hi/" . $REFSEQ->{$refdb}->{$refid}->length
	. ")  <===> Query $qryid"
	. " ($qry_coord_lo-$qry_coord_hi/" . $QRYSEQ->{$qrydb}->{$qryid}->length
	. ") (", $alno->get_query_coord_span_length, ")\n";

    # find ref genes overlapping this alignment and try to map to qry genes

    # Go through all the reference features
    foreach my $ref_fid (sort {$REFFEAT->{$refdb}->{$a}->start <=> $REFFEAT->{$refdb}->{$b}->start} keys %{$REFFEAT->{$refdb}}) {
	my $ref_feature = $REFFEAT->{$refdb}->{$ref_fid};
	# if the feature isn't on the reference sequence-of-interest, go to the next one
	if ($ref_feature->seq_id ne $refid) { next }

	# Is this feature already mapped? Skip it.
	if ($MAPREF->{$ref_fid}->{skip} == 1) { next }

	# Does this feature overlap the alignment region?
	# if it ends before this region or starts after this region, skip to the next one.
	if ($ref_feature->end < $ref_coord_lo) { next }
	if ($ref_feature->start > $ref_coord_hi) { last }

	# determine by how much this gene overlaps the alignment region
	my $ref_feat_len = $ref_feature->length;
	my $ol = &overlaps($ref_feature->start, $ref_feature->end, $ref_coord_lo, $ref_coord_hi);

	if ($ol) {
	    my $locus_tag;
	    if ($ref_feature->has_tag('locus_tag')) {
		$locus_tag = ($ref_feature->get_tag_values('locus_tag'))[0];
	    } else { $locus_tag = $ref_feature->display_name }

	    printf LOG "\treference %s (%i/%i %s) overlaps this alignment ($ol)\n", ($locus_tag, $ref_feature->start, $ref_feature->end, $ref_feature->strand);

	    # Transform the coordinates
	    my ($mapped_start,
		$mapped_end,
		$mapped_strand) = &map_coords($ref_feature->start,
					      $ref_feature->end,
					      $ref_feature->strand,
					      $remapper, $alno);
	    print LOG "\treference $locus_tag maps to query $mapped_start/$mapped_end on strand $mapped_strand\n";

	    # is this reference feature partial?
	    my $partial = 0; 
	    my $ref_partial = 0;
	    if ($ref_feature->start == 1 || $ref_feature->end == $REFSEQ->{$refdb}->{$refid}->length()) {
		print LOG "\t\treference gene $locus_tag may be partial (" . $ref_feature->start . "/" . $ref_feature->end . ")\n";
		$ref_partial = 1;
	    }
	    # is the mapping partial? if the mapped region is not the same length as the reference gene, then yes.
	    my $seq_mapping_partial = 0;
	    if (($mapped_end - $mapped_start) + 1 < $ref_feature->length) {
		print LOG "\t\tonly part of $locus_tag maps to $qryid\n";
		$seq_mapping_partial = 1;
	    }

	    # if there is a mappable feature, look for reportable differences
	    my $qry_partial = 0;
	    my ($qry_feature, $problem) = &find_feature_from_coords($QRYFEAT->{$qrydb}, $qryid, $mapped_start, $mapped_end, $mapped_strand, $ref_feature->primary_tag);
	    if (defined $qry_feature) {
		if ($qry_feature->start == 1 || $qry_feature->end == $QRYSEQ->{$qrydb}->{$qryid}) {
		    print LOG "\t\tQuery feature may be partial\n";
		    $qry_partial = 1;
		}
		$problem .= $qry_feature->start == $mapped_start ? "" : "mapped qry start=$mapped_start;";
		$problem .= $qry_feature->end   == $mapped_end   ? "" : "mapped qry end=$mapped_end;";
	    }

	    # if there is not a mappable feature, instantiate a new one.
	    else {
		my $feat_name = sprintf "NEW%05d", (++$NEW);
		print LOG "\tmaking new feature $feat_name\n";
		$QRYFEAT->{$qrydb}->{$feat_name} = Bio::SeqFeature::Generic->new( 
		    -start        => $mapped_start, 
		    -end          => $mapped_end,
		    -strand       => $mapped_strand, 
		    -primary      => $ref_feature->primary_tag,
		    -display_name => $feat_name,
		    -seq_id       => $qryid,
		    -tag          => { locus_tag        => $locus_tag } );
		$qry_feature = $QRYFEAT->{$qrydb}->{$feat_name};
	    }

	    # are there any problems with the mapping?
	    &determine_problem($ref_feature, $qry_feature, \$problem);

	    # test length of mapping
	    my $gene_mapping_partial = 0;
	    if (abs($qry_feature->length() - $ref_feat_len)/$ref_feat_len > $percent_length_difference_allowed ||
		abs($qry_feature->length() - $ref_feat_len)/$qry_feature->length() > $percent_length_difference_allowed) {
		print LOG "\tmapping doesn't meet percent length difference criteria\n";
		$gene_mapping_partial = 1;
	    }

	    # does the mapped position coincide with either end of the molecule (probable partial)
	    if ($mapped_start == 1 || $mapped_end == 1 ||
		$mapped_start == $QRYSEQ->{$qrydb}->{$qryid}->length() ||
		$mapped_end == $QRYSEQ->{$qrydb}->{$qryid}->length()) {
		$partial = 2;
	    }

	    if ($ref_partial && !$qry_partial) {
		$problem .= "partial gene mapping to complete gene;";
	    } elsif (!$ref_partial && $qry_partial) {
		$problem .= "complete gene mapping to partial gene;";
	    }
	    if ($seq_mapping_partial) {
		$problem .= "partial mapping of ref gene to contig;";
	    }
	    if ($gene_mapping_partial) {
		$problem .= "partial mapping of ref gene to qry gene;";
	    }
	    
	    # is the feature a split feature? If so, map the individual parts and drop a
	    # NCBI style coordinate string in the comment.
	    if ($ref_feature->location->isa('Bio::Location::SplitLocationI') ) {
		my $coord_string= "join(";
		my @sublocs = $ref_feature->location->sub_Location();
		foreach my $location ( sort { $a->start <=> $b->start } @sublocs ) {
		    my ($start, $end, $strand) = &map_coords($location->start,
							     $location->end,
							     $location->strand,
							     $remapper,
							     $alno);
		    $coord_string .= "$start..$end,";
		}
		chop $coord_string;
		$coord_string .= ")";
		if ($ref_feature->location->guide_strand < 0) { $coord_string = "complement($coord_string)"; }
		$problem .= "$coord_string;";
	    }

	    # Is the query gene already mapped to something else and
	    # is the previous mapping better?
	    # or is it a previously identified gene-to-sequence mapping?
	    # Log this secondary mapping, but don't report it.
	    if (defined $RESULTS{$qrydb}->{$qry_feature->display_name} &&
		($problem =~ /partial/ && $RESULTS{$qrydb}->{$qry_feature->display_name}->{problem} !~ /(potential|partial)/ ||
		 $qry_feature->display_name =~ /NEW/)) {
		print LOG "\t" . $qry_feature->display_name . " secondary mapping"
		    . " (first: " . $RESULTS{$qrydb}->{$qry_feature->display_name}->{feature} . ")\n";
	    }
	    else {
		# this frees the previous mapped reference gene to map to some other gene
		if (defined $RESULTS{$qrydb}->{$qry_feature->display_name}) {
		    $MAPREF->{$RESULTS{$qrydb}->{$qry_feature->display_name}->{feature}}->{skip} = 0 if ($qrydb eq $refdb); }
		
		# this populates the mapping in the results hash
		$RESULTS{$qrydb}->{$qry_feature->display_name} = {'refdb' => $refdb,
								  'feature' => $ref_fid,
								  'problem' => $problem,
		};
		$MAPREF->{$ref_feature->display_name}->{skip} = 1; # if ($qrydb eq $refdb); # to ensure ref genes only map once.
	    }
	}
    }
}

print OUT "toDataSource\ttoSeqAcc\tfeatType\tfeatAcc\ttoStart\ttoEnd\ttoStrand\tfromDataSource\tfromSeqAcc\tfromFeatType\tfromFeatAcc\tfromStart\tfromEnd\tfromStrand\n";
while (my ($db, $dbref) = each %$QRYFEAT) {
    while (my ($feat_name, $qfeatref) = each %$dbref) {
	my $qseqid = $qfeatref->seq_id;
	my $out = sprintf "%s\t%s\t%s\t%s\t%d\t%d\t%d",
	($db,$qseqid,$qfeatref->primary_tag,$feat_name,$qfeatref->start,$qfeatref->end,$qfeatref->strand);
	if (defined $RESULTS{$db}->{$feat_name}) {
	    my $rfeatref = $REFFEAT->{$RESULTS{$db}->{$feat_name}->{refdb}}->{$RESULTS{$db}->{$feat_name}->{feature}};
	    my $rseqid = $rfeatref->seq_id;
	    $out .= sprintf "\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\n",
	    ($refdb,$rseqid,$rfeatref->primary_tag,$RESULTS{$db}->{$feat_name}->{feature},
	     $rfeatref->start,$rfeatref->end,$rfeatref->strand, $RESULTS{$db}->{$feat_name}->{problem});
#	    $out .= "\t$refdb\t$rseqid\t$rfeatref->{feat_name}"
#		. "\t$rfeatref->{end5}\t$rfeatref->{end3}"
#		. "\t$RESULTS{$db}->{$feat_name}->{problem}\n";
	    if ($RESULTS{$db}->{$feat_name}->{problem} =~ /(potential|partial|secondary)/) { $out = "#" . $out }
#	    elsif ($qfeatref->{complete} ne "") { 
#		    chomp $out;
#		    $out = "#" . $out . "completed;\n"; 
#		}
	} else { $out .= "\n" }
	print OUT $out;
    }
}
exit();

sub get_sequences {
    my $file = shift;
    my $db = shift;
    my $ref = shift;

    print LOG "Parsing $file for sequences.\n";
    my $fileo=Bio::SeqIO->new(-file => $file);
    while (my $seqo = $fileo->next_seq) {
	print LOG "\tAdding " . $seqo->display_name . "\n";
	$ref->{$db}->{$seqo->display_name} = $seqo;
    }
}

sub get_features {
    my $file = shift;
    my $db = shift;
    my $ref = shift;

    print LOG "Parsing $file for features.\n";
    my $gffo = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);
    my $feature_count;
    while (my $featureo = $gffo->next_feature) {
	$feature_count++;
	my $identifier;
	if ($featureo->has_tag('locus_tag')) {
	    my @locus = $featureo->get_tag_values('locus_tag');
	    $identifier = $locus[0];
	    $featureo->display_name($identifier);
	} elsif ($featureo->display_name) {
	    $identifier = $featureo->display_name;
	} else {
	    $identifier = sprintf "%s_%05d", ($featureo->{_gsf_seq_id}, $feature_count);
	}
	if (defined $ref->{$db}->{$identifier}) {
	    # This is a split location feature
	    if ( $ref->{$db}->{$identifier}->location->isa('Bio::Location::SplitLocationI') ) {
		# This is already a split location object, just add the new sublocation
		$ref->{$db}->{$identifier}->location->add_sub_Location(Bio::Location::Simple->new(-start=>$featureo->start, -end=>$featureo->end, -strand=>$featureo->strand));
	    } else {
		# Let's assume we need to convert current feature location object to split.
		my $newSplitLocation = Bio::Location::Split->new();
		$newSplitLocation->splittype('join') if ($featureo->primary_tag eq "CDS");
		$newSplitLocation->guide_strand($ref->{$db}->{$identifier}->location->strand);
		$newSplitLocation->strand($ref->{$db}->{$identifier}->location->strand);
		# add the existing Location
		$newSplitLocation->add_sub_Location($ref->{$db}->{$identifier}->location);
		# add the current Location
		$newSplitLocation->add_sub_Location(Bio::Location::Simple->new(-start=>$featureo->start, -end=>$featureo->end, -strand=>$featureo->strand));
		# reset the feature Location object to the new split Location
		$ref->{$db}->{$identifier}->location($newSplitLocation);
	    }
	} else {
	    $ref->{$db}->{$identifier} = $featureo;
	}
    }
    $gffo->close;
    print LOG "$feature_count rows in file\n";
    print LOG scalar(keys %{$ref->{$db}}) . " features loaded\n";
}

sub run_alignment {
    my $ref_fa = shift;
    my $qry_fa = shift;
    my $program = shift;

    # make reference fasta file
#    my $ref_fa = $$ . "_ref_assembly.fa";
#    &write_fasta_from_structure($REFSEQ, $ref_fa);
    # make query fasta file
#    my $qry_fa = $$ . "_qry_assembly.fa";
#    &write_fasta_from_structure($QRYSEQ, $qry_fa);

    my $prefix = $$ . "_nucmer";
    my $mum_cmd = "$MUMMER -l $min_word_length -p $prefix $ref_fa $qry_fa";
    print LOG "Running mummer:\n$mum_cmd\n";
    system($mum_cmd);
    my $deltafile = $prefix . ".delta";
    if (! -e $deltafile) { die "Error running '$mum_cmd': Can't find delta file.\n";}
    return $deltafile;
}

sub write_fasta_from_structure {
    my $ref = shift;
    my $filename = shift;
    my $fileo = new TIGR::FASTAwriter ($filename);
    while (my ($db, $idref) = each %$ref) {
	while (my ($id, $fao) = each %$idref) {
	    $fileo->write($fao);
	}
    }
    $fileo->close;
}

sub find_feature_from_coords {
    # In this subroutine, we try to find a feature based on coords
    # look at each feature on molecule and determine if
    # 1) it's in the right direction
    # 2) it's end5 or end3 matches the query coords
    # 3) it overlaps the query coords (and is in the right frame)
    # The ref passed in only includes features on the aligned molecule
    my $ref = shift;
    my $seqid = shift;
    my $start = shift;
    my $end = shift;
    my $strand = shift;
    my $feat_type = shift;
    
    if (! ($start && $end && $strand)) { warn "Missing item '$start', '$end', '$strand'\n";}

    my $counter = 0;
    my $potm;
    my $potm_ol;
    foreach my $feat_name( sort keys %$ref) {
	my $feat_ref = $ref->{$feat_name};
	if ($feat_ref->seq_id ne $seqid) { next } # not on the same molecule
	if ($feat_ref->primary_tag ne $feat_type) { next } # not the same feat_type
	# is it going in the right direction?
	if ($strand * $feat_ref->strand >= 0 ) {
	    # Do either of the coords match?
	    if ($start == $feat_ref->start ||
		$end == $feat_ref->end) {
		print LOG "\tquery $feat_name (" . $feat_ref->start . "/" . $feat_ref->end . ") matches start and/or end and direction\n";
		return $feat_ref;
	    }
	    # otherwise, we need to check overlap and frame
	    else {
		my $ol = &overlaps($start, $end, $feat_ref->start, $feat_ref->end);
		my $frame = $start % 3 == $feat_ref->start % 3 || $end % 3 == $feat_ref->end % 3 ? 1 : 0;
		if ($ol > 0 && $frame) {
		    print LOG "\tquery $feat_name overlaps in same frame.\n";
		    return $feat_ref;
		}
	    }
	}
    }

    print LOG "\treference $start/$end $strand didn't map to an existing query feature\n";
    return undef;
}

sub determine_problem {
    my $ref = shift;
    my $qry = shift;
    my $problemr = shift;

    # are genes the same length?
    my $ref_len = abs($ref->{end5} - $ref->{end3}) + 1;
    my $qry_len = abs($qry->{end5} - $qry->{end3}) + 1;
    $$problemr .= "length_diff;" if ($ref_len != $qry_len);

    # are both genes divisible by 3
    if ($ref_len % 3 xor $qry_len % 3) {
	$$problemr .= "r\%3;" if ($ref_len % 3);
	$$problemr .= "q\%3;" if ($qry_len % 3);
    }
}

sub overlaps {
    my @data = @_;
    my @sorted = sort {$a<=>$b} @data;
    my $ol = abs($data[0] - $data[1]) + 1 + abs($data[2] - $data[3]) + 1 - ($sorted[3] - $sorted[0] + 1);
    return $ol;
}

if ($out_file) { close OUT }

sub map_coords {
    my ($start,
	$end,
	$strand,
	$remapper,
	$alno) = @_;

    my ($mapped_start, $mapped_end, $mapped_strand);
    # we do not assume the actual end5 and end3 will map
    # so we step in from either end of gene to find extent of mapping
    # first build array of coords from end5 to end3
    my @coords = $start .. $end;
    # step through this until end5 maps
    foreach my $coord (@coords) {
	$mapped_start = $remapper->transform_coordinate_within_single_alignment($Mummer_remap::REFERENCE_TO_QUERY, $coord, $alno);
	if ($mapped_start) { last }
    }
    # step through backwards until end3 maps
    foreach my $coord(reverse @coords) {
	$mapped_end = $remapper->transform_coordinate_within_single_alignment($Mummer_remap::REFERENCE_TO_QUERY, $coord, $alno);
	if ($mapped_end) { last }
    }
    if (! $mapped_start) { print LOG "\tNo end5 mappings\n" }
    if (! $mapped_end) { print LOG "\tNo end3 mappings\n" }
    my $mapped_strand;
    if (($strand == 1 && $mapped_start < $mapped_end) ||
	($strand == -1 && $mapped_start > $mapped_end)) { $mapped_strand = 1 }
    elsif (($strand == 1 && $mapped_start > $mapped_end) ||
	   ($strand == -1 && $mapped_start < $mapped_end)) { $mapped_strand = -1 }
    else { $mapped_strand = 0 }
    ($mapped_start, $mapped_end) = $mapped_start < $mapped_end ? ($mapped_start, $mapped_end) : ($mapped_end, $mapped_start);
    return ($mapped_start, $mapped_end, $mapped_strand);
}

exit(0);
