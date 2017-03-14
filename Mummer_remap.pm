#!/usr/local/bin/perl

package Mummer_remap;

use strict;
use warnings;
use Carp;

use Mummer_delta_parser;


our $VERBOSE = 0;
our $DEBUG = 0;

# constants
our $QUERY_TO_REFERENCE = 1;
our $REFERENCE_TO_QUERY = 2;


sub new {
    my $packagename = shift;
    my $delta_file = shift;

    my $delta_parser = new Mummer_delta_parser($delta_file);


    my $self = {
        accessions_to_alignments => {}, ## reference_acc -> query_acc -> [list of alignments]
       
        query_accession_mappings => {}, ## query_acc -> [ list of reference accs]
        reference_accession_mappings => {}, ## reference_acc -> [ list of query accs ]
        
        delta_parser => $delta_parser,
    };
    
    bless ($self, $packagename);

    
    $self->_init();

    return ($self);

}

####
sub _init {
    my $self = shift;
    
    my $delta_parser = $self->get_delta_parser();
    
    my $accessions_to_alignments_href = $self->{accessions_to_alignments};

    my @alignments = $delta_parser->get_alignments();
    
    my %ref_to_query_accs;
    my %query_to_ref_accs;

    my @alignment_lists;

    foreach my $alignment (@alignments) {
        my $reference_accession = $alignment->get_reference_accession();
        my $query_accession = $alignment->get_query_accession();

        my $alignment_list_aref = $accessions_to_alignments_href->{$reference_accession}->{$query_accession};
        unless (ref $alignment_list_aref) {
            $alignment_list_aref = $accessions_to_alignments_href->{$reference_accession}->{$query_accession} = [];
            
            ## acc accession mappings here:
            $ref_to_query_accs{$reference_accession}->{$query_accession} = 1;
            $query_to_ref_accs{$query_accession}->{$reference_accession} = 1;
        
            push (@alignment_lists, $alignment_list_aref);
        }

        push (@$alignment_list_aref, $alignment);
    }
    
    ## populate the accession mappings as key to list:
    foreach my $ref_acc (keys %ref_to_query_accs) {
        my $queries_href = $ref_to_query_accs{$ref_acc};
        my @query_accs = keys %$queries_href;
        
        $self->{reference_accession_mappings}->{$ref_acc} = [@query_accs];
    }

    foreach my $query_acc (keys %query_to_ref_accs) {
        my $ref_href = $query_to_ref_accs{$query_acc};
        my @ref_accs = keys %$ref_href;
        
        $self->{query_accession_mappings}->{$query_acc} = [@ref_accs];
    }
    
    ## sort all the alignment lists by alignment length:
    foreach my $alignment_list_aref (@alignment_lists) {
        @$alignment_list_aref = reverse sort { $a->{reference_coord_span_length} 
                                                         <=> 
                                               $b->{reference_coord_span_length} } @$alignment_list_aref;
    }
    

    return;
    
}

####
sub get_contig_mappings {
    ## returns the list of contig accessions that have alignments to the parameter contig.
    my $self = shift;
    my ($type, $contig) = @_;

    my $accs_href = undef;
    
    if ($type == $REFERENCE_TO_QUERY) {
        $accs_href = $self->{reference_accession_mappings};
    }
    elsif ($type == $QUERY_TO_REFERENCE) {
        $accs_href = $self->{query_accession_mappings};
    }
    else {
        confess "Error, don't understand type setting: $type\n";
    }

    my @contigs_with_alignments;
    
    my $list_aref = $accs_href->{$contig};
    if (ref $list_aref) {
        @contigs_with_alignments = @$list_aref;
    }

    return (@contigs_with_alignments);
}



####
sub get_delta_parser {
    my $self = shift;
    return ($self->{delta_parser});
}

####
sub get_accession_to_alignments_mappings {
    my $self = shift;
    
    return ($self->{accessions_to_alignments});
}


####
sub remap_contig_and_coordinates {
    my $self = shift;
    my ($type, $contig_id, @coordinates) = @_;

    
    my @mapped_contig_ids = $self->get_contig_mappings($type, $contig_id);

    my @results;

    foreach my $coordinate (@coordinates) {

        my $transformed_coord_obj = Mummer_remap::transformed_coordinate->new($contig_id, $coordinate);
        push (@results, $transformed_coord_obj);

        ## try each mapped contig ID
        
        foreach my $mapped_contig_id (@mapped_contig_ids) {
            
            my ($ref_contig, $query_contig) = ($type == $REFERENCE_TO_QUERY) ? ($contig_id, $mapped_contig_id) : ($mapped_contig_id, $contig_id);
            
            my $ids_href = { reference_accession => $ref_contig,
                             query_accession => $query_contig };
            
            
            my @mappings = $self->transform_coordinate_between_contigs($type, $ids_href, $coordinate);

            if (@mappings) {
                foreach my $mapping (@mappings) {
                    $transformed_coord_obj->add_transformed_location($mapping->{contig},
                                                                     $mapping->{coordinate},
                                                                     $mapping->{alignment} );
                }
            }
        }
    }
    return (@results);
}

####
sub remap_coordinate_pair_to_best_location {
    my $self = shift;
    my ($type, $contig_id, $end5, $end3, $percent_length_difference_allowed) = @_;
    my @location_pairs = $self->remap_coordinate_pair_to_all_location_pairs(@_);

    unless (@location_pairs) { return; }

    @location_pairs = sort {$b->{delta_span}<=>$a->{delta_span}
                               ||
                                   $a->{sum_alignment_lengths}<=>$b->{sum_alignment_lengths}} @location_pairs;

    
    my $best_location = pop @location_pairs;
    my ($new_end5, $new_end3) = ($best_location->{end5}, $best_location->{end3});
    my $new_contig_id = $best_location->{contig};
    
    return ({ contig_id => $new_contig_id,
              end5 => $new_end5, 
              end3 => $new_end3} );
}



####
sub remap_coordinate_pair_to_all_location_pairs {
    my $self = shift;
    my ($type, $contig_id, $end5, $end3, $percent_length_difference_allowed) = @_;
    
    my $span_length = abs ($end3 - $end5) + 1;

    unless ($percent_length_difference_allowed) { 
        $percent_length_difference_allowed = 5; 
    }

    my ($lend_mappings, $rend_mappings) = $self->remap_contig_and_coordinates($Mummer_remap::REFERENCE_TO_QUERY, $contig_id, ($end5, $end3));
    
    ## examine contig pairs and score by sum length
    my %contig_mappings;  # {lend}->[list], {rend}->[list]
    
    foreach my $mapping_set ( [ "lend", $lend_mappings],
                              [ "rend", $rend_mappings] ) {
        
        my ($map_type, $mappings_obj) = @$mapping_set;

        foreach my $location ($mappings_obj->get_transformed_locations()) {
            
            ## group by orientation too!
            my $alignment = $location->{alignment};
            my $orientation = $alignment->get_orientation();
            
            my $contig_id = $location->{contig} . "$;" . $orientation; # encode orientation into the contig identifier here.
            
            #ensure contig listing
            unless ($contig_mappings{$contig_id}->{lend}) {
                $contig_mappings{$contig_id}->{lend} = [];
            }
            unless ($contig_mappings{$contig_id}->{rend}) {
                $contig_mappings{$contig_id}->{rend} = [];
            }
            my $list_aref = $contig_mappings{$contig_id}->{$map_type};
            push (@$list_aref, $location);
        }

    }

    ## examine location pairs:
    my @location_pairs;
    foreach my $contig (keys %contig_mappings) {
        my $lend_list = $contig_mappings{$contig}->{lend};
        my $rend_list = $contig_mappings{$contig}->{rend};
        
        foreach my $lend_ele (@$lend_list) {
            my $lend_ele_coordinate = $lend_ele->{coordinate};
            my $lend_ele_alignment_length = $lend_ele->{alignment}->{refseq_match_length};
            
            foreach my $rend_ele (@$rend_list) {
                my $rend_ele_coordinate = $rend_ele->{coordinate};
                my $rend_ele_alignment_length = $rend_ele->{alignment}->{refseq_match_length};

                my $coord_span = abs ($rend_ele_coordinate - $lend_ele_coordinate) + 1;
                my $percent_span = $coord_span / $span_length * 100;
                my $delta_span = abs (100 - $percent_span);
                if ($delta_span <= $percent_length_difference_allowed) {
                    my $sum_alignment_lengths = $lend_ele_alignment_length + $rend_ele_alignment_length;

                    my ($contig_identifier, $orient) = split (/$;/, $contig);
                    push (@location_pairs, { contig => $contig_identifier,
                                             orient => $orient,
                                             end5 => $lend_ele_coordinate,
                                             end3 => $rend_ele_coordinate,
                                             sum_alignment_lengths => $sum_alignment_lengths, 
                                             delta_span => $delta_span });
                }
            }
        }
    }
    
    return (@location_pairs);
    
}
    

####
sub get_reference_accessions {
    my $self = shift;
    my @reference_accessions = keys %{$self->{reference_accession_mappings}};
    return (@reference_accessions);
}

####
sub get_query_accessions {
    my $self = shift;
    my @query_accessions = keys %{$self->{query_accession_mappings}};
    return (@query_accessions);
}


####
sub toString {
    my $self = shift;

    my $text = "";
    
    foreach my $reference_accession ($self->get_reference_accessions()) {
        
        $text .= "** Reference accession: $reference_accession\n";
        
        foreach my $query_accession ($self->get_contig_mappings($REFERENCE_TO_QUERY, $reference_accession)) {

            $text .= "// alignments between $reference_accession and $query_accession:\n";
            
            my @alignments = $self->get_alignments_between_query_and_reference ( { reference_accession => $reference_accession, 
                                                                                   query_accession => $query_accession } );
            
            foreach my $alignment (@alignments) {
                $text .= $alignment->toString() . "\n";
            }
            
            $text .= "\n";
        
        }

        
    }
    
    return ($text);

}


####
sub transform_coordinate_between_contigs {
    my $self = shift;
    my ($type, $ids_href, $coordinate) = @_;
    
    unless ($type == $REFERENCE_TO_QUERY || $type == $QUERY_TO_REFERENCE) {
        confess "fatal, don't understand 'type' function argument.  ";
    }
    
    my $target_contig = ($type == $REFERENCE_TO_QUERY) ? $ids_href->{query_accession} : $ids_href->{reference_accession};
    
    my @mappings;

    my @alignments = $self->get_alignments_between_query_and_reference($ids_href);
    unless (@alignments) { return (); }


    my @overlapping_alignments = $self->find_overlapping_coordinate_sets($type, $coordinate, \@alignments);
    
    foreach my $overlapping_alignment (@overlapping_alignments) {
        
        ## find the msp containing the coordinate;
        my $transformed_coordinate = $self->transform_coordinate_within_single_alignment($type, $coordinate, $overlapping_alignment);
        if ($transformed_coordinate) {
            push (@mappings, { contig => $target_contig,
                               coordinate => $transformed_coordinate,
                               alignment => $overlapping_alignment} );
        }
    }

    return (@mappings);
}


####
sub transform_coordinate_within_single_alignment {
    my $self = shift;
    my ($type, $coordinate, $alignment) = @_;
    
    my @msps = $alignment->get_MSPs();
    my $msp = $self->find_longest_overlapping_coordinate_set ($type, $coordinate, \@msps);
    unless ($msp) {
        return;
    }
    
    if ($Mummer_remap::VERBOSE) {
        print "Longest alignment MSP found spanning coordinate $coordinate:\n" 
            . $msp->toString() . "\n";
    }
        
    my $adj_coordinate = $self->transform_single_coordinate($type, $msp, $coordinate);
    
    return ($adj_coordinate);
}



####
sub transform_single_coordinate {
    my $self = shift;
    my ($type, $msp, $coordinate) = @_;
    
    my ($reference_lend, $reference_rend) = $msp->get_reference_coords();
    my ($query_lend, $query_rend) = $msp->get_query_coords();
    my $orientation = $msp->get_orientation();
    
    my $adj_coordinate;
    if ($type == $REFERENCE_TO_QUERY) {
        
        unless ($coordinate >= $reference_lend && $coordinate <= $reference_rend) {
            confess "Error, trying to convert coordinates from reference MSP, "
                . "but coordinate ($coordinate) is not located between ($reference_lend - $reference_rend) ";
        }
    
        my $delta = $coordinate - $reference_lend;
        
        if ($orientation eq '+') {
            $adj_coordinate = $query_lend + $delta;
            print "CoordConvert: case A\n" if $DEBUG;
        }
        else {
            # revcomp match
            $adj_coordinate = $query_rend - $delta;
            print "CoordConvert: case B\n" if $DEBUG;
        }
        
    }
    else {
        ## QUERY_TO_REFERENCE
        
        unless ($coordinate >= $query_lend && $coordinate <= $query_rend) {
            confess "Error, trying to convert coordinates from query MSP, "
                . "but coordinate ($coordinate) is not located between ($query_lend - $query_rend) ";
        }
        
        if ($orientation eq '+') {
            my $delta = $coordinate - $query_lend;
            $adj_coordinate = $reference_lend + $delta;
            print "CoordConvert: case C\n" if $DEBUG;
        }
        else {
            # revcomp'd 
            my $delta = $query_rend - $coordinate;
            $adj_coordinate = $reference_lend + $delta;
            print "CoordConvert: case D\n" if $DEBUG;
        }
    }
    
    return ($adj_coordinate);
}





###
sub find_longest_overlapping_coordinate_set {
    my $self = shift;
    my ($type, $coordinate, $coord_sets_aref) = @_;
    
    my @overlapping_MSPs = $self->find_overlapping_coordinate_sets($type, $coordinate, $coord_sets_aref);
    
    my $longest_MSP = undef;
    my $longest_length = 0;
    foreach my $msp (@overlapping_MSPs) {
        my $length = $msp->get_reference_coord_span_length();
        if ($length > $longest_length) {
            $longest_length = $length;
            $longest_MSP = $msp;
        }
    }

    return ($longest_MSP);

}


####
sub find_overlapping_coordinate_sets {
    my $self = shift;
    my ($type, $coordinate, $coord_sets_aref) = @_;
    
    my @overlapping_MSPs;
    
    foreach my $coordset (@$coord_sets_aref) {
        my $got_overlap_flag = 0;
        my $length;
        if ($type == $REFERENCE_TO_QUERY) {
            my ($reference_lend, $reference_rend) = $coordset->get_reference_coords();
            if ($reference_lend <= $coordinate && $coordinate <= $reference_rend) {
                ## overlap found
                $got_overlap_flag = 1;
            }
        }
        elsif ($type == $QUERY_TO_REFERENCE) {
            my ($query_lend, $query_rend) = $coordset->get_query_coords();
            if ($query_lend <= $coordinate && $coordinate <= $query_rend) {
                ## overlap found
                $got_overlap_flag = 1;
            }
        }
        else {
            confess "fatal bug, type not recognized.";
        }
        
        if ($got_overlap_flag) {
            push (@overlapping_MSPs, $coordset);
        }
    }
    
    return (@overlapping_MSPs);
}


####
sub get_alignments_between_query_and_reference {
    my $self = shift;
    my ($ids_href) = @_;
    
    my $reference_accession = $ids_href->{reference_accession} or confess "must specify reference_accession as key in hash_ref ";
    my $query_accession = $ids_href->{query_accession} or confess "must specify query_accession as key in hash_ref ";
    
    my $accessions_to_alignments_href = $self->get_accession_to_alignments_mappings();
    my $alignment_list_aref = $accessions_to_alignments_href->{$reference_accession}->{$query_accession};
    
    unless (ref $alignment_list_aref) {
        die "Error, couldn't find alignment based on ref: $reference_accession, query: $query_accession ";
    }
    
    return (@$alignment_list_aref);
    
    
}

################################################################

package Mummer_remap::transformed_coordinate;

use strict;
use warnings;

####
sub new {
    my $packagename = shift;
    my ($contig_id, $coordinate) = @_;

    my $self = { contig_id => $contig_id,
                 coordinate => $coordinate,

                 transformed_locations => undef,

             };
    
    bless ($self, $packagename);
    
    return ($self);
}

####
sub add_transformed_location {
    my $self = shift;
    my ($other_contig, $other_location, $alignment) = @_;
    
    my $transformed_locations_list_aref = $self->{transformed_locations};
    unless (ref $transformed_locations_list_aref) {
        $transformed_locations_list_aref = $self->{transformed_locations} = [];
    }

    push (@$transformed_locations_list_aref, { contig => $other_contig,
                                               coordinate => $other_location,
                                               alignment => $alignment } );
    
    return;
}

####
sub get_transformed_locations {
    my $self = shift;
    my @locations;
    my $locations_list_aref = $self->{transformed_locations};
    if (ref $locations_list_aref) {
        @locations = @$locations_list_aref;
    }
    return (@locations);
}


####
sub toString {
    my $self = shift;

    my $text = "Mummer_remap::transformed_coordinate contig: " . $self->{contig_id} . " coord: " . $self->{coordinate} . "\n";

    my $transformed_locations_aref = $self->{transformed_locations};
    if (ref $transformed_locations_aref) {
        foreach my $transformed_location (@$transformed_locations_aref) {
            my ($other_contig_id, $other_coordinate) = ($transformed_location->{contig},
                                                        $transformed_location->{coordinate});
            $text .= "\ttransformed to: other_contig: $other_contig_id, pos: $other_coordinate\n";
        }
    }
    else {
        $text .= "\tno transformed locations.\n";
    }
    return ($text);
}


1; #EOM
    
