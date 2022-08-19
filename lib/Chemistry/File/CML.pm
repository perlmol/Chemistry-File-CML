package Chemistry::File::CML;

# VERSION
# $Id$

use base 'Chemistry::File';
use Chemistry::Mol;
use Carp;
use List::Util;
use XML::LibXML;
use strict;
use warnings;

our $DEBUG = 0;

=head1 NAME

Chemistry::File::CML - CML reader/writer

=head1 SYNOPSIS

    use Chemistry::File::CML;

    # read a molecule
    my $mol = Chemistry::Mol->read('myfile.cml');

    # write a molecule
    $mol->write("myfile.cml");

    # use a molecule as a query for substructure matching
    use Chemistry::Pattern;
    use Chemistry::Ring;
    Chemistry::Ring::aromatize_mol($mol);

    my $patt = Chemistry::Pattern->read('query.mol');
    if ($patt->match($mol)) {
        print "it matches!\n";
    }

=cut

Chemistry::Mol->register_format(cml => __PACKAGE__);

=head1 DESCRIPTION

MDL Molfile (V2000) reader/writer.

This module automatically registers the 'mdl' format with Chemistry::Mol.

The first three lines of the molfile are stored as $mol->name, 
$mol->attr("mdlmol/line2"), and $mol->attr("mdlmol/comment").

This version only reads and writes some of the information available in a
molfile: it reads coordinates, atom and bond types, charges, isotopes,
radicals, and atom lists. It does not read other things such as
stereochemistry, 3d properties, etc.

This module is part of the PerlMol project, L<https://github.com/perlmol>.

=head2 Query properties

The MDL molfile format supports query properties such as atom lists, and
special bond types such as "single or double", "single or aromatic", "double or
aromatic", "ring bond", or "any". These properties are supported by this module
in conjunction with L<Chemistry::Pattern>. However, support for query properties
is currently read-only, and the other properties listed in the specification
are not supported yet.

So that atom and bond objects can use these special query options, the
conditions are represented as Perl subroutines. The generated code can be
read from the 'mdlmol/test_sub' attribute:

    $atom->attr('mdlmol/test_sub');
    $bond->attr('mdlmol/test_sub');
 
This may be useful for debugging, such as when an atom doesn't seem to match as
expected.

=head2 Aromatic Queries

To be able to search for aromatic substructures are represented by Kekule
structures, molfiles that are read as patterns (with
C<Chemistry::Pattern->read) are aromatized automatically by using the
L<Chemistry::Ring> module. The default bond test from Chemistry::Pattern::Bond
is overridden by one that checks the aromaticity in addition to the bond order.
The test is,

    $patt->aromatic ?  $bond->aromatic 
        : (!$bond->aromatic && $patt->order == $bond->order);

That is, aromatic pattern bonds match aromatic bonds, and aliphatic pattern
bonds match aliphatic bonds with the same bond order.
    

=cut

# some constants, based on tables from the file format specification

my %BOND_TYPE_EXPR = (
    4 => '($bond->aromatic)',
    5 => '($bond->order == 1 or $bond->order == 2)',
    6 => '($bond->order == 1 or $bond->aromatic)',
    7 => '($bond->order == 2 or $bond->aromatic)',
    8 => '(1)',                                         # any bond
);

my %BOND_TOPOLOGY_EXPR = (
    1 => '@{$bond->attr("ring/rings")||[]}',
    2 => '! @{$bond->attr("ring/rings")||[]}',
);

sub parse_string {
    my ($self, $s, %opts) = @_;

    my $mol_class  = $opts{mol_class}  || 'Chemistry::Mol';
    my $atom_class = $opts{atom_class} || $mol_class->atom_class;
    my $bond_class = $opts{bond_class} || $mol_class->bond_class;
    local $_;

    my $cml = XML::LibXML->load_xml( string => $s );
    my $xp = XML::LibXML::XPathContext->new( $cml );
    $xp->registerNs( 'cml', 'http://www.xml-cml.org/schema' );

    my @molecules;
    for my $molecule ($xp->findnodes( '/cml:cml/cml:molecule' )) {
        my $mol = $mol_class->new;
        push @molecules, $mol;

        $mol->name( $molecule->getAttribute( 'id' ) ) if $molecule->hasAttribute( 'id' );

        my ($atomArray) = $molecule->getChildrenByTagName( 'atomArray' );
        next unless $atomArray; # Skip empty molecules

        my %atom_by_name;
        my %hydrogens_by_id;

        # atomArray
        for my $element ($atomArray->getChildrenByTagName( 'atom' )) { # for each atom...
            my ($symbol, $charge, $hydrogen_count, $mass_number);
            my @coord3;

            next unless $element->hasAttribute( 'id' );
            my $id = $element->getAttribute( 'id' );
            my $atom = $atom_by_name{$id} = $mol->new_atom( name => $id );

            if( $element->hasAttribute( 'elementType' ) ) {
                $atom->symbol( $element->getAttribute( 'elementType' ) );
            }
            if( $element->hasAttribute( 'formalCharge' ) ) {
                $atom->formal_charge( int $element->getAttribute( 'formalCharge' ) );
            }
            if( $element->hasAttribute( 'hydrogenCount' ) ) {
                $hydrogens_by_id{$atom->id} = int $element->getAttribute( 'hydrogenCount' );
            }
            if( $element->hasAttribute( 'isotopeNumber' ) ) {
                $atom->mass_number( int $element->getAttribute( 'isotopeNumber' ) );
            }
            if( $element->hasAttribute( 'x3' ) &&
                $element->hasAttribute( 'y3' ) &&
                $element->hasAttribute( 'z3' ) ) {
                $atom->coords( map { $_ * 1 } $element->getAttribute( 'x3' ),
                                              $element->getAttribute( 'y3' ),
                                              $element->getAttribute( 'z3' ) );
            }
        }

        my @bonds;
        my( $bondArray ) = $molecule->getChildrenByTagName( 'bondArray' );
        if( $bondArray ) {
            @bonds = $bondArray->getChildrenByTagName( 'bond' );
        }

        # bond block
        for my $bond (@bonds) { # for each bond...
            my $order = my $type = $bond->getAttribute( 'order' );
            $order = 1 unless $order =~ /^[123]$/;
            $mol->new_bond(
                type => $type, 
                atoms => [map { $atom_by_name{$_} } split ' ', $bond->getAttribute( 'atomRefs2' )],
                order => $order,
            );
            #~ if ($mol->isa('Chemistry::Pattern')) {
                #~ $self->bond_expr($bond, $i, $type, $topology);
            #~ }
        }

        # calculate implicit hydrogens
        for my $id (sort keys %hydrogens_by_id) {
            my $atom = $mol->by_id( $id );
            my $explicit_hydrogens = scalar grep { $_->symbol eq 'H' }
                                                 $atom->neighbors;
            if( $explicit_hydrogens > $hydrogens_by_id{$id} ) {
                warn 'total number of attached hydrogen atoms is ' .
                     "less than the number of explicit hydrogen atoms\n";
                next;
            }
            next if $explicit_hydrogens == $hydrogens_by_id{$id};
            $atom->implicit_hydrogens( $hydrogens_by_id{$id} - $explicit_hydrogens );
        }

        if ($mol->isa('Chemistry::Pattern')) {
            require Chemistry::Ring;
            Chemistry::Ring::aromatize_mol($mol);
        }
    }

    return @molecules;
}

sub bond_expr {
    my ($self, $bond, $i, $type, $topology) = @_;
    my @bond_exprs;
    my $s = $BOND_TOPOLOGY_EXPR{$topology};
    push @bond_exprs, $s if $s;
    $s = $BOND_TYPE_EXPR{$type};
    push @bond_exprs, $s if $s;
    if (@bond_exprs) {
        my $expr = join " and ", @bond_exprs;
        my $sub_txt = <<SUB;
            sub {
                no warnings;
                my (\$patt, \$bond) = \@_;
                $expr;
            };
SUB
        print "MDLMol bond($i) sub: <<<<$sub_txt>>>\n" if $DEBUG;
        $bond->attr('mdlmol/test_sub' => $sub_txt);
        $bond->test_sub(eval $sub_txt);
    } else { # default bond sub
        $bond->test_sub(\&default_bond_test);
    }
}

sub default_bond_test {
    no warnings;
    my ($patt, $bond) = @_;
    $patt->aromatic ?  $bond->aromatic 
        : (!$bond->aromatic && $patt->order == $bond->order);
}

sub name_is {
    my ($self, $fname) = @_;
    $fname =~ /\.cml$/i;
}

sub file_is {
    my ($self, $fname) = @_;
    $fname =~ /\.cml$/i;
}

1;

=head1 SOURCE CODE REPOSITORY

L<https://github.com/perlmol/Chemistry-File-CML>

=head1 SEE ALSO

L<Chemistry::Mol>

=head1 AUTHOR

Andrius Merkys <merkys@cpan.org>

=head1 COPYRIGHT

Copyright (c) 2022 Andrius Merkys. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

