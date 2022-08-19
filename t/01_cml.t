use strict;
use warnings;

use Chemistry::File::CML;
use Test::More tests => 2;

my @mols = Chemistry::File::CML->parse_string( <<'END' );
<?xml version="1.0"?>
<cml xmlns="http://www.xml-cml.org/schema">
  <molecule id="name">
  </molecule>
</cml>
END

is( scalar @mols, 1 );
is( $mols[0]->name, 'name' );
