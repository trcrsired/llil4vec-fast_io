# gen-long-llil.pl
# Crude program to generate a LLiL test file with long names and counts
# perl gen-long-llil.pl long1.txt 600

use strict;
use warnings;
use autodie;

{
   my $ordmin = ord('a');
   my $ordmax = ord('z') + 1;

   # Generate a random word
   sub gen_random_word {
      my $word  = shift;    # word prefix
      my $nchar = shift;    # the number of random chars to append
      for my $i (1 .. $nchar) {
         $word .= chr( $ordmin + int( rand($ordmax - $ordmin) ) );
      }
      return $word;
   }
}

my $longworda = join '', 'a' .. 'z';
my $longwordz = join '', reverse('a' .. 'z');
my $longcount = 1_000_000;

sub create_long_test_file {
   my $fname   = shift;
   my $howmany = shift;
   open( my $fh_out, '>', $fname );

   # Some with no randomness
   for my $h ( 1 .. $howmany ) {
      for my $i ( 1 .. 8 ) {
         my $cnt   = $longcount + $i - 1;
         my $worda = $longworda x $i;
         my $wordz = $longwordz x $i;
         print {$fh_out} "$worda\t$cnt\n$wordz\t$cnt\n";
      }
   }

   # Some with randomness
   my $wordlen = 1;
   for my $h ( 1 .. $howmany ) {
      for my $i ( 1 .. 8 ) {
         my $cnt   = $longcount + $i - 1;
         my $worda = $longworda x $i;
         my $wordz = $longwordz x $i;
         for my $c ( 'a' .. 'z' ) {
            for my $z ( 1 .. 2 ) {
               print {$fh_out} $worda . gen_random_word( $c, $wordlen ) . "\t" . (1000000 + $z) . "\n";
               print {$fh_out} $wordz . gen_random_word( $c, $wordlen ) . "\t" . (1000000 + $z) . "\n";
            }
         }
      }
   }
}

my $outfile = shift;
my $count   = shift;
$outfile or die "usage: $0 outfile count\n";
$count   or die "usage: $0 outfile count\n";
$count =~ /^\d+$/ or die "error: count '$count' is not a number\n";
print "generating short long test file '$outfile' with count '$count'\n";
create_long_test_file( $outfile, $count );
print "file size=", -s $outfile, "\n";
