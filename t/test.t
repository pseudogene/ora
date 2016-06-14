#!/usr/bin/perl
use strict;
use warnings;
use Test::More tests => 1;
my $output) = q{};
system
  'wget -q ftp://ftp.ensembl.org/pub/release-84/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna.chromosome.10.fa.gz >/dev/null 2>/dev/null';
system 'gunzip Canis_familiaris.CanFam3.1.dna.chromosome.10.fa.gz';
system
  './scripts/or.pl -v -a -d --sequence Canis_familiaris.CanFam3.1.dna.chromosome.10.fa > test.output';
unlink 'Canis_familiaris.CanFam3.1.dna.chromosome.10.fa';
if (open my $in, q{<}, 'test.output')
{
    while (<$in>)
    {
        chomp;
        $output .= $_;
    }
    close $in;
}
unlink 'test.output';
ok( length $output > 0 , 'Run test (dog chromosome 10)');
done_testing();
