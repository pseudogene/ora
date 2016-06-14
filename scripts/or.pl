#!/usr/bin/perl -w
#
# Olfactory Receptor family Assigner (ORA)
# Copyright 2007-2008 Bekaert M <michael@batlab.eu>
#
# This work is licensed under the Creative Commons Attribution-
# Noncommercial-Share Alike 3.0 License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/
use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
use Term::ProgressBar;
use Bio::Seq;
use Bio::SeqIO;
use Bio::ORA;

#----------------------------------------
my $version = '1.5';
my $hmm     = getcwd() . '/or.hmm';
my $ref     = getcwd() . '/or.fasta';


#----------------------------------------------------------
sub hmm_disc {
    my ( $seq, $name, $id, $translate, $evalue, $format, $detail, $aug, $hmm, $filter, $frag, $subset ) = @_;
    my $chrom;
    if    ( defined $name )            { $chrom = $name . '_' . $id; }
    elsif ( defined $seq->display_id ) { $chrom = $seq->display_id; }
    else                               { $chrom = 'noname_' . $id; }
    $seq->display_id($chrom);
    my $ORA_obj = Bio::ORA->new( $seq, $translate, $aug, $hmm );
    $ORA_obj->{'_frag'} = $frag if ( defined $frag );
    if ( $ORA_obj->find($evalue) && ( !( defined $filter ) || ( ( defined $filter ) && ( substr( $ORA_obj->{'_hmmor'}[1], 2 ) == $filter ) ) ) ) {
        if ( defined $subset ) {
            my $out = new Bio::SeqIO( '-format' => 'fasta', '-file' => ">>$subset" );
            $out->write_seq($seq);
        }
        $ORA_obj->show($format, $detail);
    }
    return $ORA_obj->{'_verbose'};
}

sub fasta_filter {
    my ( $seq, $ref, $id, $translate, $evalue, $format, $detail, $aug, $hmm, $filter, $frag, $subset ) = @_;
    my $mess = '* FASTA search for ' . $seq->display_id . "\n";
    my @hits = Bio::ORA->getHits( $seq, 1, $ref );
    if ( $#hits >= 0 ) {
        for ( my $i = 0 ; $i <= $#hits ; $i++ ) {
            my ( $hitstrand, $hitstart, $hitend ) = split( /\|/, $hits[$i] );
            my $seqstart = $hitstart - 250;
            $seqstart = 1 if ( $seqstart < 1 );
            my $seqend = $hitend + 250;
            $seqend = $seq->length() if ( $seqend > $seq->length() );
            my $seq_or = Bio::Seq->new( -seq => $seq->subseq( $seqstart, $seqend ), -alphabet => 'dna', -id => $seq->display_id . ":$seqstart-$seqend" );
            my $ORA_obj = Bio::ORA->new( $seq_or, $translate, $aug, $hmm );
            $ORA_obj->{'_frag'} = $frag if ( defined $frag );
            if ( $ORA_obj->find( $evalue, $hitstrand, $seqstart, $seqend ) && ( !( defined $filter ) || ( ( defined $filter ) && ( substr( $ORA_obj->{'_hmmor'}[1], 2 ) == $filter ) ) ) ) {
                if ( defined $subset ) {
                    my $out = new Bio::SeqIO( '-format' => 'fasta', '-file' => ">>$subset" );
                    $out->write_seq($seq_or);
                }
                $ORA_obj->show($format, $detail);
            }
            $mess .= ' ' . $ORA_obj->{'_verbose'};
        }
    } else {
        $mess = '* No FASTA hit for ' . $seq->display_id . "\n";
    }
    return $mess;
}

#------------------------ Main ----------------------------
my ( $infile, $translate, $name, $filter, $subset, $aug, $resume, $frag, $progress, $maxseq );
my ( $evalue, $contigs, $bar, $verbose, $detail, $format, $out, $seqcount, $next_update ) = ( 1e-10, 0, 0, 0, 0, 'fasta', 1, 0 ,0 );
GetOptions( 'sequence:s' => \$infile, 'c!' => \$contigs, 'a!' => \$aug, 'format:s' => \$format, 'expect:f' => \$evalue, 'name:s' => \$name, 'table:i' => \$translate, 'filter:i' => \$filter, 'sub:s' => \$subset, 'resume:s' => \$resume, 'd!' => \$detail, 'v!' => \$verbose, 'size:i' => \$frag, 'p!' => \$bar );
$format = lc $format;
if ( ( defined $infile ) && ( -r $infile ) && ( !( defined $translate ) || ( ( defined $translate ) && ( $translate =~ /([1-9]|1[1-6]|2(1|2))/ ) ) ) && ( !( defined $filter ) || ( ( defined $filter ) && ( $filter =~ /^\d+$/ ) ) ) && ( $format =~ /^fasta|genbank|csv|r|tbl$/ ) ) {
    $translate = 1 unless ( defined $translate );
    print STDERR "\n..:: Olfactory Receptor Assigner (ORA) ::..\n> Standalone program version $version <\n\n" if ($verbose);
    if ( $contigs ) {
        $out--;
        my $id = 0;
        my @prehits = Bio::ORA->fastScan( $infile, $ref );
        if ( ($#prehits >= 0) && (my $inseq = Bio::SeqIO->new( '-file' => "<$infile", '-format' => 'fasta' ) ) ) {
            print STDERR $#prehits, " Informative contigs\n" if ($verbose); 
##    
    if ($bar) {
       $progress = Term::ProgressBar->new({name => 'Search', count => $#prehits, ETA => 'linear', term_width=>80 } );
       $progress->max_update_rate(1);
       $maxseq = $#prehits;
       $seqcount = -1;
    }
##  
            my $myresume = shift(@prehits); 
            while ( my $seq = $inseq->next_seq ) {
            	my $mess;
                next if (!(defined $myresume) || ($myresume ne $seq->display_id));
                $myresume = shift(@prehits);
                $next_update = $progress->update($seqcount) if($bar && (++$seqcount > $next_update) && $id);
                next if ( (defined $resume) && !$id && ($resume ne $seq->display_id) );
	            $id++;
                if ( length( $seq->seq() ) > 2500 ) { $mess = fasta_filter( $seq, $ref, $id, $translate, $evalue, $format, $detail, $aug, $hmm, $filter, $frag, $subset ); }
                else                                { $mess = hmm_disc( $seq, $name, $id, $translate, $evalue, $format, $detail, $aug, $hmm, $filter, $frag, $subset ); }
                $progress->message($mess) if($bar && $verbose && (defined $mess));
                print STDERR $mess if(!$bar && $verbose && (defined $mess));
            }
        }
    } elsif ( my $inseq = Bio::SeqIO->new( '-file' => "<$infile", '-format' => 'fasta' ) ) {
##    
    if ($bar) {
       open IN, $infile;
       while(<IN>) {
         $seqcount++ if(/^>/);
       }
       close IN;
       $progress = Term::ProgressBar->new({name => 'Search', count => $seqcount, ETA => 'linear', term_width=>80 } );
       $progress->max_update_rate(1);
       $maxseq = $seqcount;
       $seqcount = -1;
    }
##  
        $out--;
        my $id = 0;
        while ( my $seq = $inseq->next_seq ) {
        	my $mess;
	        $next_update = $progress->update($seqcount) if($bar && (++$seqcount > $next_update) && $id);
            next if ( (defined $resume) && !$id && ($resume ne $seq->display_id) );
            $id++;
            if ( length( $seq->seq() ) > 2500 ) { $mess = fasta_filter( $seq, $ref, $id, $translate, $evalue, $format, $detail, $aug, $hmm, $filter, $frag, $subset ); }
            else                                { $mess = hmm_disc( $seq, $name, $id, $translate, $evalue, $format, $detail, $aug, $hmm, $filter, $frag, $subset ); }
            $progress->message($mess) if($bar && $verbose && (defined $mess));
            print STDERR $mess if(!$bar && $verbose && (defined $mess));
        }
        $progress->update($maxseq) if($bar && ($maxseq >= $next_update));
    } else {
        print STDERR "FATAL: Incorrect file format.\n";
    }
    print STDERR "\n" if ($bar || $verbose);
} else {
    print STDERR "\n..:: Olfactory Receptor Assigner (ORA) ::..\n> Standalone program version $version <\n\nFATAL: Incorrect arguments.\nUsage: or.pl [-options] --sequence=<sequence file>\n\n Options\n   -a\n         Force the use of alternative start codons, according to the\n         current genetic code. Otherwise, ATG is the only initiation codon\n         allow.\n   --expect\n         Set the E-value threshold. This setting specifies the statistical\n         significance threshold for reporting matches against database\n         sequences. The default value (1e-10).\n   --format\n         Available output format codes include 'fasta' (FASTA format);\n         'csv' (Comma-separated values); 'genbank' (GenBank format);\n         'R' (Direct output for R).\n         By default the fasta format is used.\n   --filter\n         Show ONLY the selected family number.\n   --sub\n         Extact the sequences of the Fasta hits.\n   --resume\n         Resume the search from given sequence name.\n   --name\n         Overwrite the sequence name by the provided one. Otherwise the\n         program will use the sequence name from as input.\n   --table\n         Force a genetic code to be used for the translation of the query.\n   --size\n         Filter fragments over the specified size as functional.\n   -d\n         Print all sequence details.\n   -v\n         Print more possibly useful stuff, such as the individual scores\n         for each sequence.\n\n";
    $out++;
}
exit($out);
