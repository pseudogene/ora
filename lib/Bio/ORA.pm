## $Id: ORA.pm,v 1.5 2008/12/15 10:26:38 Europe/Dublin $
#
# Olfactory Receptor family Assigner (ORA)
# Copyright 2007-2008 Bekaert M <michael@batlab.eu>
#
# This work is licensed under a Creative Commons GNU General Public
# License License (GPL). To view a copy of this license,
# visit http://creativecommons.org/licenses/GPL/2.0/
# 
# POD documentation - main docs before the code

=head1 NAME

Bio::ORA - Olfactory Receptor family Assigner (ORA)

=head1 SYNOPSIS

Take a sequence object from, say, an inputstream, and find an Olfactory
Receptor gene. HMM profiles are used in order to identify location, frame
and orientation of such gene.

Creating the ORA object, eg:

  use Bio::ORA;
  use Bio::SeqIO;

  my $inputstream = Bio::SeqIO->new( -file => 'seqfile', -format => 'fasta' );
  my $seqobj = $inputstream->next_seq();
  my $ORA_obj = Bio::ORA->new( $seqobj );

Obtain an array holding the start point, the stop point, the DNA sequence
and amino-acid sequence, eg:

  my $array_ref = $ORA_obj->{'_result'} if ( $ORA_obj->find() );

Display result in genbank format, eg:

  $ORA_obj->show( 'genbank' );


=head1 DESCRIPTION

Bio::ORA is a featherweight object for identifying mammalian
olfactory receptor genes. The sequences should not be longer than 20kb. The
returned array include location, sequence and statistic for the putative
olfactory receptor gene. Fully functional with DNA and EST
sequence, no intron supported.

See Synopsis above for the object creation code.


=head1 DRIVER SCRIPT

  #!/usr/bin/perl -w
  use strict;

  use Bio::Seq;
  use Bio::ORA;

  my $inseq = Bio::SeqIO->new( '-file' => "< $yourfile", -format => 'fasta' );
  while (my $seq = $inseq->next_seq) {
    my $ORA_obj = Bio::ORA->new( $seq );
    if ( $ORA_obj->find() ) {
      $ORA_obj->show( 'genbank' );
    } else {
      print "  no hit!\n";
    }
  }


=head1 REQUIREMENTS

To use this module you may need:
 * Bioperl (L<http://www.bioperl.org/>) modules,
 * HMMER distribution (L<http://hmmer.janelia.org/>) and
 * FASTA distribution (L<ftp://ftp.ebi.ac.uk/pub/software/unix/fasta/>).


=head1 MORE
ORA home page is now L<http://ora.batlab.eu/>


=head1 FEEDBACK

User feedback is an integral part of the evolution of this modules. Send
your comments and suggestions preferably to author.


=head1 AUTHOR

B<Michael Bekaert> (michael@batlab.eu)

Address:
     School of Biology & Environmental Science
     University College Dublin
     Belfield, Dublin 4
     Dublin
     Ireland


=head1 COPYRIGHT

This work is licensed under a Creative Commons GNU General Public
License License (GPL). The full text of the license can be
found in the LICENSE file included with this module.


=head1 SEE ALSO

perl(1), fasta, hmmer and bioperl web sites


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _


=cut

package Bio::ORA;
use strict;
use File::Temp qw/tempfile/;
use Bio::Tools::HMMER::Results;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;
use Bio::PrimarySeq;
use Bio::PrimarySeqI;
use Bio::Tools::CodonTable;

BEGIN {
    use vars qw($VERSION @ISA $PATH_REF $PATH_HMM $PATH_TMP);
    $VERSION     = '1.5';
    @ISA         = qw(Bio::Root::Root Bio::Root::IO);

    # Default path
    $PATH_REF = './or.fasta';
    $PATH_HMM = './or.hmm';
    $PATH_TMP = '/tmp';
}


=head2 _findexec

 Title   : _findexec
 Usage   : my $path = $self->_findexec( $exec );
 Function: Find an executable file in the $PATH.
 Returns : The full path to the executable.
 Args    : $exec (mandatory) executable to be find.

=cut                       

sub _findexec
{
    my ($self, @args) = @_;
    my $exec = shift(@args);
    foreach my $p (split(/:/, $ENV{'PATH'})) {
        return "$p/$exec"
          if (-x "$p/$exec");
    }
}

=head2 new

 Title   : new
 Usage   : my $nb = Bio::ORA->new( $seqobj, $table, $aug, $hmm );
 Function: Initialize ORA object.
 Returns : An ORA object.
 Args    : $seqobj (mandatory) PrimarySeqI object (dna or rna),
           $table (optional) translation table/genetic code number,
              the default value is 1,
           $aug (optional) use other start codon than AUG (default 0),
           $hmm (optional) path to hmm profiles by default ORA looks at
             ./oaz.hmm.

=cut

sub new
{
    my ($self, @args) = @_;
    my ($seqobj, $table, $aug, $hmm) = @args;
    $self = {};
    $self->{'_aug'} = ((defined $aug && int($aug) == 1) ? 0 : 1);
    $self->{'_hmm'} = ((defined $hmm) ? $hmm : $PATH_HMM);
    $self->{'_pfam'} =
      ((defined $ENV{'HMMERDIR'})
        ? $ENV{'HMMERDIR'}
        : Bio::ORA->_findexec('hmmpfam'));
    $seqobj->throw(
               "die in _initialize, ORA.pm works only on PrimarySeqI objects\n")
      unless ($seqobj->isa("Bio::PrimarySeqI"));
    $seqobj->throw("die in _initialize, ORA.pm works only on DNA sequences\n")
      if ($seqobj->alphabet eq 'protein');
    $seqobj->throw(
              "die in _initialize, ORA.pm works only on DNA sequences < 40kb\n")
      if (length($seqobj->seq()) > 40000);
    $seqobj->throw("die in _initialize, hmmpfam not found\n")
      unless (-x $self->{'_pfam'});
    $seqobj->throw('die in _initialize, hmm profile not found at '
                   . $self->{'_hmm'} . "\n")
      unless (-f $self->{'_hmm'});
    my $chrom =
      ((defined $seqobj->display_id) ? $seqobj->display_id : 'noname');
    $seqobj = uc $seqobj->seq();
    $seqobj =~ tr/U/T/;
    $chrom  =~ tr/,/\_/;
    $self->{'_seqref'} =
      Bio::Seq->new(-seq => $seqobj, -alphabet => 'dna', -id => $chrom);
    $self->{'_table'} = ((defined $table) ? $table : 1);
    $self->{'_verbose'} = '';
    bless($self);
    return $self;
}

=head2 find

 Title   : find
 Usage   : my $bool = $ORI_obj->find( $evalue, $strand, $start, $end );
 Function: Identify an olfactory receptor protein.
 Returns : boolean.
 Args    : $evalue (optional) set the E-value (expected) threshold.
             Default is 1e-40,
           $strand(optional) strand where search should be done (1 direct,
             -1 reverse or 0 both). Default is 0,
           $start (optional) coordinate of the first nucleotide. Useful
             for coordinate calculation when first is not 1. Default is 1,
           $end (optional) coordinate of the last nucleotide. Default is
             the sequence length.

=cut

sub find
{
    my ($self, @args) = @_;
    my ($evalue, $strand, $start, $end) = @args;
    $self->{'_evalue'} = ((defined $evalue) && ($evalue > 0)) ? $evalue : 1e-40;
    $strand = 0 unless ((defined $strand) && ($strand == 1 || $strand == -1));
    $start = 1 unless ((defined $start) && $start > 1);
    $end = $self->{'_seqref'}->length
      unless ((defined $end) && $end > 1 && $end > $start);
    my $status = $self->_what_or($strand);
    if (($status) && ($self->{'_hmmor'}[0] < $self->{'_evalue'}))
    {
        $self->{'_verbose'} = '> '
          . $self->{'_hmmor'}[1]
          . ' found ('
          . $self->{'_hmmor'}[0]
          . ') for '
          . $self->{'_seqref'}->display_id
          . ' in frame '
          . ($self->{'_hmmor'}[6] > 0 ? '+' : '-')
          . $self->{'_hmmor'}[2] . "\n";
        return $self->_find_orf($self->{'_hmmor'}[6], $start, $end);
    }
    else
    {
        $self->{'_verbose'} =
          '> no hit for ' . $self->{'_seqref'}->display_id . "\n";
    }
    return 0;
}

=head2 _what_or

 Title   : _what_or
 Usage   : my $bool = $self->_what_oaz( $strand );
 Function: Use HMM profiles to identify an olfactory receptor gene.
 Returns : boolean.
 Args    : $strand (optional) strand where search should be done
           (1 direct, -1 reverse or 0 both). Default is 0.

=cut

sub _what_or
{
    my ($self, @args) = @_;
    my $strand = shift(@args);
    my ($best, $second);
    my $seq = $self->{'_seqref'};
    $seq = $seq->revcom if ($strand < 0);
    my ($TMP, $filename) = tempfile(DIR => $PATH_TMP, UNLINK => 1);
    for (my $i = 0 ; $i < 3 ; $i++)
    {
        print $TMP ">$i\n",
          $seq->translate(undef, undef, $i, $self->{'_table'})->seq(), "\n";
    }
    if ($strand == 0)
    {
        $seq = $seq->revcom;
        for (my $i = 0 ; $i < 3 ; $i++)
        {
            print $TMP ">$i-\n",
              $seq->translate(undef, undef, $i, $self->{'_table'})->seq(), "\n";
        }
    }
    close $TMP;
    system(  $self->{'_pfam'} . ' '
           . $self->{'_hmm'} . ' '
           . $filename . ' > '
           . $filename
           . '.report');
    eval
    {
        my $hmmer =
          new Bio::Tools::HMMER::Results(-file => $filename . '.report',
                                         -type => 'hmmpfam');
        foreach my $domain ($hmmer->each_Domain)
        {
            $self->{'_hmmor'} = (
                                 [
                                  $domain->evalue,
                                  $domain->hmmname,
                                  substr($domain->seq_id, 0, 1),
                                  $domain->bits,
                                  $domain->start,
                                  $domain->end,
                                  (
                                   ($strand != 0)
                                   ? $strand
                                   : ((length($domain->seq_id) > 1) ? -1 : 1)
                                  )
                                 ]
                                )
              if (!defined $self->{'_hmmor'}[0]
                  || $domain->evalue < $self->{'_hmmor'}[0]);
            $best = $domain->evalue
              if ((!defined $best) || $domain->evalue < $best);
            $second = $domain->evalue
              if ($domain->evalue > $best
                  && (!defined $second || $domain->evalue < $second));
        }
        $self->{'_hmmor'}[7] = $second
          if (defined $best && defined $self->{'_hmmor'}[0]);
    };
    unlink($filename);
    unlink($filename . '.report');
    return (defined $self->{'_hmmor'}[0]) ? 1 : 0;
}

=head2 _find_orf

 Title   : _find_orf
 Usage   : my $bool = $self->_find_or( $strand, $start, $end );
 Function: Retrieve the olfactory receptor ORF.
 Returns : boolean.
 Args    : $strand (mandatory) strand where ORA have been found
           (1 direct or -1 reverse),
           $start (mandatory) coordinate of the first nucleotide,
           $end (mandatory) coordinate of the last nucleotide.

=cut

sub _find_orf
{
    my ($self, @args) = @_;
    my ($strand, $start, $end) = @args;
    my ($position1, $position2);
    if ((defined $self->{'_hmmor'}) && (defined $self->{'_hmmor'}[0]))
    {
        my ($begin, $stop) = $self->_translation();
        my $seq = $self->{'_seqref'};
        $seq = $seq->revcom if ($strand < 0);
        $seq = $seq->seq();
        my ($i, $j) = 0;
        $self->{'_no5'} = 1 if ($self->{'_hmmor'}[4] < 9);
        my $mydna = substr(
                           $seq, 0,
                           (
                            $self->{'_hmmor'}[4] * 3 +
                              abs($self->{'_hmmor'}[2]) + 9
                           )
                          );
        if ($mydna =~ m/(($begin)((?!($stop|$begin))(.{3}))*?)$/o)
        {
            $position1 = length($`);
        }
        else
        {
            if ($mydna =~ m/(($stop)((?!($stop))(.{3}))*)$/o)
            {
                $position1 = length($`) + 3;
                $self->{'_noaug'} = 1;
            }
            elsif ($mydna =~ m/^(.{0,2})((.{3})*)$/o)
            {
                $position1 = length($1);
                $self->{'_noaug'} = 2;
            }
        }
        if (defined $position1)
        {
            my $mydna = substr($seq, $position1);
            my ($coord, $dna, $stopcodon);
            if ($mydna =~ m/^(.{3})(((?!($stop))(.{3})){194,})($stop)/o)
            {
                $dna       = $1 . $2;
                $stopcodon = $6;
            }
            elsif ($mydna =~ m/^(.{3})(((?!($stop))(.{3})){194,})(.{0,2})$/o)
            {
                $self->{'_no3'} = 1;
                $dna            = $2;
                $stopcodon      = '';
            }
            elsif ( $mydna =~ m/($stop)(((?!($stop))(.{3})){194,})($stop)/o ) {
                $self->{'_noaug'}   = 1;
                $position1          = length($`) + 3;
                $dna                = $2;
                $stopcodon          = $6;
            }
            else 
            {
                $dna = substr(
                           $seq,
                           ( $self->{'_hmmor'}[4] * 3 ) - 3,
                           ( ( $self->{'_hmmor'}[5] - $self->{'_hmmor'}[4] + 1 ) * 3 )
                          );
                $dna =~ m/^((.{3})*).{0,2}$/o;
                $dna       = $1;
                $position1 = $self->{'_hmmor'}[4] * 3;
                $self->{'_pseudo'} = 1;
                $stopcodon = '';
            }
            if (defined $dna)
            {
                if ($self->{'_hmmor'}[6] > 0)
                {
                    $position2 =
                      $start + $position1 + length($dna . $stopcodon) - 1;
                    $position1 = $start + $position1;
                    $coord =
                        ((defined $self->{'_no5'}) ? '<' : '')
                      . $position1 . '..'
                      . $position2
                      . ((defined $self->{'_no3'}) ? '>' : '');
                }
                else
                {
                    $position2 = $end - $position1;
                    $position1 =
                      $end - $position1 - length($dna . $stopcodon) + 1;
                    $coord =
                        'complement('
                      . ((defined $self->{'_no3'}) ? '<' : '')
                      . $position1 . '..'
                      . $position2
                      . ((defined $self->{'_no5'}) ? '>' : '') . ')';
                }
                if (defined $self->{'_frag'})
                {
                    $self->{'_pseudo'} = 0
                      if (($position2 - $position1 + 1) < $self->{'_frag'});
                }
                else
                {
                    $self->{'_pseudo'} = 0
                      if (($position2 - $position1 + 1) < 600);
                    $self->{'_pseudo'} = 0
                      if ( (($position2 - $position1 + 1) < 900)
                        && ($self->{'_hmmor'}[1] ne 'OR7')
                        && !(defined $self->{'_no5'} || defined $self->{'_no3'})
                      );
                    $self->{'_pseudo'} = 0
                      if ( (($position2 - $position1 + 1) < 800)
                        && ($self->{'_hmmor'}[1] eq 'OR7')
                        && !(defined $self->{'_no5'} || defined $self->{'_no3'})
                      );
                }
                $self->{'_result'} = (
                    [
                     $coord,
                     $position1,
                     $position2,
                     $self->{'_hmmor'}[6],
                     'Olfactory Receptor, family '
                       . substr($self->{'_hmmor'}[1], 2),
                     $self->{'_hmmor'}[1],
                     (
                      defined $self->{'_pseudo'}
                      ? 'Pseudogene (' . $self->{'_pseudo'} . '); '
                      : ''
                     )
                       . (
                         (defined $self->{'_no5'})
                         ? 'The sequence seems incomplete, 5\' of the CDS is missing; '
                         : ''
                       )
                       . (
                         (defined $self->{'_no3'})
                         ? 'The sequence seems incomplete, 3\' of the CDS is missing; '
                         : ''
                       )
                       . (
                         !(defined $self->{'_pseudo'})
                           && (defined $self->{'_noaug'})
                         ? 'The position of the initiation codon is not identified; '
                         : ''
                       )
                       . 'HMM for family '
                       . $self->{'_hmmor'}[1] . ': '
                       . $self->{'_hmmor'}[0] . ' ('
                       . ($self->{'_hmmor'}[7] - $self->{'_hmmor'}[0]) . ')',
                     $dna . $stopcodon,
                     Bio::Seq->new(-seq => $dna, alphabet => 'dna')
                       ->translate(undef, undef, undef, $self->{'_table'})
                       ->seq()
                    ]
                ) if (defined $coord);
            }
        }
    }
    return (defined $position2) ? 1 : 0;
}

=head2 getHits

 Title   : getHits
 Usage   : my @hits = Bio::ORA->getHits( $seq, $evalue, $ref );
 Function: Quick localization of ORs (use FASTA).
 Returns : Array of hits start/stop positions.
 Args    : $seq (mandatory) primarySeqI object (dna or rna),
           $evalue (mandatory) det the E-value threshold,
           $ref (optional) path to fasta reference file, by default ORA
             look at ./or.fasta.

=cut

sub getHits
{
    my ($self, @args) = @_;
    my ($seq, $evalue, $ref) = @args;
    $ref = ((defined $ref) ? $ref : $PATH_REF);
    $evalue = ((defined $evalue) && ($evalue > 0)) ? $evalue : 1;
    my ($TMP, $filename) = tempfile(DIR => $PATH_TMP, UNLINK => 1);
    my $fasta =
      ((defined $ENV{'FASTADIR'})
        ? $ENV{'FASTADIR'}
        : Bio::ORA->_findexec('tfasta34_t'));
    my @hits;
    if ((defined $seq) && (-x $fasta))
    {
        print $TMP ">query\n", $seq->seq(), "\n\n";
        close($TMP);
        system(  $fasta
               . ' -Q -b 1 -d 1 -H '
               . $ref . ' '
               . $filename . '> '
               . $filename
               . '.report');
        eval
        {
            my $in = new Bio::SearchIO(-format => 'fasta',
                                       -file => $filename . '.report');
            while (my $result = $in->next_result)
            {
                while (my $hit = $result->next_hit)
                {
                    while (my $hsp = $hit->next_hsp)
                    {
                        push(
                             @hits,
                             (
                                  $hit->strand('query') . '|'
                                . sprintf("%09d",$hsp->start('hit')) . '|'
                                . sprintf("%09d",$hsp->end('hit'))
                             )
                            ) if ($hsp->evalue <= $evalue);
                    }
                }
            }
        };
        unlink($filename . '.report');
        unlink($filename);
    }
    if ($#hits >= 0)
    {
        @hits = sort @hits;
        my @hits2 = ();
        for (my $i = 0 ; $i < $#hits ; $i++)
        {
            my ($hitstrand, $hitstart, $hitend) = split(/\|/, $hits[$i]);
            my $loop = 0;
            do
            {
                $loop = 0;
                my ($hitstrand_n, $hitstart_n, $hitend_n) =
                  split(/\|/, $hits[$i + 1])
                  if (defined $hits[$i + 1]);

                #echo ?
                if (
                    (defined $hits[$i + 1])
                    && (
                        (abs($hitstart_n - $hitend) < 1000)
                        || (abs($hitstart - $hitend_n) < 1000)
                        || ( ($hitstart_n < $hitend) && ($hitstart_n > $hitstart) )
                        || ( ($hitend_n < $hitend) && ($hitend_n > $hitstart) )
                        || ( ($hitstart < $hitend_n) && ($hitstart > $hitstart_n) )
                        || ( ($hitend < $hitend_n) && ($hitend > $hitstart_n) )
                       )
                    && ($hitstrand == $hitstrand_n)
                   )
                {
                    $i++;
                    $hitend   = $hitend_n   if ($hitend < $hitend_n);
                    $hitstart = $hitstart_n if ($hitstart > $hitstart_n);
                    $loop     = 1;
                }
            } until ($loop == 0);
            push(@hits2, ($hitstrand. '|' . sprintf("%09d",$hitstart) . '|' . sprintf("%09d",$hitend) ));
        }
        @hits = @hits2;
    }
    return @hits;
}

=head2 getHits

 Title   : getHits
 Usage   : my @hits = Bio::ORA->getHits( $seq, $evalue, $ref );
 Function: Quick localization of ORs (use FASTA).
 Returns : Array of hits start/stop positions.
 Args    : $seq (mandatory) primarySeqI object (dna or rna),
           $evalue (mandatory) det the E-value threshold,
           $ref (optional) path to fasta reference file, by default ORA
             look at ./or.fasta.

=cut

sub fastScan
{
    my ($self, @args) = @_;
    my ($seqfile, $ref) = @args;
    $ref = ((defined $ref) ? $ref : $PATH_REF);
    my ($TMP, $filename) = tempfile(DIR => $PATH_TMP, UNLINK => 1);
    close($TMP);
    my $fasta =
      ((defined $ENV{'FASTADIR'})
        ? $ENV{'FASTADIR'}
        : Bio::ORA->_findexec('fastx34_t'));
    my @hits;
    if ((defined $seqfile) && (-x $fasta))
    {
        system(  $fasta
               . ' -b 1 -d 1 -E 1 -H -Q '
               . $seqfile . ' '
               . $ref . '> '
               . $filename );
        eval
        {
            my $in = new Bio::SearchIO(-format => 'fasta',
                                       -file => $filename);
            while (my $result = $in->next_result)
            {
                push( @hits, $result->query_name() );
            }
        };
        unlink($filename);
    }
    return @hits;
}

=head2 show

 Title   : show
 Usage   : $ORA_obj->show( $outstyle );
 Function: Print result in various style.
 Returns : none.
 Args    : $outstyle (mandatory) 'fasta', 'genbank', 'cvs', 'xml' or 'R' style.


=cut

sub show
{
    my ($self, @args) = @_;
    my $out = shift(@args);
    $out = ((defined $out) ? $out : 'fasta');
    if ($out eq 'xml-begin')
    {
        print "<orml version=\"0.9\">\n";
        print " <analysis>\n";
        print "  <program>\n";
        print "   <prog-name>ORA.pm</prog-name>\n";
        print "   <prog-version>$VERSION</prog-version>\n";
        print "  </program>\n";
        print "  <date>\n";
        print '   <day>', (gmtime)[3], "</day>\n";
        print '   <month>', (gmtime)[4] + 1,    "</month>\n";
        print '   <year>',  (gmtime)[5] + 1900, "</year>\n";
        print "  </date>\n";
        print "  <parameter>\n";
        print '   <evalue>', shift(@args), "</evalue>\n";
        print '   <table>',  shift(@args), "</table>\n";
        print "  </parameter>\n";
        print " </analysis>\n";
    }
    elsif ($out eq 'xml-end') { print "</orml>\n"; }
    elsif ((defined $self->{'_result'}) && (defined $self->{'_result'}[8]))
    {
        my $detail = shift(@args);
        if ($out eq 'genbank')
        {
            my @dated = localtime(time);
            my %month = (
                         0  => 'JAN',
                         1  => 'FEB',
                         2  => 'MAR',
                         3  => 'APR',
                         4  => 'MAY',
                         5  => 'JUN',
                         6  => 'JUL',
                         7  => 'AUG',
                         8  => 'SEP',
                         9  => 'OCT',
                         10 => 'NOV',
                         11 => 'DEC'
                        );
            printf(
                "\nLOCUS       %-20s %7d bp            linear   UNA %02d-%3s-%04d\n",
                $self->{'_seqref'}->display_id,
                $self->{'_seqref'}->length(),
                $dated[3], $month{$dated[4]}, $dated[5] + 1900
            );
            print 'ACCESSION   ', $self->{'_seqref'}->display_id, "\n";
            print 'DEFINITION  ', $self->{'_result'}[5], ".\n";
            print
              "KEYWORDS    .\nSOURCE      Unknown.\n  ORGANISM  Unknown\n            Unclassified.\n";
            print 'COMMENT     Method: ORA v', $VERSION, ".\n";
            print "FEATURES             Location/Qualifiers\n";
            print '     source          1..', $self->{'_seqref'}->length(),
              "\n";
            print '     gene            ',
              (
                ($self->{'_result'}[3] < 0)
                ? ('complement(<'
                   . $self->{'_result'}[1] . '..'
                   . $self->{'_result'}[2] . '>)')
                : (  '<'
                   . $self->{'_result'}[1] . '..'
                   . $self->{'_result'}[2] . '>')
              ),
              "\n";
            print '                     /locus_tag="',
              $self->{'_seqref'}->display_id, "\"\n";
            print '                     /gene="', $self->{'_result'}[5], "\"\n";
            print "                     /pseudo\n"
              if (defined $self->{'_pseudo'});
            print '                     /inference="', $self->{'_hmmor'}[1],
              ' family: ', $self->{'_hmmor'}[0], ' (',
              ($self->{'_hmmor'}[7] - $self->{'_hmmor'}[0]), ")\"\n";

            if ($detail)
            {
                my $dna = $self->{'_result'}[7] . '"';
                print '                     /dna="', substr($dna, 0, 52), "\n";
                my $i = 0;
                while (length($dna) > (($i * 58) + 52))
                {
                    print '                     ',
                      substr($dna, (($i++ * 58) + 52), 58), "\n";
                }
            }
            print '     CDS             ', $self->{'_result'}[0], "\n";
            print '                     /locus_tag="',
              $self->{'_seqref'}->display_id, "\"\n";
            print '                     /gene="', $self->{'_result'}[5], "\"\n";
            my $note = $self->{'_result'}[6] . '"';
            print '                     /note="', substr($note, 0, 51), "\n";
            my $i = 0;
            while (length($note) > (($i * 58) + 51))
            {
                print '                     ',
                  substr($note, (($i++ * 58) + 51), 58), "\n";
            }
            print "                     /pseudo\n"
              if (defined $self->{'_pseudo'});
            print "                     /codon_start=1\n";
            print '                     /transl_table=', $self->{'_table'},
              "\n";
            print '                     /product="', $self->{'_result'}[4],
              ((defined $self->{'_pseudo'}) ? ', pseudogene' : ''), "\"\n";
            my $translation_issue = $self->{'_result'}[8] . '"';
            print '                     /translation="',
              substr($translation_issue, 0, 44), "\n";
            $i = 0;

            while (length($translation_issue) > (($i * 58) + 44))
            {
                print '                     ',
                  substr($translation_issue, (($i++ * 58) + 44), 58), "\n";
            }
            print "ORIGIN\n";
            my $dna = $self->{'_seqref'}->seq();
            my $j   = 0;
            while (length($dna) > ($j * 60))
            {
                $i = 0;
                printf('  %7d ', $j * 60 + 1);
                while (length($dna) > ($j * 60 + $i * 10) && $i < 6)
                {
                    print substr($dna, ($j * 60 + $i++ * 10), 10), ' ';
                }
                print "\n";
                $j++;
            }
            print "//\n";
        }
        elsif ( $out eq 'tbl' ) {    # 'Feature Table' output (for tbl2asn, NCBI)
            print '>Feature ', $self->{'_seqref'}->display_id, "\n";
            print '<', 1 , "\t>", $self->{'_seqref'}->length(), "\tgene\n";
            print "\t\t\tgene\t", $self->{'_result'}[5], "\n";
            print "\t\t\tlocus_tag\t", $self->{'_seqref'}->display_id, "\n";
            print "\t\t\tpseudo\n" if ( defined $self->{'_pseudo'} );
            print "\t\t\tnote\t", $self->{'_hmmor'}[1], ' family: ', $self->{'_hmmor'}[0], "\n";
            unless ( defined $self->{'_pseudo'} ) {
                print '<', 1 , "\t>", $self->{'_seqref'}->length() , "\tCDS\n";
                print "\t\t\tcodon_start\t", ((($self->{'_result'}[1])%3 == 0)? 3 :($self->{'_result'}[1])%3),"\n";
                print "\t\t\ttransl_table\t", $self->{'_table'}, "\n";
                print "\t\t\tproduct\t", $self->{'_result'}[4], "\n";
            }
            print "\n";
            print STDERR "\n>", $self->{'_seqref'}->display_id, ' [organism=unknown] ', $self->{'_result'}[4], " gene member, partial cds.\n";
            my $i = 0;
            my $dna = $self->{'_seqref'}->seq();
            while ( length( $dna ) > ( $i * 80 ) ) { print STDERR substr( $dna, ( $i++ * 80 ), 80 ), "\n"; }
        }
        elsif ($out eq 'csv')
        {    # CSV output
            print $self->{'_seqref'}->display_id, ',',
              ((defined $self->{'_pseudo'}) ? 'N' : 'Y'), ',',
              ($self->{'_result'}[2] - $self->{'_result'}[1] + 1), ',',
              substr($self->{'_hmmor'}[1], 2), ',', $self->{'_hmmor'}[0], ',',
              $self->{'_hmmor'}[7], ',', $self->{'_result'}[7], ',',
              $self->{'_result'}[8], "\n";
        }
        elsif ($out eq 'r')
        {    # R output
            print $self->{'_seqref'}->display_id, "\t",
              ((defined $self->{'_pseudo'}) ? '0' : '1'), "\t",
              ($self->{'_result'}[2] - $self->{'_result'}[1] + 1), "\tOR",
              sprintf("%02d", substr($self->{'_hmmor'}[1], 2)), "\t",
              $self->{'_hmmor'}[0], "\t",
              ((defined $self->{'_pseudo'}) ? $self->{'_pseudo'} : '0'), "\n";
        }
        elsif ($out eq 'xml')
        {
            print ' <sequence id="', $self->{'_seqref'}->display_id,
              ".seq\">\n";
            print "  <input>\n";
            print '   <seq type="dna" length="', $self->{'_seqref'}->length,
              '">', $self->{'_seqref'}->seq, "</seq>\n";
            print "  </input>\n";
            print '  <output id="', $self->{'_seqref'}->display_id, "\">\n";
            print '   <gene id="',  $self->{'_seqref'}->display_id, ".1\">\n";
            print '    <coord',
              (defined $self->{'_no5'}) ? ' 5prime="missing"' : '',
              (defined $self->{'_no3'}) ? ' 3prime="missing"' : '', '>',
              (
                ($self->{'_result'}[3] < 0)
                ? ('complement('
                   . ((defined $self->{'_no3'}) ? '<' : '')
                   . $self->{'_result'}[1] . '..'
                   . $self->{'_result'}[2]
                   . ((defined $self->{'_no5'}) ? '>' : '') . ')')
                : (  ((defined $self->{'_no5'}) ? '<' : '')
                   . $self->{'_result'}[1] . '..'
                     . $self->{'_result'}[2]
                     . ((defined $self->{'_no3'}) ? '>' : ''))
              ),
              "</coord>\n";
            print '    <name>', $self->{'_result'}[5], "</name>\n";
            print '    <seq type="dna" length="', length($self->{'_result'}[7]),
              '">', $self->{'_result'}[7], "</seq>\n";
            print "   </gene>\n";
            print '   <cds id="', $self->{'_seqref'}->display_id, ".2\">\n";
            print '    <coord',
              (defined $self->{'_noaug'}) ? ' start="unknown"'  : '',
              (defined $self->{'_no5'})   ? ' 5prime="missing"' : '',
              (defined $self->{'_no3'})   ? ' 3prime="missing"' : '', '>',
              $self->{'_result'}[0], "</coord>\n";
            print '    <name>',    $self->{'_result'}[5], "</name>\n";
            print '    <note>',    $self->{'_result'}[6], "</note>\n";
            print '    <product>', $self->{'_result'}[4], "</product>\n";
            print '    <seq type="prt" length="', length($self->{'_result'}[8]),
              '">', $self->{'_result'}[8], "</seq>\n";
            print '    <model hmm="', $self->{'_hmmor'}[1], '">',
              $self->{'_hmmor'}[0], "</model>\n";
            print "   </cds>\n";
            print "  </output>\n";
            print " </sequence>\n";
        }
        else
        {    # fasta format
            print "\n>", $self->{'_seqref'}->display_id, '|',
              $self->{'_result'}[5],
              ((defined $self->{'_pseudo'}) ? '|PSEUDOGENE' : ''), "\n";
            my $i = 0;
            while (length($self->{'_result'}[7]) > ($i * 80))
            {
                print substr($self->{'_result'}[7], ($i++ * 80), 80), "\n";
            }
        }
    }
}

=head2 _translation

 Title   : _translation
 Usage   : my ( $start, $end ) = $self->_translation();
 Function: format initiation and stop codons for regex.
 Returns : array with initiation and stop codons.
 Args    : none.

=cut

sub _translation
{
    my $self         = shift;
    my @table        = ('A', 'T', 'C', 'G');
    my @var          = ();
    my @var2         = ();
    my $var_i        = 0;
    my $var2_i       = 0;
    my $myCodonTable = Bio::Tools::CodonTable->new(-id => $self->{'_table'});
    for (my $i = 0 ; $i < 4 ; $i++)
    {

        for (my $j = 0 ; $j < 4 ; $j++)
        {
            for (my $k = 0 ; $k < 4 ; $k++)
            {
                $var[$var_i++] = $table[$i] 
                  . $table[$j]
                  . $table[$k]
                  if $myCodonTable->is_start_codon(
                                            $table[$i] . $table[$j] . $table[$k]
                  );
                $var2[$var2_i++] = $table[$i] 
                  . $table[$j]
                  . $table[$k]
                  if $myCodonTable->is_ter_codon(
                                            $table[$i] . $table[$j] . $table[$k]
                  );
            }
        }
    }
    @var = ('ATG') if ($self->{'_aug'});
    return (join('|', @var), join('|', @var2));
}

1;
