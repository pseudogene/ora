use Module::Build;

my $build = Module::Build->new
    ( module_name     => 'Bio::ORA',
      dist_name       => 'ora',
      dist_version    => '1.9.1',
      dist_author     => 'Michael Bekaert <michael.bekaert@stir.ac.uk>',
      dist_abstract   => 'Olfactory Receptor family Assigner (ORA)',
      license         => 'gpl',
      configure_requires => { 'Module::Build' => 0.40 },
      requires        => {
                            'perl'                      => '5.6.1',
                            'Bio::Root::Version'        => 0,
                            'Getopt::Long'              => 0,
                            'IO::File'                  => 0,
                            'Cwd'                       => 0
                         },
      script_files     => { 'scripts/or.pl' => 'or.pl',
                            'scripts/or.fasta' => 'or.fasta',
                            'scripts/or.hmm' => 'or.hmm',
                            'scripts/or.hmm.h3f' => 'or.hmm.h3f',
                            'scripts/or.hmm.h3i' => 'or.hmm.h3i',
                            'scripts/or.hmm.h3m' => 'or.hmm.h3m',
                            'scripts/or.hmm.h3p' => 'or.hmm.h3p'
      	                 },
      sign          => 1,
      create_readme => 1,
      create_makefile_pl => 'traditional',
      auto_features => {
        manpage_support => {
          description => "Create Unix man pages",
          requires    => {'Pod::Man' => 0},
        },
        HTML_support => {
          description => "Create HTML documentation",
          requires    => {'Pod::Html' => 0},
        },
      }
    );
$build->add_build_element('hmm');
$build->add_build_element('fasta');
$build->create_build_script;
