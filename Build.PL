use Module::Build;

my $build = Module::Build->new
    ( module_name     => 'Bio::ORA',
      dist_name       => 'ora',
      dist_version    => 1.5,
      dist_author     => 'Michael Bekaert <michael@batlab.eu>',
      dist_abstract   => 'Olfactory Receptor family Assigner (ORA)',
      license         => 'gpl',
      requires        => {
                            'perl'                      => '5.6.1',
                            'Bio::Root::Version'        => '1.5.2',
                            'Getopt::Long'              => 0,
                            'IO::File'                  => 0,
                            'Cwd'                       => 0,
                            'Term::ProgressBar'         => 0
                         },
      script_files     => { 'scripts/or.pl' => 'or.pl',
                            'scripts/or.fasta' => 'or.fasta',
                            'scripts/or.hmm' => 'or.hmm'
      	                 },
      sign          => 1,
      create_readme => 1,
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
