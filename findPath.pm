package findPath;

use strict;
use FindBin;

sub check_samtools
{
	my $samtools_bin = shift;
	
	my $samtools_version = `$samtools_bin --version` ;
	die "Can't parse samtools version output" unless($samtools_version =~ /samtools ([\d\.]+)/);
	$samtools_version = $1;
	my $samtools_version_numeric = $samtools_version;
	if($samtools_version_numeric =~ /^(\d+)\.(\d+)\.(\d+)$/)
	{
		$samtools_version_numeric = $1 . '.' . $2 . $3;
	}
	unless($samtools_version_numeric >= 1.3)
	{
		die "I need samtools >=1.3";
	}
}

sub get_working_dir
{
	my $workingDir_param = shift;
	
	my $this_bin_dir = $FindBin::RealBin;
	
	my %paths_ini;
	my $paths_ini = $this_bin_dir . '/paths.ini';
	if(-e $paths_ini)
	{
		open(INI, '<', $paths_ini) or die "Cannot open $paths_ini";
		while(<INI>)
		{
			chomp;
			next unless($_);
			$_ =~ s/[\n\r]//g;
			next if($_ =~ /^\s+$/);
			die "Invalid line format in $paths_ini -- expect all lines to be either empty or key=value pairs" unless($_ =~ /^(.+)=(.*)$/);
			my $id = $1;
			my @alts = split(/,/, $2);
			$paths_ini{$id} = \@alts;
		}
		close(INI);
	}
			
	my $working_dir;
	if($paths_ini{workingDir_HLA_ASM}[0] and not defined $workingDir_param)
	{
		$working_dir = $paths_ini{workingDir_HLA_ASM}[0];
		$working_dir =~ s/\$HLA\-LA\-DIR/$this_bin_dir/;
	}
	else
	{
		unless(defined $workingDir_param)
		{
			die "\n\nPlease specify a working directory via --workingDir.\n\nOutput for sample with ID \$sampleID will go a correspondingly named sub-directory of the working directory.\n\nFor example, if --workingDir is /path/working, and --sampleID is mySample, then the output will go into directory /path/working/mySample.\n\n";
		}
		$working_dir = $workingDir_param;
	}
	
	return $working_dir;
}

sub find_path
{
	my $id = shift;
	my $supplied_value = shift;
	my $forWhich = shift;
	
	my $this_bin_dir = $FindBin::RealBin;
	my %paths_ini;
	my $paths_ini = $this_bin_dir . '/paths.ini';
	if(-e $paths_ini)
	{
		open(INI, '<', $paths_ini) or die "Cannot open $paths_ini";
		while(<INI>)
		{
			chomp;
			next unless($_);
			$_ =~ s/[\n\r]//g;
			next if($_ =~ /^\s+$/);
			die "Invalid line format in $paths_ini -- expect all lines to be either empty or key=value pairs" unless($_ =~ /^(.+)=(.*)$/);
			my $id = $1;
			my @alts = split(/,/, $2);
			$paths_ini{$id} = \@alts;
		}
		close(INI);
	}
		
	if(defined $supplied_value)
	{
		die "Command-line supplied value/file for parameter $id not existing" unless(-e $supplied_value);
		return $supplied_value;
	}
	
	if(exists $paths_ini{$id})
	{
		foreach my $alternative (@{$paths_ini{$id}})
		{
			if(-e $alternative)
			{
				return $alternative;
			}
		}
	}	

	if($forWhich)
	{
		my $which_output = `which $forWhich`;
		$which_output =~ s/[\n\r]//g;
		if($which_output and (-e $which_output))
		{
			return $which_output;
		}
	}
	
	die "I couldn't figure out a path for ${id}. Order of precedence: check for parameter --${id}; check paths.ini in $this_bin_dir for a key named $id; parse, if command string defined, the output of the command 'which ${forWhich}' ('which ' means that the command string is not defined).";
}

1;
