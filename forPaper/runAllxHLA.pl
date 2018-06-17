#!/usr/bin/perl

# original location: /data/projects/phillippy/projects/MHC/HLA-PRG-LA/src/forPaper/runAllKourami.pl

use strict;
use warnings;
use Data::Dumper;

use FindBin;
use Cwd qw/realpath abs_path/;

my $base_dir_HLA_PRG_LA = $FindBin::Bin . '/../../';
$base_dir_HLA_PRG_LA = realpath( $base_dir_HLA_PRG_LA );
die unless(-e $base_dir_HLA_PRG_LA);

my $submission_dir = $base_dir_HLA_PRG_LA . '/src/temp/inferXHLAOnCluster';
unless(-e $submission_dir)
{
	mkdir($submission_dir) or die "Cannot mkdir $submission_dir";
}
die unless(-e $submission_dir);

my $cmds_fn = $submission_dir . '/qsub_commands.txt';
open(COMMANDS, '>', $cmds_fn) or die "Cannot open $cmds_fn";

my %sampleIDs;
open(F, '<', '../temp/validationBAMs.txt') or die;
while(<F>)
{
	my $l = $_;
	chomp($l);
	next if($l =~ /^\s+$/);
	next unless($l);
	my @fields = split(/\t/, $l, -1);
	my $cohort = $fields[0];
	my $path = $fields[1];
	my $sampleID = $fields[2];

	if($cohort eq 'Platinum')
	{
		die unless(defined $sampleID);
	}	
	elsif($cohort eq '1000G')
	{
		die unless($path =~ /.+\/(.+?)\./);
		$sampleID = '1000G_' . $1;
	}
	elsif($cohort eq 'MiSeq')
	{
		die unless($path =~ /AM_(.+?)\.bam/);
		$sampleID = $1;
	}
	elsif($cohort eq 'HapMap')
	{
		die unless($path =~ /HapMap_Exomes\/BAMs\/(.+?).bam/);
		$sampleID = $1;
	}
	elsif($cohort eq '')
	{
		#die Dumper(\@fields, $sampleID, $cohort, $path);	
		$sampleID = $fields[1];
		$path = $fields[2];
	}
	elsif($cohort =~ /WTSI/)
	{
		die unless($path =~ /.+(WTSI\d+)\.bam/);
		$sampleID = $cohort . '_' . $1;
	}		
	else
	{
		$sampleID = $cohort;
	}
	
	unless(-e $path)
	{	
		warn "Path $path not existing - skip.";
		next;
	}
	
	unless((-e $path . '.bai') or (-e $path . '.crai'))
	{
		warn "BAM/CRAM $path not indexed - skip.";
		next;
	}	
	
	if($cohort eq '1000G')
	{
		# $path = '/data/projects/phillippy/projects/MHC/1000G_highCoverage_250bp/B38_fromExtraction/' . $sampleID . '.bam';
	}
	
	unless(-e $path)
	{	
		warn "Path $path not existing - skip. (II)";
		next;
	}
	
	#else
	#{
	#	die "Unknown cohort $cohort";
	#}
	
	die if($sampleIDs{$sampleID});
	$sampleIDs{$sampleID}++;
	
	print $sampleID, "\n";
	
	my $captureOutput = $submission_dir . '/' . $sampleID . '.output';
	
	my $cmd = qq(/usr/bin/time -v perl runxHLA.pl --BAM $path --sampleID $sampleID &> $captureOutput);
	
	my $sampleID_for_jobName = 'J_' . $sampleID;
	#$sampleID_for_jobName =~ s/[^A-Za-z0-9]//g;

	my $output_fn = $submission_dir . '/' . $sampleID . '.bash';
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
print OUTPUT qq(#PBS -l select=1:ncpus=1:mem=32GB
#PBS -l walltime=23:00:00
#PBS -A 'IMMGEN'
#PBS -N xHLA${sampleID}
#PBS -r y
cd ${base_dir_HLA_PRG_LA}/src
cd forPaper
$cmd
);
	
	print COMMANDS 'qsub ' . abs_path($output_fn) . "\n";

}
close(F);


close(COMMANDS);

print "Generated file $cmds_fn - execute commands specified therein\n";



