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
	my $technology = $fields[3];
	
	next if($technology);
	next if($l =~ /d1000/);
	next if($l =~ /dPlatinum/);
	
	if($cohort eq 'Platinum')
	{
		die unless(defined $sampleID);
	}	
	elsif($cohort eq '1000G')
	{
		die unless($path =~ /.+\/(.+?)\./);
		$sampleID = $1;
		$path =~ s/1000G_highCoverage_250bp_ref1000G\//1000G_highCoverage_250bp_ref1000G\/1000G_/;		
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
		$path =~ s/BAMs\//BAMs\/HapMap_/;
	}
	elsif($cohort eq '')
	{
		#die Dumper(\@fields, $sampleID, $cohort, $path);	
		$sampleID = $fields[1];
		$path = $fields[2];
		$sampleID =~ s/_Platinum//;
		$path = '/gpfs/project/dilthey/projects/MHC/CRAMs_GRCh38_primary/Platinum_' . $sampleID . '.cram';
		$cohort = 'Platinum';
	}
	elsif($cohort =~ /WTSI/)
	{
		next;
		die unless($path =~ /.+(WTSI\d+)\.bam/);
		$sampleID = $cohort . '_' . $1;
	}		
	else
	{
		$sampleID = $cohort;
	}
	
	$path =~ s/\/data\/projects\/phillippy\/projects\/MHC\/HapMap_Exomes\/BAMs\//\/gpfs\/project\/dilthey\/projects\/MHC\/CRAMs_GRCh38_primary\//;
	$path =~ s/\/data\/projects\/phillippy\/projects\/MHC\/1000G_highCoverage_250bp_ref1000G\//\/gpfs\/project\/dilthey\/projects\/MHC\/CRAMs_GRCh38_primary\//;
	$path =~ s/\.bam$/.cram/;
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
	
	unless($sampleID)
	{
		die "Undefined sample ID - line $l";
	}
	
	die if($sampleIDs{$cohort . $sampleID});
	$sampleIDs{$cohort . $sampleID}++;
	
	print $cohort, " ", $sampleID, "\n";
	
	my $captureOutput = $submission_dir . '/' . $sampleID . '.output';
	
	my $cmd = qq(/usr/bin/time -v perl runxHLA.pl --BAM $path --sampleID ${cohort}_${sampleID} &> $captureOutput);
	
	# print $cmd, "\n";
	
	my $sampleID_for_jobName = 'J_' . $sampleID;
	#$sampleID_for_jobName =~ s/[^A-Za-z0-9]//g;

	my $output_fn = $submission_dir . '/' . $sampleID . '.bash';
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
print OUTPUT qq(#PBS -l select=1:ncpus=1:mem=100GB
#PBS -l walltime=05:59:00
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



