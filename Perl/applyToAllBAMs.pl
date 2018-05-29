#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

use FindBin;
use Cwd qw/realpath abs_path/;

my $base_dir_HLA_PRG_LA = $FindBin::Bin . '/../../';
$base_dir_HLA_PRG_LA = realpath( $base_dir_HLA_PRG_LA );
die unless(-e $base_dir_HLA_PRG_LA);

my $submission_dir = $base_dir_HLA_PRG_LA . '/src/temp/inferHLAonCluster';
die unless(-e $submission_dir);

my $cmds_fn = $submission_dir . '/qsub_commands.txt';
open(COMMANDS, '>', $cmds_fn) or die "Cannot open $cmds_fn";

# validationBAMs.txt comes from C:\Users\diltheyat\Documents\Oxford\documents\analysis\18 Januar 2017
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
	my $technology = $fields[2];
	#my $sampleID = $fields[2];

	my $sampleID;
	if($cohort eq '')
	{
		$sampleID = $fields[1];
		$path = $fields[2];
		$technology = $fields[3];
	}
	else
	{
		if($cohort eq 'Platinum')
		{
			die unless(defined $sampleID);
		}	
		elsif($cohort eq '1000G')
		{
			die unless($path =~ /.+\/(.+?)\.(wgs|mapped)/);
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
		elsif($cohort =~ /WTSI/)
		{
			die unless($path =~ /.+(WTSI\d+)\.bam/);
			$sampleID = $cohort . '_' . $1;
		}	
		else
		{
			die "Unknown cohort $cohort";
		}
	}
	print $sampleID, "\n";
	
	
	die "Path $path not existing" unless(-e $path);
	die "BAM not indexed" unless((-e $path . '.bai') or (-e $path . '.crai'));
		
	my $captureOutput = $submission_dir . '/' . $sampleID . '.output';
	
	my $longReads = (defined $technology) ? "--longReads $technology" : "";
	my $cmd = qq(/usr/bin/time -v perl HLA-PRG-LA.pl --BAM $path --graph PRG_MHC_GRCh38_withIMGT --sampleID $sampleID --maxThreads 7 ${longReads} &> $captureOutput);
	
	my $sampleID_for_jobName = 'J_' . $sampleID;
	#$sampleID_for_jobName =~ s/[^A-Za-z0-9]//g;

	my $output_fn = $submission_dir . '/' . $sampleID . '.bash';
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
print OUTPUT qq(#!/bin/bash
#\$ -q phillippy.q
#\$ -l mem_free=70G
#\$ -N J$sampleID
cd ${base_dir_HLA_PRG_LA}/src
$cmd
);
	
	print COMMANDS 'qsub ' . abs_path($output_fn) . "\n";

}
close(F);


close(COMMANDS);

print "Generated file $cmds_fn - execute commands specified therein\n";
