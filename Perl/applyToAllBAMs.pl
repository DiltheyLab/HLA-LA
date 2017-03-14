#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

use FindBin;
use Cwd 'realpath';

my $base_dir_HLA_PRG_LA = $FindBin::Bin . '/../../';
$base_dir_HLA_PRG_LA = realpath( $base_dir_HLA_PRG_LA );
die unless(-e $base_dir_HLA_PRG_LA);

my $submission_dir = $base_dir_HLA_PRG_LA . '/src/temp/inferHLAonCluster';
die unless(-e $submission_dir);

my $cmds_fn = $submission_dir . '/qsub_commands.txt';
open(COMMANDS, '>', $cmds_fn) or die "Cannot open $cmds_fn";

# validationBAMs.txt comes from C:\Users\diltheyat\Documents\Oxford\documents\analysis\18 Januar 2017
open(F, '<', '../validationBAMs.txt') or die;
while(<F>)
{
	my $l = $_;
	chomp($l);
	next if($l =~ /^\s+$/);
	next unless($l);
	my @fields = split(/\t/, $l);
	my $cohort = $fields[0];
	my $path = $fields[1];
	my $sampleID = $fields[2];
	
	die "Path $path not existing" unless(-e $path);
	die "BAM not indexed" unless(-e $path . '.bai');
	
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
		die unless($path =~ /temp_mapping\/(.+?)\/merged.bam/);
		$sampleID = $1;
	}
	else
	{
		die "Unknown cohort $cohort";
	}
	
	print $sampleID, "\n";
	
	my $captureOutput = $submission_dir . '/' . $sampleID . '.output';
	
	my $cmd = qq(/usr/bin/time -v ./inferHLATypes.pl --BAM $path --graph PRG_MHC_GRCh38_withIMGT --sampleID $sampleID --maxThreads 7 &> $captureOutput);
	
	my $sampleID_for_jobName = 'J_' . $sampleID;
	#$sampleID_for_jobName =~ s/[^A-Za-z0-9]//g;

	my $output_fn = $submission_dir . '/' . $sampleID . '.bash';
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
print OUTPUT qq(#!/bin/bash
#\$ -P mcvean.prjb -q long.qb
#\$ -N $sampleID_for_jobName
#\$ -pe shmem 7
cd ${base_dir_HLA_PRG_LA}/src
$cmd
);
	print COMMANDS 'qsub '.$sampleID . '.bash' . "\n";

}
close(F);


close(COMMANDS);

print "Generated file $cmds_fn - execute commands specified therein\n";
