#!/usr/bin/perl

use strict;
use Data::Dumper;

use FindBin;
use Cwd 'realpath';

my $base_dir_HLA_PRG_LA = $FindBin::Bin . '/../../';
$base_dir_HLA_PRG_LA = realpath( $base_dir_HLA_PRG_LA );
die unless(-e $base_dir_HLA_PRG_LA);

my $output_dir = $base_dir_HLA_PRG_LA . '/src/temp/inferHLAonCluster';
die unless(-e $output_dir);

my $cmds_fn = $output_dir . '/qsub_commands.txt';
open(COMMANDS, '>', $cmds_fn) or die "Cannot open $cmds_fn";

# WTSI only
# my @dirs = grep {$_ =~ /WTSI/} glob('/gpfs1/well/gsk_hla/PRG_Remapping/BAMs/*');
my @dirs = glob('/gpfs1/well/gsk_hla/PRG_Remapping/BAMs/*');

foreach my $dir (@dirs)
{
	die unless($dir =~ /.+[\\\/](.+?)$/);
	my $sampleID_full = $1;
	
	my $BAM = $dir . '/merged.bam';
	unless(-e $BAM)
	{
		warn "BAM for sample $sampleID_full not there [${BAM}]";
		next;
	}
	die unless(-e $BAM.'.bai');
	
	my @sampleID_parts = split(/_/, $sampleID_full);
	my $sampleID_reduced = $sampleID_parts[0];
	
	my $cmd = qq(${base_dir_HLA_PRG_LA}/bin_cluster3/MHC-PRG-2 --action HLA --sampleID $sampleID_reduced --BAM $BAM --outputDirectory ${base_dir_HLA_PRG_LA}/working/${sampleID_full});

	my $output_fn = $output_dir . '/' . $sampleID_full . '.bash';
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
print OUTPUT qq(#!/bin/bash
#\$ -P mcvean.prjb -q long.qb
#\$ -pe shmem 10
cd ${base_dir_HLA_PRG_LA}/src
$cmd
);
	close(OUTPUT);
	
	print COMMANDS 'qsub '.$sampleID_full . '.bash' . "\n";
}

close(COMMANDS);

print "Generated file $cmds_fn - execute commands specified therein\n";