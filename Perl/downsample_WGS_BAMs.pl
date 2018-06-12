use strict;
use Data::Dumper;
use List::Util qw/sum/;
use FindBin;


my $outputDir = '/data/projects/phillippy/projects/MHC/15X_BAMs/';
my $targetGigabases = (15 * 3.2e9) / 1e9;

my $qsub_temp_dir = '../temp/qsub_downsample';
mkdir($qsub_temp_dir);

my @inputs;
my @outputCRAMs;

my @G1000_BAMs = glob('/data/projects/phillippy/projects/MHC/1000G_highCoverage_250bp_ref1000G/*.cram');
foreach my $BAM (@G1000_BAMs)
{
	die unless($BAM =~ /.+\/(.+?)\.(.+?)$/);
	my $sampleID = $1;
	push(@inputs, $BAM);
	my $outputCRAM = '/data/projects/phillippy/projects/MHC/15X_BAMs/d1000G_' . $sampleID . '.cram';
	push(@outputCRAMs, $outputCRAM);
}
 

my $base_dir_HLA_PRG_LA = $FindBin::Bin . '/../'; 


push(@inputs, '/data/projects/phillippy/projects/MHC/NA12878_CRAM/NA12878.cram');
push(@inputs, '/data/projects/phillippy/projects/rDNA/alignGenome2/temp/NA12891_Platinum_1000G/merged.bam');
push(@inputs, '/data/projects/phillippy/projects/rDNA/alignGenome2/temp/NA12892_Platinum_1000G/merged.bam');
push(@outputCRAMs, $outputDir . '/dPlatinum_NA12878.cram');
push(@outputCRAMs, $outputDir . '/dPlatinum_NA12891.cram');
push(@outputCRAMs, $outputDir . '/dPlatinum_NA12892.cram');

die unless($#inputs == $#outputCRAMs);


my $samples_file = $outputDir . '/samples.txt';
my $qsub_file = $qsub_temp_dir . '/qsub_generation.txt';
open(SAMPLES, '>', $samples_file) or die "Cannot open $samples_file";
open(QSUB, '>', $qsub_file) or die "Cannot open $qsub_file";

for(my $BAMi = 0; $BAMi <= $#inputs; $BAMi++)
{
	my $input = $inputs[$BAMi];
	my $outputCRAM = $outputCRAMs[$BAMi];
	die "Weird output: $outputCRAM" unless($outputCRAM =~ /.+\/(.+?)\.cram$/);
	my $crai = $outputCRAM . '.crai';
	next if(-e $crai);
	my $outputFileName = $1;
	
	my $job_file = $qsub_temp_dir . '/' . $outputFileName . '.bash';
	print QSUB "qsub $job_file \n";
	
	open(JOB, ">", $job_file) or die "Cannot open $job_file";
	print JOB qq(#!/bin/bash
#\$ -l mem_free=90G
#\$ -N $outputFileName
cd $base_dir_HLA_PRG_LA
cd Perl
perl downsampleBAM.pl --input $input --output $outputCRAM --targetGigabases $targetGigabases
);
	
	close(JOB);

}

print "\nExecute commands in $qsub_file to trigger downsampling for all BAMs.\n";
