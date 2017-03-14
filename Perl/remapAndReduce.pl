#!/usr/bin/perl

use strict;
use Data::Dumper;
use List::MoreUtils qw/mesh/;
use Getopt::Long;   
use FindBin;
use Sys::Hostname;

# ./remapAndReduce.pl --BAM_in /gpfs1/well/gsk_hla/temp_mapping/PLATINUM_NA12892/merged.bam --BAM_id PLATINUM_NA12892

my $username =  getlogin();

my $max_queued_jobs = 10000;

my $path_to_script = $FindBin::Bin.'/'.$FindBin::Script;
my $samtools_bin = qq(/home/dilthey/samtools-0.1.18/samtools);
my $graphs_dir = qq(/gpfs1/well/gsk_hla/MHC-PRG-2/graphs/);
my $read_extraction_temp = qq(/gpfs1/well/gsk_hla/PRG_Remapping/FASTQ);
my $output_dir = qq(/gpfs1/well/gsk_hla/PRG_Remapping/BAMs);

my $alignReads = qq(/gpfs1/well/gsk_hla/alignGenome2/src/alignReads.pl);
my $reduceBAM = qq(/gpfs1/well/gsk_hla/MHC-PRG-2/src/Perl/reduceBAM2PRG.pl);

my $doExtraction;
my $BAM_in;
my $BAM_id;
my $graph = 'PRG_MHC_GRCh38_withIMGT';

GetOptions (
 'doExtraction:s' => \$doExtraction, 
 'BAM_in:s' => \$BAM_in, 
 'BAM_id:s' => \$BAM_id, 
 'graph:s' => \$graph, 
);

unless($BAM_in and (-e $BAM_in))
{ 
	die "Please provide --BAM_in";
}

unless($BAM_id)
{
	die "Please provide --BAM_id";
}

unless($graph)
{
	die "Please provide --graph";
}

die unless($BAM_id =~ /^[\w]+$/);

my $graph_dir = $graphs_dir . '/' . $graph;
die "Graph directory $graph_dir not existing" unless(-e $graph_dir);

my $extendedReferenceGenome = getOneLine($graph_dir . '/extendedReferenceGenomePath.txt');
die "Extended reference genome path $extendedReferenceGenome not existing" unless(-e $extendedReferenceGenome);

my $alignGenome_name_for_project = $BAM_id . '_RED_' . $graph;
my $query_project_name_cmd = qq($alignReads --action knowProject --name_for_project $alignGenome_name_for_project);
my $query_project_name_output = `$query_project_name_cmd`;

my $readExtractionDir = $read_extraction_temp . '/' . $BAM_id;
my $status_file_read_extraction = $readExtractionDir . '/extractedOK';

my $BAM_out_dir = $output_dir . '/' . $BAM_id;
my $BAM_out = $BAM_out_dir . '/merged.bam';

if(-e $BAM_out)
{
	print "BAM $BAM_out existing, abort.\n";
	exit;
}


if($doExtraction)
{
	my $lenientString = qq(VALIDATION_STRINGENCY=LENIENT);
	my $output_Picard = $readExtractionDir . '/PicardOut.txt';
	my $extraction_command = qq(java -Xmx20g -XX:-UseGCOverheadLimit -jar /apps/well/picard-tools/1.111/SamToFastq.jar ${lenientString} I=${BAM_in} F=${readExtractionDir}/F_1.fastq F2=${readExtractionDir}/F_2.fastq &> $output_Picard);
	execCmd($extraction_command);
	
	writeStatusFile($status_file_read_extraction, 1);
	
	exit;
}

if(not -e $readExtractionDir)
{
	if(-e $readExtractionDir)
	{
		execCmd('rm -rf '.$readExtractionDir);
	}
	mkdir($readExtractionDir) or die "Cannot mkdir $readExtractionDir";
	writeStatusFile($status_file_read_extraction, 0);
	
	invoke_self_and_wait("--doExtraction 1 --BAM_id $BAM_id --BAM_in $BAM_in", $status_file_read_extraction);
	
	die unless(readStatusFile($status_file_read_extraction));
	
	my $cmd_prepare = qq(perl $alignReads --reference_genome_FASTA_path $extendedReferenceGenome --action prepare --input_directory $readExtractionDir --paired_end 1 --gz 0 --name_for_project $alignGenome_name_for_project);
	
	my $try_to_delete_working_directory = 0;
	
	if(execCmd($cmd_prepare, 1))
	{
 
	}
	else
	{
		$try_to_delete_working_directory = 1;
		warn "Command $cmd_prepare could not be executed successfully";
	}
	
	if($try_to_delete_working_directory)
	{
		my $projectDir = getWorkingDirectoryAlignGenome($alignGenome_name_for_project);
		if(-e $projectDir)
		{
			die unless(index($projectDir, '/gpfs1/well/gsk_hla/temp_mapping_2') == 0);
			execCmd("rm -rf $projectDir");
		}
		
		die "For some reason, preparation command $cmd_prepare didn't work out, and we tried to delete $projectDir - perhaps succcessfully";
	}
}
else
{
	print "$readExtractionDir existing, assume prepared and ready!\n";
	die "Directory $readExtractionDir exists, but reads not extracted - might a background job still be running?" unless(readStatusFile($status_file_read_extraction));
}

while(1)
{
	my $nextSteps = getNextSteps($alignGenome_name_for_project);
	if($nextSteps eq 'KICK_OFF_ALIGNMENT')
	{
		my $cmd_kickoff = qq($alignReads --action kickOffAllAlignment --name_for_project $alignGenome_name_for_project);
		execCmd($cmd_kickoff);
	}
	elsif($nextSteps eq 'WAIT_QUEUED')
	{
		print "Status queued, sleep...\n";
		sleep(1200);
	}	
	elsif($nextSteps eq 'KICK_OFF_MERGE')
	{
		my $cmd_kickoff = qq($alignReads --action kickOffMerge --name_for_project $alignGenome_name_for_project);
		execCmd($cmd_kickoff);
	}
	elsif($nextSteps eq 'WAIT_RUNNING')
	{
		print "Status running, sleep...\n";
		sleep(1200);
	}
	elsif($nextSteps eq 'MERGE_RUNNING')
	{
		print "Merging running, sleep...\n";
		sleep(1200);
	}
	elsif($nextSteps eq 'DONE')
	{
		my $projectDir = getWorkingDirectoryAlignGenome($alignGenome_name_for_project);
		my $BAM = $projectDir . '/merged.bam';
		mkdir($BAM_out_dir);
		die unless(-e $BAM_out_dir);
		die "BAM $BAM not there" unless (-e $BAM);
		
		my $cmd_reduce = qq(perl $reduceBAM --BAM_in $BAM --BAM_out $BAM_out --graph $graph_dir);
		execCmd($cmd_reduce);
		
		print "\nCreated reduced BAM $BAM_out \n\n";
		print "Reduction done, now deleting $projectDir\n\n";
		die "Can't delete directory $projectDir" unless(index($projectDir, '/well/gsk_hla/temp_mapping_2') == 0);		
		execCmd("rm -rf $projectDir");
		
		print "Reduction done, now deleting $readExtractionDir (if present)\n\n";
		if(-e $readExtractionDir)
		{
			execCmd("rm -rf $readExtractionDir");		
		}
		
		print "\n\nDone - BAM path $BAM_out \n\n";
		exit 0;
	}
	else
	{
		die "Unknown status string from $alignReads for $alignGenome_name_for_project: $nextSteps";
	}
}

sub getMergedBAM
{
	my $project = shift;
	my $cmd = qq(perl $alignReads --action displayAlignmentStatus --name_for_project $project);
	my $cmd_output = `$cmd`;
	unless($cmd_output =~ /NEXTSTEPS_STATUS: (.+)/)
	{
		die "Could not extract next steps from output:\n\n$cmd_output\n\nCommand: $cmd";
	}
	my $mergedBAM = $1;
	$mergedBAM =~ s/[\n\r]//g;
	return $mergedBAM;
}

sub getWorkingDirectoryAlignGenome
{
	my $project = shift;
	die unless($project);

	my $cmd_getDir = qq($alignReads --action getProjectWorkingDirectory --name_for_project $project);
	my $projectDir = `$cmd_getDir`;
	
	return $projectDir;
}

sub getNextSteps
{
	my $project = shift;
	my $cmd = qq(perl $alignReads --action displayAlignmentStatus --name_for_project $project);
	my $cmd_output = `$cmd`;
	unless($cmd_output =~ /NEXTSTEPS_STATUS: (.+)/)
	{
		die "Could not extract next steps from output:\n\n$cmd_output\n\nCommand: $cmd";
	}
	my $nextsteps = $1;
	$nextsteps =~ s/[\n\r]//g;
	return $nextsteps;
}

sub execCmd
{
	my $cmd = shift;
	my $noDie = shift;
	print "Executing command $cmd\n";
	my $ret = system($cmd);
	unless($ret == 0)
	{
		unless($noDie)
		{
			die "Command $cmd failed";
		}
		return 0;
	}
	return 1;
}

sub getOneLine
{
	my $f = shift;
	open(F, '<', $f) or die "Cannot open $f";
	my $l = <F>;
	chomp($l);
	$l =~ s/[\n\r]//g;
	close(F);
	return $l;
}

sub invoke_self_and_wait
{
	my $arguments_string = shift;
	my $statusfile = shift;
	
	die unless($statusfile);
	
	my $host = hostname();
	my $command = 'perl '.$path_to_script.' '.$arguments_string;
	
	if(($host eq 'birch') or ($host eq 'elm') or ($host eq 'sequoia') or ($host eq 'banyan'))
	{
		print "From host $host, invoke myself with command:\n$command\n";
		print `$command`;
	}
	elsif($host =~ /rescomp/)
	{
		if(get_n_queuedjobs() > $max_queued_jobs)
		{
			die "Too many queued jobs " . get_n_queuedjobs();
		}
		
		my $key_for_dir = genkey();
		my $temp_directory_for_qsub = '/gpfs1/well/gsk_hla/MHC-PRG-2/qsub';		
		while(-e $temp_directory_for_qsub.'/'.$key_for_dir)
		{
			$key_for_dir = genkey();
		}
		mkdir($temp_directory_for_qsub.'/'.$key_for_dir) or die "Cannot create ".$temp_directory_for_qsub.'/'.$key_for_dir;

		my $file_for_submission = $temp_directory_for_qsub.'/'.$key_for_dir."/script.bash";
		open(F, '>', $file_for_submission) or die "Cannot open $file_for_submission";
		print F qq(#!/bin/bash
#\$ -P mcvean.prjb -q long.qb		
#\$ -pe shmem 3		
#\$ -N Subjob_reMapAndReduce
export PERL5LIB=/users/mcvean/dilthey/perl5/lib/perl5:\$PERL5LIB
echo Started analysis at: `date`
$command
echo Finished analysis at: `date`
);
		close(F);
		
		my $qsub_command = qq(qsub $file_for_submission);
		print $qsub_command, "\n";
		my $qsub_return = `$qsub_command`, "\n";
		print $qsub_return, "\n";
		die unless($qsub_return =~ /Your job (\d+)\s+/);
		my $jobID = $1;
		
		unless(get_n_queuedjobs() >= 1)
		{
			die "Weird - we have just submitted a job, and now didn't get a count >= 1 - current count: " . get_n_queuedjobs();
		}
				
		sleep(10);
		
		my $qstat_output = `qstat`;
		while($qstat_output =~ /^\s*$jobID/m)
		{
			print "Waiting for job $jobID ... \n";
			sleep(100);
			$qstat_output = `qstat`;			
		}
		# print "Output didn't contain $jobID\n--${qstat_output}\n--\n\n";
		
		unless(readStatusFile($statusfile))
		{
			die "Command $qsub_command / job $jobID has returned, but status file $statusfile not true!\n";
		}
		
	}
	else
	{
		die "Do not know which host I am on: $host\n";
	}
}


sub genkey
{
	my $key;
	my $num = $_[0] ? $_[0] : 25;
	for (my $i = 0; $i <= $num; $i++)
	{
		my $val =  int(rand(9)+1);
		$val = ($val > 9) ? 9 : $val;
		$key .= $val;
	}
	return $key;
}

sub writeStatusFile
{
	my $file = shift;
	my $status = shift;
	open(F, '>', $file) or die "Cannot open $file";
	print F $status;
	close(F);
}

sub readStatusFile
{
	my $file = shift;
	return unless(-e $file);
	open(F, '<', $file) or die "Cannot open $file";
	my $forReturn = <F>;
	close(F);
	return $forReturn;
}

sub get_n_queuedjobs
{
	die unless($username);
	die unless($username eq 'dilthey'); # This can later be removed, just a sanity check
	my $cmd = qq(qstat -u $username);
	my $return = `$cmd`;
	my @lines = split(/\n/, $return);
	return scalar(@lines);
}
