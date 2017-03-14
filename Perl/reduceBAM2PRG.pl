#!/usr/bin/perl

use strict;
use Data::Dumper;
use List::MoreUtils qw/mesh/;
use Getopt::Long;   

# ./reduceBAM2PRG.pl --BAM_in /gpfs1/well/gsk_hla/temp_mapping_2/NA19240_PRG/merged.bam --BAM_out /gpfs1/well/gsk_hla/temp_mapping_2/NA19240_PRG_REDUCED/merged.bam
# Command that should not work (for testing):
# ./reduceBAM2PRG.pl --BAM_in /gpfs1/well/gsk_hla/temp_mapping_2/NA12878_GRCH38_ALT/merged.bam --BAM_out bla

my $graph = qq(/gpfs1/well/gsk_hla/MHC-PRG-2/graphs/PRG_MHC_GRCh38_withIMGT);
my $linearALTs = qq(/gpfs1/well/gsk_hla/MHC-PRG-2/linearReferenceALTs/KIRGRCh38);
my $samtools_bin = findFileFromAlternatives('/home/dilthey/samtools-0.1.18/samtools', '/apps/well/samtools/0.1.19/bin/samtools');
my $BAM_in;
my $BAM_out;


GetOptions (
 'graph:s' => \$graph, 
 'linearALTs:s' => \$linearALTs, 
 'BAM_in:s' => \$BAM_in, 
 'BAM_out:s' => \$BAM_out, 
);

unless($graph)
{ 
	die "Please provide --graph";
}

unless($linearALTs)
{ 
	die "Please provide --linearALTs";
}

unless($BAM_in)
{
	die "Please provide --BAM_in";
}

unless($BAM_out)
{
	die "Please provide --BAM_out";
}


my %extract_partial_contig;
my %extract_complete_contig;

my $sequences_file = $graph . '/sequences.txt';
unless(-e $sequences_file)
{
	die "File $sequences_file is not present";
}

my $linearALTs_file = $linearALTs . '/extendedGenome_coveredRegions.txt';
unless(-e $linearALTs_file)
{
	die "File $linearALTs_file is not present";
}

my $process_linear_ALTs_file = sub {
	my $file = shift;

	open(SEQUENCES, '<', $file) or die "Cannot open $file";
	my $firstLine = <SEQUENCES>;
	$firstLine =~ s/[\n\r]//g;
	chomp($firstLine);
	my @header_fields = split(/\t/, $firstLine);
	while(<SEQUENCES>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		my @fields = split(/\t/, $line, -1);
		die Dumper($#fields, $#header_fields) unless($#fields == $#header_fields);
		my %l = (mesh @header_fields, @fields);
		die unless($l{'ContigID'});
		
		if($l{'Start'})
		{
			die unless($l{'Stop'});
			my $id = join('--', map {$l{$_}} qw/ContigID Start Stop/);
			$extract_partial_contig{$id}++;
		}
		else
		{
			$extract_complete_contig{$l{'ContigID'}}++;
		}
	}
	close(SEQUENCES);
};

my $process_sequences_file = sub {
	my $file = shift;

	open(SEQUENCES, '<', $file) or die "Cannot open $file";
	my $firstLine = <SEQUENCES>;
	$firstLine =~ s/[\n\r]//g;
	chomp($firstLine);
	my @header_fields = split(/\t/, $firstLine);
	while(<SEQUENCES>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		my @fields = split(/\t/, $line, -1);
		die Dumper($#fields, $#header_fields) unless($#fields == $#header_fields);
		my %l = (mesh @header_fields, @fields);
		if($l{'Chr'})
		{
			die Dumper("Undefied", $l{'Chr'}, $l{'Start_1based'}, $l{'Stop_1based'}, \@fields) unless($l{'Chr'} and $l{'Start_1based'} and $l{'Stop_1based'});
			my $id = join('--', map {$l{$_}} qw/Chr Start_1based Stop_1based/);
			$extract_partial_contig{$id}++;
		}
		else
		{
			die unless($l{'FASTAID'});
			$extract_complete_contig{$l{'FASTAID'}}++;
		}
	}
	close(SEQUENCES);
};

$process_sequences_file->($sequences_file);
$process_linear_ALTs_file->($linearALTs_file);

# make sure that we really have all contig IDs that we want to extract
my $samtools_idxstats_cmd = qq($samtools_bin idxstats $BAM_in);
my @samtools_idxstats_out = split(/\n/, `$samtools_idxstats_cmd` );
my %have_contigs = map {$_ => 1} map {my @f = split(/\t/, $_); die "Wrong field length $#f from line $_ "unless($#f == 3); $f[0]} @samtools_idxstats_out;

my @regions_complete = keys %extract_complete_contig;
my @regions_partial_contigIDs;
my @regions_partial = map {my @f = split(/\-\-/, $_); push(@regions_partial_contigIDs, $f[0]); die unless($#f == 2); $f[0].':'.$f[1].'-'.$f[2]} keys %extract_partial_contig;

foreach my $contigID (@regions_complete, @regions_partial_contigIDs)
{
	unless($have_contigs{$contigID})
	{
		die "Contig ID $contigID not in index of BAM $BAM_in";
	}
}

my $regions_string = join(' ', @regions_partial, @regions_complete); 

my $BAM_out_unsorted = $BAM_out . '.unsorted';
my $command = qq($samtools_bin view -bh $BAM_in $regions_string > $BAM_out_unsorted);
print "Now executing:\n",$command,"\n\n";
my $ret = system($command);
unless($ret == 0)
{
	die "Command $command failed";
}
die unless(-e $BAM_out_unsorted);

my $command_sort = qq($samtools_bin sort $BAM_out_unsorted $BAM_out);
print "Now executing:\n",$command_sort,"\n\n";
my $ret_sort = system($command_sort);
unless($ret_sort == 0)
{
	die "Command $command_sort failed";
}
die unless(-e $BAM_out_unsorted);
my $BAM_out_sorted = $BAM_out . '.bam';
die "Sorted file $BAM_out_sorted not there" unless(-e $BAM_out_sorted);

my $cmd_mv = qq(mv $BAM_out_sorted $BAM_out);
print "Now executing:\n",$cmd_mv,"\n\n";
my $ret_mv = system($cmd_mv);
unless($ret_mv == 0)
{
	die "Command $cmd_mv failed";
}

my $command_idx = qq($samtools_bin index $BAM_out);

print "Now executing:\n",$command_idx,"\n\n";
my $ret_idx = system($command_idx);
unless($ret_idx == 0)
{
	die "Command $command_idx failed";
}

unlink($BAM_out_unsorted) or die "Cannot delete $BAM_out_unsorted";

sub findFileFromAlternatives
{
	foreach my $a (@_)
	{
		if(-e $a)
		{
			return $a;
		}
	}
	die "Couldn't find from alternatives: ".join(", ", @_);
}
