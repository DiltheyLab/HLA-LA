#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::DB::HTS;
use Getopt::Long;   
$| = 1;

my $sampleID;
my $outputDirectory;
my $referenceGenome;
my $exonCoordinates;
my $FASTQ1;
my $FASTQ2;
my $bwa_bin;
my $samtools_bin;
my $maxThreads = 1;

GetOptions (
	'sampleID=s' => \$sampleID, 
	'outputDirectory=s' => \$outputDirectory, 
	'referenceGenome=s' => \$referenceGenome,
	'exonCoordinates=s' => \$exonCoordinates, 	
	'FASTQ1=s' => \$FASTQ1, 	
	'FASTQ2=s' => \$FASTQ2, 	
	'bwa_bin=s' => \$bwa_bin, 	
	'samtools_bin=s' => \$samtools_bin, 	
	'maxThreads:s' => \$maxThreads, 	
);


die "Reference genome $referenceGenome not there" unless(-e $referenceGenome);
die "FASTQ1 $FASTQ1 not there" unless(-e $FASTQ1);
die "FASTQ2 $FASTQ2 not there" unless(-e $FASTQ2);
die "File $exonCoordinates not there" unless(-e $exonCoordinates);

die "Please provide a bwa-indexed reference genome ($referenceGenome)" unless (-e $referenceGenome . '.sa');

my $output_fn = $outputDirectory . '/kMerCounts_' . $sampleID . '.txt';
open(KMERSOUT, '>', $output_fn) or die "Cannot open $output_fn"; # open file to make sure writable

# initial sanity check
my %know_reference_IDs = map {$_ => 1} getFASTAIdentifiers($referenceGenome);
open(EXONCOORDINATES, '<', $exonCoordinates) or die "Cannot open $exonCoordinates";
while(<EXONCOORDINATES>)
{
	chomp;
	next unless($_);
	my @fields = split(/\t/, $_);
	die unless(scalar(@fields) == 5);
	my $targetSequence = $fields[2];
	$targetSequence = '6' if($targetSequence eq 'pgf');
	die "Target sequence (from $exonCoordinates) not in reference genome $referenceGenome" unless($know_reference_IDs{$targetSequence});
}
close(EXONCOORDINATES);

# mapping

my $bam_unsorted = $outputDirectory . '/' . $sampleID . '_kMerExtraction.bam';
my $bam_sorted = $outputDirectory . '/' . $sampleID . '_kMerExtraction.sorted.bam';

print "Start mapping...\n";

my $cmd_map = qq($bwa_bin mem -t 32 $referenceGenome $FASTQ1 $FASTQ2 | $samtools_bin view -bo $bam_unsorted -);
system($cmd_map) and die "Command $cmd_map failed";

my $cmd_sort = qq($samtools_bin sort $bam_unsorted > $bam_sorted);
system($cmd_sort) and die "Command $cmd_sort failed";

my $cmd_index = qq($samtools_bin index $bam_sorted);
system($cmd_index) and die "Command $cmd_index failed";

# extraction

my %kmers;
my $sam = Bio::DB::HTS->new(-fasta => $referenceGenome, -bam => $bam_sorted, -expand_flags => 1);
open(EXONCOORDINATES, '<', $exonCoordinates) or die "Cannot open $exonCoordinates";
while(<EXONCOORDINATES>)
{
	chomp;
	next unless($_);
	my @fields = split(/\t/, $_);
	die unless(scalar(@fields) == 5);
	print "Extracting $fields[0] exon $fields[1]\n";
	my $targetSequence = $fields[2];
	my $targetStart = $fields[3];
	my $targetStop = $fields[4];
	
	my $alignment_iterator = $sam->features(-seq_id => $targetSequence, -start => $targetStart, -stop => $targetStop, -iterator => 1);	
	while(my $alignment = $alignment_iterator->next_seq)
	{
		next if($alignment->get_tag_values('UNMAPPED'));
		next if($alignment->get_tag_values('NOT_PRIMARY'));
		next if($alignment->get_tag_values('SUPPLEMENTARY'));	
		
		my ($ref,$matches,$query) = $alignment->padded_alignment;
		unless(length($ref) == length($query))
		{
			warn Dumper("Alignment length mismatch", length($ref), length($query), $alignment->query->name);
			next;
		}
		my $runningRefPos = $alignment->start - 1;
		my $extractedReadSequence = '';
		for(my $i = 0; $i < length($ref); $i++)
		{
			my $refC = substr($ref, $i, 1);
			my $queryC = substr($query, $i, 1);
			if($refC ne '-')
			{
				$runningRefPos++;
			}	
			if(($runningRefPos >= $targetStart) and ($runningRefPos <= $targetStop))
			{
				$extractedReadSequence .= $queryC;
			}
		}
		$extractedReadSequence =~ s/-//g;
		die Dumper("Weird extracted read sequence: $extractedReadSequence", $ref, $query, "Read ID ".$alignment->query->name, "Reference ID ".$alignment->seq_id) unless($extractedReadSequence =~ /^[ACGTNacgtn]*$/);
		
		my @kmers = kmers($extractedReadSequence, 31);
		foreach my $kmer (@kmers)
		{
			$kmers{$kmer}++;
		}
	}
}

foreach my $kmer (keys %kmers)
{
	print KMERSOUT join("\t", $kmer, $kmers{$kmer}), "\n";
}

close(KMERSOUT);
close(EXONCOORDINATES);

print "\nProduced file $output_fn\n";

#unlink($bam_unsorted);
#unlink($bam_sorted);
#unlink($bam_sorted.'.bai');


sub kmers
{
	my $sequence = shift;
	my $k = shift;
	if(length($sequence) < $k)
	{
		return ();
	}
	else
	{
		my @kMer_starts = (0 .. length($sequence) - $k);
		return map {my $kmer = substr($sequence, $_, $k); die unless(length($kmer) == $k); $kmer} @kMer_starts;
	}
}



sub getFASTAIdentifiers
{
	my $f = shift;
	my @forReturn;
	open(F, '<', $f) or die "Cannot open $f";
	while(<F>)
	{
		chomp;
		next unless($_);
		if(substr($_, 0, 1) eq '>')
		{
			my $id = substr($_, 1);
			$id =~ s/\s.+//;
			push(@forReturn, $id);
		}
	}
	close(F);
	return @forReturn;
}