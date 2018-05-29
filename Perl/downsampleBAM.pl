use strict;
use Data::Dumper;

my @applyToBAMs;
my @outputBAMs;

my @G1000_BAMs = glob('/data/projects/phillippy/projects/MHC/1000G_highCoverage_250bp/*.bam');
foreach my $BAM (@G1000_BAMs)
{
	die unless($BAM =~ /.+\/(.+?)\.(.+?)$/);
	my $sampleID = $1;
	push(@applyToBAMs, $BAM);
	my $outputBAM = '/data/projects/phillippy/projects/MHC/15X_BAMs/d1000G_' . $sampleID . '.bam';
	push(@outputBAMs, $outputBAM);
}


push(@applyToBAMs, '/data/projects/phillippy/projects/MHC/NA12878_CRAM/NA12878.cram');
push(@applyToBAMs, '/data/projects/phillippy/projects/rDNA/alignGenome2/temp/NA12891_Platinum_1000G/merged.bam');
push(@applyToBAMs, '/data/projects/phillippy/projects/rDNA/alignGenome2/temp/NA12892_Platinum_1000G/merged.bam');
push(@outputBAMs, '/data/projects/phillippy/projects/MHC/15X_BAMs/dPlatinum_NA12878.bam');
push(@outputBAMs, '/data/projects/phillippy/projects/MHC/15X_BAMs/dPlatinum_NA12891.bam');
push(@outputBAMs, '/data/projects/phillippy/projects/MHC/15X_BAMs/dPlatinum_NA12892.bam');

die unless($#applyToBAMs == $#outputBAMs);

$applyToBAMs[0] = 1;
for(my $BAMi = 0; $BAMi <= $#applyToBAMs; $BAMi++)
{
	my $BAM = $applyToBAMs[$BAMi];
	my $outputBAM = $outputBAMs[$BAMi];
	
	my $BAM = '/data/projects/phillippy/projects/MHC/1000G_highCoverage_250bp/B38_fromExtraction/1000G_NA20502.bam';
	my $outputBAM = $BAM . '.downsampled';
	if(-e $outputBAM . '.bai')
	{
		print "Detected BAM index file ${outputBAM}.bai, exit now - delete file if you want to repeat.\n";
		exit 0;
	}

	my $targetBases = 15 * 3.2e9;
	#$targetBases = 0.01 * 3.2e9;
	
	print "Downsample $BAM to ", sprintf("%.2f", $targetBases/1e9), "Gb.\n";

	# get read length
	my %readLengths;
	open(PIPE, "samtools view -F 0x100 -F 0x800 $BAM |") or die "Pipe command failed - have you loaded samtools?";
	while(<PIPE>)
	{
		chomp;
		next unless($_);
		my @line_fields = split(/\t/, $_);
		my $readLength = length($line_fields[9]);
		$readLengths{$readLength}++;
		last if($. > 10);
	}
	close(PIPE);

	my @readLengths_keys = keys %readLengths;
	if(scalar(@readLengths_keys) > 1)
	{
		die Dumper("More than one read length detected in $BAM", \%readLengths);
	}
	my $readLength = $readLengths_keys[0];

	print "Assume read length $readLength\n";
	 
	my %readIDs;
	open(PIPE, "samtools view -F 0x100 -F 0x800 $BAM |") or die "Pipe command failed - have you loaded samtools?";
	while(<PIPE>)
	{
		chomp;
		next unless($_);
		die "Weird samtools line: $_" unless($_ =~ /^(\S+)\t/);
		$readIDs{$1}++;
	}
	close(PIPE);

	my $reads_bases = 0;
	foreach my $readID (keys %readIDs)
	{
		my $n_reads = $readIDs{$readID};
		if($n_reads != 2)
		{
			warn "Read $readID has weird number of instances: $n_reads";
		}
		$reads_bases += ($n_reads * $readLength);
	}

	if($reads_bases == 0)
	{
		die "BAM doesn't contain reads.";
	}

	my $factor_for_downsampling = $targetBases / $reads_bases;

	print "BAM contains ", sprintf("%.2f", $reads_bases/1e9), "Gb, want ",  sprintf("%.2f", $targetBases/1e9), ",  i.e. now sample to factor ", sprintf("%.2f", $factor_for_downsampling), "\n";

	if($factor_for_downsampling > 1)
	{
		print "Don't have enough data in BAM for downsampling, exit now.\n";
		exit;
	}
	else
	{
		my %use_readIDs;
		foreach my $readID (keys %readIDs)
		{
			if(rand(1) <= $factor_for_downsampling)
			{
				$use_readIDs{$readID}++;
			}
		}
		
		print "Selected ", scalar(keys(%use_readIDs)), " / ", scalar(keys(%readIDs)), "\n";
		
		open(PIPEIN, "samtools view -h $BAM |") or die "Pipe in command failed - have you loaded samtools?";
		open(PIPEOUT, "| samtools view -bo $outputBAM -") or die "Pipe out command failed - have you loaded samtools?";
		while(<PIPEIN>)
		{
			chomp;
			next unless($_);
			if(substr($_, 0, 1) eq '@')
			{
				print PIPEOUT $_, "\n";
			}
			else
			{
				die "Weird samtools line: $_" unless($_ =~ /^(\S+)\t/);
				my $readID = $1;
				if($use_readIDs{$readID})
				{
					print PIPEOUT $_, "\n";				
				}
			}

		}
		close(PIPEIN);
		close(PIPEOUT);
		
		system("samtools index $outputBAM") and die "Indexing command for $outputBAM failed";
		
		print "\n\nDone. Produced $outputBAM\n";
	}
}
