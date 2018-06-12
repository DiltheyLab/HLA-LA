use strict;
use Data::Dumper;
use List::Util qw/sum/;
use Getopt::Long;   

# test command: perl downsampleBAM.pl --input /data/projects/phillippy/projects/MHC/NA12878_CRAM/NA12878.cram --output /data/projects/phillippy/projects/MHC/NA12878_CRAM/NA12878.cram.downsampled.cram --targetGigabases 1

my $reference = '/data/projects/phillippy/projects/rDNA/alignGenome2/modified_references/GRCh38_1000G/hs38DH.fa';

my $input;
my $output;
my $targetGigabases;

GetOptions (
 'input:s' => \$input, 
 'output:s' => \$output, 
 'targetGigabases:s' => \$targetGigabases, 
);

unless($targetGigabases)
{
	die "Please define --targetGigabases";
}

unless($output)
{
	die "Please define --output";
}


unless($input)
{
	die "Please define --input";
}

unless(-e $input)
{
	die "File --input ($input) not existing";
}

if(-e $output . '.crai')
{
	print "Detected CRAM index file ${output}.crai, exit now - delete file if you want to repeat.\n";
	exit 0;
}

my $BAM_contig_lengths = readBAMContigLengths($input);
my $reference_contig_lengths_href =  readReferenceContigLengths($reference);

foreach my $contigID (keys %$BAM_contig_lengths)
{
	die "Contig $contigID from input not in reference" unless(exists $reference_contig_lengths_href->{$contigID});
	die "Length mismatch for $contigID - BAM $BAM_contig_lengths->{$contigID} - reference $reference_contig_lengths_href->{$contigID}"  unless($BAM_contig_lengths->{$contigID} == $reference_contig_lengths_href->{$contigID});
}

my $targetBases = $targetGigabases * 1e9;

	
print "Downsample $input to ", sprintf("%.2f", $targetBases/1e9), "Gb.\n";

# get read length
my %readLengths;
open(PIPE, "samtools view -F 0x100 -F 0x800 $input |") or die "Pipe command failed - have you loaded samtools?";
while(<PIPE>)
{
	chomp;
	next unless($_);
	my @line_fields = split(/\t/, $_);
	my $readLength = length($line_fields[9]);
	$readLengths{$readLength}++;
	last if($. > 100);
}
close(PIPE);

my @readLengths_keys_sorted = sort {$b <=> $a}  keys %readLengths;
my $consideredReadsForLength = sum values %readLengths;
my $readLength = $readLengths_keys_sorted[0];

print "Assume read length $readLength, support from $readLengths{$readLength} / $consideredReadsForLength reads.\n";
	 
my %readIDs;
open(PIPE, "samtools view -F 0x100 -F 0x800 $input |") or die "Pipe command failed - have you loaded samtools?";
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

print "BAM/CRAM contains ", sprintf("%.2f", $reads_bases/1e9), "Gb, want ",  sprintf("%.2f", $targetBases/1e9), ",  i.e. now sample to factor ", sprintf("%.2f", $factor_for_downsampling), "\n";

if($factor_for_downsampling > 1)
{
	print "Don't have enough data in input for downsampling, exit now.\n";
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
	
	open(PIPEIN, "samtools view -h $input |") or die "Pipe in command failed - have you loaded samtools?";
	open(PIPEOUT, "| samtools view -C -T $reference -o $output -") or die "Pipe out command failed - have you loaded samtools?";
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
	
	system("samtools index $output") and die "Indexing command for $output failed";
	
	print "\n\nDone. Produced $output\n";
}

sub readReferenceContigLengths
{
	my $input = shift;
	my %forReturn;
	my $runningContigID;
	open(INPUT, '<', $input) or die "Cannot open $input";
	while(<INPUT>)
	{
		chomp;
		next unless($_);
		if(substr($_, 0, 1) eq '>')
		{
			substr($_, 0, 1) = '';
			$_ =~ s/\s.+//;
			$runningContigID = $_;
		}
		else
		{
			die unless(defined $runningContigID);
			$forReturn{$runningContigID} += length($_);
		}
	}
	close(INPUT);
	return \%forReturn;
}	

sub readBAMContigLengths
{
	my $input = shift;
	my %forReturn;
	die unless(-e $input);
	open(PIPE, "samtools idxstats $input |") or die "Can't open idxstats pipe";
	while(<PIPE>)
	{
		chomp;
		next unless($_);
		my @line_fields = split(/\s+/, $_);
		next if($line_fields[0] eq '*');
		$forReturn{$line_fields[0]} = $line_fields[1];
	}
	close(PIPE);
	return \%forReturn;
}