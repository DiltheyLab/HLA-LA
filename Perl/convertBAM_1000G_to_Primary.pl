use strict;
use Data::Dumper;
use Getopt::Long;   

# test command: perl convertBAM_1000G_to_Primary.pl --input /data/projects/phillippy/projects/MHC/15X_BAMs/d1000G_NA12878.cram --output /data/projects/phillippy/projects/MHC/CRAMs_GRCh38_primary/d1000G_NA12878.cram

my $primary_ref = '/data/projects/phillippy/projects/rDNA/alignGenome2/modified_references/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa';
my $G1000_ref_href = readFASTA('/data/projects/phillippy/projects/rDNA/alignGenome2/modified_references/GRCh38_1000G/hs38DH.fa');
my $primary_ref_href = readFASTA($primary_ref);
my $picard_sam2fastq_bin = '/data/projects/phillippy/software/picard-tools-2.7.1/picard.jar';
my $samtools_bin = '/data/projects/phillippy/software/samtools-1.4/samtools';


my $input;
my $output;
my $threads = 16;

GetOptions (
 'input:s' => \$input, 
 'output:s' => \$output, 
);

die unless($input);
die unless($output =~ /\.cram$/);
if(-e $output . '.crai')
{	
	die "Index file for $output existing, abort.";
}

my %primary_length_2_contigID;
foreach my $contigID (keys %$primary_ref_href)
{
	my $contigLength = length($primary_ref_href->{$contigID});
	die unless(defined $contigLength); 
	#die "Is there a duplicate contig length? $contigLength $contigID" if(exists $primary_length_2_contigID{$contigLength});
	push(@{$primary_length_2_contigID{$contigLength}}, $contigID);
}

warn Dumper(\%primary_length_2_contigID);

my %G1000_contigID_2_primary;
my %G1000_contigIDs;
my %_primary_used;
foreach my $contigID (keys %$G1000_ref_href)
{
	my $contigLength = length($G1000_ref_href->{$contigID});
	die unless(defined $contigLength);
	$G1000_contigIDs{$contigID}++;
	if(exists $primary_length_2_contigID{$contigLength})
	{
		my @contigIDs = @{$primary_length_2_contigID{$contigLength}};
		if(scalar(@contigIDs) != 1)
		{
			my %potential_sequences = map {$_ => $primary_ref_href->{$_}} @contigIDs;
			@contigIDs = grep {$G1000_ref_href->{$contigID} eq $primary_ref_href->{$_}} @contigIDs;
			die Dumper("Non-unambiguous contig assignment", \@contigIDs, $contigID, $contigLength, [\%potential_sequences, $G1000_ref_href->{$contigID}]) unless(scalar(@contigIDs) == 1);
		}
		my $primary_contig_ID = $contigIDs[0];
		if($_primary_used{$primary_contig_ID})
		{
			die "Mismatch I for $contigID / $primary_contig_ID" unless($G1000_ref_href->{$contigID} eq $primary_ref_href->{$primary_contig_ID});
			my $existing_1000G_for_primary = $_primary_used{$primary_contig_ID};
			die Dumper("Mismatch II for $existing_1000G_for_primary / $primary_contig_ID", $G1000_ref_href->{$existing_1000G_for_primary}, $primary_ref_href->{$primary_contig_ID}) unless($G1000_ref_href->{$existing_1000G_for_primary} eq $primary_ref_href->{$primary_contig_ID});
		}
		$_primary_used{$primary_contig_ID} = $contigID;
		$G1000_contigID_2_primary{$contigID} = $primary_contig_ID;
	}
}

print "\nContig mapping:\n";
foreach my $contigID (keys %G1000_contigID_2_primary)  
{
	print " - $contigID -> $G1000_contigID_2_primary{$contigID} \n";
}
print "\n";

my @contigIDs_nonPrimary;
{
	my $contigs_in_BAM_href = readIndexContigs($input);
	foreach my $contigID (keys %$contigs_in_BAM_href)
	{
		die unless(exists $G1000_contigIDs{$contigID});
		next if(exists $G1000_contigID_2_primary{$contigID});
		push(@contigIDs_nonPrimary, $contigID);
	}
}
print "Identitied non-primary contig IDs: ", join(', ', @contigIDs_nonPrimary[0 .. 9]), "\n";

my %reads_nonPrimaryAlignments;
{
	open(SAMTOOLS, "samtools view $input " . join(' ', @contigIDs_nonPrimary) . '|') or die "Cannot samtools view $input";
	while(<SAMTOOLS>)
	{
		chomp;
		next unless($_);
		die unless($_ =~ /^(\S+)\s/);
		my $readID = $1;
		$reads_nonPrimaryAlignments{$readID} = 1;
	}	
	close(SAMTOOLS);
}

print "Identified ", scalar(keys %reads_nonPrimaryAlignments), " reads with non-primary alignments.\n";

my $fn_alignOK = $output . '.alignOK.bam';
my $fn_alignNotOK = $output . '.alignNotOK.bam';
# open(ALIGNOK, "| $samtools_bin view -C -T $primary_ref -o $fn_alignOK -") or die "Pipe out command failed - have you loaded samtools?";
open(ALIGNOK, "| $samtools_bin view -bo $fn_alignOK -") or die "Pipe out command failed - have you loaded samtools?";
open(ALIGNNOTOK, "| $samtools_bin view -bo $fn_alignNotOK -") or die "Pipe out command failed - have you loaded samtools?";

open(ALIGNIN, "$samtools_bin view -H $input | ") or die "Cannot samtools view $input";
while(<ALIGNIN>)
{
	print ALIGNOK $_;
	print ALIGNNOTOK $_;
}
close(ALIGNIN);

open(ALIGNIN, "$samtools_bin view $input | ") or die "Cannot samtools view $input";

while(<ALIGNIN>)
{
	chomp;
	next unless($_);
	my $line = $_;
	die unless($line =~ /^(\S+)\t(\S+)\t(\S+)\t/);
	my $readID = $1;	
	my $flags = $2;	
	my $refID = $3;
	if($reads_nonPrimaryAlignments{$readID})
	{
		print ALIGNNOTOK $line, "\n";
	}
	else
	{	
		die unless($G1000_contigID_2_primary{$refID});
		my $new_line_beginning = $readID . "\t" . $flags . "\t" . $G1000_contigID_2_primary{$refID} . "\t";
		$line =~ s/^(\S+)\t(\S+)\t(\S+)\t/$new_line_beginning/;
		print ALIGNOK $line, "\n";
	}
}
close(ALIGNIN);
close(ALIGNOK);

print "\nGenerated\n\t$fn_alignOK and \n\t$fn_alignNotOK\n\n";

my $fn_alignNotOK_realigned = $fn_alignNotOK . '.realigned.cram';
{
	my $FASTQ_1 = $fn_alignNotOK . '_1.fastq'; 
	my $FASTQ_2 = $fn_alignNotOK . '_2.fastq';
	my $FASTQ_U = $fn_alignNotOK . '_U.fastq';
	my $picard_output = $fn_alignNotOK . '_picard.txt';

	my $FASTQ_extraction_command;
	if($picard_sam2fastq_bin =~ /SamToFastq\.jar$/)
	{
		die "Invalid code path - redirect STDERR and check that output gets interpreted correctly.";
		$FASTQ_extraction_command = qq(java -Xmx20g -XX:-UseGCOverheadLimit -jar $picard_sam2fastq_bin VALIDATION_STRINGENCY=LENIENT I=$fn_alignNotOK F=$FASTQ_1 F2=$FASTQ_2 FU=$FASTQ_U 2>&1);
	}
	elsif($picard_sam2fastq_bin =~ /picard.jar$/)
	{
		$FASTQ_extraction_command = qq(java -Xmx20g -jar $picard_sam2fastq_bin SamToFastq VALIDATION_STRINGENCY=LENIENT I=$fn_alignNotOK F=$FASTQ_1 F2=$FASTQ_2 FU=$FASTQ_U 2>&1 2> $picard_output);
	}
	else
	{
		die "I can't interpret the specified Picard command: $picard_sam2fastq_bin";
	}	
	
	$FASTQ_extraction_command .= ";grep -Fq 'picard.sam.SamToFastq done.' $picard_output || exit 1";
	system($FASTQ_extraction_command) and die "Extraction command $FASTQ_extraction_command failed";

	my $cmd_bwa = qq(bwa mem -t $threads $primary_ref $FASTQ_1 $FASTQ_2 | $samtools_bin view -C -T $primary_ref - | $samtools_bin sort --threads $threads -O CRAM > $fn_alignNotOK_realigned);

	print "Now realigning - with cmd\n\t$cmd_bwa\n\n";
	system($cmd_bwa) and die "Cannot execute $cmd_bwa";
	
	unlink($FASTQ_1);
	unlink($FASTQ_2);
	unlink($FASTQ_U);
	unlink($picard_output);
}

my $combined_output_bam = $output . '.bam';
my $combined_output_unsorted_cram = $output . '.unsorted.cram';

my $cmd_merge = qq($samtools_bin merge $combined_output_bam $fn_alignOK $fn_alignNotOK; samtools view -C -T $primary_ref $combined_output_bam > $combined_output_unsorted_cram);
system($cmd_merge) and die "Merge/view command $cmd_merge failed";

my $cmd_sort = qq($samtools_bin sort --threads $threads $combined_output_unsorted_cram > $output ; $samtools_bin index $output );
system($cmd_sort) and die "Sort/index command $cmd_sort failed";

print "\n\nGenerated output: $combined_output_bam \n";

unlink($fn_alignOK);
unlink($fn_alignNotOK);
unlink($combined_output_bam);
unlink($combined_output_unsorted_cram);

sub readIndexContigs
{
	my $file = shift;
	my %reference_contigs;
	open(PIPE, "samtools idxstats $file |") or die "Can't idxstat on $file";
	while(<PIPE>)
	{
		chomp;
		next unless($_);
		my @line_fields = split(/\t/, $_);
		next if ($line_fields[0] eq '*');
		$reference_contigs{$line_fields[0]}++;
	}
	close(PIPE);
	return \%reference_contigs;
}

sub readFASTA
{
	my $file = shift;	
	my $cut_sequence_ID_after_whitespace = shift;
	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			if($cut_sequence_ID_after_whitespace)
			{
				$line =~ s/\s+.+//;
			}
			$currentSequence = substr($line, 1);
			$R{$currentSequence} = '';
		}
		else
		{
			die "Weird input in $file" unless (defined $currentSequence);
			$R{$currentSequence} .= uc($line);
		}
	}	
	close(F);
		
	return \%R;
}


