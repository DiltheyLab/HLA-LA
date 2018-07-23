# original location: /data/projects/phillippy/projects/MHC/HLA-PRG-LA/src/forPaper/runKourami.pl
use strict;
use Data::Dumper;
use Getopt::Long;

my $output_folder = '/gpfs/project/dilthey/software/HLA-PRG-LA/fromXHLA';
my $G_groups_file = '/gpfs/project/dilthey/software/HLA-PRG-LA/hla_nom_g.txt';
my $target_reference = '/gpfs//project/dilthey/references/hg38.fa.forxHLA.fa';

my %G_full_alleles;
my %G_mapper_unambigious;
my %G_mapper_multiples;
my %G_to_underlyingAlleles;
my %alleles_to_fullGAmbiguity;
readG(); 

# die Dumper($alleles_to_fullGAmbiguity{'HLAA'});
# test command
# perl runxHLA.pl --BAM /software/HLA/master/tests/test.bam --sampleID testBAM
# perl runxHLA.pl --BAM /gpfs/project/dilthey/projects/MHC/CRAMs_GRCh38_primary/HapMap_SRR716646.cram --sampleID HapMap_SRR716646

my $BAM;
my $sampleID;
GetOptions (
	'BAM:s' => \$BAM,
	'sampleID:s' => \$sampleID,
);
die "Please provide --sampleID" unless($sampleID);
die "Please provide --BAM" unless($BAM);
die "--BAM not existing" unless(-e $BAM);

die "Please provide file ending in .cram or .bam" unless(($BAM =~ /\.cram$/) or ($BAM =~ /\.bam$/));
die "Please provide indexed file" unless((-e $BAM . '.crai') or (-e $BAM . '.bai'));

my $outputDir = $output_folder . '/' . $sampleID;
unless(-d $outputDir)
{
	mkdir($outputDir) or die "Cannot mkdir $outputDir";
}

my $final_output_file = $outputDir . '/xHLA.txt';

if(-e $final_output_file)
{
	unlink($final_output_file) or die "Cannot mkdir $final_output_file";
}	

# check that index is right reference
my $index_href = readIndexContigs($BAM);

my $reference_href = readFASTA($target_reference, 1);

die "Reference / BAM mismatch I" unless(scalar(keys %$index_href) == scalar(keys %$reference_href));
die Dumper("Reference / BAM mismatch II", [grep {not exists $reference_href->{$_}} keys %$index_href], [(keys %$reference_href)[0 .. 9]]) unless(scalar(grep {not exists $reference_href->{$_}} keys %$index_href) == 0);
 
my $BAM_for_xHLA = $BAM;
if($BAM_for_xHLA =~ /\.cram$/)
{ 
	$BAM_for_xHLA .= '.bam';
	my $extraction_command = qq(module load SamTools; samtools view -T $target_reference -bo $BAM_for_xHLA $BAM; samtools index $BAM_for_xHLA);
	system($extraction_command) and die "Could not execute $extraction_command";
}
else
{
	die unless($BAM_for_xHLA =~ /\.bam$/)
}

my $xHLA_raw_output_dir = $outputDir . '/xHLA/';
unless(-e $xHLA_raw_output_dir)
{
	mkdir($xHLA_raw_output_dir) or die "Cannot mkdir $xHLA_raw_output_dir";
}
system("rm ${xHLA_raw_output_dir}/*");


my $xHLA_output = $outputDir . '/output_xHLA.txt';
my $xHLA_cmd = qq(module load HLA; /software/HLA/master/bin/run.py --full --sample_id $sampleID --input_bam_path $BAM_for_xHLA --output_path $xHLA_raw_output_dir > $xHLA_output);
system($xHLA_cmd) and die "Cannot open $xHLA_cmd";

unless(-e $xHLA_raw_output_dir . '/_SUCCESS')
{
	die "File " . $xHLA_raw_output_dir . '/_SUCCESS' . " is missing, typing not susccessful?";
}	

unlink($BAM_for_xHLA);
unlink($BAM_for_xHLA . '.bai');

my $results_file = $xHLA_raw_output_dir . '/report-' . $sampleID . '-hla.json';
die "HLA results file $results_file not present" unless(-e $results_file);

my $results_string;
open(RESULTS, '<', $results_file) or die "Cannot open $results_file"; 
while(<RESULTS>)
{
	$results_string .= $_;
}
close(RESULTS);

die "Can't find alleles in file $results_file" unless($results_string =~ /"alleles": \[(.+?)\]/s);

my $alleles_string = $1;
$alleles_string =~ s/[\n\r\"\s]//g;

my @lines_for_output;
my @alleles = split(/,/, $alleles_string);
my %alleles_byLocus;
foreach my $allele (@alleles)
{ 
	die unless($allele =~ /^(\w+)\*(.+)$/);
	my $locus = $1;
	my $allele_noLocus = $2;
	push(@{$alleles_byLocus{$locus}}, $allele);
	die unless(scalar(@{$alleles_byLocus{$locus}}) <= 2);

	# die "Undefined allele $allele_noLocus at locus $locus" unless(exists $alleles_to_fullGAmbiguity{'HLA' . $locus}{$allele_noLocus});
	# push(@allAllelles_full, @{$alleles_to_fullGAmbiguity{'HLA' . $locus}{$allele_noLocus}});
	
	push(@lines_for_output, join("\t", $locus, scalar(@{$alleles_byLocus{$locus}}), $allele, 1, 1)); 

}

open(OUTPUT, '>', $final_output_file) or die "Cannot open $final_output_file";
print OUTPUT join("\t", qw/Locus Chromosome Allele Q1 Q2/), "\n";
foreach my $line (@lines_for_output)
{
	print OUTPUT $line, "\n";
}
close(OUTPUT);
print "\n\nDone. Produced $final_output_file \n\n";

sub readG
{
	open(GFILE, '<', $G_groups_file) or die "Cannot open $G_groups_file";
	while(<GFILE>)
	{
		my $line = $_;
		chomp($line);
		next if(substr($line, 0, 1) eq '#');
		die unless($line =~ /^(.+?)\*;(.+)$/);
		my $locus = $1;
		my $alleles = $2;
		
		my $G_group;
		if($alleles !~ /;$/)
		{
			die "Weird alleles for $locus: '$alleles'" unless($alleles =~ /^(.+);(.+?G)/);
			$alleles = $1;
			$G_group = $2; 
		}
		else
		{
			
			$alleles = substr($alleles, 0, length($alleles)-1); 
			die if($alleles =~ /\//);
			die if($alleles =~ /G/);
			$G_group = $alleles;
		}

		my @alleles = split(/\//, $alleles);
		
		for(@alleles)
		{
			$G_full_alleles{'HLA'.$locus}{$_}++;
		}
		
		$G_full_alleles{'HLA'.$locus}{$G_group}++;	
				
		$G_to_underlyingAlleles{'HLA'.$locus}{$G_group} = \@alleles;

		$alleles_to_fullGAmbiguity{'HLA'.$locus}{$G_group} = \@alleles;
		for(@alleles)
		{
			$alleles_to_fullGAmbiguity{'HLA'.$locus}{$_} = \@alleles;		
		}
		
		my $alleleMaxIdx = $#alleles;	
		for(my $i = 0; $i <= $alleleMaxIdx; $i++)
		{
			my $allele = $alleles[$i];
			
			die if(exists $G_mapper_unambigious{$locus . '*' . $allele});
			$G_mapper_unambigious{$locus . '*' . $allele} = $locus . '*' .$G_group;		
			
			# die Dumper(\%G_mapper_unambigious);
			
			#$G_mapper_unambigious{'HLA'.$locus}{$allele} = $G_group;
			
			# $G_mapper_multiples{'HLA'.$locus}{$allele} = [$G_group];
			
			# my @allele_parts = split(/:/, $allele);
			# die "Weird alllele II: $allele $alleles $line" unless($#allele_parts >= 1);
			
			# my $fourDigit = join(':', @allele_parts[0, 1]);
					
			# if(exists $G_mapper_unambigious{'HLA'.$locus}{$fourDigit})
			# {
				# my $existing_G_group = $G_mapper_unambigious{'HLA'.$locus}{$fourDigit};
				# if($existing_G_group ne $G_group)
				# {
					# # warn "$locus $allele $fourDigit existing: $existing_G_group - now want to set $G_group";
					# $G_mapper_unambigious{'HLA'.$locus}{$fourDigit} = undef;
				# }
				
				# my @existing_multiples = @{$G_mapper_multiples{'HLA'.$locus}{$fourDigit}};
				# my %_existing_multiples = map {$_ => 1} @existing_multiples;
				# if(not $_existing_multiples{$G_group})
				# {
					# push(@{$G_mapper_multiples{'HLA'.$locus}{$fourDigit}}, $G_group);
				# }
				
			# }
			# else
			# {
				# $G_mapper_unambigious{'HLA'.$locus}{$fourDigit} = $G_group;
				# $G_mapper_multiples{'HLA'.$locus}{$fourDigit} = [$G_group];
			# }
		}
	}
	close(GFILE);
}



sub readIndexContigs
{
	my $file = shift;
	my %reference_contigs;
	open(PIPE, "module load SamTools; samtools idxstats $file |") or die "Can't idxstat on $file";
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


