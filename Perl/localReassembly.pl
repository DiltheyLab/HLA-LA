use warnings;
use strict;
use FindBin;
use File::Spec;
use Getopt::Long;
use Data::Dumper; 
use Sys::Hostname;
use Cwd qw/getcwd abs_path/;
use List::MoreUtils qw/mesh/;
use List::Util qw/any all sum min max uniq shuffle/;
srand(1);

my $graph;
my $inputSAM;
my $bin_sam2alignment = $FindBin::Bin . '/../../bin/sam2alignment';
my $reference;
my $samtools_bin;
my $outputPrefix;
GetOptions (
	'graph:s' => \$graph,
	'inputSAM:s' => \$inputSAM, 
	'reference:s' => \$reference,
	'samtools_bin:s' => \$samtools_bin,
	'outputPrefix:s' => \$outputPrefix,
);

die "bin_sam2alignment not existing: $bin_sam2alignment" unless(-e $bin_sam2alignment);
die "Please specify --inputSAM" unless($inputSAM);
die "Please specify --reference" unless($reference);
die "Please specify --samtools_bin" unless($samtools_bin);

my $fn_out_genes = $outputPrefix . '.genes';
my $fn_out_MSA = $outputPrefix . '.MSA';
my $fn_out_activeAlleles = $outputPrefix . '.activeAlleles';
my $fn_out_readCoordinates = $outputPrefix . '.readCoordinates';
my $fn_out_readAlleles = $outputPrefix . '.readAlleles';
my $fn_out_readCoordinates_afterTrimming = $outputPrefix . '.readCoordinatesAfterTrimming';
my $fn_out_readCoordinates_afterTrimming_details = $outputPrefix . '.readCoordinatesAfterTrimming.details.txt';

my $fn_out_readCoverageByGene = $outputPrefix . '.readCoverageByGene.txt';
my $fn_out_activeAllelePolymorphic_byRead = $outputPrefix . '.activeAllelePolymorphic_byRead.txt';
my $fn_out_readCoverageByGeneDetails = $outputPrefix . '.readCoverageByGeneDetails.txt';
my $fn_out_referenceSelection = $outputPrefix . '.referenceSelectionInfo.txt';
my $fn_out_outerIterationsByGene = $outputPrefix . '.outerIterationsByGene.txt';

my $references_remapping_href = readFASTA($reference);

my %refSeqID_2_hap;
my %locus_has_chromosome_resolution;
my $fn_ref_seq2Hap = $reference . '.seq2Hap';
open(SEQ2HAP, '<', $fn_ref_seq2Hap) or die "Cannot open $fn_ref_seq2Hap";
while(<SEQ2HAP>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/[\n\r]//g;
	next unless($line);
	my @line_fields = split(/\t/, $line, -1);
	die unless(scalar(@line_fields) == 3);
	$line_fields[2] =~ s/ .+//;
	$refSeqID_2_hap{$line_fields[0]}{$line_fields[2]} = $line_fields[1];
	if($line_fields[1] ne '?')
	{
		$locus_has_chromosome_resolution{$line_fields[0]} = 1;
	}
}
close(SEQ2HAP);

foreach my $locus (keys %locus_has_chromosome_resolution)
{
	die if(any {$_ eq '?'} values %{$refSeqID_2_hap{$locus}});
}

my $alignments_file = $inputSAM . '.alignments';
my $sam2alignment_cmd = qq($bin_sam2alignment $inputSAM $reference > $alignments_file);
system($sam2alignment_cmd) and die "Could not execute: $sam2alignment_cmd";
print "Generated alignment file: $alignments_file\n";

my $full_graph_dir = $FindBin::RealBin . '/../../graphs/' . $graph;
die "Directory $full_graph_dir not existing" unless(-e $full_graph_dir);

my $raw_references_href = read_raw_references();
my $alignments_references_href = read_aligned_references();
my $alignments_translation_targets_href = compute_translation_targets($alignments_references_href);
my $raw_references_2AlignmentCoordinates_href = compute_alignment_coordinates($raw_references_href, $alignments_references_href);
my %translation_targets = map {$_ => 1} values %$alignments_translation_targets_href;
my %alleles_by_gene;
my %ref_allele_by_gene;
foreach my $gene_ref_id (keys %translation_targets)
{
	next if($gene_ref_id =~ /pgf/);
	my $gene = geneName($gene_ref_id);
	my @alleles = grep {$_ =~ /^$gene\*/} keys %$alignments_translation_targets_href;
	$alleles_by_gene{$gene} = \@alleles;
	$ref_allele_by_gene{$gene} = $gene_ref_id;
}
my %alleles_by_gene_forRemapping;
foreach my $allele_sequence_ID (keys %$references_remapping_href)
{
	next if($allele_sequence_ID =~ /pgf/);
	my $gene = geneName($allele_sequence_ID);
	$alleles_by_gene_forRemapping{$gene}{$allele_sequence_ID} = 1;
}

my $raw_references_doubleCheck_href = readFASTA($reference, 1);
foreach my $refID (keys %$raw_references_doubleCheck_href)
{
	die "Missing sequence: $refID" unless(exists $raw_references_href->{$refID});
	die "Sequence mismatch for sequence $refID" unless($raw_references_href->{$refID} eq $raw_references_doubleCheck_href->{$refID});
}

my %graphs;
foreach my $gene_referenceSequenceID (keys %translation_targets)
{
	next if($gene_referenceSequenceID =~ /pgf/);

	die unless(exists $alignments_references_href->{$gene_referenceSequenceID});	
	for(my $levelI = 0; $levelI <= length($alignments_references_href->{$gene_referenceSequenceID}); $levelI++)
	{
		my $nodeID = 'node' . $levelI;
		$graphs{$gene_referenceSequenceID}{nodes}{$nodeID}{level} = $levelI;
		$graphs{$gene_referenceSequenceID}{nodes}{$nodeID}{links} = {};
		$graphs{$gene_referenceSequenceID}{nodes}{$nodeID}{coverage} = 0;
	}
}

my %readID_2_haplotypes;
my $n_secondary = 0;
my $n_supplementary = 0;
my @runningAlignmentInfo;
my $runningAlignmentInfo_readID;
open(ALIGNMENTS, '<', $alignments_file) or die "Cannot open $alignments_file";
while(<ALIGNMENTS>)
{
	my $line = $_;
	chomp($line);
	next unless($line);

	my $headerLine = $line;
	my $alignment_reference = <ALIGNMENTS>; chomp($alignment_reference);
	my $alignment_read = <ALIGNMENTS>; chomp($alignment_read);
	my $qualities_read = <ALIGNMENTS>; chomp($qualities_read);
	
	die "File corruption" unless(length($alignment_reference) == length($alignment_read));
	
	if($headerLine =~ /secondary=1/)
	{
		$n_secondary++;
		next;
	}
	if($headerLine =~ /supplementary=1/)
	{
		$n_supplementary++;
		next;
	}
	
	my @headerLine_fields = split(/ /, $headerLine);
	my $readID = $headerLine_fields[0];
	(my $readID_noPair = $readID) =~ s/\/[12]//;
	
	# print "Now processing " . scalar(@runningAlignmentInfo) . " alignment...\n";
	
	if((defined $runningAlignmentInfo_readID) and ($readID_noPair ne $runningAlignmentInfo_readID))
	{
		processAlignments(\@runningAlignmentInfo);
		@runningAlignmentInfo = ();
	}
	$runningAlignmentInfo_readID = $readID_noPair;
	
	push(@runningAlignmentInfo, [\@headerLine_fields, $alignment_reference, $alignment_read, $qualities_read]);
	die Dumper(\@runningAlignmentInfo) if(scalar(@runningAlignmentInfo) > 2);

}
close(ALIGNMENTS);

my %genotypes_by_MSAcolumn;
my %genotypes_by_reads;
my %readID_2_start_stop_MSAcolumns;
my %readID_2_source_haplotype;
processAlignments(\@runningAlignmentInfo);

die Dumper("Missing graph for HLA-A", [keys %graphs]) unless(exists $graphs{'A*ref'});
#printGraph($graphs{'A*ref'});

print "\n";
print "n_secondary: $n_secondary\n";
print "n_supplementary: $n_supplementary\n";

my $skipped_read_because_multiple_genes = 0;
my %gene_2_reads;
foreach my $readID (keys %readID_2_start_stop_MSAcolumns)
{
	my @genes = keys %{$readID_2_start_stop_MSAcolumns{$readID}};
	if(scalar(@genes) > 1)
	{
		$skipped_read_because_multiple_genes++;
	}
	else
	{
		my $geneID = $genes[0];
		$gene_2_reads{$geneID}{$readID}++;
	}
} 

my %haveGeneInOutput;
my %selectedMSAReferenceAlleles_by_gene;
my %polymorphicPositionsByGene;
my %active_alleles_by_MSAposition_by_gene;
open(OUT_ACTIVEALLELES, '>', $fn_out_activeAlleles) or die "Cannot open $fn_out_activeAlleles";
open(OUT_REFERENCESELECTION, '>', $fn_out_referenceSelection) or die "Cannot open $fn_out_referenceSelection";
foreach my $gene_ref_id (sort keys %translation_targets)
{
	next if($gene_ref_id =~ /pgf/);
	my $gene = geneName($gene_ref_id);
	die unless(exists $alleles_by_gene{$gene});
	unless(exists $alleles_by_gene_forRemapping{$gene})
	{	
		warn "No remapping alleles for gene $gene [$gene_ref_id]" ;
		next;
	}
	
	unless(exists $genotypes_by_MSAcolumn{$gene})
	{
		warn "No reads data for gene $gene [$gene_ref_id]";
		next;
	}
	die unless(exists $alignments_references_href->{$gene_ref_id});
	
	$haveGeneInOutput{$gene}++;
	
	my $MSA_length = length($alignments_references_href->{$gene_ref_id});

	my %readCoverage_by_MSAcolumn = map {$_ => 0} (0 .. ($MSA_length - 1));
	my %alleles_present_in_reads;
	foreach my $position (sort {$a <=> $b} keys %{$genotypes_by_MSAcolumn{$gene}})
	{		
		die unless($position >= 0);
		die unless($position <= $MSA_length); 
	
		my @alleles_in_reads = keys %{$genotypes_by_MSAcolumn{$gene}{$position}};
		my $totalCoverage = sum values %{$genotypes_by_MSAcolumn{$gene}{$position}};
		$totalCoverage = 0 unless(defined $totalCoverage);
		
		$readCoverage_by_MSAcolumn{$position} = $totalCoverage,
		
		my $minimumCoverage = int($totalCoverage / 10);
		$minimumCoverage = 1 if($minimumCoverage < 1);
		
		# todo: consider increasing this to 2
		my @alleles_in_reads_sufficientCoverage = sort grep {$genotypes_by_MSAcolumn{$gene}{$position}{$_} >= $minimumCoverage} @alleles_in_reads;
		$alleles_present_in_reads{$position} = \@alleles_in_reads_sufficientCoverage;
	}
	
	my %remapping_alleles_present_at_position;
	my %remapping_alleles_by_MSAid;
	my %recoveredAlleles_by_MSAallele;
	my %nonRecoveredAlleles_by_MSAallele;
	my @alleles_used_for_remapping = keys %{$alleles_by_gene_forRemapping{$gene}};
	foreach my $MSA_allele_id (@alleles_used_for_remapping)
	{
		$nonRecoveredAlleles_by_MSAallele{$MSA_allele_id} = 0;
		$recoveredAlleles_by_MSAallele{$MSA_allele_id} = 0;
		
		my $allele_MSA_sequence = $alignments_references_href->{$MSA_allele_id};
		die unless(defined $allele_MSA_sequence);
		die unless(length($allele_MSA_sequence) == $MSA_length);
		for(my $i = 0; $i < $MSA_length; $i++)
		{
			my $allele_i_allele = substr($allele_MSA_sequence, $i, 1);
			$remapping_alleles_present_at_position{$i}{$allele_i_allele}++;
					
			my $allele_coverage = (exists $genotypes_by_MSAcolumn{$gene}{$i}{$allele_i_allele}) ? $genotypes_by_MSAcolumn{$gene}{$i}{$allele_i_allele} : 0;
			if($allele_coverage > 0)
			{
				$recoveredAlleles_by_MSAallele{$MSA_allele_id}++;
			}
			
			if($readCoverage_by_MSAcolumn{$i} >= 10)
			{
				if($allele_coverage == 0)
				{
					$nonRecoveredAlleles_by_MSAallele{$MSA_allele_id}++;
				}
			}
		}		
	}
	
	print OUT_REFERENCESELECTION "Gene $gene, " . scalar(@alleles_used_for_remapping) .  " MSA alleles used for remapping.\n";		
	my $alleles_per_haplotype = 2;
	my @alleles_used_for_remapping_filtered;
	if($locus_has_chromosome_resolution{$gene})
	{
		# todo: think about the case that allele 1 == allele 2
		my @MSA_seqIDs_h1 = grep {$refSeqID_2_hap{$gene}{$_} == 1} @alleles_used_for_remapping;
		my @MSA_seqIDs_h2 = grep {$refSeqID_2_hap{$gene}{$_} == 2} @alleles_used_for_remapping;
		 
		my @alleles_used_for_remapping_filtered_h1 = sort {$nonRecoveredAlleles_by_MSAallele{$a} <=> $nonRecoveredAlleles_by_MSAallele{$b}} @MSA_seqIDs_h1;
		if(scalar(@alleles_used_for_remapping_filtered_h1) > $alleles_per_haplotype)
		{
			@alleles_used_for_remapping_filtered_h1 = @alleles_used_for_remapping_filtered_h1[0 .. ($alleles_per_haplotype-1)];
		}

		my @alleles_used_for_remapping_filtered_h2 = sort {$nonRecoveredAlleles_by_MSAallele{$a} <=> $nonRecoveredAlleles_by_MSAallele{$b}} @MSA_seqIDs_h2;
		if(scalar(@alleles_used_for_remapping_filtered_h2) > $alleles_per_haplotype)
		{
			@alleles_used_for_remapping_filtered_h2 = @alleles_used_for_remapping_filtered_h2[0 .. ($alleles_per_haplotype-1)];
		}
		
		@alleles_used_for_remapping_filtered = (@alleles_used_for_remapping_filtered_h1, @alleles_used_for_remapping_filtered_h2);
		my %tmp_uniq = map {$_ => 1} @alleles_used_for_remapping_filtered;
		@alleles_used_for_remapping_filtered = sort keys %tmp_uniq;
	}
	else
	{
		@alleles_used_for_remapping_filtered = sort {$nonRecoveredAlleles_by_MSAallele{$a} <=> $nonRecoveredAlleles_by_MSAallele{$b}} @alleles_used_for_remapping;
		if(scalar(@alleles_used_for_remapping_filtered) > (2 * $alleles_per_haplotype))
		{
			@alleles_used_for_remapping_filtered = @alleles_used_for_remapping_filtered[0 .. ((2 * $alleles_per_haplotype)-1)];
		}
	}
	$selectedMSAReferenceAlleles_by_gene{$gene} = \@alleles_used_for_remapping_filtered;
	
	my %in_alleles_used_for_remapping_filtered = map {$_ => 1} @alleles_used_for_remapping_filtered;
	print OUT_REFERENCESELECTION "\tSelected ", scalar(@alleles_used_for_remapping_filtered), " MSA reference alleles (MSA length $MSA_length)\n";
	my @MSA_alleles_used_for_remapping_byNonRecovery = sort {$nonRecoveredAlleles_by_MSAallele{$a} <=> $nonRecoveredAlleles_by_MSAallele{$b}} @alleles_used_for_remapping;	
	foreach my $MSA_seq_id (@MSA_alleles_used_for_remapping_byNonRecovery)
	{
		my $allele_selected = (exists $in_alleles_used_for_remapping_filtered{$MSA_seq_id}) ? 1 : 0;
		print OUT_REFERENCESELECTION "\t\t", join("\t", $MSA_seq_id, $allele_selected, $nonRecoveredAlleles_by_MSAallele{$MSA_seq_id}, $recoveredAlleles_by_MSAallele{$MSA_seq_id}), "\n";
	}	

	# my %position_is_polymorphic_in_MSA = map {$_ => 1} grep {scalar(keys %{$remapping_alleles_present_at_position{$_}}) > 1} keys %remapping_alleles_present_at_position;
	# my %position_nonPolymorphic = map {$_ => 1} grep {not exists $position_is_polymorphic_in_MSA{$_}} keys %remapping_alleles_present_at_position;
	# my %position_nonPolymorphic_allele = map {my $pos = $_; my @alleles = keys %{$remapping_alleles_present_at_position{$_}}; die unless(scalar(@alleles) == 1); $pos => $alleles[0]} keys %position_nonPolymorphic;
	# my %positions_considered_alleles = map {$_ => $remapping_alleles_present_at_position{$_}} keys %position_is_polymorphic_in_MSA;
	# foreach my $position (sort {$a <=> $b} keys %positions_considered_alleles)
	# {
		# foreach my $allele (keys %{$positions_considered_alleles{$position}})
		# {
			# $positions_considered_alleles{$position}{$allele} = {fromMSA => 1};
		# }
	# }
	
	my %active_alleles_by_MSAposition;
	foreach my $position (sort {$a <=> $b} keys %{$genotypes_by_MSAcolumn{$gene}})
	{		
		die unless($position >= 0);
		die unless($position <= $MSA_length); 
	
		my @alleles_in_reads = keys %{$genotypes_by_MSAcolumn{$gene}{$position}};
		my $totalCoverage = sum values %{$genotypes_by_MSAcolumn{$gene}{$position}};
		$totalCoverage = 0 unless(defined $totalCoverage);
		my $minimumCoverage = int($totalCoverage / 10);
		$minimumCoverage = 1 if($minimumCoverage < 1);
		
		# todo: consider increasing this to 2
		my @alleles_in_reads_sufficientCoverage = sort grep {$genotypes_by_MSAcolumn{$gene}{$position}{$_} >= $minimumCoverage} @alleles_in_reads;
		my @alleles_in_MSA_selectedSequences = uniq(map {my $alleleID = $_; substr($alignments_references_href->{$alleleID}, $position, 1)} @alleles_used_for_remapping_filtered);
		
		die unless(scalar(@alleles_in_MSA_selectedSequences));
		
		# standard behaviour is to restrict set of alleles to the set of alleles found in the reads
		my %active_alleles_thisPosition = map {$_ => 1} @alleles_in_reads_sufficientCoverage;
		
		if($totalCoverage < 10)
		{
			foreach my $MSAallele (@alleles_in_MSA_selectedSequences)
			{
				$active_alleles_thisPosition{$MSAallele}++;
			}				
		}
		
		$active_alleles_by_MSAposition{$position} = \%active_alleles_thisPosition;
		
		if(scalar(keys %active_alleles_thisPosition) > 1)
		{
			$polymorphicPositionsByGene{$gene}{$position} = 1;
		}
		
		print OUT_ACTIVEALLELES join("\t",
			$gene, 
			$position,
			map {
				my $allele = $_;
				my $in_MSA = (exists $remapping_alleles_present_at_position{$position}{$allele}) ? 1 : 0;
				my $reads_coverage = (exists $genotypes_by_MSAcolumn{$gene}{$position}{$allele}) ? $genotypes_by_MSAcolumn{$gene}{$position}{$allele} : 0;
				join(';', $allele, $in_MSA, $reads_coverage)
			} sort keys %active_alleles_thisPosition
		), "\n";
	}
	
	$active_alleles_by_MSAposition_by_gene{$gene} = \%active_alleles_by_MSAposition;
}
close(OUT_ACTIVEALLELES);
close(OUT_REFERENCESELECTION);

# read trimming
open(OUT_READCOORDINATES, '>', $fn_out_readCoordinates) or die "Cannot open $fn_out_readCoordinates";
open(OUT_READCOORDINATES_AFTER_TRIMMING, '>', $fn_out_readCoordinates_afterTrimming) or die "Cannot open $fn_out_readCoordinates_afterTrimming";
open(OUT_READCOORDINATES_AFTER_TRIMMING_DETAILS, '>', $fn_out_readCoordinates_afterTrimming_details) or die "Cannot open $fn_out_readCoordinates_afterTrimming_details";
my %readID_2_start_stop_MSAcolumns_afterTrimming;
my %gene_2_reads_after_trimming;
my %readID_relevantPositions_afterTrimming;
foreach my $geneID (sort keys %gene_2_reads)
{
	my $MSA_length = length($alignments_references_href->{$ref_allele_by_gene{$geneID}});
	die Dumper("No MSA length for gene $geneID / $ref_allele_by_gene{$geneID}}", [keys %$alignments_references_href]) unless(defined $MSA_length);
	my @readIDs_sorted = sort {$readID_2_start_stop_MSAcolumns{$a}{$geneID}[0] <=> $readID_2_start_stop_MSAcolumns{$b}{$geneID}[0]} keys %{$gene_2_reads{$geneID}};
	foreach my $readID (@readIDs_sorted)
	{
		my $firstColumn = $readID_2_start_stop_MSAcolumns{$readID}{$geneID}[0];
		my $lastColumn = $readID_2_start_stop_MSAcolumns{$readID}{$geneID}[1];
		die unless((defined $firstColumn) and (defined $lastColumn));
		
		print OUT_READCOORDINATES join("\t", $geneID, $readID, $readID_2_start_stop_MSAcolumns{$readID}{$geneID}[0], $readID_2_start_stop_MSAcolumns{$readID}{$geneID}[1]), "\n";
		
		my @read_columns_polymorphic = grep {$polymorphicPositionsByGene{$geneID}{$_}} ($firstColumn .. $lastColumn); 
		
		if(scalar(@read_columns_polymorphic) == 0)
		{
			print OUT_READCOORDINATES_AFTER_TRIMMING_DETAILS join("\t", $geneID, $readID, scalar(@read_columns_polymorphic)), "\n";
			
		}
		else
		{
			my $minPos_polymorphic = min(@read_columns_polymorphic);
			my $maxPos_polymorphic = max(@read_columns_polymorphic);
			
			my @thisRead_relevantPositions_afterTrimming = grep {exists $genotypes_by_reads{$readID}{$geneID}{$_}} grep {($_ >= $minPos_polymorphic) and ($_ <= $maxPos_polymorphic)} ($firstColumn .. $lastColumn);

			if(@thisRead_relevantPositions_afterTrimming)
			{
				$readID_relevantPositions_afterTrimming{$readID}{$geneID} = \@thisRead_relevantPositions_afterTrimming;
				$readID_2_start_stop_MSAcolumns_afterTrimming{$readID}{$geneID} = [$minPos_polymorphic, $maxPos_polymorphic];
						
				if($readID eq 'HLADMB_h1_DMB*01:03:01:02_2806_3340_1:0:0_0:0:0_59')
				{
					# die Dumper($readID_2_start_stop_MSAcolumns_afterTrimming{$readID}{$geneID}, $readID_relevantPositions_afterTrimming{$readID}{$geneID}, $genotypes_by_reads{$readID}{$geneID});
				}
				
				print OUT_READCOORDINATES_AFTER_TRIMMING_DETAILS join("\t", $geneID, $readID, scalar(@read_columns_polymorphic), join(';', @read_columns_polymorphic), $firstColumn, $lastColumn, $minPos_polymorphic, $maxPos_polymorphic), "\n";
				print OUT_READCOORDINATES_AFTER_TRIMMING join("\t", $geneID, $readID, $minPos_polymorphic, $maxPos_polymorphic), "\n";
				
				$gene_2_reads_after_trimming{$geneID}{$readID}++;					
			}		
		}
		
	}
}
close(OUT_READCOORDINATES);
close(OUT_READCOORDINATES_AFTER_TRIMMING);
close(OUT_READCOORDINATES_AFTER_TRIMMING_DETAILS);


open(REEADCOVERAGEBYGENE_DETAILS, '>', $fn_out_readCoverageByGeneDetails) or die "Cannot open $fn_out_readCoverageByGeneDetails";
print REEADCOVERAGEBYGENE_DETAILS join("\t", "Iteration", "GeneID", "Position", "TotalCoverage", "CoverageH1", "CoverageH2", "EffectiveReads", "Alleles"), "\n";
my $used_reads_total = 0;
my %use_reads_by_gene;
my %coverages_by_gene;
my $target_coverage_per_haplotype = 4;
foreach my $geneID (sort keys %gene_2_reads_after_trimming)
{
	my $MSA_length = length($alignments_references_href->{$ref_allele_by_gene{$geneID}});
	die Dumper("No MSA length for gene $geneID / $ref_allele_by_gene{$geneID}}", [keys %$alignments_references_href]) unless(defined $MSA_length);
	
	my @readIDs_sorted = sort {$readID_2_start_stop_MSAcolumns_afterTrimming{$a}{$geneID}[0] <=> $readID_2_start_stop_MSAcolumns_afterTrimming{$b}{$geneID}[0]} keys %{$gene_2_reads_after_trimming{$geneID}};
	#@readIDs_sorted = shuffle @readIDs_sorted;
	my @remaining_read_IDs = @readIDs_sorted;
	my $outerIteration = 0;
	while(@remaining_read_IDs)
	{
		$outerIteration++;
		
		my @MSA_coverages;
		for(my $i = 0; $i < $MSA_length; $i++)
		{
			push(@MSA_coverages, [0, 0, 0]);
		}
		die unless(scalar(@MSA_coverages) == $MSA_length);
	
		my %use_reads_thisGene;
		foreach my $iteration (0, 1)
		{
			next unless(($iteration == 1) or $locus_has_chromosome_resolution{$geneID}); 
			# my $last_position;		
			foreach my $readID (@remaining_read_IDs)
			{
				next if($use_reads_thisGene{$readID});
				
				my $firstColumn = $readID_2_start_stop_MSAcolumns_afterTrimming{$readID}{$geneID}[0];
				my $lastColumn = $readID_2_start_stop_MSAcolumns_afterTrimming{$readID}{$geneID}[1];
				die unless((defined $firstColumn) and (defined $lastColumn));
				
				# if(defined $last_position)
				# {
					# die unless($last_position <= $firstColumn);
				# }
				# $last_position = $firstColumn;
				
				my @relevantColumns = ($firstColumn .. $lastColumn);
				my $lookup_index;
				if($iteration == 1)
				{
					$lookup_index = 0;
				}
				else
				{
					die unless($locus_has_chromosome_resolution{$geneID});
					$lookup_index = $readID_2_source_haplotype{$readID}{$geneID};
					die "Undefined source haplotype for read $readID ($geneID)" unless(defined $lookup_index);
				}
				next if(($iteration == 0) and ($lookup_index == 0));
				my $maxCoverage = ($iteration == 0) ? $target_coverage_per_haplotype : (2 * $target_coverage_per_haplotype);
				# print "$readID prior to increasing: ";
				# foreach my $column (@relevantColumns)
				# {
					# print $column, ": ", $MSA_coverages[$column][0], " ";
				# }
				# print "\n";
				if (all {$MSA_coverages[$_][$lookup_index] < $maxCoverage} @relevantColumns)
				{
					
					# die Dumper("There is a coverage issue", [map {[$_ , $MSA_coverages[$_][0]]} grep {$MSA_coverages[$_][0] > (2 * $target_coverage_per_haplotype)} @relevantColumns]) unless(all {$MSA_coverages[$_][0] < (2 * $target_coverage_per_haplotype)} @relevantColumns);
					$use_reads_thisGene{$readID}++;
					$used_reads_total++;
					foreach my $column (@relevantColumns)
					{
						# print "\tColumn $column from ", $MSA_coverages[$column][0], " -> ", ($MSA_coverages[$column][0] + 1), "\n";
						die Dumper("Column issue $column", $MSA_coverages[$column]) unless($MSA_coverages[$column][0] < (2 * $target_coverage_per_haplotype));
						# print "Before: \n" . Dumper(\@MSA_coverages), "\n";
						$MSA_coverages[$column][0]++;
						# print "\t\t";
						# foreach my $column (@relevantColumns)
						# {
							# print $column, ": ", $MSA_coverages[$column][0], " ";
						# }				
						# print "\n";
						# die "After: \n" . Dumper(\@MSA_coverages) . "\n";
						die Dumper("Coverage issue", [$geneID, $locus_has_chromosome_resolution{$geneID}], "iteration $iteration", "lookup_index $lookup_index", "column $column", $MSA_coverages[$column], $maxCoverage, (2 * $target_coverage_per_haplotype), \@relevantColumns) unless($MSA_coverages[$column][0] <= (2 * $target_coverage_per_haplotype));
						if(($iteration == 0) and $locus_has_chromosome_resolution{$geneID})
						{
							die unless(defined $lookup_index);
							die unless(($lookup_index == 1) or ($lookup_index == 2));
							$MSA_coverages[$column][$lookup_index]++;
						}
					}
				}
			}
		}
		
		# print Dumper($geneID, \@MSA_coverages), "\n";
		
		my $sum_coverage = sum (map {$_->[0]} @MSA_coverages);
		my $sum_1_coverage = sum (map {$_->[1]} @MSA_coverages);
		my $sum_2_coverage = sum (map {$_->[2]} @MSA_coverages);
		my $min_coverage = min (map {$_->[0]} @MSA_coverages);
		my $max_coverage = max (map {$_->[0]} @MSA_coverages);
		$coverages_by_gene{$geneID}{$outerIteration} = [$min_coverage, $max_coverage, sprintf("%.2f", $sum_coverage / $MSA_length), sprintf("%.2f", $sum_1_coverage / $MSA_length), sprintf("%.2f", $sum_2_coverage / $MSA_length)];
		$use_reads_by_gene{$geneID}{$outerIteration} = \%use_reads_thisGene;
		

		my %local_genotypes_by_readID;
		
		my @selected_readIDs_this_gene = grep {exists $genotypes_by_reads{$_}{$geneID}} keys %use_reads_thisGene;
		
		foreach my $readID (@selected_readIDs_this_gene)
		{
			my @readID_relevantPositions_afterTrimming = @{$readID_relevantPositions_afterTrimming{$readID}{$geneID}};
			die unless(scalar(@readID_relevantPositions_afterTrimming));
			foreach my $position (@readID_relevantPositions_afterTrimming)
			{
				# this mirrors code further down below
				my @alleles = keys %{$genotypes_by_reads{$readID}{$geneID}{$position}};
				# die unless(scalar(@alleles));
				my $allele; 
				if(scalar(@alleles) == 1)
				{
					$allele = $alleles[0];
				}
				else
				{
					$allele = 'N';
				}			
				$local_genotypes_by_readID{$position}{$allele}{$readID}++; 			
			}
		}	
		
		foreach my $column (0 .. $#MSA_coverages) 
		{
			my $effective_reads_this_column = 0;
			my $allele_read_info = (exists $local_genotypes_by_readID{$column}) ? (join(" ", map {
				my $allele = $_;
				my @readIDs = sort keys %{$local_genotypes_by_readID{$column}{$allele}};
				$effective_reads_this_column += scalar(@readIDs);
				$allele . ':' . join('/', @readIDs)
			} keys %{$local_genotypes_by_readID{$column}})) : '';
			print REEADCOVERAGEBYGENE_DETAILS join("\t", $outerIteration, $geneID, $column, $MSA_coverages[$column][0], $MSA_coverages[$column][1], $MSA_coverages[$column][2], $effective_reads_this_column, $allele_read_info), "\n";
		}
		
		@remaining_read_IDs = grep {not exists $use_reads_thisGene{$_}} @remaining_read_IDs;
		print "Size remaining_read_IDs for gene $geneID: ", scalar(@remaining_read_IDs), " (just selected: ", scalar(keys %use_reads_thisGene), " reads)\n";
	}
}
close(OUT_READALLELES);
close(REEADCOVERAGEBYGENE_DETAILS);

open(REEADCOVERAGEBYGENE, '>', $fn_out_readCoverageByGene) or die "Cannot open $fn_out_readCoverageByGene";
print REEADCOVERAGEBYGENE "Read selection summary:\n";
print REEADCOVERAGEBYGENE "\t", "skipped_read_because_multiple_genes", ": ", $skipped_read_because_multiple_genes, "\n";
print REEADCOVERAGEBYGENE "\t", "used_reads_total", ": ", $used_reads_total, "\n";
foreach my $geneID (sort keys %use_reads_by_gene)
{
	foreach my $outerIteration (sort keys %{$use_reads_by_gene{$geneID}})
	{
		my @reads = sort keys %{$use_reads_by_gene{$geneID}{$outerIteration}};
		print REEADCOVERAGEBYGENE "\t\t", $outerIteration, "\n";
		print REEADCOVERAGEBYGENE "\t\t\t", $geneID, ": ", scalar(@reads), " reads, coverage profile ", $coverages_by_gene{$geneID}{$outerIteration}[0], " (min), ", $coverages_by_gene{$geneID}{$outerIteration}[1], " (max), ", $coverages_by_gene{$geneID}{$outerIteration}[2], " (avg), ", $coverages_by_gene{$geneID}{$outerIteration}[3], " / ", $coverages_by_gene{$geneID}{$outerIteration}[4], " (avg. H1/H2)", "\n"; 
	}
}
close(REEADCOVERAGEBYGENE);

open(OUT_READALLELES, '>', $fn_out_readAlleles) or die "Cannot open $fn_out_readAlleles";
open(OUT_ACTIVEALLELESBYREAD, '>', $fn_out_activeAllelePolymorphic_byRead) or die "Cannot open $fn_out_activeAllelePolymorphic_byRead";
foreach my $gene_ref_id (sort keys %translation_targets)
{
	next if($gene_ref_id =~ /pgf/);
	my $gene = geneName($gene_ref_id);
	die unless(exists $alleles_by_gene{$gene});
	unless(exists $alleles_by_gene_forRemapping{$gene})
	{	
		warn "No remapping alleles for gene $gene [$gene_ref_id]" ;
		next;
	}
	
	unless(exists $genotypes_by_MSAcolumn{$gene})
	{
		warn "No reads data for gene $gene [$gene_ref_id]";
		next;
	}
	die unless(exists $alignments_references_href->{$gene_ref_id});
	
	foreach my $outerIteration (sort keys %{$use_reads_by_gene{$gene}})
	{
		my @selected_readIDs = grep {exists $genotypes_by_reads{$_}{$gene}} sort keys %{$use_reads_by_gene{$gene}{$outerIteration}};
		foreach my $readID (@selected_readIDs)
		{
			print OUT_READALLELES join("\t",
				$outerIteration,
				$readID,
				join(' ',
					map {
						my $position = $_;
						my @alleles = sort keys %{$genotypes_by_reads{$readID}{$gene}{$position}};
						my $allele; 
						if(scalar(@alleles) == 1)
						{
							$allele = $alleles[0];
						}
						else
						{
							# todo: check what is going on here!
							warn Dumper("Read alleleles != 1", $readID, $position, \@alleles);
							$allele = 'N';
						}
						$position . ':' . $alleles[0]}
						
					@{$readID_relevantPositions_afterTrimming{$readID}{$gene}}
				)
			), "\n";
			
			my %positions_thisRead_polymorphicActiveAlleles;
			foreach my $position (@{$readID_relevantPositions_afterTrimming{$readID}{$gene}})
			{
				my @activeAlleles = sort keys %{$active_alleles_by_MSAposition_by_gene{$gene}{$position}};
				if(scalar(@activeAlleles) > 1)
				{
					$positions_thisRead_polymorphicActiveAlleles{$position} = [\@activeAlleles, [sort keys %{$genotypes_by_reads{$readID}{$gene}{$position}}]];
				}
			}
			
			print OUT_ACTIVEALLELESBYREAD join("\t",
				$readID,
				scalar(keys %positions_thisRead_polymorphicActiveAlleles),
				join(';', map {$_ . ':[active:' . join('/', @{$positions_thisRead_polymorphicActiveAlleles{$_}[0]}) . ']:[read:' . join('/', @{$positions_thisRead_polymorphicActiveAlleles{$_}[1]}) . ']'} sort keys %positions_thisRead_polymorphicActiveAlleles)
			), "\n"; 
			
		}
	}
	
	# print "Considered alleles gene $gene: \n";
	# foreach my $position (sort {$a <=> $b} keys %positions_considered_alleles)
	# {
		# my @alleles_per_position;
		# print "\tPosition $position:\n";
		# foreach my $allele (keys %{$positions_considered_alleles{$position}})
		# {
			# my $fromMSA = (exists $positions_considered_alleles{$position}{$allele}{fromMSA}) ? 1 : 0;
			# my $fromPileUp = (exists $positions_considered_alleles{$position}{$allele}{fromPileUp}) ? 1 : 0;
			# print "\t\t", $allele, " [from MSA: $fromMSA / from pileup: $fromPileUp ]", "\n";
			# push(@alleles_per_position, [$allele, $fromMSA, $fromPileUp]);
		# }
		
		# my @overlapping_read_IDs = grep {exists $genotypes_by_reads{$_}{$gene}{$position}} @selected_readIDs;
		# foreach my $readID (sort @overlapping_read_IDs)
		# {
			# my @read_alleles = sort keys %{$genotypes_by_reads{$readID}{$gene}{$position}};
			# print "\t\t\t", $readID, " ", join(';', @read_alleles), "\n";
		# }		
	# }
	
	#exit;	
	#die Dumper(\%alleles_present_at_position);
	# $alignments_references_href->{$currentReferenceId}	
}

open(OUT_MSA, '>', $fn_out_MSA) or die "Cannot open $fn_out_MSA";
foreach my $gene (keys %selectedMSAReferenceAlleles_by_gene)
{
	foreach my $seqID (@{$selectedMSAReferenceAlleles_by_gene{$gene}})
	{
		my $seq = $alignments_references_href->{$seqID};
		die "Undefined sequence ID $seqID in \$alignments_references_href" unless(defined $seq);
		print OUT_MSA join("\t", $gene, $seqID, $refSeqID_2_hap{$gene}{$seqID}, $seq), "\n";
	}
}
close(OUT_MSA);

open(OUT_GENES, '>', $fn_out_genes) or die "Cannot open $fn_out_genes";
foreach my $gene_referenceSequenceID (sort keys %translation_targets)
{
	next if($gene_referenceSequenceID =~ /pgf/);
	my $gene = geneName($gene_referenceSequenceID);
	next unless($haveGeneInOutput{$gene});
	my $MSA_length = length($alignments_references_href->{$gene_referenceSequenceID});
	print OUT_GENES join("\t", $gene, $MSA_length), "\n";
}
close(OUT_GENES);

open(OUT_OUTERITERATIONSBYGENE, '>', $fn_out_outerIterationsByGene) or die "Cannot open $fn_out_outerIterationsByGene";
foreach my $gene (sort keys %use_reads_by_gene)
{
	print OUT_OUTERITERATIONSBYGENE join("\t", $gene, (max keys %{$use_reads_by_gene{$gene}})), "\n";
}
close(OUT_OUTERITERATIONSBYGENE);
sub processAlignments
{
	my $alignments_aref = shift;
	
	my @SAM_entries;
	ALIGNMENT: foreach my $read_alignment (@$alignments_aref)
	{
		(my $readID_no12 = $read_alignment->[0][0]) =~ s/\/[12]//;
		my $readID = $read_alignment->[0][0];
		
		(my $currentReferenceId, my $currentReferenceFrom_1based, my $currentReferenceTo_1based) = parseRefInfo($read_alignment->[0][1]);
		die unless($currentReferenceFrom_1based < $currentReferenceTo_1based);
		next if($currentReferenceId eq 'pgf_Ns');
		
		(my $isReverseComplement, my $readFrom, my $readTo) = parseReadInfo($read_alignment->[0][2]);
		die unless($readFrom < $readTo);
		# print $read_alignment->[0][2], " ", $isReverseComplement, "\n";
		
		die unless($read_alignment->[0][3] =~ /\((\d+)\/(\d+)\)/);
		my $nm_before = $1;
		
		my $FLAGs = parseFlags($read_alignment->[0][5]);
		
		my $aligned_reference = $read_alignment->[1];
		my $aligned_read = $read_alignment->[2];
		my $aligned_read_qual = $read_alignment->[3];
		
		(my $aligned_reference_noGaps = $aligned_reference) =~ s/-//g;
		(my $aligned_read_noGaps = $aligned_read) =~ s/-//g;
		
		my @aligned_read_noGaps_indices = grep {substr($aligned_read, $_, 1) ne '-'} (0 .. (length($aligned_read) - 1));
		my $aligned_read_qual_noGaps = join('', map {substr($aligned_read_qual, $_, 1)} @aligned_read_noGaps_indices);
		die unless(length($aligned_read_qual_noGaps) == length($aligned_read_noGaps));
		
		die unless($raw_references_href->{$currentReferenceId} =~ /^[ACGTN]+$/);
		die unless(exists $raw_references_href->{$currentReferenceId});
		my $supposed_aligned_reference_noGaps = substr($raw_references_href->{$currentReferenceId}, $currentReferenceFrom_1based - 1, $currentReferenceTo_1based - $currentReferenceFrom_1based + 1);
		
		die Dumper("Unexpected sequence mismatch (I)", [$currentReferenceId, $currentReferenceFrom_1based, $currentReferenceTo_1based, $supposed_aligned_reference_noGaps], $aligned_reference_noGaps, $aligned_reference) unless($aligned_reference_noGaps eq $supposed_aligned_reference_noGaps);
		
		next unless(exists $alignments_translation_targets_href->{$currentReferenceId});
		my $projectOnto_sequenceId = $alignments_translation_targets_href->{$currentReferenceId}; 
		
		die unless(exists $alignments_references_href->{$currentReferenceId}); 
		die unless(exists $alignments_references_href->{$projectOnto_sequenceId});
		
		my $ref2ref_alignment_start_0based = $raw_references_2AlignmentCoordinates_href->{$currentReferenceId}[$currentReferenceFrom_1based - 1];
		my $ref2ref_alignment_stop_0based = $raw_references_2AlignmentCoordinates_href->{$currentReferenceId}[$currentReferenceTo_1based - 1];
		
		my $ref2ref_alignment_currentRef = substr($alignments_references_href->{$currentReferenceId}, $ref2ref_alignment_start_0based, $ref2ref_alignment_stop_0based - $ref2ref_alignment_start_0based + 1);
		(my $ref2ref_alignment_currentRef_noGaps = $ref2ref_alignment_currentRef) =~ s/_//g;
		die unless($ref2ref_alignment_currentRef_noGaps =~ /^[ACGTN]+$/);
		die unless($ref2ref_alignment_currentRef_noGaps eq $supposed_aligned_reference_noGaps);
		
		my $ref2ref_alignment_newRef = substr($alignments_references_href->{$projectOnto_sequenceId}, $ref2ref_alignment_start_0based, $ref2ref_alignment_stop_0based - $ref2ref_alignment_start_0based + 1);
		
		(my $new_alignment_ref, my $new_alignment_read) = translateAlignment($ref2ref_alignment_newRef, $ref2ref_alignment_currentRef, $aligned_reference, $aligned_read);
		die unless($new_alignment_ref =~ /^[ACGTN\-]+$/);
		die unless($new_alignment_read =~ /^[ACGTN\-]+$/);
		
		die unless(exists $graphs{$projectOnto_sequenceId});
		
		my $running_ref_character = 0;
		
		# (my $new_alignment_ref_noGaps = $new_alignment_ref) =~ s/[_\-]//g;
		# die "Weird characters new_alignment_ref_noGaps: \n $new_alignment_ref_noGaps" unless($new_alignment_ref_noGaps =~ /^[ACGTN]+$/);
		# die Dumper("Sequence mismatch", $new_alignment_ref_noGaps, $ref2ref_alignment_currentRef_noGaps) unless($new_alignment_ref_noGaps eq $ref2ref_alignment_currentRef_noGaps);

		if($projectOnto_sequenceId !~ /pgf/)
		{		
			my $gene = geneName($projectOnto_sequenceId);	

			# todo read print
			# print $readID, " ", $gene, "\n";

			my @MSA_columns;
			my @read_characters;
			
			my $emitEmptyColumn = sub {
				my $thisEmissionColumn = shift;
				my $lastEmissionColumn = shift;
				die unless(defined $thisEmissionColumn);
				die unless(defined $lastEmissionColumn);
				for(my $i = $lastEmissionColumn; $i < ($thisEmissionColumn - 1); $i++)
				{
					push(@MSA_columns, $i+1);
					push(@read_characters, '_');
					$genotypes_by_MSAcolumn{$gene}{$i+1}{'_'}++;
					$genotypes_by_reads{$readID_no12}{$gene}{$i+1}{'_'}++;
					# print $i, " ", "-", "\n";
				}
			};

						
			my $running_readAllele_origin_MSAcolumn;
			my $running_readAllele;
			my $emitted_readAllles_all;
			my $running_ref_characters = 0;	
			my $last_emitted_readAllele_column = undef;
			for(my $i = 0; $i < length($aligned_read); $i++)
			{
				my $c_ref = substr($aligned_reference, $i, 1);	
				my $c_read = substr($aligned_read, $i, 1); 
				
				my $c_ref_isGap = (($c_ref eq '-') or ($c_ref eq '_'));
				die if(($i == 0) and ($c_ref_isGap));
				
				my $c_ref_alignmentColumn = ($c_ref_isGap) ? undef : $raw_references_2AlignmentCoordinates_href->{$currentReferenceId}[$currentReferenceFrom_1based - 1 + $running_ref_characters];
				if(defined $c_ref_alignmentColumn)
				{
					my $alignmentColumn_fromReferenceMSA = substr($alignments_references_href->{$currentReferenceId}, $c_ref_alignmentColumn, 1);
					die unless($c_ref eq $alignmentColumn_fromReferenceMSA);
				}
				die if(($i > 0) and $running_readAllele_origin_MSAcolumn and (not defined $running_readAllele_origin_MSAcolumn));
				
				if(not $c_ref_isGap)
				{
					if($running_readAllele)
					{
						if(defined $last_emitted_readAllele_column)
						{
							$emitEmptyColumn->($running_readAllele_origin_MSAcolumn, $last_emitted_readAllele_column);						
						}
						push(@MSA_columns, $running_readAllele_origin_MSAcolumn);
						push(@read_characters, $running_readAllele);		
						$genotypes_by_MSAcolumn{$gene}{$running_readAllele_origin_MSAcolumn}{$running_readAllele}++;	
						$genotypes_by_reads{$readID_no12}{$gene}{$running_readAllele_origin_MSAcolumn}{$running_readAllele}++;								
						# print $running_readAllele_origin_MSAcolumn, " ", $running_readAllele, "\n";
						$emitted_readAllles_all .= $running_readAllele;
						$last_emitted_readAllele_column = $running_readAllele_origin_MSAcolumn;
					}
					
					$running_readAllele = '';
					$running_readAllele_origin_MSAcolumn = $c_ref_alignmentColumn;
				}
				
				$running_readAllele .= $c_read;
				$running_ref_characters++ if(not $c_ref_isGap);				
			}
			if($running_readAllele)
			{
				if(defined $last_emitted_readAllele_column)
				{
					$emitEmptyColumn->($running_readAllele_origin_MSAcolumn, $last_emitted_readAllele_column);				
				}
				$genotypes_by_MSAcolumn{$gene}{$running_readAllele_origin_MSAcolumn}{$running_readAllele}++;
				$genotypes_by_reads{$readID_no12}{$gene}{$running_readAllele_origin_MSAcolumn}{$running_readAllele}++;
				push(@MSA_columns, $running_readAllele_origin_MSAcolumn);
				push(@read_characters, $running_readAllele);						
				# print $running_readAllele_origin_MSAcolumn, " ", $running_readAllele, "\n";
				$emitted_readAllles_all .= $running_readAllele;
			}
			die unless($emitted_readAllles_all eq $aligned_read);
			
			die unless(scalar(@MSA_columns) == scalar(@read_characters));
			
			# check that columns are contiguous
			my $last_column = $MSA_columns[$#MSA_columns];
			my $first_column = $MSA_columns[0];
			die unless($first_column < $last_column);
			
			my $covered_columns = $last_column - $first_column + 1;
			die Dumper("Non-contiguous MSA columns [1]", \@MSA_columns) unless(scalar(@MSA_columns) == $covered_columns);
			for(my $i = 0; $i < $#MSA_columns; $i++)
			{
				die Dumper("Non-contiguous MSA columns [2]", $i, $MSA_columns[$i], $MSA_columns[$i+1], $aligned_reference, $aligned_read)  unless($MSA_columns[$i] == ($MSA_columns[$i+1] - 1));
			}
			
			my $min_MSA_column_wholeRead = (exists $readID_2_start_stop_MSAcolumns{$readID_no12}{$gene}) ? min($readID_2_start_stop_MSAcolumns{$readID_no12}{$gene}[0], $first_column) : $first_column;
			my $max_MSA_column_wholeRead = (exists $readID_2_start_stop_MSAcolumns{$readID_no12}{$gene}) ? max($readID_2_start_stop_MSAcolumns{$readID_no12}{$gene}[1], $last_column) : $last_column;

			$readID_2_start_stop_MSAcolumns{$readID_no12}{$gene} = [$min_MSA_column_wholeRead, $max_MSA_column_wholeRead];
			
			die Dumper("No reference haplotype data for $gene / $currentReferenceId", [keys %refSeqID_2_hap], [keys %{$refSeqID_2_hap{$gene}}]) unless(exists $refSeqID_2_hap{$gene}{$currentReferenceId});
			if((exists $readID_2_source_haplotype{$readID_no12}{$gene}) and ($readID_2_source_haplotype{$readID_no12}{$gene} ne $refSeqID_2_hap{$gene}{$currentReferenceId}))
			{
				$readID_2_source_haplotype{$readID_no12}{$gene} = 0;
			}
			else
			{
				$readID_2_source_haplotype{$readID_no12}{$gene} = $refSeqID_2_hap{$gene}{$currentReferenceId};
			}
			# todo read print
			# print "\t", join(' ', @MSA_columns), "\n";
			# print "\t", join(' ', @read_characters), "\n";
					
		}
		
		{
			my $running_readAllele;
			my $running_readAllele_origin_nodeID;
			my $running_ref_characters = 0;
			for(my $i = 0; $i < length($aligned_read); $i++)
			{
				my $c_ref = substr($aligned_reference, $i, 1);	
				my $c_read = substr($aligned_read, $i, 1); 
				
				my $c_ref_isGap = (($c_ref eq '-') or ($c_ref eq '_'));
				my $c_ref_alignmentColumn = ($c_ref_isGap) ? undef : $raw_references_2AlignmentCoordinates_href->{$currentReferenceId}[$currentReferenceFrom_1based - 1 + $running_ref_characters];
				
				# print join("\t", $i, $c_ref_alignmentColumn, $c_ref, $c_read, $running_ref_characters, $c_ref_alignmentColumn), "\n";

				my $node_before_refC = ($c_ref_isGap) ? undef : 'node' . $c_ref_alignmentColumn;
				my $node_after_refC = ($c_ref_isGap) ? undef : 'node' . ($c_ref_alignmentColumn + 1);
				die Dumper("Missing node (I) $node_before_refC", $c_ref_alignmentColumn, $currentReferenceFrom_1based - 1 + $running_ref_characters, $currentReferenceFrom_1based, $running_ref_characters) unless((not defined $node_before_refC) or (exists $graphs{$projectOnto_sequenceId}{nodes}{$node_before_refC}));
				die Dumper("Missing node (II) $node_after_refC", $currentReferenceFrom_1based - 1 + $running_ref_characters, $currentReferenceFrom_1based, $running_ref_characters) unless((not defined $node_after_refC) or (exists $graphs{$projectOnto_sequenceId}{nodes}{$node_after_refC}));
				
				if((not defined $running_readAllele_origin_nodeID) and (defined $node_before_refC))
				{
					die Dumper("Gap processing issue", $i, [$c_ref, $c_read], [$aligned_reference, $aligned_read], [$new_alignment_ref, $new_alignment_read])	if($c_ref_isGap);
					$running_readAllele_origin_nodeID = $node_before_refC;
				}
				
				if((not $c_ref_isGap) and $running_readAllele)
				{
					# print "\t", "Add ", $running_readAllele_origin_nodeID, " -> ", $running_readAllele_origin_nodeID, "\n";
					
					die unless((defined $running_readAllele_origin_nodeID) and (defined $node_before_refC));
					die Dumper("Something weird is going on here", $running_readAllele_origin_nodeID, $node_before_refC, $i, [$aligned_read, $aligned_reference], $running_readAllele)  if($running_readAllele_origin_nodeID eq $node_before_refC);
					# my $novel = exists $graphs{$projectOnto_sequenceId}{nodes}{$running_readAllele_origin_nodeID}{links}{$node_after_refC}{$running_Allele};
					# print "Adding link $nodeI_before -> $nodeI_after from $currentReferenceId with allele $running_read_allele (novel: $novel)\n";
					
					$graphs{$projectOnto_sequenceId}{nodes}{$running_readAllele_origin_nodeID}{links}{$node_before_refC}{$running_readAllele}++;
					$graphs{$projectOnto_sequenceId}{nodes}{$running_readAllele_origin_nodeID}{coverage}++;
					
					$running_readAllele = '';
					$running_readAllele_origin_nodeID = undef;
				}
				
				if((not defined $running_readAllele_origin_nodeID) and (defined $node_before_refC))
				{
					die Dumper("Gap processing issue", $i, [$c_ref, $c_read], [$aligned_reference, $aligned_read], [$new_alignment_ref, $new_alignment_read])	if($c_ref_isGap);
					$running_readAllele_origin_nodeID = $node_before_refC;
				}
				
				$running_readAllele .= $c_read;
				$running_ref_characters++ if(not $c_ref_isGap);
						
				
				# if(($i == (length($aligned_read) - 1)) and not $c_ref_isGap)
				# {
					# die unless((defined $node_before_refC) and (defined $node_after_refC));				
					# $graphs{$projectOnto_sequenceId}{nodes}{$running_readAllele_origin_nodeID}{links}{$node_after_refC}{$running_readAllele}++;
					# $graphs{$projectOnto_sequenceId}{nodes}{$running_readAllele_origin_nodeID}{coverage}++;
					# $graphs{$projectOnto_sequenceId}{nodes}{$node_after_refC}{coverage}++;
									 
				# }
			}
		}
	}
}

print "\n\nDone. Produced files ${outputPrefix}.*\n\n";

# my $fn_output_BAM = 'bla.bam'; 
# my $cmd_samtools_sort = qq(module load SamTools; samtools sort -o $fn_output_BAM $outputSAM; samtools index $fn_output_BAM);
# system($cmd_samtools_sort) and die "Could not execute: $cmd_samtools_sort\n";

# print "\n\nDone. Produced $fn_output_BAM\n\n";

sub printGraph
{
	my $graph = shift;
	die unless(exists $graph->{nodes});
	
	my @nodeIDs = keys %{$graph->{nodes}};
	my %knowNode = map {$_ => 1} @nodeIDs;
	
	die unless(all {exists $graph->{nodes}{$_}{level}} @nodeIDs);
	@nodeIDs = sort {$graph->{nodes}{$a}{level} <=> $graph->{nodes}{$b}{level}} @nodeIDs;
	
	print "Graph printout\n";
	foreach my $nodeID (@nodeIDs)
	{
		print "Node $nodeID - level $graph->{nodes}{$nodeID}{level} - coverage $graph->{nodes}{$nodeID}{coverage}\n";
		die unless(exists $graph->{nodes}{$nodeID}{links});
		die unless(exists $graph->{nodes}{$nodeID}{coverage});
		my @outlinks = keys %{$graph->{nodes}{$nodeID}{links}};
		foreach my $outlink_nodeID (@outlinks)
		{
			print "\t--> $outlink_nodeID\n";
			die unless(exists $knowNode{$outlink_nodeID});
			my @edge_labels = keys %{$graph->{nodes}{$nodeID}{links}{$outlink_nodeID}};
			die unless(scalar(@edge_labels));
			foreach my $edge_label (@edge_labels)
			{
				my $edge_coverage = $graph->{nodes}{$nodeID}{links}{$outlink_nodeID}{$edge_label};
				print "\t\t\t $edge_label [$edge_coverage]", "\n";
			}
						
		}
	}
}
sub compressCIGAR
{
	my $CIGAR_in = shift;
	my $runningSymbol;
	my $runningSymbol_n;
	
	my $CIGAR_out = '';
	foreach my $c (split(//, $CIGAR_in))
	{
		if((not defined $runningSymbol) or ($runningSymbol ne $c))
		{
			if(defined $runningSymbol)
			{
				die unless(defined $runningSymbol_n);
				$CIGAR_out .= $runningSymbol_n . $runningSymbol;
			}
			
			$runningSymbol = $c;
			$runningSymbol_n = 0;
		}
		$runningSymbol_n++;
	}
	if(defined $runningSymbol)  
	{
		die unless(defined $runningSymbol_n);
		$CIGAR_out .= $runningSymbol_n . $runningSymbol;
	}
	return $CIGAR_out;
}

sub translateAlignment
{
	my $refref_new = shift;
	my $refref_old = shift;
	my $oldref_ref = shift;
	my $oldref_read = shift;
	
	die unless(length($refref_new) == length($refref_old));
	die unless(length($oldref_ref) == length($oldref_read));
	
	(my $refref_old_noGaps = $refref_old) =~ s/[_\-]//g;
	(my $oldref_ref_noGaps = $oldref_ref) =~ s/[_\-]//g;
	
	die Dumper("Unexpected sequence mismatch (II)", $refref_old_noGaps, $oldref_ref_noGaps)  unless($refref_old_noGaps eq $oldref_ref_noGaps);
	
	my @oldref_ref_noGaps_alignedCharacters = ('') x (length($oldref_ref_noGaps) + 1);
	die unless($#oldref_ref_noGaps_alignedCharacters == length($oldref_ref_noGaps));
	
	my $running_ref_character = 0;
	for(my $i = 0; $i < length($oldref_ref); $i++)
	{
		my $c_ref = substr($oldref_ref, $i, 1);
		my $c_read = substr($oldref_read, $i, 1);
		
		if(($c_ref ne '-') and ($c_ref ne '_'))
		{
			$running_ref_character++;
		}
		
		$oldref_ref_noGaps_alignedCharacters[$running_ref_character] .= $c_read;
	}
	
	die unless(join('', @oldref_ref_noGaps_alignedCharacters) eq  $oldref_read);
	
	my $forReturn_ref_alignment = getEquivalentGapCharacters($oldref_ref_noGaps_alignedCharacters[0]);
	my $forReturn_read_alignment = $oldref_ref_noGaps_alignedCharacters[0];
	
	my $runningCharacter_oldRef = -1;
	for(my $i = 0; $i < length($refref_new); $i++)
	{
		my $c_ref_new = substr($refref_new, $i, 1);
		my $c_ref_old = substr($refref_old, $i, 1);
		
		my $new_ref_gap = (($c_ref_new eq '-') or ($c_ref_new eq '_'));
		my $old_ref_gap = (($c_ref_old eq '-') or ($c_ref_old eq '_'));
		
		next if($new_ref_gap and $old_ref_gap);
		
		$runningCharacter_oldRef++ if(! $old_ref_gap);
		
		my $toAdd_ref;
		my $toAdd_read;
		if($new_ref_gap and !$old_ref_gap)
		{
			$toAdd_read = $oldref_ref_noGaps_alignedCharacters[$runningCharacter_oldRef + 1];		
			$toAdd_ref = getEquivalentGapCharacters($toAdd_read); 
		}
		elsif(!$new_ref_gap and $old_ref_gap)
		{
			$toAdd_ref = $c_ref_new;
			$toAdd_read = '-';
		}
		else
		{
			die unless((!$new_ref_gap) and (!$old_ref_gap));
			$toAdd_ref = $c_ref_new;
			$toAdd_read = $oldref_ref_noGaps_alignedCharacters[$runningCharacter_oldRef + 1];
			my $missing_gaps = length($toAdd_read) - length($toAdd_ref);
			die unless($missing_gaps >= 0);
			$toAdd_ref .= getNGaps($missing_gaps);
		}
		
		die unless(length($toAdd_ref) == length($toAdd_read));
		$forReturn_ref_alignment .= $toAdd_ref;
		$forReturn_read_alignment .= $toAdd_read;
	}
	
	(my $oldref_read_noGaps = $oldref_read) =~ s/[\-_]//g;
	(my $forReturn_read_alignment_noGaps = $forReturn_read_alignment) =~ s/\-//g;
	die unless($forReturn_read_alignment_noGaps eq $oldref_read_noGaps);
	
	(my $forReturn_ref_alignment_noGaps = $forReturn_ref_alignment) =~ s/\-//g;
	(my $refref_new_noGaps = $refref_new) =~ s/[\-_]//g;
	die unless($forReturn_ref_alignment_noGaps eq $refref_new_noGaps);
	
	return ($forReturn_ref_alignment, $forReturn_read_alignment);
}	

sub getNGaps
{
	my $n = shift;
	my $forReturn = ('-' x $n);	
	die unless(length($forReturn) == $n);
	return $forReturn;
}

sub getEquivalentGapCharacters
{
	my $input = shift;
	return getNGaps(length($input));
}
	
sub parseReadInfo
{
	my $input = shift;
	die "Weird read info: $input" unless($input =~ /^read:(\d+)-(\d+)$/);
	if($1 < $2)
	{
		return (0, $1, $2);
	}
	else
	{
		return (1, $2, $1);
	}
}	

sub parseRefInfo
{
	my $input = shift;
	die "Weird reference info: $input" unless($input =~ /^(\S+):(\d+)-(\d+)$/);
	return ($1, $2, $3);
}

sub parseFlags
{
	my $input = shift;
	die "Weird reference info: $input" unless($input =~ /^flags=(\d+)$/);
	return $1;
}

sub read_raw_references
{
	my @filesToRead = glob($full_graph_dir . '/pseudoGenomic_fullLengthMapping/raw_*');
	my $combined_raw_references_href = {};
	foreach my $file (@filesToRead)
	{
		my $thisFile_href = readFASTA($file, 1);
		foreach my $seqID (keys %$thisFile_href)
		{
			next unless($seqID =~ /\*/);
			die "Duplicate sequence ID $seqID (now in file $file" if(exists $combined_raw_references_href->{$seqID});
			$combined_raw_references_href->{$seqID} = $thisFile_href->{$seqID};
		}
	}
	
	my $fn_pgf_n_masked = $full_graph_dir . '/pseudoGenomic_fullLengthMapping/PGF_with_Ns.fa';
	my $pgf_href = readFASTA($fn_pgf_n_masked);
	die unless(scalar(keys %$pgf_href) == 1);
	foreach my $seqID (keys %$pgf_href)
	{
		die "Duplicate sequence ID $seqID (now in file $fn_pgf_n_masked" if(exists $combined_raw_references_href->{$seqID});
		$combined_raw_references_href->{$seqID} = $pgf_href->{$seqID};
	}
	
	return $combined_raw_references_href;
}	

sub read_aligned_references
{
	my @filesToRead = glob($full_graph_dir . '/pseudoGenomic_fullLengthMapping/alignments_*');
	my $combined_aligned_references_href = {};
	foreach my $file (@filesToRead)
	{
		next if($file =~ /controlUnmodified/);		
		my $thisFile_href = readFASTA($file, 1);
		foreach my $seqID (keys %$thisFile_href)
		{
			next unless($seqID =~ /\*/);		
			die if(exists $combined_aligned_references_href->{$seqID});
			$combined_aligned_references_href->{$seqID} = $thisFile_href->{$seqID};
		}
	}
	
	my $fn_pgf_n_masked = $full_graph_dir . '/pseudoGenomic_fullLengthMapping/PGF_with_Ns.fa';
	my $pgf_href = readFASTA($fn_pgf_n_masked);
	die unless(scalar(keys %$pgf_href) == 1);
	foreach my $seqID (keys %$pgf_href)
	{
		die "Duplicate sequence ID $seqID (now in file $fn_pgf_n_masked" if(exists $combined_aligned_references_href->{$seqID});
		$combined_aligned_references_href->{$seqID} = $pgf_href->{$seqID};
	}
	
	return $combined_aligned_references_href;
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
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}

sub compute_translation_targets
{
	my $alignments_references_href = shift;
	my $translation_targets_href = {};
	foreach my $key (keys %$alignments_references_href)
	{
		if($key eq 'pgf_Ns')
		{
			$translation_targets_href->{$key} = $key;
		}
		else
		{
			die unless($key =~ /^(\S+?)\*(\S+)$/);
			my $designated_target = $1 . '*ref';
			if(exists $alignments_references_href->{$designated_target})
			{
				$translation_targets_href->{$key} = $designated_target;
			}
			else
			{
				warn "No target for $key";
			}
		}
	}
	return $translation_targets_href;
}


sub compute_alignment_coordinates
{
	my $raw_references_href = shift;
	my $alignments_references_href = shift;
	
	die unless(scalar(keys %$raw_references_href) == scalar(keys %$alignments_references_href));
	die unless(all {exists $alignments_references_href->{$_}} keys %$raw_references_href);
	
	my %raw_2_alignment_allSequences;
	
	foreach my $seqID (keys %$alignments_references_href)
	{
		my @raw_2_alignment;
		$#raw_2_alignment = length($raw_references_href->{$seqID}) - 1;
		my $alignment_seq = $alignments_references_href->{$seqID};
		die "Weird characters in alignment: $alignment_seq" unless($alignment_seq =~ /^[ACGTN_]+$/);
		my $currentRawCoordinate = -1;
		for(my $i = 0; $i < length($alignment_seq); $i++)
		{
			if(substr($alignment_seq, $i, 1) ne '_')
			{
				$currentRawCoordinate++;
				$raw_2_alignment[$currentRawCoordinate] = $i;
			}
		}
		die unless(all {defined $_} @raw_2_alignment);
		$raw_2_alignment_allSequences{$seqID} = \@raw_2_alignment;
	}
	
	return \%raw_2_alignment_allSequences;
}


sub geneName
{
	my $alleleName = shift;
	die "Weird gene ID $alleleName" unless($alleleName =~ /^(\w+)\*/); 
	my $gene = $1;	
	return $gene;
}