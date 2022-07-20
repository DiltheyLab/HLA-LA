use warnings;
use strict;
use FindBin;
use Getopt::Long;
use Data::Dumper; 
use List::Util qw/min max/;
use List::MoreUtils qw/all mesh/;
use Cwd qw/abs_path/;
use lib '.';

use Util;
use VCFFunctions;
$| = 1;

my $this_bin_dir = $FindBin::RealBin;
my $graphs_dir = $this_bin_dir . '/../../graphs/';

my $g_file = $this_bin_dir . '/../hla_nom_g.txt';
die "Missing file: $g_file" unless(-e $g_file);

my $wgsim_bin = 'wgsim';
my $bwa_bin = 'bwa';
my $samtools_bin = '/usr/bin/samtools';
my $minimap2_bin = '/usr/bin/minimap2';


my $HLA_LA_bin = qq($this_bin_dir/../HLA-LA.pl);
die unless(-e $HLA_LA_bin);
$HLA_LA_bin = abs_path($HLA_LA_bin);

my $working_dir = qq($this_bin_dir/../../working/);
die unless(-d $working_dir);
$working_dir = abs_path($working_dir);

my $globalAlignment_bin = qq($this_bin_dir/../globalAlignment.pl);
die unless(-e $globalAlignment_bin);
$globalAlignment_bin = abs_path($globalAlignment_bin);

my $includes=" -I" . join(" -I", (@INC, abs_path($this_bin_dir . '/..')));

my $GRCh38_refGenome = '/home/dilthey/GRCh38_full_analysis_set_plus_decoy_hla.fa';
my $GRCh38_refGenome_primary = '/home/dilthey/Homo_sapiens.GRCh38.dna.primary_assembly.fa';
die unless(-e $GRCh38_refGenome);
die unless(-e $GRCh38_refGenome_primary);

my $hla_gen = '../hla_gen.fasta';
my $graph = 'PRG_MHC_GRCh38_withIMGT';
my $outputPrefix = '../../haplotypeSimulations';
my $n = "5,5,5";
my @mutationRates = (0, 1/5000, 1/500);
my $targetCoverage = 15;
my $readLength = 100;
my $action = 'evaluate';
my $action2 = 'do';
my $individualPrefix = 'flS_Indiv';
my $evaluateMergeMode = 2;
my $evaluateDirectory = 'remap';
GetOptions (
	'hla_gen:s' => \$hla_gen,
	'individualPrefix:s' => \$individualPrefix,
	'graph:s' => \$graph,
	'action:s' => \$action,
	'action2:s' => \$action2,
	'outputPrefix:s' => \$outputPrefix,
	'targetCoverage:s' => \$targetCoverage,
	'n:s' => \$n,
	'evaluateMergeMode:s' => \$evaluateMergeMode,
	'evaluateDirectory:s' => \$evaluateDirectory,
);
	
die "Please specify --action (simulate, applyAll)" unless($action);

my $graph_dir = $graphs_dir . '/' . $graph;

unless(-d $outputPrefix)
{
	die "Please provide existing directory as output prefix (--outputPrefix)";
}

my @mutationRates_n = split(/,/, $n);
die unless(scalar(@mutationRates) == scalar(@mutationRates_n));

my $outputFn_samples = $outputPrefix . '/' . $individualPrefix . '_samples.txt';
my $genomic_alignments_dir = $this_bin_dir . '/../IMGT_genomic_alignments';
die "Missing directory §genomic_alignments_dir" unless(-d $genomic_alignments_dir);

print "Action: $action\n";
if($action eq 'evaluate')
{
	print "--evaluateMergeMode $evaluateMergeMode / --evaluateDirectory $evaluateDirectory\n";
	
	my $sampleID = 'simWithNovel7';
	
	my $sampleID_workingDir = $working_dir . '/simWithNovel7';
	die "Directory $sampleID_workingDir not existing - has inference been carried out?" unless(-d $sampleID_workingDir);
	
	my $prefix_inference = $sampleID_workingDir . '/' . $evaluateDirectory . '/haplotypeHMMInput.fullLengthInference.mm' . $evaluateMergeMode;
	
	my $fasta_haplotypes = $prefix_inference . '.fasta';
	my $fasta_haplotypes_graphLevels = $prefix_inference . '.fasta.graphLevels';
	die "Missing file '$fasta_haplotypes' - has the second inference step been carried out?" unless(-e $fasta_haplotypes);
	die "Missing file '$fasta_haplotypes_graphLevels' - has the second inference step been carried out?" unless(-e $fasta_haplotypes_graphLevels);
	print "\tReading file: $fasta_haplotypes\n";
	
	my $prefix_truth = '/home/dilthey/HLA-LA-devel/haplotypeSimulations/simWithNovel';
	my $fn_truth_details = $prefix_truth . '_HLAHaplotypeDetails.txt';
	die "Missing truth file: '$fn_truth_details'" unless(-e $fn_truth_details);

	my %truth;
	open(TRUTH, '<', $fn_truth_details) or die "Cannot open $fn_truth_details";
	my $truth_headerLine = <TRUTH>;
	chomp($truth_headerLine);
	my @truth_header_fields = split(/\t/, $truth_headerLine);
	while(<TRUTH>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless(scalar(@line_fields) == scalar(@truth_header_fields));
		my %line_hash = (mesh @truth_header_fields, @line_fields);
		my $individualID = $line_hash{IndividualID};
		my $locus = $line_hash{Locus};
		my $haplotype = $line_hash{Haplotype};
		$truth{$individualID}{$locus}{$haplotype} = \%line_hash;
	}
	close(TRUTH);
	my $inferred_haplotypes_graphLevels_href = Util::readFASTA($fasta_haplotypes_graphLevels);
	my %inferred_haplotypes_withGraphLevelSplit;
	my %inferred_haplotypes_graphLevels;
	#my %inferred_haplotypes_graphLevels;
	foreach my $fastaID (keys %$inferred_haplotypes_graphLevels_href)
	{
		if($fastaID =~ /Levels/)
		{
			die unless($fastaID =~ /^(\w+)-Levels$/);	
			$inferred_haplotypes_withGraphLevelSplit{$1}{Levels} = $inferred_haplotypes_graphLevels_href->{$fastaID};
		}
		else
		{
			die unless($fastaID =~ /^(\w+)-(H(1|2))$/);
			die unless(($3 eq '1') or ($3 eq '2'));
			$inferred_haplotypes_withGraphLevelSplit{$1}{$3} = $inferred_haplotypes_graphLevels_href->{$fastaID};
		}
	}
	die "Issue with files $fasta_haplotypes / $fasta_haplotypes_graphLevels" unless(all {(exists $inferred_haplotypes_withGraphLevelSplit{$_}{1}) and (exists $inferred_haplotypes_withGraphLevelSplit{$_}{2}) and (exists $inferred_haplotypes_withGraphLevelSplit{$_}{Levels})} keys %inferred_haplotypes_withGraphLevelSplit);
	print "Found " . scalar(keys %inferred_haplotypes_withGraphLevelSplit) . " genes with inference.\n";
	
	my %inferred_haplotypes;
	my %inferred_haplotypes_byGraphLevel;
	my %inferred_haplotypes_byGraphLevel_asHashRef;
	my %inferred_genotypes_byGraphLevel;
	foreach my $locus (sort keys %inferred_haplotypes_withGraphLevelSplit)
	{
		die unless(defined $inferred_haplotypes_withGraphLevelSplit{$locus}{1});
		my @h1_split = split(/;/, $inferred_haplotypes_withGraphLevelSplit{$locus}{1});
		my @h2_split = split(/;/, $inferred_haplotypes_withGraphLevelSplit{$locus}{2});
		my @levels_split = split(/;/, $inferred_haplotypes_withGraphLevelSplit{$locus}{Levels});
		die unless(scalar(@h1_split) == scalar(@h2_split));
		die unless(scalar(@h1_split) == scalar(@levels_split));
		my $h1_raw = join('', @h1_split);
		my $h2_raw = join('', @h2_split);
		$inferred_haplotypes{$locus}{1} = $h1_raw;
		$inferred_haplotypes{$locus}{2} = $h2_raw;
		for(my $i = 0; $i <= $#h1_split; $i++)
		{
			$inferred_haplotypes_byGraphLevel{$locus}{1}[$i] = $h1_split[$i];
			$inferred_haplotypes_byGraphLevel{$locus}{2}[$i] = $h2_split[$i];
			$inferred_haplotypes_byGraphLevel_asHashRef{$locus}{1}{$i} = $h1_split[$i];
			$inferred_haplotypes_byGraphLevel_asHashRef{$locus}{2}{$i} = $h2_split[$i];
			$inferred_genotypes_byGraphLevel{$locus}[$i] = [$h1_split[$i], $h2_split[$i]];
		} 
	}
	
	
	foreach my $locus (sort keys %inferred_haplotypes)
	{
		unless((exists $truth{$sampleID}) and (exists $truth{$sampleID}{$locus}))
		{
			warn "Combination $sampleID / $locus not found in truth data"; 
			next;
		}
		
		print "\t", $locus, "\n";
		my $tempFn_prefix = $fasta_haplotypes . '.comparison.tmpAln.' . $locus;		
		
		
		my @inferred_haplotypes = map {$inferred_haplotypes{$locus}{$_}} (1, 2);
		my @truth_haplotypes = map {$truth{$sampleID}{$locus}{$_}{underlyingSequence_aligned_postMutationSeq}} (1, 2);
		my @truth_haplotypes_preMutation = map {$truth{$sampleID}{$locus}{$_}{underlyingSequence_aligned_preMutationSeq}} (1, 2);
		my @truth_haplotypes_alignmentCoordinates = map {[split(/;/, $truth{$sampleID}{$locus}{$_}{underlyingSequence_aligned_alignmentCoordinates})]} (1, 2);
		my @truth_haplotypes_noGaps = map {my $h = $_; $h =~ s/[\-_]//g; $h} @truth_haplotypes;
		
		my @inferred_haplotypes_byGraphLevel = map {$inferred_haplotypes_byGraphLevel{$locus}{$_}} (1, 2);
		my @thisGene_inferred_haplotypes_byGraphLevel_asHashRef = map {$inferred_haplotypes_byGraphLevel_asHashRef{$locus}{$_}} (1, 2);
		my $inferred_genotypes_byGraphLevel_aref = $inferred_genotypes_byGraphLevel{$locus};

		my @truth_haplotypes_byGraphLevel;
		my @truth_haplotypes_byGraphLevel_preMutation;
		my $minGraphLevel_bothHaplotypes;
		my $maxGraphLevel_bothHaplotypes;
		foreach my $h (1, 2)
		{
			die unless(length($truth_haplotypes[$h-1]) == scalar(@{$truth_haplotypes_alignmentCoordinates[$h-1]}));
			die unless(length($truth_haplotypes_preMutation[$h-1]) == scalar(@{$truth_haplotypes_alignmentCoordinates[$h-1]}));
			my $runningGraphLevel;
			for(my $i = 0; $i < length($truth_haplotypes[$h-1]); $i++)
			{
				if($truth_haplotypes_alignmentCoordinates[$h-1][$i] != -1)
				{
					if(defined $runningGraphLevel)
					{
						die Dumper("Graph level weirdness", $h, $i, [@{$truth_haplotypes_alignmentCoordinates[$h-1]}[$i - 3 .. $i + 3]]) unless($truth_haplotypes_alignmentCoordinates[$h-1][$i] > $runningGraphLevel);
					}
					$runningGraphLevel = $truth_haplotypes_alignmentCoordinates[$h-1][$i];
				}
				die unless(defined $runningGraphLevel);
				$truth_haplotypes_byGraphLevel[$h-1]{$runningGraphLevel} = '' unless(exists $truth_haplotypes_byGraphLevel[$h-1]{$runningGraphLevel});
				$truth_haplotypes_byGraphLevel[$h-1]{$runningGraphLevel} .= substr($truth_haplotypes[$h-1], $i, 1);
				
				$truth_haplotypes_byGraphLevel_preMutation[$h-1]{$runningGraphLevel} = '' unless(exists $truth_haplotypes_byGraphLevel_preMutation[$h-1]{$runningGraphLevel});
				$truth_haplotypes_byGraphLevel_preMutation[$h-1]{$runningGraphLevel} .= substr($truth_haplotypes_preMutation[$h-1], $i, 1);				
			}
			my $min_graphLevel = min(keys %{$truth_haplotypes_byGraphLevel[$h-1]});
			my $max_graphLevel = max(keys %{$truth_haplotypes_byGraphLevel[$h-1]});
			
			die Dumper("Missing some graph levels", $min_graphLevel, $max_graphLevel) unless(all {exists $truth_haplotypes_byGraphLevel[$h-1]{$_}} ($min_graphLevel .. $max_graphLevel));
			
			if(defined $minGraphLevel_bothHaplotypes)
			{
				die Dumper("Graph Level mismatch (min)", $minGraphLevel_bothHaplotypes, $min_graphLevel, [@{$truth_haplotypes_alignmentCoordinates[0]}[0..5]], [@{$truth_haplotypes_alignmentCoordinates[1]}[0..5]]) unless($minGraphLevel_bothHaplotypes == $min_graphLevel);
				die Dumper("Graph Level mismatch (max)", $maxGraphLevel_bothHaplotypes, $max_graphLevel) unless($maxGraphLevel_bothHaplotypes == $max_graphLevel);
			}
			else
			{
				$minGraphLevel_bothHaplotypes = $min_graphLevel;
				$maxGraphLevel_bothHaplotypes = $max_graphLevel;				
			}
		}
		
		my %graphLevels_with_mutation;
		foreach my $h (1, 2)
		{
			foreach my $graphLevel (sort keys %{$truth_haplotypes_byGraphLevel[$h-1]})
			{
				if(($graphLevel >= $minGraphLevel_bothHaplotypes) and ($graphLevel <= $maxGraphLevel_bothHaplotypes))
				{
					if($truth_haplotypes_byGraphLevel[$h-1]{$graphLevel} ne  $truth_haplotypes_byGraphLevel_preMutation[$h-1]{$graphLevel})
					{
						my $mutationID = 'h' . $h . ':' . $graphLevel . ':' . $truth_haplotypes_byGraphLevel_preMutation[$h-1]{$graphLevel} . '->' . $truth_haplotypes_byGraphLevel[$h-1]{$graphLevel};
						$graphLevels_with_mutation{$graphLevel}{$mutationID}++;
					}
				}
			}
		}
		
		my $truth_genotypes_byGraphLevel_href = {};
		foreach my $h (1, 2)
		{
			foreach my $graphLevel (sort keys %{$truth_haplotypes_byGraphLevel[$h-1]})
			{
				push(@{$truth_genotypes_byGraphLevel_href->{$graphLevel}}, $truth_haplotypes_byGraphLevel[$h-1]{$graphLevel});
			}
		}
		
		my $haplotypeComparison = sub {
			my $ref_href = shift;
			my $inference_href = shift;
			die unless(defined $ref_href);
			die unless(defined $inference_href);
			my @differences;
			for(my $graphLevel = $minGraphLevel_bothHaplotypes; $graphLevel <= $maxGraphLevel_bothHaplotypes; $graphLevel++)
			{
				die unless(defined $ref_href->{$graphLevel});
				die unless(defined $inference_href->[$graphLevel]);
				if($ref_href->{$graphLevel} ne $inference_href->[$graphLevel])
				{
					my $surroundingLower = $graphLevel - 10; $surroundingLower = $minGraphLevel_bothHaplotypes if($surroundingLower < $minGraphLevel_bothHaplotypes);
					my $surroundingUpper = $graphLevel + 10; $surroundingUpper = $maxGraphLevel_bothHaplotypes if($surroundingUpper > $maxGraphLevel_bothHaplotypes);
					my $surrounding_ref = join('', map {$ref_href->{$_}} ($surroundingLower .. $surroundingUpper));
					my $surrounding_inference = join('', map {$inference_href->[$_]} ($surroundingLower .. $surroundingUpper));
					
					$surrounding_ref =~ s/[\-_]//g;
					$surrounding_inference =~ s/[\-_]//g;
					if($surrounding_ref ne $surrounding_inference)
					{
						push(@differences, $graphLevel . ':' . $ref_href->{$graphLevel} . '!' . $inference_href->[$graphLevel]);
					}
				}
			}
			return \@differences;
		};

		my %haplotype_errors_byGraphLevel;
		my $minEditDistance;
		my $minEditDistance_ordering;
		my $minEditDistance_diffs;
		foreach my $ordering ([1, 2], [2, 1])
		{
			my $comparison_h1 = $haplotypeComparison->($truth_haplotypes_byGraphLevel[0], $inferred_haplotypes_byGraphLevel[$ordering->[0]-1]);
			my $comparison_h2 = $haplotypeComparison->($truth_haplotypes_byGraphLevel[1], $inferred_haplotypes_byGraphLevel[$ordering->[1]-1]);
			
			my $editDistance_cum = scalar(@$comparison_h1) + scalar(@$comparison_h2); 
			
			if(1 == 0)
			{
				print "\t\t\t\tHaplotype comparison ordering $ordering->[0], $ordering->[1]\n";				
				print "\tEdit distance  h1: ", scalar(@$comparison_h1), "\n";
				print "\tEdit distance  h2: ", scalar(@$comparison_h2), "\n";
				print "\tEdit distance cum.: ", $editDistance_cum, "\n";
			}
			if((not defined $minEditDistance) or ($editDistance_cum < $minEditDistance))
			{
				$minEditDistance = $editDistance_cum;
				$minEditDistance_ordering = [@$ordering];
				$minEditDistance_diffs = [(map {['h1', $_]} @$comparison_h1), (map {['h2', $_]} @$comparison_h2)];
			}
		}
		
		print "\t\tHaplotype comparison\n";
		print "\t\t\t", $minEditDistance, " alleles difference (in terms of graph levels)\n";
		if($minEditDistance)
		{
			print "\t\t\t", join(' - ', @$minEditDistance_ordering), "\n";
			print "\t\t\t", join('; ', map {join(';' , @$_)} @$minEditDistance_diffs), "\n";
			foreach my $errorInfo (@$minEditDistance_diffs)
			{
				die unless($errorInfo->[1] =~ /^(\d+):/);
				$haplotype_errors_byGraphLevel{$1}++;
			}
		}
		
		my %genotype_errors_byGraphLevel;
		my $genotype_differences = 0;
		my @genotype_differences_which;
		for(my $graphLevel = $minGraphLevel_bothHaplotypes; $graphLevel <= $maxGraphLevel_bothHaplotypes; $graphLevel++)
		{
			die unless(defined $truth_genotypes_byGraphLevel_href->{$graphLevel});
			die unless(defined $inferred_genotypes_byGraphLevel_aref->[$graphLevel]);
			my @alleles_truth = @{$truth_genotypes_byGraphLevel_href->{$graphLevel}};
			my @alleles_inference = @{$inferred_genotypes_byGraphLevel_aref->[$graphLevel]};
			die unless(scalar(@alleles_truth) == 2);
			die unless(scalar(@alleles_inference) == 2);
	
			my $surroundingLower = $graphLevel - 10; $surroundingLower = $minGraphLevel_bothHaplotypes if($surroundingLower < $minGraphLevel_bothHaplotypes);
			my $surroundingUpper = $graphLevel + 10; $surroundingUpper = $maxGraphLevel_bothHaplotypes if($surroundingUpper > $maxGraphLevel_bothHaplotypes);
			
			my @alleles_truth_surrounding = map {my $h = $_; join('', map {my $gL = $_; die $gL unless(defined $truth_genotypes_byGraphLevel_href->{$gL}); $truth_genotypes_byGraphLevel_href->{$gL}[$h]} ($surroundingLower .. $surroundingUpper))} (0, 1);
			my @alleles_inference_surrounding = map {my $h = $_; join('', map {$inferred_genotypes_byGraphLevel_aref->[$_][$h]} ($surroundingLower .. $surroundingUpper))} (0, 1);
			die unless(scalar(@alleles_truth_surrounding) == 2);
			die unless(scalar(@alleles_inference_surrounding) == 2);
			@alleles_truth_surrounding = map {$_ =~ s/[\-_]//g; $_} @alleles_truth_surrounding;
			@alleles_inference_surrounding = map {$_ =~ s/[\-_]//g; $_} @alleles_inference_surrounding;
	
			# if($graphLevel > ($minGraphLevel_bothHaplotypes + 20))
			# {
				# my $h = 0;
				# die join('', map {my $gL = $_; die $gL unless(defined $truth_genotypes_byGraphLevel_href->{$gL}); $truth_genotypes_byGraphLevel_href->{$gL}[$h]} ($surroundingLower, $surroundingUpper));
				# die Dumper($surroundingLower, $surroundingUpper, \@alleles_truth_surrounding, \@alleles_inference_surrounding);
			# }
			my $minDifference;
			foreach my $ordering ([1, 2], [2, 1])
			{
				my $diff_h1 = (($alleles_truth[0] eq $alleles_inference[$ordering->[0]-1]) ? 0 : 1);
				my $diff_h2 = (($alleles_truth[1] eq $alleles_inference[$ordering->[1]-1]) ? 0 : 1);
				my $diff_cum = $diff_h1 + $diff_h2;
				my $diff_h1_surrounding = (($alleles_truth_surrounding[0] eq $alleles_inference_surrounding[$ordering->[0]-1]) ? 0 : 1);
				my $diff_h2_surrounding = (($alleles_truth_surrounding[1] eq $alleles_inference_surrounding[$ordering->[1]-1]) ? 0 : 1);
				my $diff_cum_surrounding = $diff_h1_surrounding + $diff_h2_surrounding;
				if($diff_cum)
				{
					# warn Dumper($diff_cum, $ordering, [\@alleles_truth, \@alleles_inference], [\@alleles_truth_surrounding, \@alleles_inference_surrounding]);
				}
				$diff_cum = 0 if($diff_cum_surrounding == 0);
				if((not defined $minDifference) or ($diff_cum < $minDifference))
				{
					$minDifference = $diff_cum;
				}
			}
			$genotype_differences += $minDifference;
			if($minDifference != 0)
			{
				push(@genotype_differences_which, $graphLevel . ':(' . join('/', @alleles_truth) . ')!(' . join('/', @alleles_inference) . ')'); 
				$genotype_errors_byGraphLevel{$graphLevel}++;
			}
		}

		print "\t\tGenotype comparison:\n";			
		print "\t\t\t", $genotype_differences, " genotypic difference (in terms of graph levels)\n";
		if($genotype_differences)
		{
			print "\t\t\t", join('; ', @genotype_differences_which), "\n";
		}
		print "\t\tMutations:\n";				
		print "\t\t\t", join('; ', map {sort keys %{$graphLevels_with_mutation{$_}}} sort keys %graphLevels_with_mutation), "\n";		
		
		if(1 == 0)
		{
			# my @aligned_references = map {$truth{$sampleID}{$locus}{$_}{underlyingSequence_aligned_ref}} (1, 2);
			# my @aligned_preMutationSeq = map {$truth{$sampleID}{$locus}{$_}{underlyingSequence_aligned_preMutationSeq}} (1, 2);
			# my @aligned_postMutationSeq = map {$truth{$sampleID}{$locus}{$_}{underlyingSequence_aligned_postMutationSeq}} (1, 2);
			
			# print "aligned_references\n";
			# print "\t1: ", length($aligned_references[0]), "\n";
			# print "\t2: ", length($aligned_references[1]), "\n";
			# print "aligned_preMutationSeq\n";
			# print "\t1: ", length($aligned_preMutationSeq[0]), "\n";
			# print "\t2: ", length($aligned_preMutationSeq[1]), "\n";
			# print "aligned_postMutationSeq\n";
			# print "\t1: ", length($aligned_postMutationSeq[0]), "\n";
			# print "\t2: ", length($aligned_postMutationSeq[1]), "\n";
			
			my @inferred_haplotypes_noGaps = map {my $h = $_; $h =~ s/[\-_]//g; $h} @inferred_haplotypes;
			
			if(1 == 0)
			{
				print "\t$locus Now compare sequences of length "  .
					length($inferred_haplotypes_noGaps[0]) . ' / ' . length($inferred_haplotypes_noGaps[1]) . " (inference) and " . 
					length($truth_haplotypes_noGaps[0]) . ' / ' . length($truth_haplotypes_noGaps[1]) . " (truth).\n";
			}
						
			my $align_sequences =  sub {
				my $seq_ref = shift;
				my $seq_query = shift;
							
				die Dumper("Query shorter than reference, this is unexpected in this context", length($seq_ref), length($seq_query)) unless(length($seq_ref) <= length($seq_query));
				
				system("rm -rf ${tempFn_prefix}*");
				
				my $fn_ref =  $tempFn_prefix . '.ref';
				my $fn_query =  $tempFn_prefix . '.query';
				my $fn_output =  $tempFn_prefix . '.output';
				my $fn_stdout =  $tempFn_prefix . '.stdout';
			
				Util::writeFASTA($fn_ref, {contig => $seq_ref});									
				Util::writeFASTA($fn_query, {ref => $seq_query});
						
				my $cmd_align = qq(bash -c "perl $includes $globalAlignment_bin --use_minimap2 1 --reference $fn_query --query $fn_ref --output $fn_output --samtools_bin $samtools_bin --minimap2_bin $minimap2_bin 2>&1 > $fn_stdout");
				system($cmd_align) and die "Could not execute; '$cmd_align'";
													
				open(OUTPUT, '<', $fn_output) or die "Cannot open $fn_output";
				my $output_headerLine = <OUTPUT>;
				chomp($output_headerLine);
				die "Unexpected header line (I) in $fn_output: '$output_headerLine'" unless($output_headerLine =~ /^\d+ (\d+)-(\d+) \+(\d+)-(\d+)$/);
				
				my $alignment_firstQueryPos_0based = $1;
				my $alignment_lastQueryPos_0based = $2;
				my $alignment_firstRefPos_0based = $3;			
				my $alignment_lastRefPos_0based = $4;
				
				die "Unexpected header line (II) in $fn_output: '$output_headerLine'" unless($alignment_firstRefPos_0based == 0);
				#die "Unexpected header line (III) in $fn_output: '$output_headerLine'" unless($3 == 0);
				#die Dumper("Unexpected header line (IV) in $fn_output: '$output_headerLine'", $2, length($seq_ref)) unless($2 == (length($seq_ref)-1));
				die "Unexpected header line (V) in $fn_output: '$output_headerLine'" unless($alignment_lastRefPos_0based == (length($seq_ref)-1));

				my $alignment_ref = <OUTPUT>; chomp($alignment_ref);
				my $alignment_query = <OUTPUT>; chomp($alignment_query);
				close(OUTPUT);

				die unless(length($alignment_ref) == length($alignment_query));
				
				
				die "Cannot parse alignment_ref output '$alignment_ref'" unless($alignment_ref =~ /^(-*).+?(-*)$/);
				my $alignment_ref_gaps_beginning = $1;
				my $alignment_ref_gaps_end = $2;
				# print "Trimming " . length($alignment_gaps_beginning) . " gaps at the beginning and " . length($alignment_gaps_end) . " at the end.\n_mutations";
				
				$alignment_ref = substr($alignment_ref, length($alignment_ref_gaps_beginning), length($alignment_ref) - length($alignment_ref_gaps_beginning) - length($alignment_ref_gaps_end));
				$alignment_query = substr($alignment_query, length($alignment_ref_gaps_beginning), length($alignment_query) - length($alignment_ref_gaps_beginning) - length($alignment_ref_gaps_end));
				die unless(length($alignment_ref) == length($alignment_query));
				
				system("rm -rf ${tempFn_prefix}*");
				
				return [$alignment_ref, $alignment_query];
				
			};

			my $editDistance = sub {
				my $ref = shift;
				my $query = shift;
				die unless(length($ref) == length($query));
				my $editDistance = 0;
				my @diff;
				for(my $i = 0; $i < length($ref); $i++)
				{
					my $c_1 = substr($ref, $i, 1);
					my $c_2 = substr($query, $i, 1);
					if($c_1 ne $c_2)
					{
						$editDistance++;
						push(@diff, [$i, $c_1, $c_2]);
					}
				}
				return ($editDistance, \@diff);
			};
						
			my $minEditDistance;
			my $minEditDistance_ordering;
			my $minEditDistance_diffs;
			foreach my $ordering ([1, 2], [2, 1])
			{
				my $alignment_h1 = $align_sequences->($truth_haplotypes_noGaps[0], $inferred_haplotypes_noGaps[$ordering->[0]-1]);
				my $alignment_h2 = $align_sequences->($truth_haplotypes_noGaps[1], $inferred_haplotypes_noGaps[$ordering->[1]-1]);
				my ($editDistance_h1, $diff_h1_aref) = $editDistance->(@$alignment_h1);
				my ($editDistance_h2, $diff_h2_aref) = $editDistance->(@$alignment_h2);
				my $editDistance_cum = $editDistance_h1 + $editDistance_h2;
				
				if(1 == 0)
				{
					print "Haplotype comparison ordering $ordering->[0], $ordering->[1]\n";				
					print "\tEdit distance  h1: ", $editDistance_h1, "\n";
					print "\tEdit distance  h2: ", $editDistance_h2, "\n";
					print "\tEdit distance cum.: ", $editDistance_cum, "\n";
				}
				if((not defined $minEditDistance) or ($editDistance_cum < $minEditDistance))
				{
					$minEditDistance = $editDistance_cum;
					$minEditDistance_ordering = [@$ordering];
					$minEditDistance_diffs = [(map {['h1', @$_]} @$diff_h1_aref), (map {['h2', @$_]} @$diff_h2_aref)];
				}
			}
			
			my $fn_msa_in = $tempFn_prefix . '.forMSA.fa';
			my $fn_msa_out = $tempFn_prefix . '.fromrMSA.fa';
			Util::writeFASTA($fn_msa_in, {t1 => $truth_haplotypes_noGaps[0], t2 => $truth_haplotypes_noGaps[1], i1 => $inferred_haplotypes_noGaps[0], i2 => $inferred_haplotypes_noGaps[1]});
			my $mafft_cmd = qq(mafft $fn_msa_in > $fn_msa_out);
			system($mafft_cmd) and die "Could not execute: $mafft_cmd\n";
			
			print "\t\t Haplotype edit distance: ", $minEditDistance, "\n";
			print "\t\t\t", join(' - ', @$minEditDistance_ordering), "\n";
			print "\t\t\t", join('; ', map {join(';' , @$_)} @$minEditDistance_diffs), "\n";
			print "\t\t\t", "MSA: $fn_msa_out", "\n";
		}
		
		my @printSequences = (
			[truth_h1_preMutation => $truth_haplotypes_byGraphLevel_preMutation[0]],
			[truth_h2_preMutation => $truth_haplotypes_byGraphLevel_preMutation[1]],
			[truth_h1_postMutation => $truth_haplotypes_byGraphLevel[0]],
			[truth_h2_postMutation => $truth_haplotypes_byGraphLevel[1]],
			[inference_h1 => $thisGene_inferred_haplotypes_byGraphLevel_asHashRef[0]],
			[inference_h2 => $thisGene_inferred_haplotypes_byGraphLevel_asHashRef[1]],
		);
		
		my %maxLength_byGraphLevel;
		foreach my $printSeqTuple (@printSequences)
		{
			my $printSeqId = $printSeqTuple->[0];
			for(my $graphLevel = $minGraphLevel_bothHaplotypes; $graphLevel <= $maxGraphLevel_bothHaplotypes; $graphLevel++)
			{
				die $printSeqId unless(ref($printSeqTuple->[1]) eq 'HASH');
				die unless(defined $printSeqTuple->[1]{$graphLevel});
				my $thisLevel_length = length($printSeqTuple->[1]{$graphLevel});
				if((not exists $maxLength_byGraphLevel{$graphLevel}) or ($thisLevel_length > $maxLength_byGraphLevel{$graphLevel}))
				{
					$maxLength_byGraphLevel{$graphLevel} = $thisLevel_length;
				}
			}
		}
		
		my $lengthAlign = sub {
			my $s = shift;
			my $graphLevel = shift;
			my $fillCharacter = shift;
			
			die unless(defined $maxLength_byGraphLevel{$graphLevel});
			my $missing_characters = $maxLength_byGraphLevel{$graphLevel} - length($s);
			die unless($missing_characters >= 0);
			
			$fillCharacter = '-' if(not defined $fillCharacter);
			my $gapC = $fillCharacter x $missing_characters;
			die unless(length($gapC) == $missing_characters);
			return $s . $gapC;
		};
		
		open(OUTALN, '>', $tempFn_prefix) or die "Cannot open $tempFn_prefix";
		print OUTALN $minGraphLevel_bothHaplotypes, "-", $maxGraphLevel_bothHaplotypes, "\n";
		
		for(my $graphLevel = $minGraphLevel_bothHaplotypes; $graphLevel <= $maxGraphLevel_bothHaplotypes; $graphLevel++)
		{
			my $hasMutation = (exists $graphLevels_with_mutation{$graphLevel});
			my $printString = ($hasMutation) ? 'M' : ' ';
			print OUTALN $lengthAlign->($printString, $graphLevel);			
		}
		print OUTALN "\n";
		for(my $graphLevel = $minGraphLevel_bothHaplotypes; $graphLevel <= $maxGraphLevel_bothHaplotypes; $graphLevel++)
		{
			my $printString = ' ';
			$printString = 'h' if($haplotype_errors_byGraphLevel{$graphLevel});
			$printString = 'g' if($genotype_errors_byGraphLevel{$graphLevel});
			print OUTALN $lengthAlign->($printString, $graphLevel);			
		}	
		print OUTALN "\n";		
		foreach my $printSeqTuple (@printSequences)
		{
			my $printSeqId = $printSeqTuple->[0];
			print OUTALN '>', $printSeqId, "\n";
			for(my $graphLevel = $minGraphLevel_bothHaplotypes; $graphLevel <= $maxGraphLevel_bothHaplotypes; $graphLevel++)
			{
				print OUTALN $lengthAlign->($printSeqTuple->[1]{$graphLevel}, $graphLevel);
			}
			print OUTALN "\n";
		}		
		close(OUTALN);
		
		
		# $align_sequences->($truth_haplotypes_noGaps[0], $inferred_haplotypes_noGaps[0]);
	}
}
elsif($action eq 'applyAll')
{
	my @sampleIDs_and_BAMs;
	open(SAMPLES, '<', $outputFn_samples) or die "Cannot open $outputFn_samples";
	while(<SAMPLES>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split("\t", $line);
		die unless(scalar(@line_fields) == 3);
		push(@sampleIDs_and_BAMs, [$line_fields[0], $line_fields[1]]);
	}
	close(SAMPLES);
	
	for(my $sampleI = 0; $sampleI <= $#sampleIDs_and_BAMs; $sampleI++)
	{
		my $sampleID = $sampleIDs_and_BAMs[$sampleI][0];
		my $BAM = $sampleIDs_and_BAMs[$sampleI][1];
		my $results_file = $working_dir . '/' . $sampleID . '/hla/R1_bestguess.txt';
		die unless(-e $BAM);
		if(-e $results_file)
		{
			print "File $results_file existing, skip step 1..\n";
		}
		else
		{
			my $HLA_LA_cmd = qq(perl $HLA_LA_bin --BAM $BAM --graph PRG_MHC_GRCh38_withIMGT --sampleID $sampleID --maxThreads 7);
			if($action2 eq 'do')
			{
				print "Making inference for sample $sampleID ($sampleI / $#sampleIDs_and_BAMs)\n";
				print "\t", $HLA_LA_cmd, "\n";
				system($HLA_LA_cmd) and die "HLA*LA command '$HLA_LA_cmd' failed";
			}
			else
			{
				print $HLA_LA_cmd, "\n";
			}
		}
			
		my $HLA_LA_cmd_II = qq(perl $HLA_LA_bin --graph PRG_MHC_GRCh38_withIMGT --sampleID $sampleID --maxThreads 7 --action call2);
		if($action2 eq 'do')
		{
			print "Making inference (II) for sample $HLA_LA_cmd_II ($sampleI / $#sampleIDs_and_BAMs)\n";
			print "\t", $HLA_LA_cmd_II, "\n";
			system($HLA_LA_cmd_II) and warn "HLA*LA command (II) '$HLA_LA_cmd_II' failed";
		}
		else
		{
			print $HLA_LA_cmd_II, "\n";
		}			
		
	}
	
	print "\n\n.Action '$action' done.\n\n";	
}
elsif($action eq 'simulate')
{
	unless($n =~ /^\d+,\d+,\d+$/)
	{
		die "Please provide a simple comma-separated list for parameter --n (individuals with no mutations, with rate 1:5000, with rate 1:500";
	}
	 
	my $IMGT_genomic_alignments_href = Util::read_genomic_alignments($genomic_alignments_dir, 0);
	my $gen_fasta_href = {};
	foreach my $locus (sort keys %$IMGT_genomic_alignments_href)
	{
		next if($locus eq 'HLA-DOA'); # the sequence that remains after N- and gap-trimming is too short: "HLA-DOA extract from 1601 to 1700 (full length 3653)"
		my $extract_from;
		my $extract_to;
		foreach my $allele (sort keys %{$IMGT_genomic_alignments_href->{$locus}})
		{
			my $allele_sequence = uc($IMGT_genomic_alignments_href->{$locus}{$allele});
			unless($allele_sequence =~ /^[ACGTN_]+$/)
			{
				(my $allele_sequence_illegalCharacters = $allele_sequence) =~ s/[^ACGTN_]//g;
				die "Illegal characters in IMGT alignment: $allele_sequence_illegalCharacters";
			}
			die unless($allele_sequence =~ /^([N_]*)(.+?)([N_]*)$/);
			my $gaps_front = $1;
			my $gaps_end = $3;
			my $thisAllele_from = length($gaps_front);
			my $thisAllele_to = length($allele_sequence) - length($gaps_end) - length($gaps_front) - 1;
			if((not defined $extract_from) or ($thisAllele_from > $extract_from))
			{
				$extract_from = $thisAllele_from;
			}
			if((not defined $extract_to) or ($thisAllele_to < $extract_to))
			{
				$extract_to = $thisAllele_to;
			}
			# print join("\t", $allele, $gaps_front, $gaps_end), "\n";
		}
		die unless((defined $extract_from) and (defined $extract_to));
		
		if($locus eq 'HLA-DRB1')
		{
			$extract_to -= 160;
			warn "Descreasing DRB1 - extract_to by 160 (reverse-complement!)";
		}
		print "$locus extract from $extract_from to $extract_to (full length " . length((values %{$IMGT_genomic_alignments_href->{$locus}})[0]) . ")\n";
		
		my %illegalCharacters;
		foreach my $allele (sort keys %{$IMGT_genomic_alignments_href->{$locus}})
		{
			my $allele_sequence = uc($IMGT_genomic_alignments_href->{$locus}{$allele});
			my $allele_sequence_extracted = substr($allele_sequence, $extract_from, $extract_to - $extract_from + 1);
			$allele_sequence_extracted =~ s/_//g;
			unless($allele_sequence_extracted =~ /^[ACGT]+$/)
			{
				warn "Sequence of allele $allele contains illegal characters";	
				$illegalCharacters{$allele} = 1;
			}
			$gen_fasta_href->{$allele} = $allele_sequence_extracted;
		}	
		
		if(scalar(keys %illegalCharacters))
		{
			print "Illegal characters detected in " . scalar(keys %illegalCharacters) . " sequences of $locus, delete alleles...\n";
			foreach my $allele (keys %illegalCharacters)
			{
				delete ($gen_fasta_href->{$allele});
			}
		}
	}

	my $fn_sequences = $graph_dir . '/sequences.txt';
	die "Missing file $fn_sequences" unless(-e $fn_sequences); 

	my %knownAlleles_in_PRG;
	my %knownAlleles_2_locus;
	loadAllelesFromGraph($graph_dir, \%knownAlleles_in_PRG, \%knownAlleles_2_locus);

	my %pseudoGenomic_alignments;
	my %pseudoGenomic_refPos;
	my %pseudoGenomic_refAlleles;
	VCFFunctions::readPseudoGenomicSequences($graph_dir . '/pseudoGenomic_fullLengthMapping', \%pseudoGenomic_alignments, \%pseudoGenomic_refPos);
	
	foreach my $locus (keys %pseudoGenomic_alignments)
	{
		my @alleles = keys %{$pseudoGenomic_alignments{$locus}};
		my @refAlleles = grep {$_ =~ /\*ref/} @alleles;
		die unless(scalar(@refAlleles) == 1);
		$pseudoGenomic_refAlleles{$locus} = $refAlleles[0];
	}

	my %HLAtypes_PGF;
	VCFFunctions::readPGFAlleles('PGF_loci_and_alleles.txt', \%HLAtypes_PGF);
	# foreach my $locus (keys %HLAtypes_PGF)
	# {
		# print $locus, "\t", $HLAtypes_PGF{$locus}{PGFAllele}, "\n";
	# }

	my $PGF_sequence;
	if(-e '../pgf_GRCh38.fa')
	{
		my $pgf_href = Util::readFASTA('../pgf_GRCh38.fa', 1);
		die unless(scalar(values %$pgf_href) == 1);
		$PGF_sequence = (values %$pgf_href)[0];
	}
	else
	{
		die "No pgf_GRCh38.fa found - slow";
		$PGF_sequence = VCFFunctions::getPGFSequence($graph_dir);
	}

	my $PGF_sequence_rc = Util::reverseComplement($PGF_sequence);

	my %overlapStatsByGene;
	my %useAllelesForSimulations;
	my %useAllelesForSimulations_positionInPseudoGenomicAllele;
	# my $gen_fasta_href = Util::readFASTA($hla_gen, 0);

	my %IMGT_strand_relative_2_PGF;
	my $fn_strand = '../_fullLengthSimulations_strand.txt';
	if(-e $fn_strand)
	{
		print "Reading strand info from file $fn_strand\n";
		open(STRAND, '<', $fn_strand) or die "Cannot open $fn_strand";
		while(<STRAND>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line);
			die unless(scalar(@line_fields) == 2);
			$IMGT_strand_relative_2_PGF{$line_fields[0]} = $line_fields[1];
		}
		close(STRAND);	
	}
	else
	{
		print "Re-computing strand info...\n";
		my $alleleID_i = 0;
		foreach my $alleleID (sort keys %$gen_fasta_href)
		{
			$alleleID_i++;
			print "\r    ", $alleleID_i, " / ", scalar(keys %$gen_fasta_href);
			my @alleleID_parts = split(/ /, $alleleID);
			die unless(scalar(@alleleID_parts) == 4);
			my $coreAlleleID = $alleleID_parts[1];
			die unless($coreAlleleID =~ /^(\w+?)\*(.+)$/);
			my $locus = $1;
			my $allele_numeric = $2;
			
			$overlapStatsByGene{$locus} = [0, 0] unless(defined $overlapStatsByGene{$locus});
			
			if(not exists $IMGT_strand_relative_2_PGF{$locus})
			{
				my $strand;
				if(index($PGF_sequence, $gen_fasta_href->{$alleleID}) != -1)
				{
					$strand = '+';
				}
				elsif(index($PGF_sequence_rc, $gen_fasta_href->{$alleleID}) != -1)
				{
					die if(defined $strand);
					$strand = '-';
				}
				if(defined $strand)
				{
					if(exists $IMGT_strand_relative_2_PGF{$locus})
					{
						die unless($IMGT_strand_relative_2_PGF{$locus} eq $strand);
					}
					$IMGT_strand_relative_2_PGF{$locus} = $strand;
				}
			}
		}	
		
		open(STRAND, '>', $fn_strand) or die "Cannot open $fn_strand";
		foreach my $gene (sort keys %IMGT_strand_relative_2_PGF)
		{
			print STRAND join("\t", $gene, $IMGT_strand_relative_2_PGF{$gene}), "\n";
		}
		close(STRAND);
	}

	my $alleleID_i = 0;
	foreach my $alleleID (sort keys %$gen_fasta_href)
	{
		$alleleID_i++;
		my $coreAlleleID = $alleleID;
		die unless($coreAlleleID =~ /^(\w+?)\*(.+)$/);
		my $locus = $1;
		my $allele_numeric = $2;
		
		$overlapStatsByGene{$locus} = [0, 0] unless(defined $overlapStatsByGene{$locus});
		
		if(exists $knownAlleles_in_PRG{$coreAlleleID})
		{
			$overlapStatsByGene{$locus}[0]++;
			$useAllelesForSimulations{$locus}{$locus . '*' . $allele_numeric} = $gen_fasta_href->{$alleleID};
		}
		else
		{
			$overlapStatsByGene{$locus}[1]++;
		}
	}

	print "\n";

	foreach my $locus (keys %useAllelesForSimulations)
	{
		my $needRC = 0;
		if(exists $IMGT_strand_relative_2_PGF{$locus})
		{
			if($IMGT_strand_relative_2_PGF{$locus} eq '-')
			{
				$needRC = 1;
			}
		}
		else
		{
			warn "No strand info for gene $locus\n";
		}
		if($needRC)
		{
			my @alleles = sort keys %{$useAllelesForSimulations{$locus}};
			foreach my $allele (@alleles)
			{
				$useAllelesForSimulations{$locus}{$allele} = Util::reverseComplement($useAllelesForSimulations{$locus}{$allele});
			}
		}
	}

	my @simulateLoci;
	print "Overlap statistics (GeneID, Shared, NonShared):\n";
	my %locus_2_locusHla;
	my %locus_alignment_start;
	foreach my $locus (sort keys %overlapStatsByGene)
	{
		if(exists $pseudoGenomic_alignments{$locus})
		{
			$locus_2_locusHla{$locus} = $locus;
		}
		else
		{
			if(exists $pseudoGenomic_alignments{'HLA-' . $locus})
			{
				$locus_2_locusHla{$locus} = 'HLA-' . $locus;
			}		
			else
			{
				# warn "Missing locus '$locus' in pseudo-genomic sequences - skip";
				delete $useAllelesForSimulations{$locus};
				next;
				die Dumper("Missing locus '$locus' in pseudo-genomic sequences", [keys %pseudoGenomic_alignments]) 			
			}
		}
		
		print "\t - $locus $overlapStatsByGene{$locus}[0] $overlapStatsByGene{$locus}[1] ";
		if($overlapStatsByGene{$locus}[0])
		{
			my @availableAlleles = sort keys %{$useAllelesForSimulations{$locus}};
			my @alleleLengths = map {length($useAllelesForSimulations{$locus}{$_})} @availableAlleles;
			print "[" . min(@alleleLengths) . " - " . max(@alleleLengths) . " bp] ";
										
			my $translate_coordinate_in_gap_sequenced = sub {
				my $gapped_sequence = shift;
				my $position = shift;
				my $runningPos_noGap = -1;
				for(my $i = 0; $i < length($gapped_sequence); $i++)
				{
					my $c = substr($gapped_sequence, $i, 1);
					if(($c ne '_') and ($c ne '-'))
					{
						$runningPos_noGap++;
					}
					if($runningPos_noGap == $position)
					{
						return $i;
					}
				}
				die;
			};
				
			my $minAlleleLength = 0.75 * max(@alleleLengths);
			
			my %removeAlleles;
			my %alignmentStarts;
			foreach my $iteration (0, 1)
			{
				my $bestAlignmentStart;
				if($iteration == 1)
				{
					my @starts = sort {$alignmentStarts{$b} <=> $alignmentStarts{$a}} keys %alignmentStarts;
					$bestAlignmentStart = $starts[0];
				}
					
				foreach my $allele (@availableAlleles)
				{
					my $alleleLength = length($useAllelesForSimulations{$locus}{$allele});
					if($alleleLength >= $minAlleleLength)
					{
						my $pseudoGenomic_withGaps = $pseudoGenomic_alignments{$locus_2_locusHla{$locus}}{$allele};
						(my $pseudoGenomic_noGaps = $pseudoGenomic_withGaps) =~ s/[\-_]//g;
						my $IMGT_genomic_start_in_noGaps = index($pseudoGenomic_noGaps, $useAllelesForSimulations{$locus}{$allele});
						next if($IMGT_genomic_start_in_noGaps == -1);							
						my $IMGT_genomic_start_in_Gaps = $translate_coordinate_in_gap_sequenced->($pseudoGenomic_withGaps, $IMGT_genomic_start_in_noGaps);	
						if($iteration == 0)
						{
							$alignmentStarts{$IMGT_genomic_start_in_Gaps}++;
						}
						else
						{
							if($IMGT_genomic_start_in_Gaps != $bestAlignmentStart)
							{
								$removeAlleles{$allele} = 1;
							}							
						}
					}
				}				
			}
						
			my $survivingAlleles = 0;			
			foreach my $allele (@availableAlleles)
			{
				my $alleleLength = length($useAllelesForSimulations{$locus}{$allele});
				if(($alleleLength >= $minAlleleLength) and (not $removeAlleles{$allele}))
				{
					my $pseudoGenomic_withGaps = $pseudoGenomic_alignments{$locus_2_locusHla{$locus}}{$allele};
					(my $pseudoGenomic_noGaps = $pseudoGenomic_withGaps) =~ s/[\-_]//g;
					
					my $IMGT_genomic_start_in_noGaps = index($pseudoGenomic_noGaps, $useAllelesForSimulations{$locus}{$allele});
					unless($IMGT_genomic_start_in_noGaps != -1)
					{
						if($locus eq 'DPA1')
						{
							# warn "Cannot find the start position of $allele the IMGT sequence in the pseudo-genomic sequence - remove";
						} 					
						delete $useAllelesForSimulations{$locus}{$allele};
						next;
					}
					my $IMGT_genomic_lastPos_in_noGaps = $IMGT_genomic_start_in_noGaps + length($useAllelesForSimulations{$locus}{$allele}) - 1;
					
					$survivingAlleles++;
					die "Missing pseudo-genomic sequence for allele $allele" unless(exists $pseudoGenomic_alignments{$locus_2_locusHla{$locus}}{$allele});
					
					my $IMGT_genomic_start_in_Gaps = $translate_coordinate_in_gap_sequenced->($pseudoGenomic_withGaps, $IMGT_genomic_start_in_noGaps);
					my $IMGT_genomic_lastPos_in_Gaps = $translate_coordinate_in_gap_sequenced->($pseudoGenomic_withGaps, $IMGT_genomic_lastPos_in_noGaps);
					{
						# sanity check
						my $_prefix = substr($pseudoGenomic_withGaps, 0, $IMGT_genomic_start_in_Gaps);
						$_prefix =~ s/[\-_]//g;
						die unless($_prefix eq substr($pseudoGenomic_noGaps, 0, $IMGT_genomic_start_in_noGaps));
						
						my $_middle = substr($pseudoGenomic_withGaps, $IMGT_genomic_start_in_Gaps, $IMGT_genomic_lastPos_in_Gaps - $IMGT_genomic_start_in_Gaps + 1);
						$_middle =~ s/[\-_]//g;		
						die unless($_middle eq substr($pseudoGenomic_noGaps, $IMGT_genomic_start_in_noGaps, $IMGT_genomic_lastPos_in_noGaps - $IMGT_genomic_start_in_noGaps + 1));
						die unless($_middle eq $useAllelesForSimulations{$locus}{$allele});
					}
					
					$useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$allele}{raw} = [$IMGT_genomic_start_in_noGaps, $IMGT_genomic_lastPos_in_noGaps];
					$useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$allele}{alignment} = [$IMGT_genomic_start_in_Gaps, $IMGT_genomic_lastPos_in_Gaps];	
					 
					if($locus eq 'DRB1')
					{
						# my $alignment_extract = substr($pseudoGenomic_alignments{$locus_2_locusHla{$locus}}{$allele}, $IMGT_genomic_start_in_Gaps - 10, 20);
						# my $alignment_start = substr($pseudoGenomic_alignments{$locus_2_locusHla{$locus}}{$allele}, $IMGT_genomic_start_in_Gaps, 10);
						# print join("\t", $locus, $allele, $IMGT_genomic_start_in_Gaps, $IMGT_genomic_start_in_noGaps, $alignment_extract, $alignment_start, length($useAllelesForSimulations{$locus}{$allele})), "\n";
					}
					
					if(exists $locus_alignment_start{$locus})
					{
						unless($useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$allele}{alignment}[0] == $locus_alignment_start{$locus})
						{
							die Dumper("Locus alignment start mismatch", $locus, $allele, $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$allele}{alignment}[0], $locus_alignment_start{$locus});
						}
					}
					else
					{
						$locus_alignment_start{$locus} = $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$allele}{alignment}[0];
					}
					

				}
				else
				{
					delete $useAllelesForSimulations{$locus}{$allele};
				}
			}
			print " -> " . scalar(keys %{$useAllelesForSimulations{$locus}}) . " [>= $minAlleleLength bp - " . scalar(keys %removeAlleles) . " removed because of alignment start mismatch] ";
			if(scalar(keys %{$useAllelesForSimulations{$locus}}) > 0)
			{
				push(@simulateLoci, $locus);
			}
		}
		print "\n";
	}

	print "... simulate " . scalar(@simulateLoci) . " loci.\n\n";
	@simulateLoci = sort @simulateLoci;
	
	my $outputFn_HLATypes = $outputPrefix . '/' . $individualPrefix . '_HLATypes.txt';
	my $outputFn_HLADetails = $outputPrefix . '/' . $individualPrefix . '_HLAHaplotypeDetails.txt';
	foreach my $f ($outputFn_HLATypes, $outputFn_HLADetails, $outputFn_samples)
	{
		if(-e $f)
		{
			unlink($f) or die "Cannot unlink $f";
		}				
	}

	print "Simulating...\n";
	my @sampleIDs_and_BAMs;
	my %HLA_truth_types;
	my $globalIndivI = 0;
	my %locus_start_PGF;
	for(my $mutationRateI = 0; $mutationRateI <= $#mutationRates; $mutationRateI++)
	{
		my $indiv_n = $mutationRates_n[$mutationRateI];
		my $rate = $mutationRates[$mutationRateI];
		print "Generating $indiv_n individuals with mutation rate $rate .. \n";
		for(my $indivI = 0; $indivI < $indiv_n; $indivI++)
		{
			$globalIndivI++;
			my $sampleID = $individualPrefix . $globalIndivI;
			print "\t", $sampleID, "\n";
			my $outputDir = $outputPrefix . '/' . $sampleID;
			unless(-e $outputDir)
			{
				mkdir($outputDir) or die "Cannot mkdir $outputDir";
			}
			 
			my $fastq_combined_1 = $outputDir . '/fastq_1.fastq';
			my $fastq_combined_2 = $outputDir . '/fastq_2.fastq';
			my $BAM_combined = $outputDir . '/reads.bam';
			my $BAM_combined_primary = $outputDir . '/reads_primary.bam';
			foreach my $f ($fastq_combined_1, $fastq_combined_2, $BAM_combined)
			{
				if(-e $f)
				{
					unlink($f) or die "Cannot unlink $f";
				}				
			}
			
			my %fasta_for_simulation;
			foreach my $locus (@simulateLoci)
			{
				my @alleles = sort keys %{$useAllelesForSimulations{$locus}};
				for(my $hap = 0; $hap <= 1; $hap++)
				{
					my $selectedAllele = $alleles[int(rand(scalar(@alleles)))];
					$HLA_truth_types{$sampleID}{$locus}{'alleles'}{'h' . $hap} = $selectedAllele;
					$HLA_truth_types{$sampleID}{$locus}{'mutations'}{'h' . $hap} = {};
					
					my $alleleSequence = $useAllelesForSimulations{$locus}{$selectedAllele};
					die unless($alleleSequence);
					
					my $alignmentCoordinate_start = $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[0];
					
					die unless(exists $pseudoGenomic_alignments{$locus_2_locusHla{$locus}}{$selectedAllele});
					my $alleleSequence_aligned = 
						substr($pseudoGenomic_alignments{$locus_2_locusHla{$locus}}{$selectedAllele},
							   $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[0],
							   $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[1] - $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[0] + 1);

					(my $alleleSequence_aligned_noGaps = $alleleSequence_aligned) =~ s/[\-\_]//g;
					die unless($alleleSequence_aligned_noGaps eq $alleleSequence);
					
					my $refAllele = $pseudoGenomic_refAlleles{$locus_2_locusHla{$locus}}; die unless($refAllele);
					my $alleleRef_aligned = 
						substr($pseudoGenomic_alignments{$locus_2_locusHla{$locus}}{$refAllele},
							   $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[0],
							   $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[1] - $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[0] + 1);
					die unless(length($alleleRef_aligned) == length($alleleSequence_aligned));

					my $refAllele_offset_sequence = substr($pseudoGenomic_alignments{$locus_2_locusHla{$locus}}{$refAllele}, 0, $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[0]);
					(my $refAllele_offset_sequence_noGaps = $refAllele_offset_sequence) =~ s/[\-_]//g;
		
					my $alleleSequence_ref_aligned;
					my $alleleSequence_mutated_aligned;
					my $alleleSequence_original_aligned;
					my @alignedAlignmentCoordinates;
					
					if($rate > 0)
					{
						my $runningAllelePos_nonAligned = -1;
						for(my $i = 0; $i < length($alleleSequence_aligned); $i++)
						{
							my $oldAllele = substr($alleleSequence_aligned, $i, 1);
							my $oldAllele_ref = substr($alleleRef_aligned, $i, 1);
							my $isGap = (($oldAllele eq '-') or ($oldAllele eq '_'));
							die unless(length($oldAllele) == length($oldAllele_ref));
							
							$runningAllelePos_nonAligned++ if(not $isGap);
							
							if((rand(1) <= $rate) and not $isGap)
							{
								my $newAllele = generateMutation($oldAllele);
								die if($newAllele eq $oldAllele);
								if(length($newAllele) == length($oldAllele))
								{
									$alleleSequence_ref_aligned .= $oldAllele_ref;
									$alleleSequence_mutated_aligned .= $newAllele;
									$alleleSequence_original_aligned .= $oldAllele;	
									push(@alignedAlignmentCoordinates, $i);									
								}
								else
								{
									die unless(length($newAllele) == 2);
									die unless(substr($newAllele, 0, 1) eq $oldAllele);
									
									$alleleSequence_ref_aligned .= ($oldAllele_ref . '-');		
									$alleleSequence_mutated_aligned .= $newAllele;
									$alleleSequence_original_aligned .= ($oldAllele . '-');		
									push(@alignedAlignmentCoordinates, $i);									
									push(@alignedAlignmentCoordinates, -1);									

								}
								$HLA_truth_types{$sampleID}{$locus}{'mutations'}{'h' . $hap}{$runningAllelePos_nonAligned . ':' . $oldAllele . '->' . $newAllele}++;
							}
							else
							{
								$alleleSequence_ref_aligned .= $oldAllele_ref;
								$alleleSequence_mutated_aligned .= $oldAllele;
								$alleleSequence_original_aligned .= $oldAllele;
								push(@alignedAlignmentCoordinates, $i);
							}
						}
					}
					else
					{
						$alleleSequence_ref_aligned = $alleleRef_aligned;
						$alleleSequence_mutated_aligned = $alleleSequence_aligned;
						$alleleSequence_original_aligned = $alleleSequence_aligned;	
						@alignedAlignmentCoordinates = (0 .. (length($alleleSequence_aligned)-1));
					}
					@alignedAlignmentCoordinates = map {my $v = $_; $v += $alignmentCoordinate_start if($v != -1); $v} @alignedAlignmentCoordinates;
					
					die unless(length($alleleSequence_mutated_aligned) == length($alleleSequence_original_aligned));
					die unless(length($alleleSequence_ref_aligned) == length($alleleSequence_original_aligned));
					(my $alleleSequence_mutated_raw = $alleleSequence_mutated_aligned) =~ s/[-_]//g;
					
					$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{raw}{'h' . $hap} = $alleleSequence_mutated_raw;
					
					$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{reference}{'h' . $hap} = $alleleSequence_ref_aligned;
					$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{afterMutation}{'h' . $hap} = $alleleSequence_mutated_aligned;
					$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{beforeMutation}{'h' . $hap} = $alleleSequence_original_aligned;
					$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentCoordinates}{'h' . $hap} = \@alignedAlignmentCoordinates;

					$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInAlignedPseudogenomic}{'h' . $hap} = $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[0];
					$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInRawPseudogenomic}{'h' . $hap} = $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{raw}[0];
					$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInPGF}{'h' . $hap} = 'NA';

					
					die unless(exists $pseudoGenomic_refPos{$locus_2_locusHla{$locus}});
					
					if($pseudoGenomic_refPos{$locus_2_locusHla{$locus}}[0] eq 'pgf')
					{
						$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInPGF}{'h' . $hap} = 
							$pseudoGenomic_refPos{$locus_2_locusHla{$locus}}[1]+ 
							length($refAllele_offset_sequence_noGaps);
							
						if(exists $locus_start_PGF{$locus})
						{
							unless($HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInPGF}{'h' . $hap} == $locus_start_PGF{$locus})
							{
								die Dumper(
									"PGF start mismatch",
									$locus_start_PGF{$locus},
									$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInPGF}{'h' . $hap}
								); 
							}
							
						}
						else
						{
							$locus_start_PGF{$locus} = $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInPGF}{'h' . $hap};
						}
					}
					
	#				$useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$allele}{alignment} = $IMGT_genomic_start_in_Gaps;
					
					my $fasta_for_wgsim = $outputDir . '/' . $locus . '_h' . $hap . '.input.fa';
					my $fastq_from_wgsim_1 = $outputDir . '/' . $locus . '_h' . $hap . '.output_1.fastq';
					my $fastq_from_wgsim_2 = $outputDir . '/' . $locus . '_h' . $hap . '.output_2.fastq';
					
					die unless($alleleSequence_mutated_raw =~ /^[ACGTN]+$/i);
					Util::writeFASTA($fasta_for_wgsim, {'HLA' . $locus . '_h' . $hap . '_' . $selectedAllele => $alleleSequence_mutated_raw});
					
					my $want_read_pairs = (length($alleleSequence_mutated_raw) * $targetCoverage) / (2 * $readLength);
					my $cmd_wgsim = qq($wgsim_bin -e 0.01 -d 500 -N $want_read_pairs -1 $readLength -2 $readLength -r 0 -R 0 -X 0 $fasta_for_wgsim $fastq_from_wgsim_1 $fastq_from_wgsim_2);
					system($cmd_wgsim) and die "Could not execute: $cmd_wgsim";
					
					my $bases_read_1 = countBasesInFASTQ($fastq_from_wgsim_1);
					my $bases_read_2 = countBasesInFASTQ($fastq_from_wgsim_2);
					
					# die Dumper($locus, length($alleleSequence_mutated_raw), $bases_read_1, $bases_read_2, ($bases_read_1+$bases_read_2)/length($alleleSequence_mutated_raw));
					
					my $cmd_cat = qq(cat $fastq_from_wgsim_1 >> $fastq_combined_1 && cat $fastq_from_wgsim_2 >> $fastq_combined_2);
					system($cmd_cat) and die "Cat command '$cmd_cat' failed";
					
					# unlink($fasta_for_wgsim);
					unlink($fastq_from_wgsim_1);
					unlink($fastq_from_wgsim_2);
				}
				my $n_mutations = scalar(keys %{$HLA_truth_types{$sampleID}{$locus}{'mutations'}{'h1'}}) + scalar(keys %{$HLA_truth_types{$sampleID}{$locus}{'mutations'}{'h1'}});
				print "\t\t$locus $n_mutations mutations\n"
			}
		
			my $cmd_bwa = qq($bwa_bin mem $GRCh38_refGenome $fastq_combined_1 $fastq_combined_2 | $samtools_bin sort --reference $GRCh38_refGenome -O BAM - > $BAM_combined && samtools index $BAM_combined);
			system($cmd_bwa) and die "BWA mapping failed - command '$cmd_bwa'";
			
			my $cmd_bwa_II = qq($bwa_bin mem $GRCh38_refGenome_primary $fastq_combined_1 $fastq_combined_2 | $samtools_bin sort --reference $GRCh38_refGenome_primary -O BAM - > $BAM_combined_primary && samtools index $BAM_combined_primary);
			system($cmd_bwa_II) and die "BWA mapping failed - command '$cmd_bwa_II'";
			
			push(@sampleIDs_and_BAMs, [$sampleID, abs_path($BAM_combined), $rate]);
		}
	}


	open(HLATYPES, '>', $outputFn_HLATypes) or die "Cannot open $outputFn_HLATypes";
	open(HLADETAILS, '>', $outputFn_HLADetails) or die "Cannot open $outputFn_HLADetails";
	open(SAMPLES, '>', $outputFn_samples) or die "Cannot open $outputFn_samples";
	foreach my $sample (@sampleIDs_and_BAMs)
	{
		print SAMPLES join("\t", @$sample), "\n";
	}
	close(SAMPLES);
	
	print HLATYPES join("\t", 'IndividualID', @simulateLoci), "\n";
	print HLADETAILS join("\t",
		'IndividualID', 'Locus', 'Haplotype', 'Allele', 'Mutations',
		'rawUnderlyingSequence', 'underlyingSequence_aligned_ref', 'underlyingSequence_aligned_postMutationSeq', 'underlyingSequence_aligned_preMutationSeq', 'underlyingSequence_aligned_alignmentCoordinates',
		'alignmentStart_pseudoGenomic_aligned', 'alignmentStart_pseudoGenomic_raw', 'alignmentStart_PGF',
	), "\n";
	foreach my $sampleID (sort keys %HLA_truth_types)
	{
		my @outputFields_HLATypes = ($sampleID);
		foreach my $locus (@simulateLoci)
		{
			my @alleles = map {$HLA_truth_types{$sampleID}{$locus}{'alleles'}{'h' . $_}} (0, 1);
			push(@outputFields_HLATypes, join('/', @alleles));
			foreach my $haplotype (0, 1)
			{ 
				my @outputFields_Details = ($sampleID, $locus, $haplotype + 1);
				my $haplotype_id = 'h' . $haplotype;
				push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleles'}{$haplotype_id});
				push(@outputFields_Details, join('; ', sort keys %{$HLA_truth_types{$sampleID}{$locus}{'mutations'}{$haplotype_id}}));
				push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{'raw'}{$haplotype_id});
				push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'reference'}{$haplotype_id});
				push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'afterMutation'}{$haplotype_id});
				push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'beforeMutation'}{$haplotype_id});
				die unless(scalar(@{$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'alignmentCoordinates'}{$haplotype_id}}) == length( $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'beforeMutation'}{$haplotype_id}));
				push(@outputFields_Details, join(';', @{$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'alignmentCoordinates'}{$haplotype_id}}));
				push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'alignmentStartInAlignedPseudogenomic'}{$haplotype_id});
				push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'alignmentStartInRawPseudogenomic'}{$haplotype_id});
				push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'alignmentStartInPGF'}{$haplotype_id});
				print HLADETAILS join("\t", @outputFields_Details), "\n";
			}
		}
		print HLATYPES join("\t", @outputFields_HLATypes), "\n";
	}
	close(HLATYPES);
}
else
{
	die "Unknown --action: $action";
}

sub produceVCFFromTruth
{
	my $individualID = shift;
	my $locus = shift;
	my $truth_href = shift;
	my $outputFn = shift;
	
	foreach my $h (1, 2)
	{
		die unless(exists $truth_href->{$h});
		die unless(exists $truth_href->{$h}{underlyingSequence_aligned_ref});
		die unless(exists $truth_href->{$h}{underlyingSequence_aligned_preMutationSeq});
		die unless(exists $truth_href->{$h}{underlyingSequence_aligned_postMutationSeq});
	}

	my @aligned_references = map {$truth_href->{$_}{underlyingSequence_aligned_ref}} (1, 2);
	my @aligned_preMutationSeq = map {$truth_href->{$_}{underlyingSequence_aligned_preMutationSeq}} (1, 2);
	my @aligned_postMutationSeq = map {$truth_href->{$_}{underlyingSequence_aligned_postMutationSeq}} (1, 2);
	
	# print "aligned_references\n";
	# print "\t1: ", length($aligned_references[0]), "\n";
	# print "\t2: ", length($aligned_references[1]), "\n";
	# print "aligned_preMutationSeq\n";
	# print "\t1: ", length($aligned_preMutationSeq[0]), "\n";
	# print "\t2: ", length($aligned_preMutationSeq[1]), "\n";
	# print "aligned_postMutationSeq\n";
	# print "\t1: ", length($aligned_postMutationSeq[0]), "\n";
	# print "\t2: ", length($aligned_postMutationSeq[1]), "\n";
	
	my @aligned_references_noGaps = map {my $h = $_; $h =~ s/[\-_]//g; $h} @aligned_references;
	
	open(VCF, '>', $outputFn) or die "Cannot open $outputFn";
	print VCF "##fileformat=VCFv4.2", "\n";
	print VCF '#', join("\t", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", $individualID), "\n";
	
	my $runningAllele_h1;
	my $runningAllele_h2;
	my $runningAllele_h1_preMutation;
	my $runningAllele_h2_preMutation;
	my $runninReference_h1;
	my $runninReference_h2;
	
	my $runningAlleles_ref_start;
	my $runningAllele_h1_start;
	my $runningAllele_h2_start;
	my $runningPosition_h1_noGaps;
	my $runningPosition_h2_noGaps;
	my $runningPosition_h1_withGaps;
	my $runningPosition_h2_withGaps;

	$runninReference_h1 = substr($aligned_references[0], 0, 1);
	$runninReference_h2 = substr($aligned_references[1], 0, 1);		
	$runningAllele_h1_preMutation = substr($aligned_preMutationSeq[0], 0, 1);
	$runningAllele_h2_preMutation = substr($aligned_preMutationSeq[1], 0, 1);
	$runningAllele_h1 = substr($aligned_postMutationSeq[0], 0, 1);
	$runningAllele_h2 = substr($aligned_postMutationSeq[1], 0, 1);

	die unless($runninReference_h1 !~ /[\-_]/);
	die unless($runninReference_h2 !~ /[\-_]/);
	die unless($runninReference_h1 eq $runninReference_h2);
	die unless($runningAllele_h1 eq $runninReference_h1);
	die unless($runningAllele_h2 eq $runninReference_h2);
	
	$runningAlleles_ref_start = 0;
	$runningAllele_h1_start = 0;
	$runningAllele_h2_start = 0;
	$runningPosition_h1_noGaps = 0;
	$runningPosition_h2_noGaps = 0;
	$runningPosition_h1_withGaps = 0;
	$runningPosition_h2_withGaps = 0;
			
	while(($runningPosition_h1_withGaps < (length($aligned_postMutationSeq[0])-1)) and ($runningPosition_h2_withGaps < (length($aligned_postMutationSeq[1]))-1))
	{	
		my $add_h1 = 0;
		my $add_h2 = 0;
		
		if($runningPosition_h1_noGaps == $runningPosition_h2_noGaps)
		{
			$add_h1 = 1;
			$add_h2 = 1;
		}
		elsif($runningPosition_h1_noGaps < $runningPosition_h2_noGaps)
		{
			$add_h1 = 1;
		}
		else
		{
			die unless($runningPosition_h1_noGaps > $runningPosition_h2_noGaps);
			$add_h2 = 1;
		}
		
		$runningPosition_h1_withGaps += $add_h1;
		$runningPosition_h2_withGaps += $add_h2;
		
		die Dumper("Error", $runningPosition_h1_withGaps, length($aligned_references[0]), length($aligned_postMutationSeq[0])) unless($runningPosition_h1_withGaps < length($aligned_references[0]));
		die unless($runningPosition_h2_withGaps < length($aligned_references[1]));
		
		my $refC_h1 = substr($aligned_references[0], $runningPosition_h1_withGaps, 1);
		my $hapC_h1 = substr($aligned_postMutationSeq[0], $runningPosition_h1_withGaps, 1);
		my $hapC_h1_preMutation = substr($aligned_preMutationSeq[0], $runningPosition_h1_withGaps, 1);
		
		my $isGap_h1 = (($refC_h1 eq '-') or ($refC_h1 eq '_') or ($hapC_h1 eq '-') or ($hapC_h1 eq '_'));
		my $isGap_h1_preMutation = (($refC_h1 eq '-') or ($refC_h1 eq '_') or ($hapC_h1_preMutation eq '-') or ($hapC_h1_preMutation eq '_'));
		
		if(!(($refC_h1 eq '-') or ($refC_h1 eq '_')))
		{
			$runningPosition_h1_noGaps++ if($add_h1);
		}
		
		if($add_h1)
		{
			$runninReference_h1 .= $refC_h1;
			$runningAllele_h1 .= $hapC_h1;
			$runningAllele_h1_preMutation .= $hapC_h1_preMutation;
		}
		
		
		my $refC_h2 = substr($aligned_references[1], $runningPosition_h2_withGaps, 1);
		my $hapC_h2 = substr($aligned_postMutationSeq[1], $runningPosition_h2_withGaps, 1);
		my $hapC_h2_preMutation = substr($aligned_preMutationSeq[1], $runningPosition_h2_withGaps, 1);
		
		my $isGap_h2 = (($refC_h2 eq '-') or ($refC_h2 eq '_') or ($hapC_h2 eq '-') or ($hapC_h2 eq '_'));	
		my $isGap_h2_preMutation = (($refC_h2 eq '-') or ($refC_h2 eq '_') or ($hapC_h2_preMutation eq '-') or ($hapC_h2_preMutation eq '_'));	

		if(!(($refC_h2 eq '-') or ($refC_h2 eq '_')))
		{
			$runningPosition_h2_noGaps++ if($add_h2);
		}
		
		if($add_h2)
		{
			$runninReference_h2 .= $refC_h2;
			$runningAllele_h2 .= $hapC_h2;
			$runningAllele_h2_preMutation .= $hapC_h2_preMutation;
		}
		
		die unless($runningPosition_h1_noGaps < length($aligned_references_noGaps[0]));
		die unless($runningPosition_h2_noGaps < length($aligned_references_noGaps[0]));
		
		if(($runningPosition_h1_withGaps > 1) and ($hapC_h1 eq $refC_h1) and ($hapC_h2 eq $refC_h2) and ($hapC_h1_preMutation eq $refC_h1) and ($hapC_h2_preMutation eq $refC_h2) and (! $isGap_h1) and (! $isGap_h2) and (! $isGap_h1_preMutation) and (! $isGap_h2_preMutation) and ($runningPosition_h1_noGaps == $runningPosition_h2_noGaps))
		{
			(my $runninReference_h1_noGaps = $runninReference_h1) =~ s/[\-_]//g;
			(my $runninReference_h2_noGaps = $runninReference_h2) =~ s/[\-_]//g;
			die unless($runninReference_h1_noGaps eq $runninReference_h2_noGaps);
			
			# die unless(substr($runninReference_h1_noGaps, length($runninReference_h1_noGaps) - 1, 1) eq substr($runninReference_h2_noGaps, length($runninReference_h2_noGaps) - 1, 1));
			
			(my $runningAllele_h1_noGaps = $runningAllele_h1) =~ s/[\-_]//g;
			(my $runningAllele_h1_preMutation_noGaps = $runningAllele_h1_preMutation) =~ s/[\-_]//g;
			(my $runningAllele_h2_noGaps = $runningAllele_h2) =~ s/[\-_]//g;				
			(my $runningAllele_h2_preMutation_noGaps = $runningAllele_h2_preMutation) =~ s/[\-_]//g;				
			
			# print join("\t", $runningPosition_h1_withGaps, $runningPosition_h2_withGaps, $runningAllele_h1_noGaps, $runningAllele_h2_noGaps, $runningAllele_h1_preMutation_noGaps, $runningAllele_h2_preMutation_noGaps), "\n";
			die unless(substr($runningAllele_h1_noGaps, 0, 1) eq substr($runningAllele_h2_noGaps, 0, 1));
			die unless(substr($runningAllele_h1_preMutation_noGaps, 0, 1) eq substr($runningAllele_h2_preMutation_noGaps, 0, 1));
			die unless(substr($runninReference_h1_noGaps, 0, 1) eq substr($runningAllele_h1_noGaps, 0, 1));
			die unless(substr($runninReference_h1_noGaps, 0, 1) eq substr($runningAllele_h1_preMutation_noGaps, 0, 1));
			
			die unless(length($runningAllele_h1_noGaps) >= 2);
			die unless(length($runningAllele_h2_noGaps) >= 2);
			
			my $runningRef_variantPart = substr($runninReference_h1_noGaps, 0, length($runninReference_h1_noGaps) - 1);
			my $runningAllele_h1_noGaps_variantPart = substr($runningAllele_h1_noGaps, 0, length($runningAllele_h1_noGaps) - 1);
			my $runningAllele_h2_noGaps_variantPart = substr($runningAllele_h2_noGaps, 0, length($runningAllele_h2_noGaps) - 1);
			
			my $runningAllele_h1_preMutation_noGaps_variantPart = substr($runningAllele_h1_preMutation_noGaps, 0, length($runningAllele_h1_noGaps) - 1);
			my $runningAllele_h2_preMutation_noGaps_variantPart = substr($runningAllele_h2_preMutation_noGaps, 0, length($runningAllele_h2_noGaps) - 1);

			if(($runningRef_variantPart ne $runningAllele_h1_noGaps_variantPart) or ($runningRef_variantPart ne $runningAllele_h2_noGaps_variantPart))
			{
				my $novel_h1 = (($runningAllele_h1_noGaps_variantPart eq $runningAllele_h1_preMutation_noGaps_variantPart) ? 0 : 1);
				my $novel_h2 = (($runningAllele_h2_noGaps_variantPart eq $runningAllele_h2_preMutation_noGaps_variantPart) ? 0 : 1);
				
				my $positionForOutput = $runningAlleles_ref_start;
				
				if(length($runningAllele_h1_noGaps_variantPart) == length($runningAllele_h2_noGaps_variantPart))
				{
					$runningRef_variantPart = substr($runningRef_variantPart, 1);
					$runningAllele_h1_noGaps_variantPart = substr($runningAllele_h1_noGaps_variantPart, 1);
					$runningAllele_h2_noGaps_variantPart = substr($runningAllele_h2_noGaps_variantPart, 1);
					
					$runningAllele_h1_preMutation_noGaps_variantPart = substr($runningAllele_h1_preMutation_noGaps_variantPart, 1);
					$runningAllele_h2_preMutation_noGaps_variantPart = substr($runningAllele_h2_preMutation_noGaps_variantPart, 1);
					$positionForOutput++;
				}
				
				my %alleles = ($runningRef_variantPart => 0);
				my $allele_2_code = sub {
					my $allele = shift;							
					unless(exists $alleles{$allele}) 
					{
						$alleles{$allele} = (scalar(keys(%alleles))+1);
					}
					return $alleles{$allele};
				};
				$allele_2_code->($runningAllele_h1_noGaps_variantPart);
				$allele_2_code->($runningAllele_h2_noGaps_variantPart);
				
				my @alleles_sorted = sort {$alleles{$a} <=> $alleles{$b}} keys %alleles;
				die unless($alleles_sorted[0] eq $runningRef_variantPart);
				shift(@alleles_sorted);
				
				# print VCF '#', join("\t", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", $individualID), "\n";

				print VCF join("\t",
					$locus,
					$positionForOutput,
					".",
					$runningRef_variantPart,
					join(',', @alleles_sorted),
					".",
					".",
					"Novel:" . join('/', $novel_h1, $novel_h2),
					"GT",
					join('/', map {$allele_2_code->($_)} ($runningAllele_h1_noGaps_variantPart, $runningAllele_h2_noGaps_variantPart)),
				), "\n";
				
				# print join("\t", $positionForOutput, $runningRef_variantPart, $runningAllele_h1_noGaps_variantPart . '/' . $runningAllele_h2_noGaps_variantPart, $runningAllele_h1_preMutation_noGaps_variantPart . '/' . $runningAllele_h2_preMutation_noGaps_variantPart, $novel_h1 . '/' . $novel_h2), "\n";
			}
			
			# print Dumper($runninReference_h1, $runninReference_h2, $runningAllele_h1, $runningAllele_h2, $runningAllele_h1_preMutation, $runningAllele_h2_preMutation);
			

			$runninReference_h1 = substr($runninReference_h1, length($runninReference_h1) - 1, 1);
			$runninReference_h2 = substr($runninReference_h2, length($runninReference_h2) - 1, 1);
			$runningAllele_h1 = substr($runningAllele_h1, length($runningAllele_h1) - 1, 1);
			$runningAllele_h2 = substr($runningAllele_h2, length($runningAllele_h2) - 1, 1);
			$runningAllele_h1_preMutation = substr($runningAllele_h1_preMutation, length($runningAllele_h1_preMutation) - 1, 1);
			$runningAllele_h2_preMutation = substr($runningAllele_h2_preMutation, length($runningAllele_h2_preMutation) - 1, 1);
			
			$runningAlleles_ref_start = $runningPosition_h1_noGaps;
			
			die unless($runningAllele_h1 eq $runningAllele_h2);
			die unless($runningAllele_h1_preMutation eq $runningAllele_h2_preMutation);
			die unless($runninReference_h1 eq $runningAllele_h1);
			die Dumper($runninReference_h1, $runninReference_h2, $runningAllele_h1, $runningAllele_h2, $runningAllele_h1_preMutation, $runningAllele_h2_preMutation) unless($runninReference_h1 eq $runningAllele_h1_preMutation);
		}
	}
	
	close(VCF);
	print "Generated $outputFn", "\n";		
}

sub generateMutation
{
	my $nuc_in = shift;
	die unless($nuc_in =~ /^[ACGT]$/);
	my @available = qw/A C G T - I/;
	my $selectedAllele;
	do {
		$selectedAllele = $available[int(rand(scalar(@available)))];
		if($selectedAllele eq 'I')
		{
			my $selectedAllele_insertion = $available[int(rand(scalar(@available) - 2))];
			die unless($selectedAllele_insertion =~ /^[ACGT]$/);
			$selectedAllele = $nuc_in . $selectedAllele_insertion;
		}
	} while($selectedAllele eq $nuc_in);
	return $selectedAllele;
}

sub loadAllelesFromGraph
{
	my $dir_graph = shift;
	my $knownAlleles_href = shift;
	my $alleles_2_locus_href = shift;
		
	my @prg_files = grep {$_ =~ /(exon_\d+\.txt)|(intron\d+\.txt)/} glob("${dir_graph}/PRG/*.txt");
	foreach my $prg_file (@prg_files)
	{
		open(PRG, '<', $prg_file) or die "Cannot open $prg_file";
		my $prg_headerLine = <PRG>;
		while(<PRG>)
		{
			my $line = $_;
			if($line =~ /^(\S+) /)
			{
				my $allele = $1;
				if($allele =~ /^(\w+?)\*(.+)$/)
				{
					my $locus = $1;
					my $allele_numeric = $2;
					$knownAlleles_href->{$allele}++;
					$alleles_2_locus_href->{$allele} = $locus;
				}				
			}
		}
		close(PRG);
	}
}

sub countBasesInFASTQ
{
	my $file = shift;
	my $totalSeq = 0;
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		die unless(substr($line, 0, 1) eq '@');
		my $seq = <F>;
		chomp($seq);
		my $plus = <F>;
		chomp($plus);
		die unless(substr($plus, 0, 1) eq '+');
		my $qual = <F>;
		chomp($qual);
		die unless(length($seq) == length($qual));
		$totalSeq += length($seq);
	}
	close(F);
	return $totalSeq;
}


