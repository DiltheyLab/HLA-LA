#!/usr/bin/perl -w

# ./HLAtypeinference_validation.pl --sampleIDs SRR702070_HapmapExomes_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/SRR1000GGSK_HapMap_combined.txt
# ./HLAtypeinference_validation.pl --sampleIDs NA12878_Platinum_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
# ./HLAtypeinference_validation.pl --sampleIDs NA12891_Platinum_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
# ./HLAtypeinference_validation.pl --sampleIDs NA12892_Platinum_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
# ./HLAtypeinference_validation.pl --sampleIDs all_1000G --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
# ./HLAtypeinference_validation.pl --sampleIDs withA_NA18939_1000G_RED_FASTQ --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended

# ./HLAtypeinference_validation.pl --sampleIDs NA12878_Platinum_RED,NA12891_Platinum_RED,NA12892_Platinum_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended

BEGIN {
	use FindBin;
	push(@INC, $FindBin::Bin);
	push(@INC, $FindBin::RealBin);
}

use strict;
use List::MoreUtils qw/all mesh any /;
use List::Util qw/sum/;
use Data::Dumper;
use Getopt::Long;   
use Sys::Hostname;
use File::Copy;
use File::Basename;
use Storable;
use simpleHLA;
use File::Basename;

my $kMer_size = 55;  

							

# my @testCases = (
	# [[qw/A A/], [qw/A A/]],
	# [[qw/? A/], [qw/A A/]],
	# [[qw/A ?/], [qw/A A/]],
	# [[qw/A T/], [qw/A A/]],
	# [[qw/A A/], [qw/T A/]],
	# [[qw/A C/], [qw/G T/]],
	# [[qw/A C/], [qw/? T/]],
	# [[qw/? C/], [qw/G T/]],
	# [[qw/? ?/], [qw/G T/]],
	# [[qw/? T/], [qw/? T/]],
	# [[qw/? C/], [qw/? T/]],
	
# );
# foreach my $testCase (@testCases)
# {
	# print join(' vs ', map {join('/', @$_)} @$testCase), "   ", join(' ', compatibleStringAlleles($testCase->[0], $testCase->[1])), "\n";
# }
# exit;

# input parameters

my $sampleIDs = '';
my $trueHLA;
my $trueHaplotypes;
my $T = 0;
my $minCoverage = 0;
my $all_2_dig = 0;
my $only_4_dig = 1;
my $reduce_to_4_dig = 0;

my $fromPHLAT = 0;
my $fromHLAreporter = 0;
my $fromKourami = 0;


my @loci_for_check = qw/A B C DQA1 DQB1 DRB1/;

my $exon_folder = qq(../graphs/PRG_MHC_GRCh38_withIMGT/PRG/);
unless(-e $exon_folder)
{
	die "Please provide a kMerified graph -- exon folder not there!";
}

GetOptions (
 'sampleIDs:s' => \$sampleIDs, 
 'trueHLA:s' => \$trueHLA,
 'trueHaplotypes:s' => \$trueHaplotypes, 
 #'validation_round:s' => \$validation_round,
 'T:s' => \$T,
 'minCoverage:s' => \$minCoverage,
 'all_2_dig:s' => \$all_2_dig,
 'only_4_dig:s' => \$only_4_dig,
 'fromPHLAT:s' => \$fromPHLAT,
 'fromHLAreporter:s' => \$fromHLAreporter,
 'fromKourami:s' => \$fromKourami,
 'reduce_to_4_dig:s' => \$reduce_to_4_dig,
);         

die if($fromPHLAT and $fromHLAreporter);
die if($fromPHLAT and $fromKourami);
die if($fromHLAreporter and $fromKourami);

if($minCoverage)
{
	print "Minimum coverage threshold in place: $minCoverage\n";
}

if($sampleIDs =~ /all_(.+)/)
{
	my $searchString = $1;
	my @alternatives = glob('../working/*');
	@alternatives = grep {(-d $_) and ($_ =~ /$searchString/)} @alternatives;
	@alternatives = map {die unless($_ =~ /.+[\\\/](.+)/); $1} @alternatives;
	$sampleIDs = join(',', @alternatives);
}

my @sampleIDs = split(/,/, $sampleIDs);
if(@sampleIDs)
{
	foreach my $sampleID (@sampleIDs)
	{
		unless($sampleID =~ /^[\w]+$/)
		{
			die "Please provide only sample IDs with normal characters.";
		}
	}
}	

my $sample_IDs_abbr;
if($sampleIDs =~ /^allSimulations(_\w+)?/)
{
	my $addFilter = $1;
	my @dirs;
	if($addFilter)
	{
		@dirs = grep {$_ =~ /I\d+_simulations${addFilter}/} grep {-d $_} glob('../tmp/hla/*');
	}
	else
	{
		@dirs = grep {$_ =~ /I\d+_simulations/} grep {-d $_} glob('../tmp/hla/*');
	}
	
	@sampleIDs = map {die "Can't parse $_" unless($_ =~ /tmp\/hla\/(.+)/); $1} @dirs;
	
	if($sampleIDs =~ /^all_simulations_I(\d+)/i)
	{
		my $iteration = $1;
		@sampleIDs = grep {$_ =~ /^I${iteration}_/i} @sampleIDs;
	}
	
	my $debug = 1;
	if($debug)
	{
		@sampleIDs = grep {die unless($_ =~ /sample(\d+)$/); ($1 < 30)} @sampleIDs;
	}
	
	$sample_IDs_abbr = $sampleIDs;
}
elsif($sampleIDs =~ /^all/)
{
	my @dirs = grep {$_ !~ /simulations/} grep {-d $_} glob('../tmp/hla/*');
	@sampleIDs = map {die "Can't parse $_" unless($_ =~ /tmp\/hla\/(.+)/); $1} @dirs;
	
	if($sampleIDs =~ /^all_I(\d+)/i)
	{
		my $iteration = $1;
		@sampleIDs = grep {$_ =~ /^I${iteration}_/i} @sampleIDs;
	}
	else
	{
		die "Does this make sense?";
		@sampleIDs = grep {$_ =~ /^$sampleIDs/i} @sampleIDs;	
	}
	
	$sample_IDs_abbr = $sampleIDs;
}
else
{
	$sample_IDs_abbr = join('_', @sampleIDs);
	if(length($sample_IDs_abbr) > 50)
	{
		$sample_IDs_abbr = substr($sample_IDs_abbr, 0, 50);
	}
}

my $validation_round = 'R1';
die "Please specify --trueHLA for validation" unless($trueHLA);
		
# read reference dataset
my %reference_data;
open(REFERENCE, "<", $trueHLA) or die "Cannot open $trueHLA";
my $headerLine = <REFERENCE>;
chomp($headerLine);
$headerLine =~ s/\n//g;
$headerLine =~ s/\r//g;
my @header_fields = split(/[\t ]/, $headerLine);
@header_fields = map {if($_ =~ /HLAD((QA)|(QB)|(RB))$/){$_ .= '1';} $_} @header_fields;	
while(<REFERENCE>)
{
	my $line = $_;
	chomp($line);
	
	$line =~ s/\n//g;
	$line =~ s/\r//g;
	
	next unless($line);
	
	my @fields = split(/[\t ]/, $line);
	my %line = (mesh @header_fields, @fields);
	
	my $primary_key = $line{'IndividualID'};
	$reference_data{$primary_key} = \%line;
}
close(REFERENCE);

my %imputed_HLA;
my %imputed_HLA_Q;
my %imputed_HLA_avgCoverage;
my %imputed_HLA_lowCoverage;
my %imputed_HLA_minCoverage;

my %sample_noI_toI;

my $total_imputations = 0;

my %missing_reference_data;

mkdir('temp');
mkdir('temp/hla_validation');

my $summary_file = 'temp/summary_' . $sample_IDs_abbr . '.txt';	
open(SUMMARY, '>', $summary_file) or die "Cannot open $summary_file";
print SUMMARY $sample_IDs_abbr, "\n";
print SUMMARY "\t", join("\t", qw/Locus N CallRate Accuracy/), "\n";
foreach my $sampleID (@sampleIDs)
{
	my $sampleID_noI = $sampleID;
	$sampleID_noI =~ s/^I\d+_//g;
	
	
	my $bestGuess_file;
	
	if($fromPHLAT)
	{
		$bestGuess_file = '/gpfs1/well/gsk_hla/PHLAT/'.$sampleID.'/'.$validation_round.'_bestguess.txt';	
		unless(-e $bestGuess_file)
		{
			warn "Best-guess file $bestGuess_file not existing";
			next;
		}		
	}
	elsif($fromHLAreporter)
	{
		$bestGuess_file = '/gpfs1/well/gsk_hla/HLAreporter/results/'.$sampleID.'/'.$validation_round.'_bestguess.txt';	
		unless(-e $bestGuess_file)
		{
			warn "Best-guess file $bestGuess_file not existing";
			next;
		}			
	}
	elsif($fromKourami) 
	{
		$bestGuess_file = '/data/projects/phillippy/projects/MHC/HLA-PRG-LA/fromKourami/' . $sampleID . '/KouramiHLA.txt';
		unless(-e $bestGuess_file)
		{
			warn "Best-guess file $bestGuess_file not existing";
			next;
		}		 
	}
	else
	{
		$bestGuess_file = '../working/'.$sampleID.'/hla/'.$validation_round.'_bestguess.txt';
		unless(-e $bestGuess_file)
		{
			warn "Best-guess file $bestGuess_file not existing";
			next;
		}		
	}
		  
	open(BESTGUESS, '<', $bestGuess_file) or die "Cannot open $bestGuess_file";
	my $bestguess_header_line = <BESTGUESS>;
	chomp($bestguess_header_line);
	my @bestguess_header_files = split(/\t/, $bestguess_header_line);
	while(<BESTGUESS>)
	{
		my $line = $_;
		chomp($line);
		my @line_fields = split(/\t/, $line);
		my %line_hash = (mesh @bestguess_header_files, @line_fields);
		
		my $Q = $line_hash{'Q1'};
		$imputed_HLA_Q{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = $Q;
		if($Q < $T)
		{
			$imputed_HLA{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = '??:??';			
		}
		else
		{
			die unless(defined $line_hash{'Allele'});
			$imputed_HLA{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = $line_hash{'Allele'};			
		}
		
		$total_imputations++;
		
		
		if($line_hash{'Chromosome'} eq '1')
		{
			$imputed_HLA_avgCoverage{$line_hash{'Locus'}}{$sampleID_noI} = $line_hash{'AverageCoverage'};
			$imputed_HLA_lowCoverage{$line_hash{'Locus'}}{$sampleID_noI} = $line_hash{'CoverageFirstDecile'};
			$imputed_HLA_minCoverage{$line_hash{'Locus'}}{$sampleID_noI} = $line_hash{'MinimumCoverage'};
			
			if($minCoverage and ($line_hash{'MinimumCoverage'} < $minCoverage))
			{
				$imputed_HLA{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = '??:??';								
			}
		}
	}	
	close(BESTGUESS);
	
	die if(exists $sample_noI_toI{$sampleID_noI});
	$sample_noI_toI{$sampleID_noI} = $sampleID;
}

print "\nTotal imputations (for comparisons): ", $total_imputations, "\n";
	
my $debug = 0;
my $comparisons = 0;
my $compare_problems = 0;
my %locus_avgCoverages;
my %locus_lowCoverages;
my %locus_minCoverages;
my @allLoci_allIndivs_avgCoverage;

my %problem_locus_detail;
my %problem_locus_examined;
my %problem_haplo_counter;
my %problem_haplo_detail;
my %imputations_predictions;
my %reference_predictions;
my %imputed_HLA_Calls;
my %quality_measures; # not used
my $pileup_href = {};

# die Dumper(\%imputed_HLA);

my @loci = sort keys %imputed_HLA;

my $process_quality_measures = sub {};

my $PP_to_basket = sub {
	my $PP = shift;
	die unless(($PP >= 0) && ($PP <= 1));
	my $basket = int($PP * 10);
	$basket = 9 if($basket == 10);
	return $basket;			
};

my %types_as_validated;
foreach my $locus (sort @loci)
{
	my $arbitraty_indiv = (keys %reference_data)[0];
	my $reference_data_prefix;
	if(defined $reference_data{$arbitraty_indiv}{'HLA'.$locus})
	{
		$reference_data_prefix = 'HLA';
	}
	elsif($reference_data{$arbitraty_indiv}{$locus})
	{
		$reference_data_prefix = '';
	}
	else
	{
		next;
	}
	
	
	my %calibration_baskets;
	my %coverage_over_samples;
	my %coverage_over_samples_individualValues;
	my $coverage_over_samples_nSamples = 0;
	
	my $add_to_calibration_basket = sub {
		my $str_correct = shift;
		my $PP = shift;
		my $weight = shift;
		
		die unless(($str_correct eq 'correct') or ($str_correct eq 'incorrect'));
		die unless(defined $PP);
		die unless(defined $weight);
		
		push(@{$calibration_baskets{$PP_to_basket->($PP)}{$str_correct}}, {PP => $PP, weight => $weight});
	};	
	
	$problem_locus_examined{$locus} = 0;
	$problem_locus_detail{$locus} = 0;
	my @indivIDs = sort keys %{$imputed_HLA{$locus}};
	
	my $indivI_processedForEvaluation = -1;
	INDIV: foreach my $indivID (@indivIDs)
	{	
		$| = 1;
		
		# print "\t", $indivID, "\t", $locus, "\n";
					
		$debug = 0;
		
		my @imputed_hla_values = map { $imputed_HLA{$locus}{$indivID}{$_} } keys %{$imputed_HLA{$locus}{$indivID}};
		
		if(grep {simpleHLA::autoHLA_is2digit($_)} @imputed_hla_values)
		{
			die "Warning: 2-digit alleles detected in the inference set\n" . Dumper(\@imputed_hla_values); # 2-digit test
		}
		
		my @imputed_hla_values_q = map { my $r = $imputed_HLA_Q{$locus}{$indivID}{$_}; die unless(defined $r); $r } keys %{$imputed_HLA{$locus}{$indivID}};
		
		die "Undefined HLA ".join(', ', @imputed_hla_values) unless(scalar(grep {defined $_} @imputed_hla_values) == scalar(@imputed_hla_values));
				
		my @reference_hla_values;
		
		next INDIV unless($#imputed_hla_values == 1);
		
		my $reference_lookup_ID = $indivID;
		if($reference_lookup_ID =~ /^withA_/)
		{
			$reference_lookup_ID =~ s/withA_//;			
		}					
		if($reference_lookup_ID =~ /_FASTQ$/)
		{
			$reference_lookup_ID =~ s/_FASTQ$//;			
		}					
		if($reference_lookup_ID =~ /^downsample_/)
		{
			$reference_lookup_ID =~ s/^downsample_(I\d+_)?//;
			$reference_lookup_ID =~ s/_DSC\d+_\d+//;				
		}
		if($reference_lookup_ID =~ /^C_Platinum_/)
		{
			$reference_lookup_ID =~ s/C_Platinum_//;			
		}		
		if($reference_lookup_ID =~ /^Platinum_/)
		{
			$reference_lookup_ID =~ s/Platinum_//;			
		}			
		if($reference_lookup_ID =~ /_1000G/i)
		{
			$reference_lookup_ID =~ s/_1000G//i;			
		}				
		if($reference_lookup_ID =~ /1000G_/i)
		{
			$reference_lookup_ID =~ s/1000G_//i;			
		}
		if($reference_lookup_ID =~ /_PLATINUM/i)
		{
			$reference_lookup_ID =~ s/_PLATINUM//i;			
		}					
		
		if($reference_lookup_ID =~ /_HapmapExomes/i)
		{
			$reference_lookup_ID =~ s/_HapmapExomes//i;			
		}		
		
		if($reference_lookup_ID =~ /_RED/i)
		{
			$reference_lookup_ID =~ s/_RED//i;			
		}					

		if($reference_lookup_ID =~ /_PacBio/i)
		{
			$reference_lookup_ID =~ s/_PacBio//i;			
		}			
		if($reference_lookup_ID =~ /_Nanopore/i)
		{
			$reference_lookup_ID =~ s/_Nanopore//i;			
		}		
		unless(exists $reference_data{$reference_lookup_ID})
		{
			$missing_reference_data{$reference_lookup_ID}{$locus}++;
			# warn "No reference data for $locus $indivID";
		}
		next INDIV unless(exists $reference_data{$reference_lookup_ID});
		
		$reference_data{$reference_lookup_ID}{$reference_data_prefix.$locus} or die;
		
		if($reference_data{$reference_lookup_ID}{$reference_data_prefix.$locus})
		{
			@reference_hla_values = split(/\//, $reference_data{$reference_lookup_ID}{$reference_data_prefix.$locus});
		}
		
		$types_as_validated{$reference_lookup_ID}{$locus} = \@reference_hla_values;

		die Dumper($reference_data{$reference_lookup_ID}, \@reference_hla_values) unless($#reference_hla_values == 1);				
					
		die "Undefined HLA ".join(', ', @reference_hla_values) unless(scalar(grep {defined $_} @reference_hla_values) == scalar(@reference_hla_values));			
		@reference_hla_values = grep {! &simpleHLA::modernHLA_is_missing($_)} @reference_hla_values;

		if($#reference_hla_values == -1)
		{
			$missing_reference_data{$reference_lookup_ID}{$locus}++;
			next;
		}
					
		if($only_4_dig)
		{
			next unless (&simpleHLA::HLA_is4digit($reference_hla_values[0]) and (($#reference_hla_values == 0) || (&simpleHLA::HLA_is4digit($reference_hla_values[1]))));
		}
		
		if($reduce_to_4_dig)
		{
			my @reference_hla_values_before = @reference_hla_values;
			@reference_hla_values = map {
				my @inputAlleles = split(/;/, $_);
				join(';', map {simpleHLA::HLA_4digit($_)} @inputAlleles)		
			} @reference_hla_values;
			
			print "Before:\n\t".join(' / ', @reference_hla_values_before), "\n";
			print "After:\n\t".join(' / ', @reference_hla_values), "\n\n";
		}
				
		$imputed_HLA_Calls{$locus}{sum} += scalar(@imputed_hla_values);		
		$indivI_processedForEvaluation++;
		
		my @imputed_present = map {(! &simpleHLA::is_missing($_)) ? 1 : 0} @imputed_hla_values;
		my @imputed_hla_values_q_new;
		for(my $i = 0; $i <= $#imputed_hla_values; $i++)
		{
			my $Q = $imputed_hla_values_q[$i];
			die unless(defined $Q);
			die unless(defined $imputed_present[$i]);
			if($imputed_present[$i])
			{
				push(@imputed_hla_values_q_new, $Q);
			}
		}
		
		@imputed_hla_values = grep {! &simpleHLA::is_missing($_)} @imputed_hla_values;
		@imputed_hla_values_q = @imputed_hla_values_q_new;
		die unless($#imputed_hla_values == $#imputed_hla_values_q);
		
		$imputed_HLA_Calls{$locus}{called} += scalar(@imputed_hla_values);
		
		if($indivID =~ /wtsi/i)
		{	
			my @imp_before = @imputed_hla_values;
			@imputed_hla_values = map {my @alleles = split(/;/, $_); join(';', map {($_ =~ /[GNLQ]$/) ? ($_) : ($_, $_ . ':01', $_ . ':01:01', )} @alleles)} @imputed_hla_values;
			#die Dumper(\@imp_before, \@imputed_hla_values);
		}
		if($all_2_dig)
		{
			@reference_hla_values = map {join(';', map {&simpleHLA::autoHLA_2digit($_)} split(/;/, $_))} @reference_hla_values;
			@imputed_hla_values = map {join(';', map {&simpleHLA::autoHLA_2digit($_)} split(/;/, $_))} @imputed_hla_values;
		}
				
		if($locus eq 'B')
		{
		#	print Dumper($indivID, \@reference_hla_values, \@imputed_hla_values), "\n";
		}
		my $comparisons_before = $comparisons;
		my $problem_locus_detail_before = $problem_locus_detail{$locus};
		
		if($#imputed_hla_values == 1)
		{
			if($#reference_hla_values > -1)
			{
				if($#reference_hla_values == 0)
				{
					if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]) or &compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]))
					{
						$comparisons++;
						$problem_locus_examined{$locus}++;
						
						if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]))
						{
							$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;
							
							$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);														
						}
						elsif(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]))
						{
							$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;

							$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);								
						}
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 1, 1);
						
					}	
					else
					{
						$comparisons++;
						$compare_problems++;
						$problem_haplo_counter{$indivID}++;
						$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (1 problem of 1)';
						$problem_locus_detail{$locus}++;
						$problem_locus_examined{$locus}++;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect} += 0.5;
						$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect} += 0.5;
						
						$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 0.5);
						$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[1], 0.5);
						
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 1);
					}
				}
				elsif($#reference_hla_values == 1)
				{
					if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]) and &compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[1]))
					{
						$comparisons += 2;
						$problem_locus_examined{$locus} += 2;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;	
						$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;	
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 2, 2);	
						
						$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
						$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);
						
					}		
					elsif(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]) and &compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
					{
						$comparisons += 2;
						$problem_locus_examined{$locus} += 2;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
						$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;			
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 2, 2);		

						$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
						$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);
						
					}
					else
					{
						if(
							&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]) or &compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[1]) or
							&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]) or &compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1])
						)
						{
							$comparisons += 2;
							$compare_problems++;
							$problem_haplo_counter{$indivID}++;
							$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (1 problem of 2)';						
							$problem_locus_detail{$locus}++;
							$problem_locus_examined{$locus} += 2;
							
							if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]))
							{
								$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
								$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect}++;
								$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect}++;		

								$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
								$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[1], 1);									
							}
							elsif(&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[1]))
							{
								$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;					
								$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;		
								
								$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);
								$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);	
							}
							elsif(&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
							{
								$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
								$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
								$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect}++;	
								
								$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
								$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[1], 1);	
							}
							elsif(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]))
							{
								$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;					
								$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;		

								$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);
								$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);										
							}
							
							$process_quality_measures->($locus, $quality_measures{$indivID}, 1, 2);
						}
						else
						{
							$comparisons += 2;
							$compare_problems += 2;
							$problem_haplo_counter{$indivID} += 2;
							$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (2 problems of 2)';						
							$problem_locus_detail{$locus} += 2;
							$problem_locus_examined{$locus} += 2;
							
							$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;					
							$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect}++;
							$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect}++;	
							
							$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);
							$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[1], 1);									
							
							$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 2);						
						}
					}																
				}
				else
				{
					die;
				}
			}
		}
		elsif($#imputed_hla_values == 0)
		{
			if($#reference_hla_values > -1)
			{
				if($#reference_hla_values == 0)
				{
					if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]))
					{
						$comparisons++;
						$problem_locus_examined{$locus}++;
						$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;
						
						$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
							
						$process_quality_measures->($locus, $quality_measures{$indivID}, 1, 1);
						
					}	
					else
					{
						$comparisons++;
						$compare_problems++;
						$problem_haplo_counter{$indivID}++;
						$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (1 problem of 1)';
						$problem_locus_detail{$locus}++;
						$problem_locus_examined{$locus}++;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;
						
						$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 1);									
					}					
				}
				elsif($#reference_hla_values == 1)
				{
					if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]) or &compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
					{
						$comparisons += 1;
						$problem_locus_examined{$locus} += 1;
						
						if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]))
						{
							$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;

							$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
							
						}
						elsif(&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
						{
							$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;	

							$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);								
						}
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 1, 1);

					}		
					else
					{
						$comparisons++;
						$compare_problems++;
						$problem_haplo_counter{$indivID}++;
						$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (1 problem of 1)';
						$problem_locus_detail{$locus}++;
						$problem_locus_examined{$locus}++;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect} += 0.5;
						$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect} += 0.5;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;	
						
						$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);

						$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 1);								
					}
				}
				else
				{
					die;
				}
			}		
		}
		
		my $thisIndiv_problems = $problem_locus_detail{$locus} - $problem_locus_detail_before;
			
		my $avgCoverage = $imputed_HLA_avgCoverage{$locus}{$indivID};
		my $lowCoverage = $imputed_HLA_lowCoverage{$locus}{$indivID};
		my $minCoverage = $imputed_HLA_minCoverage{$locus}{$indivID};

		# average coverages
		if(($thisIndiv_problems > 0))
		{
			push(@{$locus_avgCoverages{$locus}{problems}}, $avgCoverage);
			push(@{$locus_lowCoverages{$locus}{problems}}, $lowCoverage);
			push(@{$locus_minCoverages{$locus}{problems}}, $minCoverage);
		}
		else
		{
			push(@{$locus_avgCoverages{$locus}{ok}}, $avgCoverage);
			push(@{$locus_lowCoverages{$locus}{ok}}, $lowCoverage);
			push(@{$locus_minCoverages{$locus}{ok}}, $minCoverage);	
		}
		
		push (@allLoci_allIndivs_avgCoverage, $avgCoverage);
		
		# print "\t", $thisIndiv_problems, "\n";
		
		# just for debugging - deactivated
		if($thisIndiv_problems == 0)
		{
			if($locus eq 'A')
			{
				# print join(' vs ', join('/', @reference_hla_values), join('/', @imputed_hla_values)), "\n";
			}	
		}
		
		my $thisIndiv_comparions = $comparisons - $comparisons_before;
		my $thisIndiv_OK = $thisIndiv_comparions - $thisIndiv_problems;
		
		my $indivID_withI = $sample_noI_toI{$indivID};
		die unless(defined $indivID_withI);				
		my $pileup_file = qq(../working/$indivID_withI/hla/${validation_round}_pileup_${locus}.txt);
			
		my $coverages_href = load_coverages_from_pileup($pileup_file);

		my @k_coverages_existing;
		if($indivI_processedForEvaluation > 0)
		{
			foreach my $exon (keys %coverage_over_samples)
			{
				foreach my $exonPos (keys %{$coverage_over_samples{$exon}})
				{
					push(@k_coverages_existing, $exon . '-/-' . $exonPos);
				}
			}					
		}
		

		my @k_coverages_thisSample;
		foreach my $exon (keys %$coverages_href)
		{
			foreach my $exonPos (keys %{$coverages_href->{$exon}})
			{
				$coverage_over_samples{$exon}{$exonPos} += $coverages_href->{$exon}{$exonPos};
				push(@{$coverage_over_samples_individualValues{$exon}{$exonPos}}, $coverages_href->{$exon}{$exonPos});
				push(@k_coverages_thisSample, $exon . '-/-' . $exonPos);
			}
		}	
		$coverage_over_samples_nSamples++;
		
		if($indivI_processedForEvaluation > 0)
		{
			my ($n_shared, $aref_l1_excl, $aref_l2_excl) = list_comparison(\@k_coverages_existing, \@k_coverages_thisSample); 
			if(($#{$aref_l1_excl} == -1) and ($#{$aref_l2_excl} == -1))
			{	
				die unless($n_shared == scalar(@k_coverages_existing));
				die unless($n_shared == scalar(@k_coverages_thisSample));					
			}
			else
			{
				warn Dumper("There is a problem with per-exon coverage numbers.", $indivI_processedForEvaluation, scalar(@k_coverages_existing), scalar(@k_coverages_thisSample), scalar(@$aref_l1_excl), scalar(@$aref_l2_excl));
			}
		}

		if($thisIndiv_problems)
		{
			print join("\t", $indivID, $locus, $thisIndiv_problems . ' problems',  "Reference " . $reference_data{$reference_lookup_ID}{$reference_data_prefix.$locus}, "Inference " . join('/', @imputed_hla_values)), "\n";
		}
		
		if(($thisIndiv_problems > 0) and (not $all_2_dig) and not ($fromPHLAT) and not($fromHLAreporter) and not($fromKourami) and (not ($indivID) =~ /wtsi/i))
		{
			my %readIDs;
			
			unless(-e $pileup_file)
			{
				die "Can't find pileup file $pileup_file";
			}
			
			load_pileup($pileup_href, $pileup_file, $indivID_withI);

		
			my $output_fn = 'temp/hla_validation/pileup_'.$validation_round.'_'.$indivID_withI.'_'.$locus.'.txt';
			open(my $output_fh, '>', $output_fn) or die "Cannot open $output_fn";
			print $output_fh join("\t", $indivID_withI, $locus, $thisIndiv_OK), "\n";
			
			unless(scalar(@imputed_hla_values) == 2)
			{
				warn "Can't produce pileup for $locus / $indivID";
				next;
			}
			my $inferred = \@imputed_hla_values;
			
			my $truth = \@reference_hla_values;
			
			if(scalar( grep {$_ =~ /XX/} @$truth) > 0)
			{
				warn "Can't produce pileup for $locus / $indivID";			
				next;
			}
			print {$output_fh} $inferred->[0], "\t\t\t", $truth->[0], "\n";
			print {$output_fh} $inferred->[1], "\t\t\t", $truth->[1], "\n\n";
			
			print {$output_fh} join("\t", "Inferred", "", "", "Apparently true", ""), "\n";
			
			my @exons = print_which_exons($locus);				
			my @inferred_trimmed = twoClusterAlleles($inferred);
			
			foreach my $exon (@exons)
			{
				my $file = find_exon_file($locus, $exon);
				my $sequences = read_exon_sequences($file);
									
				my @validated_extended = twoValidationAlleles_2_proper_names($truth, $locus, $sequences);
					
				my $oneAllele = (keys %$sequences)[0];
				my $length = length($sequences->{$oneAllele});
									
				# print "-", $sequences->{$oneAllele}, "-\n\n";
				
				die Dumper("Truth alleles not equal to 2?", \@validated_extended) unless(scalar(grep {$sequences->{$_}} @validated_extended) == 2);
				die unless(scalar(grep {$sequences->{$_}} @inferred_trimmed) == 2);
				
				print {$output_fh} join("\t",
					$inferred_trimmed[0].' Exon '.$exon,
					$inferred_trimmed[1].' Exon '.$exon,
					"",
					$validated_extended[0].' Exon '.$exon,
					$validated_extended[1].' Exon '.$exon,
					), "\n";
				
				# print Dumper($locus, [keys %{$pileup_href->{$indivID_withI}}]);

				for(my $i = 0; $i < $length; $i++)
				{
					my @chars_for_print;
					
					
					push(@chars_for_print, map {substr($sequences->{$_}, $i, 1)} @inferred_trimmed);
					push(@chars_for_print, '');
					push(@chars_for_print, map {substr($sequences->{$_}, $i, 1)} @validated_extended);
					
					if($chars_for_print[0] ne $chars_for_print[1])
					{
						$chars_for_print[2] = '!';
					}
					
					if($chars_for_print[0] ne $chars_for_print[3])
					{
						$chars_for_print[5] = '!';
					}
					else
					{
						$chars_for_print[5] = '';
					}
					
					if($chars_for_print[1] ne $chars_for_print[4])
					{
						$chars_for_print[6] = '!';
					}			
					else
					{
						$chars_for_print[6] = '';				
					}
					
			
					if(exists $pileup_href->{$indivID_withI}{'HLA'.$locus})
					{
						my $pileUpString = $pileup_href->{$indivID_withI}{'HLA'.$locus}[$exon-2][$i];
						die unless(defined $pileup_href->{$indivID_withI});
						die unless(defined $pileup_href->{$indivID_withI}{'HLA'.$locus});
						die unless(defined $pileup_href->{$indivID_withI}{'HLA'.$locus}[$exon-2]);
						# die "Problem with pileup for $locus / $indivID / $exon / $i " unless(defined $pileUpString);
						# next unless(defined $pileUpString);
						
						
						my %idt_per_allele;
						my $pileUpString_forSplit = $pileUpString;
						#my $pileUpString_forSplit = 'GGGGAAATACCTCAGGGAGT (35, 35, 35, 61, 62, 57, 52, 41, 56, 49, 46, 52, 42, 45, 36, 52, 51, 51, 56, 56) ';
						$pileUpString_forSplit =~ s/\([\d\s,]+\)/()/g;
						#die $pileUpString_forSplit;
						my @pileUpString_parts = split(/,/, $pileUpString_forSplit);
						
						my @all_alleles;
						my @all_identities;
						foreach my $element (@pileUpString_parts)
						{
							die Dumper("Weird element", $element, \@pileUpString_parts, $pileUpString_forSplit) unless($element =~ /(\S+) \(\d*\) \[(.+)\]/);
							my $allele = $1;
							my $details = $2;
							# e.g. 198 | 1 | 1 1 | 1 0.983762 | 
													 
							die "Weird details string $details" unless($details =~ /pairsDistance ([\-\d]+) \| alignmentLength ([\d]+) \| ([\d\.]+) \| ([\d\.]+) ([\d\.]+) \| ([\d\.]+) ([\-\d\.]+) \| /);
							my $idt1 = $6;
							my $idt2 = $7; 
							my $meanIdt;
							if($idt2 == -1)
							{ 
								$idt2 = $idt1;							
								$meanIdt = ($idt1 + $idt2) / 2;							
								push(@{$idt_per_allele{$allele}}, $idt1);							
							}
							else
							{				
								$meanIdt = ($idt1 + $idt2) / 2;							
								push(@{$idt_per_allele{$allele}}, $idt1, $idt2);														
							}
							
							push(@all_alleles, $allele);
							push(@all_identities, $meanIdt);
						}
						
						my @all_indices = (0 .. $#all_alleles);
						@all_indices = reverse sort {$all_identities[$a] <=> $all_identities[$b]} @all_indices;
						# die Dumper("Pileup too short", $length, scalar(@{$pileup_href->{$indivID_withI}{'HLA'.$locus}[$exon-2]})) unless(defined $pileUpString);
						unless(($chars_for_print[5] eq '!') or ($chars_for_print[6] eq '!'))
						{
							# $pileUpString =~ s/\[.+?\]//g;
						}
						
						my $printPileUpDetails = (1 or (($chars_for_print[5] eq '!') or ($chars_for_print[6] eq '!')));
		
						my $getShortRead = sub {
							my $read = shift;
							my $rE = shift;
							
							my @readIDs = split(/ /, $read);
							die unless($#readIDs == 1);
							
							
							
							my $rI; 
							if(exists $readIDs{$readIDs[0]})
							{
								$rI = $readIDs{$readIDs[0]};
							}
							else
							{
								my $nR = scalar(keys %readIDs);
								$nR++;
								$rI = "Read${nR}X";
								$readIDs{$readIDs[0]} = $rI;
								$readIDs{$readIDs[1]} = $rI;
							}
							
							if($printPileUpDetails)
							{
								return $rI.$rE;
							}
							else
							{
								return $rE;
							}
						};
						
						$pileUpString =~ s/(\@\@.+?)(\](\,|$))/$getShortRead->($1, $2);/ge;
						
						if(($chars_for_print[5] eq '!') or ($chars_for_print[6] eq '!'))
						{
							push(@chars_for_print, $pileUpString);
						}
						
						my @str_mean_idt;
						foreach my $key (keys %idt_per_allele)
						{
							my $m = mean($idt_per_allele{$key});
							push(@str_mean_idt, $key . '=' . scalar(@{$idt_per_allele{$key}}) . ";" . $m);
						}
						
						push(@chars_for_print, join(" ", @str_mean_idt));
						
						my $max_index_for_print = ($#all_indices >= 20) ? 19 : $#all_indices;
						push(@chars_for_print, join("", map {$all_alleles[$_]} @all_indices[0..$max_index_for_print]));
					}
					print {$output_fh} join("\t", @chars_for_print), "\n";
				}	
			}
			
			foreach my $readID (keys %readIDs)
			{
				print {$output_fh} $readIDs{$readID}, "\t", $readID, "\n";
			}
			close($output_fh);					
		}
	}
	
	if($T == 0)
	{
		my $totalAlleles = 0;
		
		my $calibration_file = 'temp/calibration_' . $locus . '_' . $sample_IDs_abbr . '.txt';	
		open(CALIBRATION, ">", $calibration_file) or die "Cannot open $calibration_file";
		print CALIBRATION join("\t", qw/Bin MeanPP PercCorrect NCorrect NIncorrect/), "\n";
		for(my $i = 0; $i <= 9; $i++)
		{
			my $meanPP = 0;
			my $percCorrect = 0;
			
			my $nCorrect = 0;
			my $nIncorrect = 0; 
			
			my $mean_normalize = 0;
			foreach my $elem (@{$calibration_baskets{$i}{correct}})			
			{
				$nCorrect += $elem->{weight};
				$meanPP += $elem->{PP}*$elem->{weight};
				$mean_normalize += $elem->{weight};
			}
			foreach my $elem (@{$calibration_baskets{$i}{incorrect}})			
			{
				$nIncorrect += $elem->{weight};
				$meanPP += $elem->{PP}*$elem->{weight};
				$mean_normalize += $elem->{weight};				
			}			
			
			if(($nCorrect+$nIncorrect) != 0)
			{
				$percCorrect = $nCorrect/($nCorrect+$nIncorrect);
			}
			
			if($mean_normalize != 0)
			{
				$meanPP = $meanPP / $mean_normalize;
			}
			
			print CALIBRATION join("\t",
				$i,
				sprintf("%.2f", $meanPP),
				sprintf("%.2f", $percCorrect),
				sprintf("%.2f", $nCorrect),
				sprintf("%.2f", $nIncorrect),
			), "\n";
			
			$totalAlleles += ($nCorrect + $nIncorrect);
		}
		close(CALIBRATION);		
		
		if(($totalAlleles > 0) or (defined $imputed_HLA_Calls{$locus}{called}))
		{
			die "Problem $totalAlleles vs $imputed_HLA_Calls{$locus}{called}" unless($totalAlleles == $imputed_HLA_Calls{$locus}{called});
		}		
		my $spatial_coverage_file = 'temp/spatialCoverage_' . $locus . '_' . $sample_IDs_abbr . '.txt';	
		open(SPATIALCOVERAGE, ">", $spatial_coverage_file) or die "Cannot open $spatial_coverage_file";			
		foreach my $exon (sort {$a <=> $b} keys %coverage_over_samples)
		{
			foreach my $exonPos (sort {$a <=> $b} keys %{$coverage_over_samples{$exon}})
			{
				my @individualValues = @{$coverage_over_samples_individualValues{$exon}{$exonPos}};
				@individualValues = sort {$a <=> $b} @individualValues;
				my $idx_10 = int($#individualValues * 0.1 + 0.5);
				my $idx_90 = int($#individualValues * 0.9 + 0.5);
				print SPATIALCOVERAGE join("\t", $exon, $exonPos, $exon . '-' . $exonPos, $coverage_over_samples{$exon}{$exonPos} / $coverage_over_samples_nSamples, $individualValues[$idx_10], $individualValues[$idx_90]), "\n";
			}
		}				
		close(SPATIALCOVERAGE)		
	}
			
	my $CR = sprintf("%.2f", $imputed_HLA_Calls{$locus}{sum} ? ($imputed_HLA_Calls{$locus}{called}/$imputed_HLA_Calls{$locus}{sum}) : 0);
	my $accuracy = sprintf("%.2f", 1 - (($problem_locus_examined{$locus} != 0) ? $problem_locus_detail{$locus}/$problem_locus_examined{$locus} : 0));
			
	print SUMMARY "\t", join("\t", $locus, $problem_locus_examined{$locus}, $CR, $accuracy), "\n";
}

my $comparions_OK = $comparisons - $compare_problems;
print "\nComparisons: $comparisons -- OK: $comparions_OK\n";
	
open(TMP_OUTPUT, '>', 'temp/hla_validation/validation_summary.txt') or die;
print "\nPER-LOCUS SUMMARY:\n";
foreach my $key (sort keys %problem_locus_detail)
{
	my $OK_locus = $problem_locus_examined{$key} - $problem_locus_detail{$key};
	my $accuracy = sprintf("%.2f", 1 - (($problem_locus_examined{$key} != 0) ? $problem_locus_detail{$key}/$problem_locus_examined{$key} : 0));
	my $CR = sprintf("%.2f", $imputed_HLA_Calls{$key}{sum} ? ($imputed_HLA_Calls{$key}{called}/$imputed_HLA_Calls{$key}{sum}) : 0);
	print "\t", $key, ": ", $OK_locus, " of ", $problem_locus_examined{$key}, ",\tAccuracy ", $accuracy, " ";
	print "\tCall rate: ", $CR,  "\n";
	
	my @fields = ($key, $problem_locus_examined{$key}, $CR, $accuracy);
	print TMP_OUTPUT join("\t", @fields), "\n";

}
close(TMP_OUTPUT);	

close(SUMMARY);

if(not $fromKourami)
{
	
	print "\nCorrect vs incorrect coverages per locus:\n";
	foreach my $locus (sort keys %problem_locus_detail)
	{
		my @avg_minMax_ok = min_avg_max(@{$locus_avgCoverages{$locus}{ok}});
		my @low_minMax_ok = min_avg_max(@{$locus_lowCoverages{$locus}{ok}});
		my @min_minMax_ok = min_avg_max(@{$locus_minCoverages{$locus}{ok}});
			
		my @avg_minMax_problems = min_avg_max(@{$locus_avgCoverages{$locus}{problems}});
		my @low_minMax_problems = min_avg_max(@{$locus_lowCoverages{$locus}{problems}});
		my @min_minMax_problems = min_avg_max(@{$locus_minCoverages{$locus}{problems}});
		
		print "\t", $locus, "\n";

		print "\t\tAverage ", join(' / ', @avg_minMax_ok), " vs ", join(' / ', @avg_minMax_problems), " [problems]", "\n";
		print "\t\tLow ", join(' / ', @low_minMax_ok), " vs ", join(' / ', @low_minMax_problems), " [problems]", "\n";
		print "\t\tMin ", join(' / ', @min_minMax_ok), " vs ", join(' / ', @min_minMax_problems), " [problems]", "\n";
	}


	print "Over all loci, all individuals:\n";
	my @avg_avg_minMax = min_avg_max(@allLoci_allIndivs_avgCoverage);
	print "\tAverage ", join(' / ', @avg_avg_minMax), "\n";
	print "\tLow ", join(' / ', @avg_avg_minMax), "\n";
	print "\tMin ", join(' / ', @avg_avg_minMax), "\n";

	if(scalar(keys %missing_reference_data))
	{
		print "Missing reference data for individuals:\n";
		foreach my $indivID (keys %missing_reference_data)
		{
			print " - ", $indivID, "\n";
		}
	}


	print "\n";
}

open(TYPES, '>', '_types_as_validated.txt') or die;
my @loci_for_print = sort {$a cmp $b} @loci;
print TYPES join(' ', 'IndividualID', map {'HLA' . $_ } @loci_for_print), "\n";
foreach my $indivID (sort keys %types_as_validated)
{
	print TYPES join(' ', $indivID, map {((exists $types_as_validated{$indivID}{$_}) and (scalar(@{$types_as_validated{$indivID}{$_}}) > 0) ) ? join('/', @{$types_as_validated{$indivID}{$_}}) : '????/????'} @loci_for_print), "\n";
}
close(TYPES);


sub compatibleStringAlleles
{
	my $alleles_validation = shift;
	my $alleles_inference = shift;

	die unless($#{$alleles_inference} == 1);
	die unless($#{$alleles_validation} == 1);

	my ($OK1, $NOTOK1, $MISSING1) = compatibleStringAlleles_noFlip($alleles_validation, $alleles_inference);
	
	my $alleles_validation_flipped = [reverse(@$alleles_validation)];
	my ($OK2, $NOTOK2, $MISSING2) = compatibleStringAlleles_noFlip($alleles_validation_flipped, $alleles_inference);
	
	# die unless($MISSING1 == $MISSING2);
	
	if(($NOTOK1 < $NOTOK2))
	{
		return ($OK1, $NOTOK1, $MISSING1);
	}
	elsif($NOTOK1 == $NOTOK2)
	{
		return (($OK1 > $OK2) ? ($OK1, $NOTOK1, $MISSING1) : ($OK2, $NOTOK2, $MISSING2));		
	}
	else
	{
		return ($OK2, $NOTOK2, $MISSING2);
	}
	
	
}


sub compatibleStringAlleles_noFlip
{
	my $alleles_validation = shift;
	my $alleles_inference = shift;

	die unless($#{$alleles_inference} == 1);
	die unless($#{$alleles_validation} == 1);

	my $alleles_OK = 0;
	my $alleles_NOTOK = 0;
	my $alleles_MISSING = 0;
	
	for(my $aI = 0; $aI < 2; $aI++)
	{	
		my ($OK, $NOTOK, $MISSING) = (compatibleStringAlleles_individual(undef, $alleles_validation->[$aI], $alleles_inference->[$aI]));
		$alleles_OK += $OK;
		$alleles_NOTOK += $NOTOK;
		$alleles_MISSING += $MISSING;
	}
		
	return ($alleles_OK, $alleles_NOTOK, $alleles_MISSING);
}


sub compatibleStringAlleles_individual
{
	my $locus = shift;
	my $allele_validation = shift;	
	my $allele_inference = shift;

	if(($allele_validation =~ /\?/) || ($allele_inference =~ /\?/))
	{
		return (0, 0, 1);
	}
	else
	{
		if($allele_validation eq $allele_inference)
		{
			return (1, 0, 0);
		}
		else
		{
			return (0, 1, 0);
		}
	}
}



# sub compatibleAlleles
# {

	# my $alleles_validation = shift;
	# my $alleles_inference = shift;

	# die unless($#{$alleles_inference} == 1);
	# die unless($#{$alleles_validation} == 1);

	# my $r1 = compatibleAlleles_noFlip($alleles_validation, $alleles_inference);
	
	# my $alleles_validation_flipped = [reverse(@$alleles_validation)];
	# my $r2 = compatibleAlleles_noFlip($alleles_validation_flipped, $alleles_inference);
	
	# return (($r1 > $r2) ? $r1 : $r2);
# }


# sub compatibleAlleles_noFlip
# {

	# my $alleles_validation = shift;
	# my $alleles_inference = shift;

	# die unless($#{$alleles_inference} == 1);
	# die unless($#{$alleles_validation} == 1);

	
	# my $alleles_compatible = 0;
	
	# for(my $aI = 0; $aI < 2; $aI++)
	# {	

		# $alleles_compatible += (compatibleAlleles_individual(undef, $alleles_validation->[$aI], $alleles_inference->[$aI]));
	# }

		
	# return $alleles_compatible;
# }


sub compatibleAlleles_individual
{
	my $locus = shift;
	my $allele_validation = shift;
	my $allele_inference = shift;
	
	# die Dumper($allele_validation, $allele_inference);
	
	$allele_validation =~ s/^.+\*//;
	
	my @components_allele_validation = split(/;/, $allele_validation);
	if(scalar(@components_allele_validation) > 1)
	{
		my $found_compatible = 0;
		foreach my $validation_allele (@components_allele_validation)
		{
			$found_compatible += compatibleAlleles_individual($locus, $validation_allele, $allele_inference);
			last if($found_compatible);
		}
		if($found_compatible)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}	
	
	my $allele_validation_original = $allele_validation;
		
	
	die "Weird allele $locus $allele_validation" unless(length($allele_validation) >= 4);
	die unless($allele_inference =~ /:/);

	
	#die unless($allele_validation =~ /^\d\d/);

	my $components_validation;
	if($allele_validation !~ /\:/)
	{
		die if($allele_validation =~ /\//);
		$allele_validation = substr($allele_validation, 0, 2).':'.substr($allele_validation, 2);
		$components_validation = 2;
		$allele_validation =~ s/g//i;
		
	}
	else
	{
		my @_components = split(/\:/, $allele_validation);
		$components_validation = scalar(@_components);
		#die "Weird allele code $allele_validation" if($allele_validation =~ /g/i);
		# die Dumper("Weird allele code (II) $allele_validation", $allele_validation, $allele_inference)  if($allele_validation =~ /D/i);
			
	}
	die unless(defined $components_validation);
	die unless($components_validation >= 2);

	
	my $true_allele = $allele_validation;

		
	my $inferred_alleles = $allele_inference;
	my @inferred_alleles = split(/;/, $inferred_alleles);
	
	# if inferred allele has fewer components than validation allele: add 0s, mismatch if validation allele is not 0 at this position
	# => conservative
	# if inferred allele has more components than validation allele, truncate inferred allele
	# e.g. inferred: 02:03:01, validation allele: 03:01, truncated inferred allele: 02:03
	# if inferred and validation allele have different 4-digit groups, we always generate a mismatch
	# if inferred and validation allele have the same 4-digit group, we generate a match
	# e.g. inferred 02:03:01, validation: 02:03
	# is this a problem? No. Consider two scenarios:
	# - validation allele typed to complete 4-digit (amino acid sequence, all exons) resolution: this is what we want
	# - validation allele typed to SBT G-groups: this case is not a problem.
	# 		Case 1: specified as G, i.e. 02:03G
	#			This is a bad example, because 02:03G is not a valid G allele - needs to be six digits, e.g.
	#			02:03:01G. This is fine as long as our output is also at G group level (i.e. no
	#			higher resolution than exon 2 / exons 2,3), as 02:03:01 will be one of our output alleles
	#			if and only if we infer 02:03:01G as a result.
	# 			Gs are confusing and we will deactivate them now.
	#			(Or enable proper translation)
	#  		Case 2: specified as explicit allele list, separated with slashes - we now assume that 02:03 is 
	#   			one such explicit member of a list.
	#			This case does not exist. If 02:03:01 exists as a named allele, other alleles like 02:03:02 must also exist,
	#			and 02:03 would not be an individual member of any G group (because all members with identical
	#			4-digit group would be differentiated by their 6-digit groups).

	
	@inferred_alleles = map {
		# die "Can't parse allele $_" unless($_ =~ /^\s*\w+\*(\d\d\d?)\:(\d\d\d?N?).*$/);
		# my $allele = $1.':'.$2;
		# $allele
		
		die "Can't parse allele $_" unless($_ =~ /^\s*(\w+)\*([\d\:N]+Q?)L?S?$/);
		
		my $allele = $2;
		
		my @components_allele = split(/:/, $allele);
		my @components_allele_rightLength;
		
		for(my $i = 0; $i < $components_validation; $i++)
		{
			if($i <= $#components_allele)
			{
				push(@components_allele_rightLength, $components_allele[$i]);
			}
			else
			{
				push(@components_allele_rightLength, '00');
			}
		}
		$allele = join(':', @components_allele_rightLength);
		
	} @inferred_alleles;
	
	my %inferred = map {$_ => 1} @inferred_alleles;
	
	if($inferred{$true_allele})
	{
		# print Dumper($allele_validation_original, $allele_inference, \@inferred_alleles, 1), "\n" if ($locus eq 'DQB1');	
		
		return 1;  
	}		
	else
	{
		# print Dumper($allele_validation_original, $allele_inference, \@inferred_alleles, 0), "\n"  if ($locus eq 'DQB1');	
	
		return 0;
	}
}

sub intersection
{
	my $aref1 = shift;
	my $aref2 = shift;

	my %h1 = map {$_ => 1} @$aref1;
	my %h2 = map {$_ => 1} @$aref2;
	
	my %c = map {$_ => 1} ((keys %h1), (keys %h2));
	
	my @ret = grep {$h1{$_} and $h2{$_}} keys %c;
	
	return @ret;
}

sub print_which_exons
{
	my $locus = shift;
	if((length($locus) == 1) or ($locus =~ /HLA\w/))
	{
		return (2, 3);
	}
	else
	{
		return (2);
	}
}	


sub load_pileup
{
	my $r_href = shift;
	my $file = shift;
	my $indivID = shift;
	
	die unless($file =~ /R\d+_pileup_(\w+).txt$/);
	my $locus = $1;
	
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line, -1);
		die unless(($#f >= 3) or ($#f == 2));
		if($#f >= 3)
		{
			my $exon = $f[0];
			my $exonPos = $f[1];
			my $pileUp = $f[3];		
			my $pileUp_post = $f[4];		
			$r_href->{$indivID}{'HLA'.$locus}[$exon][$exonPos] = $pileUp_post . ' ' . $pileUp;
		}
		else
		{
			my $exon = $f[0];
			my $exonPos = $f[1];
			my $pileUp = '';		
			$r_href->{$indivID}{'HLA'.$locus}[$exon][$exonPos] = $pileUp;
		}
	}	
	close(F);
}

sub load_coverages_from_pileup
{
	my $file = shift;
	
	die unless($file =~ /R\d+_pileup_(\w+).txt$/);
	my $locus = $1;
	
	my %forReturn;
	
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		die "Cannot parse pileup coverage line" unless($line =~ /^(\d+)\s(\d+)\s(\d+)/);
		my $exon = $1;
		my $exonPos = $2;
		my $coverage = $3;				
		
		die "Pos twice in $file ? $exon $exonPos" if(defined $forReturn{$exon}{$exonPos});
		$forReturn{$exon}{$exonPos} = $coverage;
	}	
	close(F);
	
	return \%forReturn;
}


sub twoValidationAlleles_2_proper_names
{
	my $alleles_validation = shift;
	my $locus = shift;
	my $exon_sequences = shift;
	
	die unless($#{$alleles_validation} == 1);
	
	my @forReturn;
	
	$locus =~ s/HLA//;
	for(my $aI = 0; $aI < 2; $aI++)
	{
		my @vA = split(/;/, $alleles_validation->[$aI]);
		REFALLELE: foreach my $validation_allele (@vA)
		{
			my $original_validation_allele = $validation_allele;
			
			if($validation_allele =~ /\:/)
			{
				$validation_allele = $locus . '*' . $validation_allele;
				# $validation_allele =~ s/g//;
				# my @components = split(/\:/, $validation_alleles),
				
				# die unless(length($validation_allele) >= 4);
				
				# $validation_allele = substr($validation_allele, 0, 2) . ':' . substr($validation_allele, 2);
				# $validation_allele = $locus . '*' . $validation_allele;
			
			}
			else
			{
				$validation_allele =~ s/g//;
				die unless(length($validation_allele) >= 4);
				
				$validation_allele = substr($validation_allele, 0, 2) . ':' . substr($validation_allele, 2);
				$validation_allele = $locus . '*' . $validation_allele;
			}
			my @validation_alleles = ($validation_allele);
			
			my @history_extensions = ($validation_allele);
			my $extensions = 0;
			my $notFoundFirstRound = 0;
			while( not exists $exon_sequences->{$validation_allele} )
			{
				$extensions++;
			
				$validation_allele .= ':01';
				push(@history_extensions, $validation_allele);
				if($extensions > 3)
				{
					$notFoundFirstRound = 1;
					last;
				}
				
			}
			
			if($notFoundFirstRound)
			{
				my $extensions = 0;
				while(scalar(grep {exists $exon_sequences->{$_}} @validation_alleles) == 0)
				{
					$extensions++;
					my $viMax = $#validation_alleles;
					for(my $vI = 0; $vI <= $viMax; $vI++)
					{
						my $a1 = $validation_alleles[$vI];
						$validation_alleles[$vI] .= ':01';
						push(@validation_alleles, $a1.':02');
					}
					
					$validation_allele .= ':01';
					if($extensions > 3)
					{
						# print Dumper([grep {$_ =~ /DQB1/} keys %$exon_sequences]);
						if($validation_allele eq $vA[$#vA])
						{
							#die Dumper("Can't identify (II) exon alleles for $alleles_validation->[$aI]", \@validation_alleles, $locus, $alleles_validation, $extensions, [(keys %$exon_sequences)[0 .. 10]], $validation_allele, $original_validation_allele, \@history_extensions);
							warn Dumper("Can't identify (II) exon alleles for $alleles_validation->[$aI]", \@validation_alleles, $locus, $alleles_validation, $extensions, [(keys %$exon_sequences)[0 .. 10]], $validation_allele, $original_validation_allele, \@history_extensions);
							push(@forReturn, $alleles_validation->[$aI]);
							last REFALLELE;
						}
						else
						{
							next REFALLELE;
						}
					}
				}		
				
				my @foundAlleles = grep {exists $exon_sequences->{$_}} @validation_alleles;
				die unless(scalar(@foundAlleles) > 0);
				push(@forReturn, $foundAlleles[0]);	
				last REFALLELE;				
				
			}
			else
			{
				push(@forReturn, $validation_allele);
				last REFALLELE;
			}
		}
	}
			
	return @forReturn;
}

sub list_comparison
{
	my $list1_aref = shift;
	my $list2_aref = shift;
	
	my %l1 = map {$_ => 1} @$list1_aref;
	my %l2 = map {$_ => 1} @$list2_aref;
	
	unless(scalar(keys %l1) == scalar(@$list1_aref))
	{
		die "List 1 non-unique elements";
	}
	
	unless(scalar(keys %l2) == scalar(@$list2_aref))
	{
		die "List 2 non-unique elements";
	}	
	
	my %combined = map {$_ => 1} (@$list1_aref, @$list2_aref);
	
	my @l1_exclusive;
	my @l2_exclusive;
	my $n_shared = 0;
	foreach my $e (keys %combined)
	{
		if($l1{$e} and $l2{$e})
		{
			$n_shared++;
		}
		elsif($l1{$e})
		{
			push(@l1_exclusive, $e);
		}
		elsif($l2{$e})
		{
			push(@l2_exclusive, $e);
		}
		else
		{
			die;
		}
				
	}
	
	die unless(($n_shared + scalar(@l1_exclusive) + scalar(@l2_exclusive)) == scalar(keys %combined));
	
	return ($n_shared, \@l1_exclusive, \@l2_exclusive);
}

sub twoClusterAlleles
{
	my $alleles_inference = shift;
	die Dumper("Don't have two inferred alleles", $alleles_inference) unless($#{$alleles_inference} == 1);
	
	my @forReturn;
	
	for(my $aI = 0; $aI < 2; $aI++)
	{
		my $inferred_alleles = $alleles_inference->[$aI];
		my @inferred_alleles = split(/;/, $inferred_alleles);
	
		$inferred_alleles[0] =~ s/^\s+//;
		
		push(@forReturn, $inferred_alleles[0]);
	}
	
	return @forReturn;
}


sub read_exon_sequences
{
	my $file = shift;
	my %r;
	open(F, '<', $file) or die "Cannot open $file";
	my $firstLine = <F>;
	chomp($firstLine);
	my @header_fields = split(/ /, $firstLine);
	die unless($header_fields[0] eq 'IndividualID');
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\r//g;
		$line =~ s/\n//g;
		my @line_fields = split(/ /, $line);
		my $allele = shift(@line_fields);
		$r{$allele} = join('', @line_fields);
	}
	close(F);
	
	return \%r;
}


sub find_exon_file
{
	my $locus = shift;
	my $exon = shift;
	$locus =~ s/HLA//;
	
	opendir(my $dh, $exon_folder) or die;
	my @exon_folder_files = readdir($dh);
	closedir($dh);
	
	@exon_folder_files = grep {$_ =~ /HLA-${locus}_\d+_exon_${exon}\.txt$/} @exon_folder_files;
	if($#exon_folder_files != 0)
	{
		die Dumper("Can't find exon file for $locus // $exon", @exon_folder_files, $exon_folder);
	}
	else
	{
		return $exon_folder.'/'.$exon_folder_files[0];
	}	
}


sub min_avg_max
{
	my @v = @_;
	@v = sort {$a <=> $b} @v;
	if($#v >= 1)
	{
		die unless($v[0] <= $v[1]);
	}
	
	if(scalar(@v) == 0)
	{
		return ('', '', '');
	}
	if(scalar(@v) == 1)
	{
		return($v[0], $v[0], '');
	}
	
	my $min = $v[0];
	my $max = $v[$#v];
	  
	my $sum = 0;
	foreach my $vE (@v)
	{
		$sum += $vE;
	}
	
	my $avg = '';
	if(scalar(@v) > 0)
	{
		$avg = $sum / scalar(@v);
	}
	
	return ($min, $avg, $max);
}

sub mean
{
	my $aref = shift;
	my @a = @$aref;
	if(scalar(@a))
	{
		return sum(@a)/scalar(@a);
	}
	else
	{
		return 0;
	}
}
