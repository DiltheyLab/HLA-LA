use warnings;
use strict;
use FindBin;
use Getopt::Long;
use Data::Dumper; 
use List::Util qw/min max/;
use List::MoreUtils qw/all mesh/;
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
my $samtools_bin = 'samtools';

my $GRCh38_refGenome = '/home/dilthey/GRCh38_full_analysis_set_plus_decoy_hla.fa';
die unless(-e $GRCh38_refGenome);

my $hla_gen = '../hla_gen.fasta';
my $graph = 'PRG_MHC_GRCh38_withIMGT';
my $outputPrefix = '../../haplotypeSimulations';
my $n = "0,0,1";
my $targetCoverage = 15;
my $readLength = 100;
GetOptions (
	'hla_gen:s' => \$hla_gen,
	'graph:s' => \$graph,
	'outputPrefix:s' => \$outputPrefix,
	'n:s' => \$n,
);
	
my $graph_dir = $graphs_dir . '/' . $graph;

unless($hla_gen and (-e $hla_gen))
{
	die "No IMGT/HLA genomic file specified (--hla_gen), or specified file '$hla_gen' not present";
}
unless(-d $outputPrefix)
{
	die "Please provide existing directory as output prefix (--outputPrefix)";
}
unless($n =~ /^\d+,\d+,\d+$/)
{
	die "Please provide a simple comma-separated list for parameter --n (individuals with no mutations, with rate 1:5000, with rate 1:500";
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
my $gen_fasta_href = Util::readFASTA($hla_gen, 0);

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
	my @alleleID_parts = split(/ /, $alleleID);
	die unless(scalar(@alleleID_parts) == 4);
	my $coreAlleleID = $alleleID_parts[1];
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
		my $minAlleleLength = 0.75 * max(@alleleLengths);
		my $survivingAlleles = 0;
		
		foreach my $allele (@availableAlleles)
		{
			my $alleleLength = length($useAllelesForSimulations{$locus}{$allele});
			if($alleleLength >= $minAlleleLength)
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
			}
			else
			{
				delete $useAllelesForSimulations{$locus}{$allele};
			}
		}
		print " -> " . scalar(keys %{$useAllelesForSimulations{$locus}}) . " [>= $minAlleleLength bp] ";
		if(scalar(keys %{$useAllelesForSimulations{$locus}}) > 0)
		{
			push(@simulateLoci, $locus);
		}
	}
	print "\n";
}
print "... simulate " . scalar(@simulateLoci) . " loci.\n\n";
@simulateLoci = sort @simulateLoci;

my @mutationRates = (0, 1/5000, 1/500);
my @mutationRates_n = split(/,/, $n);
die unless(scalar(@mutationRates) == scalar(@mutationRates_n));

my $outputFn_HLATypes = $outputPrefix . '/HLATypes.txt';
my $outputFn_HLADetails = $outputPrefix . '/HLAHaplotypeDetails.txt';
foreach my $f ($outputFn_HLATypes, $outputFn_HLADetails)
{
	if(-e $f)
	{
		unlink($f) or die "Cannot unlink $f";
	}				
}

print "Simulating...\n";
my %HLA_truth_types;
my $globalIndivI = 0;
for(my $mutationRateI = 0; $mutationRateI <= $#mutationRates; $mutationRateI++)
{
	my $indiv_n = $mutationRates_n[$mutationRateI];
	my $rate = $mutationRates[$mutationRateI];
	print "Generating $indiv_n individuals with mutation rate $rate .. \n";
	for(my $indivI = 0; $indivI < $indiv_n; $indivI++)
	{
		$globalIndivI++;
		my $sampleID = 'Indiv' . $globalIndivI;
		print "\t", $sampleID, "\n";
		my $outputDir = $outputPrefix . '/' . $sampleID;
		unless(-e $outputDir)
		{
			mkdir($outputDir) or die "Cannot mkdir $outputDir";
		}
		 
		my $fastq_combined_1 = $outputDir . '/fastq_1.fastq';
		my $fastq_combined_2 = $outputDir . '/fastq_2.fastq';
		my $BAM_combined = $outputDir . '/reads.bam';
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

				
	
				my $alleleSequence_ref_aligned = $alleleRef_aligned;
				my $alleleSequence_mutated_aligned = $alleleSequence_aligned;
				my $alleleSequence_original_aligned = $alleleSequence_aligned;
				
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
							}
							else
							{
								die unless(length($newAllele) == 2);
								die unless(substr($newAllele, 0, 1) eq $oldAllele);
								$alleleSequence_ref_aligned .= ($oldAllele_ref . '-');		
								$alleleSequence_mutated_aligned .= $newAllele;
								$alleleSequence_original_aligned .= ($oldAllele . '-');																
							}
							$HLA_truth_types{$sampleID}{$locus}{'mutations'}{'h' . $hap}{$runningAllelePos_nonAligned . ':' . $oldAllele . '->' . $newAllele}++;
						}
						else
						{
							$alleleSequence_ref_aligned .= $oldAllele_ref;
							$alleleSequence_mutated_aligned .= $oldAllele;
							$alleleSequence_original_aligned .= $oldAllele;
						}
					}
				}
				else
				{
					$alleleSequence_ref_aligned = $alleleRef_aligned;
					$alleleSequence_mutated_aligned = $alleleSequence_aligned;
					$alleleSequence_original_aligned = $alleleSequence_aligned;					
				}
				die unless(length($alleleSequence_mutated_aligned) == length($alleleSequence_original_aligned));
				die unless(length($alleleSequence_ref_aligned) == length($alleleSequence_original_aligned));
				(my $alleleSequence_mutated_raw = $alleleSequence_mutated_aligned) =~ s/-//g;
				
				$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{raw}{'h' . $hap} = $alleleSequence_mutated_raw;
				
				$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{reference}{'h' . $hap} = $alleleSequence_ref_aligned;
				$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{afterMutation}{'h' . $hap} = $alleleSequence_mutated_aligned;
				$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{beforeMutation}{'h' . $hap} = $alleleSequence_original_aligned;

				$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInAlignedPseudogenomic}{'h' . $hap} = $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{alignment}[0];
				$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInRawPseudogenomic}{'h' . $hap} = $useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{raw}[0];
				$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInPGF}{'h' . $hap} = 'NA';


				if(exists $pseudoGenomic_refPos{$locus_2_locusHla{$locus}})
				{
					$HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{alignmentStartInPGF}{'h' . $hap} = 
						$pseudoGenomic_refPos{$locus_2_locusHla{$locus}}+ 
						$useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$selectedAllele}{raw}[0];
				}
				
#				$useAllelesForSimulations_positionInPseudoGenomicAllele{$locus}{$allele}{alignment} = $IMGT_genomic_start_in_Gaps;
				
				my $fasta_for_wgsim = $outputDir . '/' . $locus . '_h' . $hap . '.input.fa';
				my $fastq_from_wgsim_1 = $outputDir . '/' . $locus . '_h' . $hap . '.output_1.fastq';
				my $fastq_from_wgsim_2 = $outputDir . '/' . $locus . '_h' . $hap . '.output_2.fastq';
				
				Util::writeFASTA($fasta_for_wgsim, {'HLA' . $locus . '_h' . $hap => $alleleSequence_mutated_raw});
				
				my $want_read_pairs = (length($alleleSequence_mutated_raw) * $targetCoverage) / (2 * $readLength);
				my $cmd_wgsim = qq($wgsim_bin -e 0.01 -d 400 -N $want_read_pairs -1 $readLength -2 $readLength -r 0 -R 0 -X 0 $fasta_for_wgsim $fastq_from_wgsim_1 $fastq_from_wgsim_2);
				system($cmd_wgsim) and die "Could not execute: $cmd_wgsim";
				
				my $cmd_cat = qq(cat $fastq_from_wgsim_1 >> $fastq_combined_1 && cat $fastq_from_wgsim_2 >> $fastq_combined_2);
				system($cmd_cat) and die "Cat command '$cmd_cat' failed";
				
				unlink($fasta_for_wgsim);
				unlink($fastq_from_wgsim_1);
				unlink($fastq_from_wgsim_2);
			}
			my $n_mutations = scalar(keys %{$HLA_truth_types{$sampleID}{$locus}{'mutations'}{'h1'}}) + scalar(keys %{$HLA_truth_types{$sampleID}{$locus}{'mutations'}{'h1'}});
			print "\t\t$locus $n_mutations mutations\n"
		}
	
		my $cmd_bwa = qq($bwa_bin mem $GRCh38_refGenome $fastq_combined_1 $fastq_combined_2 | $samtools_bin sort --reference $GRCh38_refGenome -O BAM - > $BAM_combined && samtools index $BAM_combined);
		#system($cmd_bwa) and die "BWA mapping failed - command '$cmd_bwa'";
	}
}


open(HLATYPES, '>', $outputFn_HLATypes) or die "Cannot open $outputFn_HLATypes";
open(HLADETAILS, '>', $outputFn_HLADetails) or die "Cannot open $outputFn_HLADetails";
print HLATYPES join("\t", 'IndividualID', @simulateLoci), "\n";
print HLADETAILS join("\t",
	'IndividualID', 'Locus', 'Haplotype', 'Mutations',
	'rawUnderlyingSequence', 'underlyingSequence_aligned_ref', 'underlyingSequence_aligned_postMutationSeq', 'underlyingSequence_aligned_preMutationSeq',
	'alignmentStart_pseudoGenomic_aligned', 'alignmentStart_pseudoGenomic_raw', 'alignmentStart_PGF',
), "\n";
foreach my $sampleID (sort keys %HLA_truth_types)
{
	my @outputFields_HLATypes = ($sampleID);
	foreach my $locus (@simulateLoci)
	{
		my @alleles = map {$HLA_truth_types{$sampleID}{$locus}{'alleles'}{'h' . 1}} (0, 1);
		push(@outputFields_HLATypes, join('/', @alleles));
		foreach my $haplotype (0, 1)
		{
			my @outputFields_Details = ($sampleID, $locus, $haplotype + 1);
			my $haplotype_id = 'h' . $haplotype;
			push(@outputFields_Details, join('; ', sort keys %{$HLA_truth_types{$sampleID}{$locus}{'mutations'}{$haplotype_id}}));
			push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{'raw'}{$haplotype_id});
			push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'reference'}{$haplotype_id});
			push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'afterMutation'}{$haplotype_id});
			push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'beforeMutation'}{$haplotype_id});
			push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'alignmentStartInAlignedPseudogenomic'}{$haplotype_id});
			push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'alignmentStartInRawPseudogenomic'}{$haplotype_id});
			push(@outputFields_Details, $HLA_truth_types{$sampleID}{$locus}{'alleleHaplotypes'}{aligned}{'alignmentStartInPGF'}{$haplotype_id});
			print HLADETAILS join("\t", @outputFields_Details), "\n";
		}
	}
	print HLATYPES join("\t", @outputFields_HLATypes), "\n";
}
close(HLATYPES);

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


