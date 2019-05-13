#!/usr/bin/env perl

BEGIN {
	use FindBin;
	push(@INC, $FindBin::Bin);
	push(@INC, $FindBin::RealBin);
}

use strict;
use List::MoreUtils qw/mesh all/;
use List::Util qw/max min/;
use Getopt::Long;
use Data::Dumper;
use Storable qw/store retrieve dclone/;

my $HLA;
my $B37 = 1;
my $chr6_B37 = '/gpfs/project/dilthey/projects/JasonMHC/chr6.GRCh37.fasta';
my $chr6_B38 = '/gpfs/project/dilthey/projects/JasonMHC/chr6.GRCh38.fasta';
my $prefix = 'HLAtypes';
GetOptions (
	'HLA:s' => \$HLA,
	'B37:s' => \$B37,
	'chr6_B37:s' => \$chr6_B37,
	'chr6_B38:s' => \$chr6_B38,
	'prefix:s' => \$prefix,
);

die "Please specify parameter --HLA (specifying the HLA types to be converted into VCF format; HLA*ASM format)" unless($HLA);
die "Specified file $HLA (via --HLA) not existing" unless(-e $HLA);
my $ref_sequence;
my $ref_sequence_startCoordinateField;
my $ref_sequence_stopCoordinateField;
if($B37)
{
	die "Please provide a one-contig reference file of chromosome 6, B37 (parameter --chr6_B37)" unless($chr6_B37);
	die "Input file for one-contig reference file of chromosome 6, B37 ($chr6_B37) not existing. Please modify via parameter --chr6_B37." unless(-e $chr6_B37);
	my $chr6_href = readFASTA($chr6_B37);
	unless(scalar((values %$chr6_href)) == 1)
	{
		die "Reference file for chromosome 6, B37 (parameter --chr6_B37) is not one-contig";
	}
	$ref_sequence = uc((values %$chr6_href)[0]);
	$ref_sequence_startCoordinateField = 'FirstBase_B37_0based';
	$ref_sequence_stopCoordinateField = 'LastBase_B37_0based';
}
else
{ 
	die "Please provide a one-contig reference file of chromosome 6, B38 (parameter --chr6_B38)" unless($chr6_B38);
	die "Input file for one-contig reference file of chromosome 6, B38 ($chr6_B38) not existing. Please modify via parameter --chr6_B38." unless(-e $chr6_B38);
	my $chr6_href = readFASTA($chr6_B38);
	unless(scalar((values %$chr6_href)) == 1)
	{
		die "Reference file for chromosome 6, B37 (parameter --chr6_B38) is not one-contig";
	}
	$ref_sequence = uc((values %$chr6_href)[0]);
	$ref_sequence_startCoordinateField = 'FirstBase_B38_0based';
	$ref_sequence_stopCoordinateField = 'LastBase_B38_0based';	
}

# HLA gene classes
my %classI = map {$_ => 1} qw/HLA-A HLA-B HLA-C/;
my %classII = map {$_ => 1} qw/HLA-DQA1 HLA-DQB1 HLA-DRB1 HLA-DRB3 HLA-DRB4/;

# read G group info
my %G_full_alleles;
my %G_to_underlyingAlleles;
my %alleles_to_fullGAmbiguity;
my %G_mapper_unambigious;
readG();

# read specified HLA types
my %HLAtypes_by_locus;
open(HLA, '<', $HLA) or die "Cannot open $HLA";
while(<HLA>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/[\n\r]//g;
	next unless($line);
	my @line_fields = split(/\t/, $line);
	
	die "Invalid line in $HLA - expect tab-separated: HLA-LOCUS allele1 allele2. Example: HLA-A\t01:01:01:01\t11:01:01:01" unless(scalar(@line_fields) == 3);
	die "Invalid line in $HLA - expect tab-separated: HLA-LOCUS allele1 allele2. Example: HLA-A\t01:01:01:01\t11:01:01:01" unless($line_fields[0] =~ /^HLA-/);
	die "Duplicate locus in $HLA for $line_fields[0]?" if(exists $HLAtypes_by_locus{$line_fields[0]});
	
	my @alleles_tr;
	foreach my $allele ($line_fields[1], $line_fields[2])
	{
		die "Invalid line in $HLA - expect tab-separated: HLA-LOCUS allele1 allele2. Example: HLA-A\t01:01:01G\t11:01:01:01.\n\nConversion is ALWAYS carried out a G group resolution." unless($line_fields[1] =~ /\:/);
		my $locus_for_gGroupLookup = $line_fields[0];
		$locus_for_gGroupLookup =~ s/HLA-/HLA/;
		die "Unknown locus $line_fields[0]" unless(exists $G_to_underlyingAlleles{$locus_for_gGroupLookup});
		die "Can't resolve allele $line_fields[0] $allele at G-group resolution" unless(defined $G_to_underlyingAlleles{$locus_for_gGroupLookup}{$allele});		
		my $allele_tr = $G_to_underlyingAlleles{$locus_for_gGroupLookup}{$allele}[0];
		$allele_tr =~ s/^\S+\*//;
		push(@alleles_tr, $allele_tr);
	}
			
	$HLAtypes_by_locus{$line_fields[0]} = \@alleles_tr;
	
}
close(HLA);

# read reference (PGF) HLA types
my %HLAtypes_PGF;
open(PGF, '<', 'Perl/PGF_loci_and_alleles.txt') or die "Cannot open Perl/PGF_loci_and_alleles.txt";
my $PGFalleles_header_line = <PGF>;
$PGFalleles_header_line =~ s/[\n\r]//g;
my @PGF_header_fields = split(/\t/, $PGFalleles_header_line);
while(<PGF>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/[\n\r]//g;
	next unless($line);
	my @line_fields = split(/\t/, $line, -1);
	my %line = (mesh @PGF_header_fields, @line_fields);
	
	die unless(exists $line{Locus});
	
	die unless(exists $line{PGFAllele});
	die unless($line{PGFAllele} =~ /^(\w+)\*(.+)$/);
	my $locus = $1;
	my $allele_numeric = $2;
	my $locus_with_HLA = 'HLA' . $locus;
	
	next unless((exists $classI{'HLA-' . $locus}) or (exists $classII{'HLA-' . $locus}));

	$HLAtypes_PGF{$line{Locus}} = \%line;

	unless((exists $alleles_to_fullGAmbiguity{$locus_with_HLA}) and (exists $alleles_to_fullGAmbiguity{$locus_with_HLA}{$allele_numeric}))
	{
		$allele_numeric .= 'G';
		die Dumper("Undefined allele $line{PGFAllele}", [(keys %{$alleles_to_fullGAmbiguity{$locus_with_HLA}})[0 .. 10]]) unless((exists $alleles_to_fullGAmbiguity{$locus_with_HLA}) and (exists $alleles_to_fullGAmbiguity{$locus_with_HLA}{$allele_numeric}));
	}
}
close(PGF);

my $exon_allele_alignments_href;
{
	unless(-d 'temp')
	{
		mkdir('temp') or die "Cannot mkdir temp"; 
	}

	my $exon_allele_alignments_href_file = 'temp/_exon_allele_alignments_href';
	if(-e $exon_allele_alignments_href_file)
	{
		$exon_allele_alignments_href = retrieve $exon_allele_alignments_href_file;
	}
	else
	{
		$exon_allele_alignments_href = read_exon_alignments(1);
		store($exon_allele_alignments_href, $exon_allele_alignments_href_file);
	}
}
		
my %alleles_2_features;
foreach my $locus (keys %$exon_allele_alignments_href)
{
	foreach my $allele (keys %{$exon_allele_alignments_href->{$locus}})
	{
		my @exons = split(/\|/, $exon_allele_alignments_href->{$locus}{$allele});
		die unless(scalar(@exons) > 1);
		for(my $exonI = 1; $exonI <= scalar(@exons); $exonI++)
		{
			my $exonSeq = $exons[$exonI-1];
			$exonSeq =~ s/[\-_]//g;
			$exonSeq = uc($exonSeq);
			die Dumper("Weird exon sequence", $exonSeq) unless($exonSeq =~ /^[ACGTN]*$/);
			$alleles_2_features{$locus}{$allele}{'exon'.$exonI} = $exonSeq;
		}
	}
}		

my $fn_BED = $prefix . '.bed';
my $fn_VCF = $prefix . '.VCF';
my $fh_BED;
my $fh_VCF;
open($fh_BED, '>', $fn_BED) or die "Cannot open $fn_BED";
open($fh_VCF, '>', $fn_VCF) or die "Cannot open $fn_VCF";
foreach my $locus (sort {$HLAtypes_PGF{$a}{$ref_sequence_startCoordinateField} <=> $HLAtypes_PGF{$b}{$ref_sequence_startCoordinateField}} keys %HLAtypes_by_locus)
{
	die unless($locus =~ /HLA\-/);
	(my $locus_noHLA = $locus) =~ s/HLA-//;
	die Dumper("Unknown locus $locus", [keys %HLAtypes_PGF]) unless(exists $HLAtypes_PGF{$locus});
	
	my $locus_key_for_exon_alignments = $locus;
	$locus_key_for_exon_alignments =~ s/DRB\d/DRB/; 
	die Dumper("Unknown locus $locus II", [keys %alleles_2_features]) unless(exists $alleles_2_features{$locus_key_for_exon_alignments});
	
	print $locus, "\n";
	my $PGF_reference_allele = $HLAtypes_PGF{$locus}{PGFAllele};
	
	die "Weird reference allele '$PGF_reference_allele'" unless($PGF_reference_allele =~ /^(\w+)\*(.+)$/);
	$PGF_reference_allele =~ s/^(\w+)\*(.+)$//;
	my $locus_noHLA = $1;
	die 'HLA-' . $locus_noHLA . " ne $locus" unless('HLA-' . $locus_noHLA eq $locus);
	my $PGF_reference_allele_numeric = $2;
	die unless(exists $alleles_to_fullGAmbiguity{'HLA' . $locus_noHLA});
	unless(exists $alleles_to_fullGAmbiguity{'HLA' . $locus_noHLA}{$PGF_reference_allele_numeric})
	{
		$PGF_reference_allele_numeric .= 'G';
	}
	die unless(exists $alleles_to_fullGAmbiguity{'HLA' . $locus_noHLA}{$PGF_reference_allele_numeric});
	my @G_group_reference_alleles = @{$alleles_to_fullGAmbiguity{'HLA' . $locus_noHLA}{$PGF_reference_allele_numeric}};
	
	my $PGF_reference_allele_with_locus = $locus_noHLA . '*' . $G_group_reference_alleles[0];

	die Dumper("Unknown reference allele $PGF_reference_allele_with_locus ($locus)", [keys %{$alleles_2_features{$locus_key_for_exon_alignments}}]) unless(exists $alleles_2_features{$locus_key_for_exon_alignments}{$PGF_reference_allele_with_locus});

	my $gene_start = $HLAtypes_PGF{$locus}{$ref_sequence_startCoordinateField};
	my $gene_stop = $HLAtypes_PGF{$locus}{$ref_sequence_stopCoordinateField};
	die unless((defined $gene_start) and (defined $gene_stop));
	
	my $gene_sequence = uc(substr($ref_sequence, $gene_start, $gene_stop - $gene_start + 1));
	die unless($gene_sequence);
	my $gene_strand = $HLAtypes_PGF{$locus}{Strand};
	die unless(($gene_strand eq '+') or ($gene_strand eq '-'));
	
	my @relevantExons = ($classI{$locus}) ? ('exon2', 'exon3') : ('exon2');
	if($gene_strand eq '-')
	{
		@relevantExons = reverse @relevantExons;
	}
	foreach my $exonID (@relevantExons)
	{
		my $reference_alignment = uc($alleles_2_features{$locus_key_for_exon_alignments}{$PGF_reference_allele_with_locus}{$exonID});
		die unless($reference_alignment);
		if($gene_strand eq '-')
		{
			$reference_alignment = reverseComplement($reference_alignment);
		}
		
		(my $reference_sequence_noGaps = $reference_alignment) =~ s/[\-_]//g;
		die unless($reference_sequence_noGaps =~ /^[ACGTN]+$/);
		
		my $exon_start_within_gene = index($gene_sequence, $reference_sequence_noGaps);
		die if($exon_start_within_gene == -1);
		die unless(index($gene_sequence, $reference_sequence_noGaps, $exon_start_within_gene + 1) == -1);
		
		#print "\t", $reference_sequence_noGaps, "\n";
		#print "\t", $gene_sequence, "\n";
		#warn if($exon_start_within_gene == -1);
		
		my @allele_alignments;
		my @sample_alleles = @{$HLAtypes_by_locus{$locus}};
		foreach my $allele (@sample_alleles)
		{
			my $allele_with_locus = $locus_noHLA . '*' . $allele;
			die Dumper("Unknown allele $allele_with_locus", [keys %{$alleles_2_features{$locus_key_for_exon_alignments}}]) unless(exists $alleles_2_features{$locus_key_for_exon_alignments}{$allele_with_locus});
			my $exonSequence_aligned = $alleles_2_features{$locus_key_for_exon_alignments}{$allele_with_locus}{$exonID};
			if($gene_strand eq '-')
			{
				$exonSequence_aligned = reverseComplement($exonSequence_aligned);
			}			
			push(@allele_alignments, $exonSequence_aligned);
			die unless(length($exonSequence_aligned) == length($reference_alignment));
			
			die unless(substr($exonSequence_aligned, 0, 1) eq substr($reference_sequence_noGaps, 0, 1));
			die unless(substr($exonSequence_aligned, length($exonSequence_aligned) - 1, 1) eq substr($reference_sequence_noGaps, length($reference_sequence_noGaps) - 1, 1));			
		}
		
		my $exon_first_coordinate = $gene_start + $exon_start_within_gene;
		my $exon_last_coordinate = $exon_first_coordinate + length($reference_sequence_noGaps) - 1;
		toVCF(
			$exon_first_coordinate,
			$reference_alignment,
			$reference_sequence_noGaps,
			\@allele_alignments,
			$fh_VCF,
			"chr6",
			\$ref_sequence,
		);
		
		print {$fh_BED} join("\t", $locus . '-' . $exonID, $exon_first_coordinate + 1, $exon_last_coordinate + 1), "\n";
		# die Dumper($reference_alignment, @allele_alignments);
	}
}
	
print "\n\nProduced files:\n\t$fn_VCF\n\t$fn_BED\n\n";
close($fh_BED);
close($fh_VCF);


sub toVCF
{
	my $alignments_start_0based = shift;
	my $reference_aligned = shift;
	my $reference_unaligned = shift;
	my $sequences_aligned_aref = shift;
	my $fh_out = shift;
	my $referenceSequenceID = shift;
	my $completeReferenceSequence_sref = shift;
	
	die unless(all {length($_) == length($reference_aligned)} @$sequences_aligned_aref);
	
	
	my @startingHaplotypes_pos_2_alleles;
	
	foreach my $alignment (@$sequences_aligned_aref)
	{
		my $alignment_asHash_href = {};

		my $start_pos = 0 - 1;
		my $ref_pos = $start_pos;
		my $running_gaps = 0;

		for(my $i = 0; $i < length($alignment); $i++)
		{
			my $c_ref = substr($reference_aligned, $i, 1);
			my $c_query = substr($alignment, $i, 1);
			if(($c_ref eq '-') or ($c_ref eq '*'))
			{
				die unless($i != 0);
			}
			else
			{
				$ref_pos++;
				die Dumper("Mismatch", $c_ref, $ref_pos, substr($reference_unaligned, $ref_pos, 1)) unless($c_ref eq substr($reference_unaligned, $ref_pos, 1));
			}
			
			$alignment_asHash_href->{$ref_pos} .= $c_query;
		}
		
		die Dumper("End pos mismatch", $ref_pos, length($reference_unaligned)) unless($ref_pos == (length($reference_unaligned)-1));	
		push(@startingHaplotypes_pos_2_alleles, $alignment_asHash_href);
	}
	
	my $N_startingHaplotypes = scalar(@$sequences_aligned_aref);
	
	my $posI_startedOpenAlleles = -1;
	my @runningOpenAlleles;
	$#runningOpenAlleles = $N_startingHaplotypes - 1;
	for(my $haplotypeI = 0; $haplotypeI < $N_startingHaplotypes; $haplotypeI++)
	{
		$runningOpenAlleles[$haplotypeI] = '';
	}
	
	my @runningOpenAlleles_which;
	$#runningOpenAlleles_which = $N_startingHaplotypes - 1;
	
	my $inVariableStretch = 0;
	
	for(my $posI = 0; $posI < length($reference_unaligned); $posI++)
	{
		my $refC = substr($reference_unaligned, $posI, 1);
		my $allReferenceIdentical = 1;
		my @runningOpenAlleles_thisPos;
		$#runningOpenAlleles_thisPos = $N_startingHaplotypes - 1;
		for(my $haplotypeI = 0; $haplotypeI < $N_startingHaplotypes; $haplotypeI++)
		{
			if(exists $startingHaplotypes_pos_2_alleles[$haplotypeI]{$posI})
			{
				my $allele = $startingHaplotypes_pos_2_alleles[$haplotypeI]{$posI};
				if($allele ne $refC)
				{
					$allReferenceIdentical = 0;
				}
				$runningOpenAlleles[$haplotypeI] .= $allele;
				$runningOpenAlleles_thisPos[$haplotypeI] = 1;
			}
		}
		
		my $printedWarning = 0;
		
		for(my $haplotypeI = 0; $haplotypeI < $N_startingHaplotypes; $haplotypeI++)
		{
			if($runningOpenAlleles_which[$haplotypeI] != $runningOpenAlleles_thisPos[$haplotypeI])
			{
				if($runningOpenAlleles_which[$haplotypeI] && (not $runningOpenAlleles_thisPos[$haplotypeI]))
				{			
					if(not $printedWarning)
					{
						warn "(I) $posI @runningOpenAlleles_which @runningOpenAlleles_thisPos";
					}
					$printedWarning = 1;
					$posI_startedOpenAlleles = -1;
				}
				elsif((not $runningOpenAlleles_which[$haplotypeI]) && $runningOpenAlleles_thisPos[$haplotypeI])
				{
					if(! $allReferenceIdentical)
					{
						if(! $printedWarning)
						{
							warn "(II) $posI @runningOpenAlleles_which @runningOpenAlleles_thisPos";
						}
						
						$printedWarning = 1;
						$posI_startedOpenAlleles = -1,
					}					
				}
				else
				{
					die;
				}
			}
		}
		

		if($allReferenceIdentical)
		{
			if($posI_startedOpenAlleles != -1)
			{
				my %runningAlleles;
				my @runningAlleles_vec;

				for(my $haplotypeI = 0; $haplotypeI < $N_startingHaplotypes; $haplotypeI++)
				{
					my $runningAllele_noGaps = substr($runningOpenAlleles[$haplotypeI], 0, length($runningOpenAlleles[$haplotypeI]) - 1);
					$runningAllele_noGaps =~ s/[\-_]//g;
					
					if($runningOpenAlleles_which[$haplotypeI])
					{
						die unless(length($runningAllele_noGaps) > 0);
						die Dumper("Some mismatch", substr($runningOpenAlleles[$haplotypeI], length($runningOpenAlleles[$haplotypeI])-1, 1), $refC) unless(substr($runningOpenAlleles[$haplotypeI], length($runningOpenAlleles[$haplotypeI])-1, 1) eq $refC);
						$runningAlleles{$runningAllele_noGaps}++;
						push(@runningAlleles_vec, $runningAllele_noGaps);
					}
					else
					{
						die unless(length($runningAllele_noGaps) == 0);
					}
				}
				
				my $runningAllele_lastRefPos = $posI - 1;
				my $runningAlleles_refSpan = substr($reference_unaligned, $posI_startedOpenAlleles, $runningAllele_lastRefPos - $posI_startedOpenAlleles + 1);

				if((scalar(keys %runningAlleles) > 1) or ((scalar(keys %runningAlleles) > 0) and ($runningAlleles_vec[0] ne $runningAlleles_refSpan)))
				{
					my @alleles_ordered;
					my %allele_2_i;
					$allele_2_i{$runningAlleles_refSpan} = 0;
					push(@alleles_ordered, $runningAlleles_refSpan);
					
					foreach my $oneAllele (keys %runningAlleles)
					{
						if(not exists $allele_2_i{$oneAllele})
						{
							$allele_2_i{$oneAllele} = scalar(@alleles_ordered);
							push(@alleles_ordered, $oneAllele);
						}
					}
					
					my $removeInitialXBases = 0;
					{
						my $allAllelesSameLength = 1;						
						my $alleleLength = -1;
						my $minAlleleLength = -1;
						foreach my $oneAllele (@alleles_ordered)
						{
							if($alleleLength == -1)
							{
								$alleleLength = length($oneAllele);
								$minAlleleLength = $alleleLength;
							}	
							if($alleleLength != length($oneAllele))
							{
								$allAllelesSameLength = 0;
							}
							if(length($oneAllele) < $minAlleleLength)
							{
								$minAlleleLength = length($oneAllele);
							}
						}
						
						if($minAlleleLength > 1)
						{
							for(my $allelePos = 0; $allelePos < ($minAlleleLength-1); $allelePos++)
							{	
								my $allAllelesIdenticalAtPos = 1;
								my $consensusBase;						
								foreach my $oneAllele (@alleles_ordered)
								{
									die unless($allelePos < length($oneAllele));
									my $thisAlleleBase = substr($oneAllele, $allelePos, 1);
									
									$consensusBase = $thisAlleleBase if(not defined $consensusBase);
									$allAllelesIdenticalAtPos = 0 if($thisAlleleBase ne $consensusBase);
								}
								if($allAllelesIdenticalAtPos)
								{
									$removeInitialXBases = $allelePos+1;
								}
								else
								{
									last;
								}
							}
						}						
					}
					
					my @alleles_ordered_forPrint = map {
						die unless($removeInitialXBases < length($_));
						substr($_, $removeInitialXBases);
					} @alleles_ordered;
					
					my @genotypes = map {
						if(defined $_)
						{
							die unless(defined $allele_2_i{$_});
							$allele_2_i{$_}
						}
						else
						{
							"?"
						}
					} @runningAlleles_vec;
					
					my $startCoordinate_1based = $alignments_start_0based + $posI_startedOpenAlleles+1+$removeInitialXBases;
					
					print {$fh_out} join("\t",
						$referenceSequenceID,
						$startCoordinate_1based,
						".",
						$alleles_ordered_forPrint[0],
						join(",", @alleles_ordered_forPrint[1 .. $#alleles_ordered_forPrint]),
						".",
						"PASS",
						".",
						"GT",
						join("|", @genotypes),
					), "\n";	

					die unless(substr($$completeReferenceSequence_sref, $startCoordinate_1based - 1, length($alleles_ordered_forPrint[0])) eq $alleles_ordered_forPrint[0]);					
				}
			}
			
			for(my $haplotypeI = 0; $haplotypeI < $N_startingHaplotypes; $haplotypeI++)
			{
				if($runningOpenAlleles_thisPos[$haplotypeI])
				{	
					die unless(substr($runningOpenAlleles[$haplotypeI], length($runningOpenAlleles[$haplotypeI]) - 1, 1) eq $refC);
					$runningOpenAlleles[$haplotypeI] = $refC;
				}
				else
				{
					$runningOpenAlleles[$haplotypeI] = "";
				}
			}
			
			$posI_startedOpenAlleles = $posI;
			$inVariableStretch = 0;
		}
		else
		{
			if($posI_startedOpenAlleles != -1)
			{
				$inVariableStretch = 1;
			}
		}
		
		@runningOpenAlleles_which = @runningOpenAlleles_thisPos;
				
	}
}

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	return reverse($kMer);
	return $kMer;
}
	
	
sub readG
{
	my $g_file = qq(reference_HLA_ASM/hla_nom_g.txt);
	open(GFILE, '<', $g_file) or die "Cannot open G groups file $g_file - either you don't have installed the HLA*ASM data package, or you are not calling me from the correct directory (src within the HLA*LA install)";
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

sub read_exon_alignments
{
	my $preserveParts = shift;
	
	my %alignments;
	my @files =  glob ('reference_HLA_ASM/exonAlignments/*_nuc.txt'); 
	foreach my $file (@files)
	{
		die unless($file =~ /.+(\\|\/)(\w+)_nuc.txt$/);
		my $locus = $2;
		unless(($locus =~ /TAP/) or ($locus =~ /MIC/))
		{
			$locus = 'HLA-' . $locus;
		}
		
		print "Reading $locus ...\n";
		my $alignment_href = read_HLA_alignment($file, $preserveParts);
		
		$alignments{$locus} = $alignment_href;
		
	}
	return \%alignments;
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



