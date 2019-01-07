#!/usr/bin/perl

BEGIN {
	use FindBin;
	push(@INC, $FindBin::Bin);
}

use warnings;
use strict;
use FindBin;
use Data::Dumper;
use Storable qw/store retrieve dclone/;
use List::MoreUtils qw/mesh all/;
use List::Util qw/max min/;
use Text::LevenshteinXS qw(distance);
use Getopt::Long;
use Cwd qw/getcwd abs_path/;
use findPath;
use Bio::DB::HTS;

my $scriptPath = $FindBin::Bin;
my $includes=" -I" . join(" -I", @INC);

$| = 1;
my $this_bin_dir = $FindBin::RealBin;

# Test command
# perl HLA-ASM.pl --assembly_fasta /data/projects/phillippy/scratch/trioCanu/NA12878/mhc/nanopolish/dad_methylation/dad_nm_mhc.fasta --sampleID NA12878_nanopolish_dad_methylation --truth reference_HLA_ASM/truth_NA12878.txt

my $sampleID;
my $assembly_fasta;
my $tolerantMissing  = 0;
my $truthFile;
my $samtools_bin;
my $bwa_bin;
my $nucmer_bin;
my $dnadiff_bin;
my $workingDir_param;
GetOptions (
	'assembly_fasta:s' => \$assembly_fasta,
	'sampleID:s' => \$sampleID,
	'tolerantMissing:s' => \$tolerantMissing, 	
	'truthFile:s' => \$truthFile, 	
	'samtools_bin:s' => \$samtools_bin, 	
	'bwa_bin:s' => \$bwa_bin, 	
	'nucmer_bin:s' => \$nucmer_bin, 	
	'dnadiff_bin:s' => \$dnadiff_bin, 	
	'workingDir:s' => \$workingDir_param,	
);

unless($assembly_fasta and $sampleID)
{
	die "Please specify parameters --assembly_fasta and --sampleID. --assembly_fasta should specify a path to a unified FASTA of your assembly with unique contig IDs, and sampleID should be an alphanumeric sample ID";
}

die "File --assembly_fasta not existing" unless(-e $assembly_fasta);
die "Please specify an alphanumeric --sampleID" unless($sampleID =~ /^\w+$/);

if($truthFile)
{
	die "--truthFile $truthFile not existing" unless(-e $truthFile);
	$truthFile = abs_path($truthFile);
}	

$assembly_fasta = abs_path($assembly_fasta);
chdir($this_bin_dir) or die "Cannot chdir into $this_bin_dir";

my $GRCh38_primary = 'reference_HLA_ASM/hg38.primary.fna';
die "Reference reference_HLA_ASM/hg38.primary.fna not existing -- have you downloaded and installed the HLA-ASM data package?" unless(-e $GRCh38_primary);

$samtools_bin = findPath::find_path('samtools_bin', $samtools_bin, 'samtools');
$bwa_bin = findPath::find_path('bwa_bin', $bwa_bin, 'bwa');
$nucmer_bin = findPath::find_path('nucmer_bin', $nucmer_bin, 'nucmer');
$dnadiff_bin = findPath::find_path('dnadiff_bin', $dnadiff_bin, 'dnadiff');

findPath::check_samtools($samtools_bin);



my %classI = map {$_ => 1} qw/HLA-A HLA-B HLA-C/;
my %classII = map {$_ => 1} qw/HLA-DQA1 HLA-DQB1 HLA-DRB1 HLA-DRB3 HLA-DRB4/;


my $working_dir = findPath::get_working_dir($workingDir_param);

unless(-e $working_dir)
{
	mkdir($working_dir) or die "Cannot mkdir $working_dir";
}

my $temp_alignments_dir = $working_dir . '/' . $sampleID . '/temp_alignments_separate/';
unless(-e  $working_dir . '/' . $sampleID)
{
	mkdir( $working_dir . '/' . $sampleID) or die "Cannot mkdir " .  $working_dir . '/' . $sampleID;
}
unless(-e $temp_alignments_dir)
{
	mkdir($temp_alignments_dir) or die "Cannot mkdir $temp_alignments_dir";
}		

mkdir('temp') unless(-d 'temp');
die unless(-d 'temp');

my $fn_assembly_path = $temp_alignments_dir . '/assemblyPath.txt';
if(not -e $fn_assembly_path)
{
	open(ASSEMBLYPATH, '>', $fn_assembly_path) or die "Cannot open $fn_assembly_path";
	print ASSEMBLYPATH $assembly_fasta, "\n";
	close(ASSEMBLYPATH);
}
else
{
	open(ASSEMBLYPATH, '<', $fn_assembly_path) or die "Cannot open $fn_assembly_path";
	my $firstLine = <ASSEMBLYPATH>;
	chomp($firstLine);
	close(ASSEMBLYPATH);
	unless($firstLine eq $assembly_fasta)
	{
		die "The file $fn_assembly_path indicates that you already used this sample ID, but with a different assembly file. As I am not smart enough to figure out which files need to be re-generated under these circumstances, it's best to delete the whole directory $temp_alignments_dir and start over. Alternatively, just use a different sample ID.";
	}
}


my $verbose_reading = 0;

my %relevantRegions = (
	'MHC' => ['CM000668.2', 28510120, 33480577, 'reference_HLA_ASM/MHC8.mfa'],
	# 'KIR' => ['CM000681.2', 54025634, 55084318, 'KIR36_plus.mfa'],
);

my %G_full_alleles;
my %G_mapper_unambigious;
my %G_mapper_multiples;
my %G_to_underlyingAlleles;
my %alleles_to_fullGAmbiguity;
readG();

my %sample_truth;
if($truthFile)
{
	open(TRUTH, '<', $truthFile) or die "Cannot open $truthFile";
	while(<TRUTH>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		next unless($line);
		my @line_fields = split(/\t/, $line);
		
		die "Invalid line in $truthFile - expect tab-separated: HLA-LOCUS allele1 allele2. Example: HLA-A\t01:01:01:01\t11:01:01:01" unless(scalar(@line_fields) == 3);
		die "Invalid line in $truthFile - expect tab-separated: HLA-LOCUS allele1 allele2. Example: HLA-A\t01:01:01:01\t11:01:01:01" unless($line_fields[0] =~ /^HLA-/);
		die "Duplicate locus in $truthFile for $line_fields[0]?" if(exists $sample_truth{$line_fields[0]});
		
		my @alleles_tr;
		foreach my $allele ($line_fields[1], $line_fields[2])
		{
			die "Invalid line in $truthFile - expect tab-separated: HLA-LOCUS allele1 allele2. Example: HLA-A\t01:01:01G\t11:01:01:01.\n\nTyping is ALWAYS carried out a G group resolution." unless($line_fields[1] =~ /\:/);
			my $locus_for_gGroupLookup = $line_fields[0];
			$locus_for_gGroupLookup =~ s/HLA-/HLA/;
			die "Unknown locus $line_fields[0]" unless(exists $G_to_underlyingAlleles{$locus_for_gGroupLookup});
			die "Can't resolve allele $line_fields[0] $allele at G-group resolution" unless(defined $G_to_underlyingAlleles{$locus_for_gGroupLookup}{$allele});		
			my $allele_tr = $G_to_underlyingAlleles{$locus_for_gGroupLookup}{$allele}[0];
			$allele_tr =~ s/^\S+\*//;
			push(@alleles_tr, $allele_tr);
		}
				
		$sample_truth{$line_fields[0]} = \@alleles_tr;
		
	}
	close(TRUTH);
}

my $fn_contigs = $assembly_fasta;
my $contigs_href = readFASTA($fn_contigs, 1);

validate_assembly($assembly_fasta);
my $mcoords = produce_nucmer($assembly_fasta);
my $fn_nucmer_wcoords = $mcoords;

open(DOTPLOTDATA, '>', $temp_alignments_dir . '/dotPlot.tab') or die;
open(EDITDISTANCE, '>', $temp_alignments_dir . '/editDistances.tab') or die;
open(CONTIGLENGTHS, '>', $temp_alignments_dir . '/refContigLengths.tab') or die;
open(GENEPOSITIONS, '>', $temp_alignments_dir . '/genePositions.tab') or die;
open(GENESEQUENCES, '>', $temp_alignments_dir . '/geneSequences.txt') or die;
print DOTPLOTDATA join("\t", "region", "contigID", "refID", "refCoordinate", "contigCoordinate"), "\n";
print EDITDISTANCE join("\t", "region", "contigID", "refID", "editDistance"), "\n";
print CONTIGLENGTHS join("\t", "region", "refID", "length"), "\n";
print GENEPOSITIONS join("\t", "contigID", "gene", "definedFromRefID", "startInContig_0Based", "stopInContig_0Based", "capturedParts"), "\n";

my %results;
foreach my $region (sort keys %relevantRegions)
{
	my $regionInfo = $relevantRegions{$region};
	my $references_href = readFASTA($regionInfo->[3]);
	
	my $exon_allele_alignments_href;
	#my $genomic_allele_alignments_withPARTS_href;
	my %sequenceExtractionCoordinates_start;
	my %sequenceExtractionCoordinates_stop;
	my %sequenceExtractionCoordinates_refDBStrand;	
	my %alleles_2_features;
	
	if($region eq 'MHC')
	{
		my $exon_allele_alignments_href_file = 'temp/_exon_allele_alignments_href';
		#my $genomic_allele_alignments_withPARTS_href_file = 'temp/_genomic_allele_alignments_withPARTS_href';
		
		if(-e $exon_allele_alignments_href_file)
		{
			$exon_allele_alignments_href = retrieve $exon_allele_alignments_href_file;
			#$genomic_allele_alignments_withPARTS_href = retrieve $genomic_allele_alignments_withPARTS_href_file;
		}
		else
		{
			$exon_allele_alignments_href = read_exon_alignments(1);
			store($exon_allele_alignments_href, $exon_allele_alignments_href_file);
			
			#$genomic_allele_alignments_withPARTS_href = read_genomic_alignments(1);
			#store($genomic_allele_alignments_withPARTS_href, $genomic_allele_alignments_withPARTS_href_file);
		}
		
		open(EXONS, '<', $regionInfo->[3] . '.exons') or die;
		my $exons_headerLine = <EXONS>;
		chomp($exons_headerLine);
		my @exons_header_fields = split(/\t/, $exons_headerLine);
		while(<EXONS>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line);
			die unless($#line_fields == $#exons_header_fields);
			my %line = (mesh @exons_header_fields, @line_fields);
						
			push(@{$sequenceExtractionCoordinates_start{$line{refID}}{$line{Start_0based}}}, [$line{geneName}, $line{subComponent}]);
			push(@{$sequenceExtractionCoordinates_stop{$line{refID}}{$line{Stop_0based}}}, [$line{geneName}, $line{subComponent}]);
			
			if(not defined $sequenceExtractionCoordinates_refDBStrand{$line{refID}}{$line{geneName}})
			{
				$sequenceExtractionCoordinates_refDBStrand{$line{refID}}{$line{geneName}} = $line{alleleDatabase_strand};
			}
			die unless($sequenceExtractionCoordinates_refDBStrand{$line{refID}}{$line{geneName}} eq $line{alleleDatabase_strand});
		}	
		close(EXONS);
		
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
		
		foreach my $DRBallele (keys %{$alleles_2_features{'HLA-DRB'}})
		{
			die unless($DRBallele =~ /(DRB\d)/);
			my $realLocus = $1;
			$alleles_2_features{'HLA-'.$realLocus}{$DRBallele} = $alleles_2_features{'HLA-DRB'}{$DRBallele};
		}
		delete $alleles_2_features{'HLA-DRB'};
		
		foreach my $locus (keys %sample_truth)
		{
			foreach my $truthAllele (@{$sample_truth{$locus}})
			{
				my $truth_key = $locus . '*' . $truthAllele;
				$truth_key =~ s/HLA\-//;
				die Dumper("Can't find sequence for $locus $truth_key (specified in file $truthFile)", ["Examples", (keys %{$alleles_2_features{$locus}})[0 .. 9]], [keys %alleles_2_features]) unless(exists $alleles_2_features{$locus}{$truth_key});
			}					
		}
		
		#open(DRB3, '>', 'NA12878_DRB3_sequence.fa') or die;
		#my $DRB3_sequence = $exon_allele_alignments_href->{'HLA-DRB'}{'DRB3*01:01:02:01'};
		#$DRB3_sequence =~ s/-//g;
		#$DRB3_sequence =~ s/_//g;
		#$DRB3_sequence =~ s/\|//g;
		#print DRB3 '>DRB3*01:01:02:01', "\n";
		#print DRB3 $DRB3_sequence, "\n";
		#close(DRB3);
		#exit;
		# die Dumper(keys %{$alleles_2_features{'HLA-DRB3'}});
		# $alleles_2_features{'HLA-DRB1'} = $alleles_2_features{'HLA-DRB'};
	}
	my %overlappingContigs;
	my %overlappingContigs_idty;
	my %overlappingContig_strands;
	my %overlappingContig_strand_called;
	my %strand_per_originalContig;
	my %strand_per_haplotig;
	
	my %alphanumeric_keys;
	foreach my $refContig (keys %$references_href)
	{		
		my $refID_alphanumeric = $refContig;
		$refID_alphanumeric =~ s/\W//g;
		die Dumper("Duplicate ID?", $refID_alphanumeric, \%alphanumeric_keys) if($alphanumeric_keys{$refID_alphanumeric});
		$alphanumeric_keys{$refID_alphanumeric}++;
	}
	
	my %alignments_per_refContig;
	foreach my $refContigID (sort keys %$references_href)
	{
		next if($refContigID eq 'gi|568335997|gb|KI270758.1| Homo sapiens chromosome 6 genomic contig, GRCh38 reference assembly alternate locus group ALTREFLOCI8');
		my $refContigSequence = $references_href->{$refContigID};		
		print CONTIGLENGTHS join("\t", $region, $refContigID, length($refContigSequence)), "\n";
	}
	
		
	if($fn_nucmer_wcoords)
	{
		open(COORDS, '<', $fn_nucmer_wcoords) or die "Cannot open $fn_nucmer_wcoords";
		while(<COORDS>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line);
			
			my $refID = $line_fields[11];
			my $contigID = $line_fields[12];
			my $refStart = $line_fields[0];
			my $refStop = $line_fields[1];
			my $contigStart = $line_fields[2];
			my $contigStop = $line_fields[3];
			my $strand = ($contigStart < $contigStop) ? '+' : '-';
			my $identity = $line_fields[6];
			next unless($refID eq $regionInfo->[0]);
			
			if(overlap($regionInfo->[1], $regionInfo->[2], $refStart, $refStop))
			{
				my $useKey = $contigID;
				if(exists $contigs_href->{$useKey . '_pilon'})
				{
					$useKey = $useKey . '_pilon';
				}
				elsif($contigs_href->{$useKey . '_pilon_pilon'})
				{
					$useKey = $useKey . '_pilon_pilon';
				}
				$overlappingContigs{$useKey}++;
				$overlappingContigs_idty{$useKey} = $identity;
				$overlappingContig_strands{$useKey}{$strand}++;
			} 
		}
		close(COORDS);
		
		print "${region}-overlapping contigs:\n";
		print join("\n", map {' - ' . $_ . ' [' . $overlappingContigs_idty{$_} . '] - ' . (exists $overlappingContig_strands{$_}{'+'} ? $overlappingContig_strands{$_}{'+'} : 0) . '+ / ' . (exists $overlappingContig_strands{$_}{'-'} ? $overlappingContig_strands{$_}{'-'} : 0) . '-'} keys %overlappingContigs), "\n";
			
		foreach my $contigID (keys %overlappingContigs)
		{
			my @strands = keys(%{$overlappingContig_strands{$contigID}});
			@strands = reverse sort {$overlappingContig_strands{$contigID}{$a} <=> $overlappingContig_strands{$contigID}{$b}} @strands;
			my $strand = $strands[0];			
			
			unless(scalar(@strands) == 1)
			{
				warn join("\t", $contigID, @strands), "\n";
				warn "No enquivocal strand call for $contigID - now $strand - , check manually!";
			}
			$strand_per_originalContig{$contigID} = $strand;	
		}
		
		# foreach my $contigID (keys %overlappingContigs)
		# {
			# foreach my $haplotigID (keys %haplotigs)
			# {
				# if(substr($haplotigID, 0, length($contigID)) eq $contigID)
				# {
					# $overlappingContigs{$haplotigID} = 1;
					
					# my $haplotigSequence = $haplotigs{$haplotigID};
					
					# my $jaccard_plus = kMerSharingSubset($contigs_href->{$contigID}, $haplotigSequence);
					# my $jaccard_minus = kMerSharingSubset($contigs_href->{$contigID}, reverseComplement($haplotigSequence));
					
					# print "Contig $contigID [", length($contigs_href->{$contigID}), "], haplotig $haplotigID [", length($haplotigSequence), "], plus $jaccard_plus, minus $jaccard_minus \n";
					
					# my $haplotig_strand_relative_to_originalContig = ($jaccard_plus > $jaccard_minus) ? '+' : '-';
					
					# my $haplotig_strand_absolute;
					# if($haplotig_strand_relative_to_originalContig eq '+')
					# {
						# $haplotig_strand_absolute = $strand_per_originalContig{$contigID};
					# }
					# else
					# {
						# $haplotig_strand_absolute = ($strand_per_originalContig{$contigID} eq '-') ? '+' : '-';
					# }
					# $strand_per_haplotig{$haplotigID} = $haplotig_strand_absolute;
					
					# $contigs_href->{$haplotigID} = $haplotigSequence;
				# }
			# }
		# }
	}
	else
	{
		%overlappingContigs = map {$_ => 1} keys %$contigs_href;
		%strand_per_originalContig = map {$_ => '+'} keys %$contigs_href;
	}


   # todo get rid of!
   # $strand_per_originalContig{'1_tig00000010_B_pilon'} = '-';
   # $strand_per_originalContig{'0_tig00000049_A_pilon'} = '-';
   # $strand_per_originalContig{'1_tig00000003_B_pilon'} = '-';
   # $strand_per_originalContig{'0_tig00000001_A_pilon'} = '-';
   # $strand_per_originalContig{'0_tig00000022_A_pilon'} = '-';
   # $strand_per_originalContig{'1_tig00000012_B_pilon'} = '-';
   # $strand_per_originalContig{'0_tig00000008_A_pilon'} = '-';
		   
	# die Dumper((sort keys %overlappingContigs), \%strand_per_originalContig);
	
	foreach my $contigID (sort keys %overlappingContigs)
	{
		# next unless($contigID eq 'tig00019339_pilon');
		# next unless($contigID eq 'tig00007299_pilon');
		
		unless(exists $contigs_href->{$contigID})
		{
			if($tolerantMissing)
			{
				warn "Contig $contigID undefined?";
				next;
			}
			else
			{
				die "Contig $contigID undefined?";
			}
		}	
		my $contigID_alphanumeric = $contigID;
		$contigID_alphanumeric =~ s/\W//g;
		die if($alphanumeric_keys{$contigID});
		$alphanumeric_keys{$contigID}++;
			
		my $strand;
		# if(not exists $haplotigs{$contigID})
		{
			$strand = $strand_per_originalContig{$contigID};
		}
		# else
		# {
			# $strand = $strand_per_haplotig{$contigID};		
		# }
		die unless(defined $strand);
		
		$overlappingContig_strand_called{$contigID} = $strand;
		
		my $contigSequence = $contigs_href->{$contigID};
		if(defined $strand and ($strand eq '-'))
		{
			$contigSequence = reverseComplement($contigSequence);
		}
		
		my %alignments_per_refContig;
		foreach my $refContigID (sort keys %$references_href)
		{
			next if($refContigID eq 'gi|568335997|gb|KI270758.1| Homo sapiens chromosome 6 genomic contig, GRCh38 reference assembly alternate locus group ALTREFLOCI8');
			my $refContigSequence = $references_href->{$refContigID};
			
			my $refID_alphanumeric = $refContigID;
			$refID_alphanumeric =~ s/\W//g;
		
			my $fn_output =  $temp_alignments_dir . '/' . $refID_alphanumeric . '_' . $contigID_alphanumeric;
			
			my $T = '';
			my $Q = '';
			
			my $TStart_0based;
			my $TStop_0based;
			my $QStart_0based;
			my $QStop_0based;
			my $score;
			my $alignmentStrand;

			my $aligner = 'bwa';
			if($aligner eq 'edlib')
			{
				unless(-e $fn_output)
				{
					writeFASTA('_contig', {contig => $contigSequence});
					writeFASTA('_ref', {ref => $refContigSequence});
					print "Aligning $contigID [${strand}]vs $refContigID \n";
					my $edlib_cmd = qq(/data/projects/phillippy/software/edlib/build/bin/aligner -n 1 -m HW -p _contig _ref > $fn_output);
					system($edlib_cmd);
				}
				
				{
					open(F, '<', $fn_output) or die "Cannot open $fn_output";
					while(<F>)
					{
						my $line = $_;
						if($line =~ /score = (\d+)/)
						{
							die "More than one alignment?" if(defined $score);
							$score = $1;
						}
						else
						{
							if($line =~ /^(T|Q)\: (\w+) \((\d+) - (\d+)\)/)
							{
								my $TQ = $1;
								my $alignment = $2;
								my $start = $3;
								my $stop = $4;
								
								if($TQ eq 'T')
								{
									$T .= $alignment;
									$TStart_0based = $start if(not defined $TStart_0based);
									$TStop_0based = $stop;
								}
								else
								{
									die unless($TQ eq 'Q');
									$Q .= $alignment;
									$QStart_0based = $start if(not defined $QStart_0based);
									$QStop_0based = $stop;							
								}
							}				
						}
					}
					close(F);
				}
			}
			elsif($aligner eq 'nucmer')
			{
				die "Not fully implemented";
				my $alignment_f = $fn_output . '.alignment';

				print "Alignment file $alignment_f\n";
				unless(-e $alignment_f)
				{
					print "Re-computing...\n";
					writeFASTA('_contig', {contig => $contigSequence});
					writeFASTA('_ref', {ref => $refContigSequence});
					print "Aligning $contigID [${strand}]vs $refContigID \n";
					my $mummer_cmd = qq(/data/projects/phillippy/software/MUMmerSensitive/nucmer -p $fn_output _ref _contig);
					system($mummer_cmd) and die "Nucmer command $mummer_cmd failed";
					
					my $delta_file = $fn_output . '.delta';
					my $delta_file_g = $fn_output . '.delta.g';
					my $filter_command = "delta-filter -g $delta_file > $delta_file_g";
					system($filter_command) and die "Nucmer command $mummer_cmd failed";
					
					my $alignment_command = "show-aligns -r $delta_file_g ref contig > $alignment_f";
					system($alignment_command) and die "Nucmer command $alignment_command failed";
				}		
				
					open(ALIGNMENT, '<', $alignment_f) or die "Cannot open $alignment_f";
					while(<ALIGNMENT>)
					{
						my $line = $_;
						chomp($line);
						next unless($line);
						
						my $inAlignment = 0;
						if($line =~ /BEGIN alignment \[ ((\+1)|(-1)) (\d+) - (\d+) \| ((\+1)|(-1)) (\d+) - (\d+) \]/)
						{
							my $plusMinus_1 = $1;
							my $start_ref = $4;
							my $stop_ref = $5;
							
							my $plusMinus_2 = $6;
							my $start_seq = $9;
							my $stop_seq = $10;
							
							# warn Dumper($plusMinus_1, $start_ref, $stop_ref, $plusMinus_2, $start_seq, $stop_seq, $line);
							
							die unless($plusMinus_1 eq '+1');
							die unless($plusMinus_2 eq '+1');
							
							die unless($start_ref < $stop_ref);
							die unless($start_seq < $stop_seq);
							
							$inAlignment = 1;
							
							if(not defined $TStart_0based)
							{
								$TStart_0based = $start_ref - 1;
								$QStart_0based = 0;
								if($start_seq != 1)
								{
									my $missing_gaps = $start_seq - 1;
									$T .= ('-' x $missing_gaps);
									die unless(length($T) == $missing_gaps);
									$Q .= substr($contigSequence, 0, $missing_gaps);
								}
							}
							
							$TStop_0based = $stop_ref - 1;
							
							if(defined $QStop_0based)
							{
								my $noGap_query = $Q;
								$noGap_query =~ s/-//g;
								die unless(length($noGap_query) == ($QStop_0based+1));
								die Dumper($start_seq, $QStop_0based) unless(($start_seq-1) > $QStop_0based);
								my $missing_query_gaps_in_between = ($start_seq-1) - $QStop_0based - 1;
								$Q .= substr($contigSequence, $QStop_0based+1, $missing_query_gaps_in_between);
								$T .= ('-' x $missing_query_gaps_in_between);
								
								my $missing_ref_sequence = ($start_ref-1) - $TStop_0based - 1;
								$T .= substr($refContigSequence, $TStop_0based + 1, $missing_ref_sequence);
								$Q .= ('-' x $missing_ref_sequence);								
							}
							$QStop_0based = $stop_seq - 1;
						}
						elsif($line =~ /END alignment/)
						{
							$inAlignment = 0;
						}
						else
						{
							if($line =~ /^\d+\s+([\w\.]+)/)
							{
								my $ref_seq = $1;
								die "Weird sequence: '$ref_seq'" unless($ref_seq =~ /^[ACGTN\.]+$/i);
								$ref_seq =~ s/\./-/g;
								
								my $nextLine = <ALIGNMENT>;
								chomp($nextLine);
								die "Weird nextLine $nextLine" unless($nextLine =~ /^\d+\s+([\w\.]+)/);
								
								my $query_seq = $1;
								die unless($query_seq =~ /^[ACGTN\.]+$/i);
								$query_seq =~ s/\./-/g;				

								$Q .= $query_seq;
								$T .= $ref_seq;
							}
						}		
					}
					
					if($QStop_0based != (length($contigSequence)-1))
					{
						my $missing_gaps_end = (length($contigSequence) - 1) - $QStop_0based;
						die unless($missing_gaps_end > 0);
						$Q .= substr($contigSequence, $QStop_0based+1);
						$T .= ('-' x $missing_gaps_end);						
					}
					
					my $Q_noGap = $Q;
					my $T_noGap = $T;
					$Q_noGap =~ s/-//g;
					$T_noGap =~ s/-//g;
					die unless(uc($Q_noGap) eq uc($contigSequence));
					die unless(length($T_noGap) == ($TStop_0based - $TStart_0based + 1));
					die unless(uc($T_noGap) eq uc(substr($refContigSequence, $TStart_0based, $TStop_0based - $TStart_0based + 1)));
				}
				elsif($aligner eq 'bwa')
				{

					my $fn_output_refined = $fn_output . '.2';
					unless(-e $fn_output)				
					{
						writeFASTA($temp_alignments_dir . '/_ref', {ref => $refContigSequence});									
						writeFASTA($temp_alignments_dir . '/_contig', {contig => $contigSequence});
											
						my $cmd_align = qq(perl $includes $scriptPath/globalAlignment.pl --reference ${temp_alignments_dir}/_ref --query ${temp_alignments_dir}/_contig --output $fn_output --samtools_bin $samtools_bin --bwa_bin $bwa_bin);
						print "Command: $cmd_align\n";
						system($cmd_align) and die "Command $cmd_align failed";
					}
					
					unless(1 or (-e $fn_output_refined))
					{
						
						open(ALIGNMENT, '<', $fn_output) or die "Cannot open $fn_output";
						my $l1 = <ALIGNMENT>; chomp($l1);
						my $l2 = <ALIGNMENT>; chomp($l2);
						my $l3 = <ALIGNMENT>; chomp($l3);
						close(ALIGNMENT);
						
						my $fn_forMAFFT = $fn_output_refined . '.forMAFFT';
						my $fn_fromMAFFT = $fn_output_refined . '.fromMAFFT';
						
						my $alignmentLength_before = length($l2);
						
						$l2 =~ s/-//g;
						$l3 =~ s/-//g;
						
						writeFASTA($fn_forMAFFT, {'ref'=> $l2, 'contig' => $l3});
						
						my $mafft_bin = '/data/projects/phillippy/projects/rDNA/mafft/mafft-7.273-with-extensions/install/bin/mafft';
						my $cmd_mafft = qq($mafft_bin --retree 2 --maxiterate 2 $fn_forMAFFT > $fn_fromMAFFT);
						if(system($cmd_mafft))
						{
							die "Mafft command $cmd_mafft failed";
						}	
						
						my $fromMAFFT_href = readFASTA($fn_fromMAFFT);
						die unless($fromMAFFT_href->{ref});
						die unless($fromMAFFT_href->{contig});
						
						print "MAFFT: alignment length went from ", $alignmentLength_before, " to ", length($fromMAFFT_href->{ref}), "\n";	

						open(OUTPUT, '>', $fn_output_refined) or die;
						print OUTPUT $l1, "\n";
						print OUTPUT $fromMAFFT_href->{ref}, "\n";
						print OUTPUT $fromMAFFT_href->{contig}, "\n";
						close(OUTPUT);
					}  
					
					open(ALIGNMENT, '<', $fn_output) or die "Cannot open $fn_output";
					my $l1 = <ALIGNMENT>; chomp($l1);
					my $l2 = <ALIGNMENT>; chomp($l2);
					my $l3 = <ALIGNMENT>; chomp($l3);
					die "Weird first line $l1 in $fn_output" unless($l1 =~ /(\d+) (\d+)-([\-\d]+) (\+|-)(\d+)-(\d+)/);
					$score = $1;
					$TStart_0based = $2;
					$TStop_0based = $3; 
					$alignmentStrand = $4;
					$QStart_0based = $5;
					$QStop_0based = $6;
					close(ALIGNMENT);
					# warn Dumper("Weird strand $alignmentStrand from file $fn_output", $contigID, $overlappingContig_strands{$contigID}, $refContigID) unless($alignmentStrand eq '+');
					
					$T = $l2;
					$Q = $l3;
					my $Q_noGap = $Q;
					my $T_noGap = $T;
					$Q_noGap =~ s/-//g;
					$T_noGap =~ s/-//g;
					
					if($alignmentStrand eq '+')
					{
						die Dumper("Alignment contig sequence not equal to original contig sequence", length($Q_noGap), length($contigSequence)) unless(uc($Q_noGap) eq uc($contigSequence));
					}
					else
					{
						die Dumper("Alignment contig sequence not equal to original contig sequence", length($Q_noGap), length($contigSequence)) unless(uc($Q_noGap) eq reverseComplement(uc($contigSequence)));
					}
					die unless(length($T_noGap) == ($TStop_0based - $TStart_0based + 1));
					
					$T_noGap = uc($T_noGap);
					my $T_noGaps_equiv = uc(substr($refContigSequence, $TStart_0based, $TStop_0based - $TStart_0based + 1));
					for(my $posI = 0; $posI < length($T_noGap); $posI++)
					{
						my $c1 = substr($T_noGap, $posI, 1);
						my $c2 = substr($T_noGaps_equiv, $posI, 1);
						next if(($c1 eq 'N') or ($c2 eq 'N'));
						die unless($c1 eq $c2);
						
					}
					# die Dumper($TStart_0based, $TStop_0based, length($T_noGap), ()) unless(uc($T_noGap) eq uc();					
				}
								
				die unless(length($T) == length($Q));
				die unless((defined $TStart_0based) and (defined $QStart_0based));
				die unless($QStart_0based == 0);
				die unless($QStop_0based == (length($contigSequence)-1));
				
				$alignments_per_refContig{$refContigID} = [$score, $T, $Q, [$TStart_0based, $TStop_0based], [$QStart_0based, $QStop_0based], $alignmentStrand];
				
				print EDITDISTANCE join("\t", $region, $contigID, $refContigID, $score), "\n";
				
			}

			my %alignment_plots_per_refContig;
			print "Alignment statistics for $contigID - length ", length($contigSequence), "\n";
			foreach my $refContigID (keys %alignments_per_refContig)
			{
				print " - ", $refContigID, ": strand $alignments_per_refContig{$refContigID}[5], edit distance ", $alignments_per_refContig{$refContigID}[0], "\n";
				
				my $runningPos_T = $alignments_per_refContig{$refContigID}[3][0];
				my $runningPos_Q = $alignments_per_refContig{$refContigID}[4][0];
				my $L = length($alignments_per_refContig{$refContigID}[1]);
				for(my $i = 0; $i < $L; $i++)
				{ 
					my $T_c = substr($alignments_per_refContig{$refContigID}[1], $i, 1);
					my $Q_c = substr($alignments_per_refContig{$refContigID}[2], $i, 1);
					
					if(($i % 100) == 0)
					{
						push(@{$alignment_plots_per_refContig{$refContigID}}, [$runningPos_T, $runningPos_Q]);
						
						print DOTPLOTDATA join("\t", $region, $contigID, $refContigID, $runningPos_T, $runningPos_Q), "\n";
					}
					
					if(($T_c ne '-'))
					{
						$runningPos_T++;
					}
					if(($Q_c ne '-'))
					{
						$runningPos_Q++;
					}				
				}
				
			}
			
			my @sorted_refContigIDs = sort {$alignments_per_refContig{$a}[0] <=> $alignments_per_refContig{$b}[0]} keys %alignments_per_refContig;
			if(scalar(@sorted_refContigIDs) > 1)
			{
				# die Dumper($sorted_refContigIDs[0], $alignments_per_refContig{$sorted_refContigIDs[0]}[0]);
				die unless($alignments_per_refContig{$sorted_refContigIDs[0]}[0] <= $alignments_per_refContig{$sorted_refContigIDs[1]}[0]);
			}	

			
			foreach my $useIndexForTyping (0)
			{
				my %extractedSequences;
				
				my $refContigID = $sorted_refContigIDs[$useIndexForTyping];
				
				print "Extract sequences from alignment to $refContigID\n";
				
				next unless($alignments_per_refContig{$refContigID}[5] eq '+');
				my $alignment_running_reference = $alignments_per_refContig{$refContigID}[3][0] - 1;
				my $L = length($alignments_per_refContig{$refContigID}[1]);
				my $alignment_running_query = -1;
				my %found_starts_alignmentCoordinates;
				my %found_stops_alignmentCoordinates;			
				for(my $i = 0; $i < $L; $i++)
				{ 
					my $T_c = substr($alignments_per_refContig{$refContigID}[1], $i, 1);
					my $Q_c = substr($alignments_per_refContig{$refContigID}[2], $i, 1);
					
					if($T_c ne '-')
					{
						$alignment_running_reference++;
					}

					# if($Q_c ne '-')
					# {
						# die unless($Q_c =~ /^[ACGTN]$/);
						# $alignment_running_query++;
						# die unless(substr($contigSequence, $alignment_running_query, 1) eq $Q_c);
					# }
					
					if($T_c ne '-')
					{
						if(exists $sequenceExtractionCoordinates_start{$refContigID}{$alignment_running_reference})
						{
							my @ks = map {join(';-;', @$_)} @{$sequenceExtractionCoordinates_start{$refContigID}{$alignment_running_reference}};
							foreach my $k (@ks)
							{
								die Dumper("Alread found starting coordinates", $k) if(defined $found_starts_alignmentCoordinates{$k});
								$found_starts_alignmentCoordinates{$k} = $i;
								# die if(defined $found_starts_queryCoordinates{$k});						
								# $found_starts_queryCoordinates{$k} = $alignment_running_query;
							}
						}
					}
					
					if($T_c ne '-')
					{					
						if(exists $sequenceExtractionCoordinates_stop{$refContigID}{$alignment_running_reference})
						{
							my @ks = map {join(';-;', @$_)} @{$sequenceExtractionCoordinates_stop{$refContigID}{$alignment_running_reference}};
							foreach my $k (@ks)
							{
								die if(defined $found_stops_alignmentCoordinates{$k});
								$found_stops_alignmentCoordinates{$k} = $i;
								# die if(defined $found_stops_queryCoordinates{$k});
								# $found_stops_queryCoordinates{$k} = $alignment_running_query;
								
								# if(exists $found_starts_alignmentCoordinates{$k})
								# {
									# my $alignmentSeq = substr($alignments_per_refContig{$refContigID}[2], $found_starts_alignmentCoordinates{$k}, $i - $found_starts_alignmentCoordinates{$k} + 1);
									# my $contigSeq = substr($contigSequence, $found_starts_queryCoordinates{$k}, $alignment_running_query - $found_starts_queryCoordinates{$k} + 1);
									# print $found_starts_alignmentCoordinates{$k}, " - ", $i, "\n", $found_starts_queryCoordinates{$k}, " - ", $alignment_running_query, "\n", $alignmentSeq, "\n", $contigSeq, "\n\n";
								# }
							}
						}	
					}					
				}
								
				my %joint_keys = map {$_ => 1} (keys %found_starts_alignmentCoordinates, keys %found_stops_alignmentCoordinates);
				my %positions_per_gene;
				my %positions_per_gene_inAlignment;
				my %components_per_gene;
				foreach my $k (keys %joint_keys)
				{
					if((exists $found_starts_alignmentCoordinates{$k}) and (exists $found_stops_alignmentCoordinates{$k}))
					{
						my $l_extract_contig_alignment = $found_stops_alignmentCoordinates{$k} - $found_starts_alignmentCoordinates{$k} + 1;
						my $extracted_contig_alignment = substr($alignments_per_refContig{$refContigID}[2], $found_starts_alignmentCoordinates{$k}, $l_extract_contig_alignment);
						$extractedSequences{$k} = $extracted_contig_alignment;
						

						
						my $extracted_contig_alignment_noGap = $extracted_contig_alignment;
						$extracted_contig_alignment_noGap =~ s/-//g;
						
						
						die unless($k =~ /^(.+);-;(.+?)$/);
						my $gene = $1;
						my $component = $2;
						die unless(defined $component);
						push(@{$components_per_gene{$gene}}, $component);
						
						#print "EXTRACT $gene $component from $found_starts_alignmentCoordinates{$k} length $l_extract_contig_alignment:\n";
						#print substr($alignments_per_refContig{$refContigID}[1], $found_starts_alignmentCoordinates{$k}, $l_extract_contig_alignment), "\n";
						#print substr($alignments_per_refContig{$refContigID}[2], $found_starts_alignmentCoordinates{$k}, $l_extract_contig_alignment), "\n";
						#print "\n";
						
						if(length($extracted_contig_alignment_noGap))
						{
							my $skipped_query_alignment = '';
							if($found_starts_alignmentCoordinates{$k})
							{
								$skipped_query_alignment = substr($alignments_per_refContig{$refContigID}[2], 0, $found_starts_alignmentCoordinates{$k});
							}
							$skipped_query_alignment =~ s/-//g;

							my $first_queryCoordinate = length($skipped_query_alignment);
							my $last_queryCoordinate = $first_queryCoordinate + length($extracted_contig_alignment_noGap) - 1;
							
							my $corresponding_raw_contig_sequence = substr($contigSequence, $first_queryCoordinate, $last_queryCoordinate - $first_queryCoordinate + 1);
							die Dumper("Mismatch") unless($corresponding_raw_contig_sequence eq $extracted_contig_alignment_noGap);
							
							push(@{$positions_per_gene{$gene}}, $first_queryCoordinate, $last_queryCoordinate);
							push(@{$positions_per_gene_inAlignment{$gene}}, $found_starts_alignmentCoordinates{$k}, $found_stops_alignmentCoordinates{$k});
						}
					}
				}	
			   # if(1 == 0)	
# {
				# foreach my $gene ('HLA-DQB1')
				# {
					# my @gene_positions = sort @{$positions_per_gene_inAlignment{$gene}};
					# my $firstGenePos = $gene_positions[0];
					# my $lastGenePos = $gene_positions[$#gene_positions];
					
					# print $gene, " from $firstGenePos length ", $lastGenePos - $firstGenePos + 1, "\n";
					# print substr($alignments_per_refContig{$refContigID}[1], $firstGenePos, $lastGenePos - $firstGenePos + 1), "\n";
					# print substr($alignments_per_refContig{$refContigID}[2], $firstGenePos, $lastGenePos - $firstGenePos + 1), "\n";
					
					# print substr($alignments_per_refContig{$refContigID}[1], $firstGenePos - 4000, 4000), "\n";
					# print substr($alignments_per_refContig{$refContigID}[2], $firstGenePos - 4000, 4000), "\n";
					

					# print substr($alignments_per_refContig{$refContigID}[1], $lastGenePos + 4000, 4000), "\n";
					# print substr($alignments_per_refContig{$refContigID}[2], $lastGenePos + 4000, 4000), "\n";
					
					# exit;
				# }
# }				
				
				foreach my $gene (sort {min($positions_per_gene{$a}) <=> min($positions_per_gene{$b})} keys %positions_per_gene)
				{
					my @positions = @{$positions_per_gene{$gene}};
					my @positions_alignment = @{$positions_per_gene_inAlignment{$gene}};
					if(@positions)
					{
						die unless(defined $components_per_gene{$gene});
						print GENEPOSITIONS join("\t", 
							$contigID,
							$gene,
							$refContigID, 
							min(@positions),
							max(@positions),
							join(';', sort @{$components_per_gene{$gene}})
						), "\n";
						
						my $minPos_onContig = min(@positions);
						my $maxPos_onContig = max(@positions);
						
						my $minPos_alignment = min(@positions_alignment);
						my $maxPos_alignment = max(@positions_alignment);			

						if(0 and $gene eq 'HLA-DRB3')
						{
							my $alignmentL = $maxPos_alignment - $minPos_alignment + 1;
							print $gene, " complete sequence from $minPos_onContig - $maxPos_onContig on contig, $minPos_alignment - $maxPos_alignment in alignment:\n";
							print substr($alignments_per_refContig{$refContigID}[1], $minPos_alignment, $alignmentL), "\n";
							print substr($alignments_per_refContig{$refContigID}[2], $minPos_alignment, $alignmentL), "\n";
							
							print "Pre:\n";
							print substr($alignments_per_refContig{$refContigID}[1], $minPos_alignment-5000, 5000), "\n";
							print substr($alignments_per_refContig{$refContigID}[2], $minPos_alignment-5000, 5000), "\n";

							print "Pre 2:\n";
							print substr($alignments_per_refContig{$refContigID}[1], $minPos_alignment-10000, 5000), "\n";
							print substr($alignments_per_refContig{$refContigID}[2], $minPos_alignment-10000, 5000), "\n";
														
							
							print "Pre 3:\n";
							print substr($alignments_per_refContig{$refContigID}[1], $minPos_alignment-15000, 5000), "\n";
							print substr($alignments_per_refContig{$refContigID}[2], $minPos_alignment-15000, 5000), "\n";																			
						}
					}
				}
				
				
				if($region eq 'MHC')
				{
					print "MHC analysis $contigID\n";
					# 6-digit G typing
					LOCUS: foreach my $locus (keys %classI, keys %classII)
					{
						die Dumper("Undefined locus $locus", [keys %alleles_2_features]) unless(defined $alleles_2_features{$locus});
						unless(defined $sequenceExtractionCoordinates_refDBStrand{$refContigID}{$locus})
						{
							# warn "No strand information for $refContigID, $locus";
							next LOCUS;
						}
						
						my $strand_relative_to_refContigID = $sequenceExtractionCoordinates_refDBStrand{$refContigID}{$locus};
						
					
						my %present_components = map {my @c = split(/;-;/, $_); die unless(scalar(@c) == 2); if($c[0] eq $locus){($c[1] => 1)}else{()}} keys %extractedSequences;
						# print Dumper($locus, \%present_components, \%extractedSequences);

						my @compatible_alleles;
						my @compare = ($classI{$locus}) ? ('exon2', 'exon3') : ('exon2');
						#print "Locus $locus, comparing with ", scalar(keys %{$alleles_2_features{$locus}}), " alleles.\n";
						
						my %equivalence_classes;
						
						ALLELE: foreach my $allele (keys %{$alleles_2_features{$locus}})
						{
							my $combined_allele_seq = join(";", map {$alleles_2_features{$locus}{$allele}{$_}} @compare);
							# die $combined_allele_seq;
							push(@{$equivalence_classes{$combined_allele_seq}}, $allele);
						}
						
						my %editDistance_perClass;
						EQUIVALENCECLASS: foreach my $class (keys %equivalence_classes)
						{
							my $representativeAllele = $equivalence_classes{$class}[0];
							my $compatible = 1;
							foreach my $compare (@compare)
							{
								die unless(defined $alleles_2_features{$locus}{$representativeAllele}{$compare});
								my $extractedSequenceKey = join(';-;', $locus, $compare);
								unless(exists $extractedSequences{$extractedSequenceKey})
								{
									# print "Gene $locus missing sequence $extractedSequenceKey\n";
									next LOCUS;
								}
															
								my $sequenceForComparison = ($strand_relative_to_refContigID eq '+') ? $extractedSequences{$extractedSequenceKey} : reverseComplement($extractedSequences{$extractedSequenceKey});
								my $sequenceForComparison_withGaps = $extractedSequences{$extractedSequenceKey};
								$sequenceForComparison =~ s/-//g;
								$sequenceForComparison = uc($sequenceForComparison);
								die Dumper("Weird sequence for comparison", $sequenceForComparison, $extractedSequences{$extractedSequenceKey}) unless($sequenceForComparison =~ /^[ACGTN]*$/);
								
								my $editDistance = distance($sequenceForComparison, $alleles_2_features{$locus}{$representativeAllele}{$compare});
								if(1 or length($sequenceForComparison) == length($alleles_2_features{$locus}{$representativeAllele}{$compare}))
								{
									#print $locus, " ", $representativeAllele, " ", $compare, " ", $editDistance, "\n";
									#print $sequenceForComparison, "\n", $alleles_2_features{$locus}{$representativeAllele}{$compare}, "\n\n";
								}
								
								if(0 and $locus eq 'HLA-DRB1')
								{
									print $locus, " ", $representativeAllele, " ", $compare, " ", $editDistance, "\n";
									print "Assembly: ", $sequenceForComparison, "\n", "Ass. with gaps, no revCmp: ", $sequenceForComparison_withGaps, "\n", $representativeAllele, ';', $compare, ": ", $alleles_2_features{$locus}{$representativeAllele}{$compare}, "\n\n";									
								}
								
								$editDistance_perClass{$class} += $editDistance;
								unless($sequenceForComparison eq $alleles_2_features{$locus}{$representativeAllele}{$compare})
								{
									#print "Mismatch ", length($sequenceForComparison), " vs ", length($alleles_2_features{$locus}{$allele}{$compare}), "\n";
									$compatible = 0;
									# next ALLELE;
								}
							}
							if($compatible)
							{
								push(@compatible_alleles, $class);
							}
						}
						
						my @classKeys_sorted = sort {$editDistance_perClass{$a} <=> $editDistance_perClass{$b}} keys %editDistance_perClass;
						my $minimum_edit_distance = $editDistance_perClass{$classKeys_sorted[0]};
						my @classKeys_best = grep {$editDistance_perClass{$_} == $minimum_edit_distance} keys %editDistance_perClass;
						my @classKeys_best_alleles = map {@{$equivalence_classes{$_}}} @classKeys_best;
						
						my $min_distance_call_to_truth = 'NA';
						my $min_distance_assembly_truth = 'NA';
						my $min_distance_call_to_truth_which = 'NA';
						my $min_distance_assembly_truth_which = 'NA';
												
						my $presentComponents = join(', ', sort keys %present_components);
						my $useComponents = join(', ', sort @compare);
						print "\tEdit distances for $locus (from $refContigID - present components: ", $presentComponents , ")\n";
						
						if(exists $sample_truth{$locus})
						{		
							my @truth_alleles = @{$sample_truth{$locus}};

							my %distances_to_truth_calls_perAllele;
							my %distances_to_truth_assembly_perAllele;
							
							die "No sample truth data for $locus" unless(exists $sample_truth{$locus});
							
							(my $locus_noHLA = $locus) =~ s/HLA-//;
							
							foreach my $truthAllele (@truth_alleles)
							{
								my $truthAllele_G = (allele_list_to_G([$truthAllele], $locus))[0];
							
								my $truth_key = $locus . '*' . $truthAllele;
								$truth_key =~ s/HLA\-//;
								die "Can't find sequence for $truth_key" unless(exists $alleles_2_features{$locus}{$truth_key});
															
								foreach my $classKey (@classKeys_best)
								{
									my @alleles = @{$equivalence_classes{$classKey}};

									# distance to call
									{
										my $distance = 0;
										foreach my $compare (@compare)
										{
											my $truth_sequence = $alleles_2_features{$locus}{$truth_key}{$compare};
											my $called_sequence = $alleles_2_features{$locus}{$alleles[0]}{$compare};
											die unless($truth_sequence and $called_sequence);
											$distance += distance($truth_sequence, $called_sequence);
										}
										
										my $comparisonAllele_G = (allele_list_to_G([$alleles[0]], $locus))[0];
										
										my $key = $comparisonAllele_G . '_v/s_' . $truthAllele_G;
										$distances_to_truth_calls_perAllele{$key} = $distance;
									}
								}
								
								# distance to assembly
								{
									my $distance = 0;
									foreach my $compare (@compare)
									{
										my $truth_sequence = $alleles_2_features{$locus}{$truth_key}{$compare};
										die unless($truth_sequence);
										
										my $extractedSequenceKey = join(';-;', $locus, $compare);
										die unless(exists $extractedSequences{$extractedSequenceKey});
																	
										my $sequenceForComparison = ($strand_relative_to_refContigID eq '+') ? $extractedSequences{$extractedSequenceKey} : reverseComplement($extractedSequences{$extractedSequenceKey});
										$sequenceForComparison =~ s/-//g; 
										$sequenceForComparison = uc($sequenceForComparison);
										die Dumper("Weird sequence for comparison", $sequenceForComparison, $extractedSequences{$extractedSequenceKey}) unless($sequenceForComparison =~ /^[ACGTN]*$/);
									
										# warn "No sequenceForComparison for $locus" unless($truth_sequence and $sequenceForComparison);
										$distance += distance($truth_sequence, $sequenceForComparison);
									} 
									
									$distances_to_truth_assembly_perAllele{$truthAllele} = $distance;
								}

								
								my @distances_to_truth_calls_perAllele_keys_sorted = sort {$distances_to_truth_calls_perAllele{$a} <=> $distances_to_truth_calls_perAllele{$b}} keys %distances_to_truth_calls_perAllele;
								my @distances_to_truth_assembly_perAllele_keys_sorted = sort {$distances_to_truth_assembly_perAllele{$a} <=> $distances_to_truth_assembly_perAllele{$b}} keys %distances_to_truth_assembly_perAllele;
							
								my $distances_to_truth_calls_perAllele_keys_best = $distances_to_truth_calls_perAllele_keys_sorted[0];
								my $distances_to_truth_assembly_perAllele_keys_best = $distances_to_truth_assembly_perAllele_keys_sorted[0];
								
								my @distances_to_truth_calls_perAllele_keys_allBest = grep {$distances_to_truth_calls_perAllele{$_} == $distances_to_truth_calls_perAllele{$distances_to_truth_calls_perAllele_keys_best}} keys %distances_to_truth_calls_perAllele; 
								my @distances_to_truth_assembly_perAllele_keys_allBest = grep {$distances_to_truth_assembly_perAllele{$_} == $distances_to_truth_assembly_perAllele{$distances_to_truth_assembly_perAllele_keys_best}} keys %distances_to_truth_assembly_perAllele;
								
								my %_classKeys_best_alleles = map {$_ => 1} @classKeys_best_alleles;
								$min_distance_call_to_truth = $distances_to_truth_calls_perAllele{$distances_to_truth_calls_perAllele_keys_best};
								$min_distance_assembly_truth = $distances_to_truth_assembly_perAllele{$distances_to_truth_assembly_perAllele_keys_best};
								$min_distance_call_to_truth_which = join(' ', @distances_to_truth_calls_perAllele_keys_allBest, $locus);
								$min_distance_assembly_truth_which = join(' ', allele_list_to_G(\@distances_to_truth_assembly_perAllele_keys_allBest, $locus));
								# die unless(all {$_classKeys_best_alleles{$_}} @distances_to_truth_calls_perAllele_keys_allBest);
								# die unless(all {$_classKeys_best_alleles{$_}} @distances_to_truth_assembly_perAllele_keys_allBest);
							}
							
							print "\t\t", "Lowest edit distances", "\n", "\t\t\t", join(' ', @classKeys_best_alleles), "\n\t\t\tassembly to call: ", $editDistance_perClass{$classKeys_sorted[0]}, "; call to truth: $min_distance_call_to_truth; assembly to truth: $min_distance_assembly_truth\n";
							
							die unless($editDistance_perClass{$classKeys_sorted[0]} == $minimum_edit_distance);
							
							foreach my $compare (@compare)
							{
								my $extractedSequenceKey = join(';-;', $locus, $compare);
								print GENESEQUENCES join("\t", $locus, 'extracted', $extractedSequenceKey, $extractedSequences{$extractedSequenceKey}), "\n";
							}
							
							# print "Compatible alleles for $locus: ", join(", ", @compatible_alleles), "\n";
						}
						
						$results{$contigID}{$locus} = [join(' ', allele_list_to_G(\@classKeys_best_alleles)), $useComponents, $minimum_edit_distance, $min_distance_assembly_truth, $min_distance_call_to_truth, $min_distance_assembly_truth_which, $min_distance_call_to_truth_which];
						
					}	
				}
			}
		}
		
}

close(DOTPLOTDATA);


my %results_by_locus;
foreach my $contigID (sort {$a cmp $b} keys %results)
{
	foreach my $locus (sort {$a cmp $b} keys %{$results{$contigID}})
	{
		$results_by_locus{$locus}{$contigID} = $results{$contigID}{$locus};
	}
}

my $summaryFn = $working_dir . '/' . $sampleID . '/summary.txt';

open(SUMMARY, '>', $summaryFn) or die "Cannot open $summaryFn";
print SUMMARY join("\t", qw/contigID locus calledGenotypes components editDistance_calledGenotypes_assembly minEditDistance_assembly_truth minEditDistance_calledGenotype_truth minEditDistance_assembly_truth_whichAlleles minEditDistance_calledGenotype_truth_whichAlleles/), "\n";
foreach my $locus (sort {$a cmp $b} keys %results_by_locus)
{
	foreach my $contigID (sort {$a cmp $b} keys %{$results_by_locus{$locus}})
	{
		print SUMMARY join("\t", $contigID, $locus, @{$results_by_locus{$locus}{$contigID}}), "\n";
	}
}

close(SUMMARY);

print "\n\nDone. Summary output is in $summaryFn\n\n";

		
		
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

sub writeFASTA
{
	my $file = shift;
	# print "Writing $file\n";
	my $href = shift;
	open(F, '>', $file) or die "Cannot open $file";
	foreach my $key (keys %$href)
	{
		my $seq = $href->{$key};
		print F '>', $key, "\n";
		# print "\t", $key, "\t", length($seq), "\n";
		while($seq)
		{
			my $toPrint;
			if(length($seq) > 50)
			{
				$toPrint = substr($seq, 0, 50);
				substr($seq, 0, 50) = '';
			}
			else
			{
				$toPrint = $seq;
				$seq = '';
			}	
			print F $toPrint, "\n";
		}
	}
	close(F);	
}

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	return reverse($kMer);
	return $kMer;
}

sub overlap
{
	my $eins_left = shift;
	my $eins_right = shift;
	my $zwei_left = shift;
	my $zwei_right = shift;
	
	die unless($eins_left < $eins_right);
	die unless($zwei_left < $zwei_right);
	my $eins_length = $eins_right - $eins_left + 1;
	my $zwei_length = $zwei_right - $zwei_left + 1;   
	die unless($eins_length > 0);
	die unless($zwei_length > 0);
		
	my $O;
	if($eins_length > $zwei_length)
	{
		$O = overlap_eins_larger($eins_left, $eins_right, $zwei_left, $zwei_right);
		die unless($O >= 0);
		
	}
	else
	{
		$O = overlap_eins_larger($zwei_left, $zwei_right, $eins_left, $eins_right);
		die unless($O >= 0);	
	}
	
	return $O;
}

sub overlap_eins_larger
{
	my $eins_left = shift;
	my $eins_right = shift;
	my $zwei_left = shift;
	my $zwei_right = shift;
		
	die unless($eins_left < $eins_right);
	die unless($zwei_left < $zwei_right);
	my $eins_length = $eins_right - $eins_left + 1;
	my $zwei_length = $zwei_right - $zwei_left + 1;
	die unless($eins_length > 0);
	die unless($zwei_length > 0);
	die Dumper("L error", $eins_left, $eins_right, $zwei_left, $zwei_right, $eins_length, $zwei_length) unless($eins_length >= $zwei_length);
	
	if(($eins_left <= $zwei_left) and ($eins_right >= $zwei_right))
	{
		return $zwei_length;
	}
	elsif(($zwei_left >= $eins_left) and ($zwei_left <= $eins_right))
	{
		die unless($zwei_right >= $eins_right);
		return ($eins_right - $zwei_left + 1);
	}
	elsif(($zwei_right >= $eins_left) and ($zwei_right <= $eins_right))
	{
		die unless($zwei_left <= $eins_left);
		return ($zwei_right - $eins_left + 1);
	}
	else
	{
		return 0;
	}
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

sub read_genomic_alignments
{
	my $preserveParts = shift;
	
	my %alignments;
	my @files =  glob ('reference_HLA_ASM/genomicAlignments/*_gen.txt'); 
	foreach my $file (@files)
	{
		die unless($file =~ /.+(\\|\/)(\w+)_gen.txt$/);
		my $locus = $2;
		unless(($locus =~ /TAP/) or ($locus =~ /MIC/))
		{
			$locus = 'HLA-' . $locus;
		}
		
		print "Reading $locus ...\n" if ($verbose_reading);
		
		my $alignment_href = read_HLA_alignment($file, $preserveParts);
		
		$alignments{$locus} = $alignment_href;
		
	}
	return \%alignments;
}


sub read_HLA_alignment
{
	my $file = shift;
	my $preserveParts = shift;
	
	my %sequences;
	
	my $inContent = 0;
	my $justSawgDNA = 0;
	my $firstSequenceAftergDNA;
	
	open(ALIGNMENT, '<', $file) or die "Cannot open $file";
	while(<ALIGNMENT>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		if($line =~ /[cg]DNA/)
		{
			$inContent = 1;
			die if($justSawgDNA);
			$justSawgDNA = 1;
			next;
		}
		next if($line =~ /^[\s\|]+$/);
		next if($line =~ /Please see/);
		next if($line =~ /AA codon/);
		
		if($inContent)
		{
			die Dumper("Invalid line", $file, $., $line) unless($line =~ /^\s*(\w+?)\*(\S+?) (.+)$/);
			my $locus = $1;
			my $allele = $2;
			my $sequence = $3;
			
			$sequence =~ s/\s//g;
			unless($preserveParts)
			{
				$sequence =~ s/\|//g;
			}
			
			if($sequence =~ /\-/)
			{
				die unless($firstSequenceAftergDNA);
				die unless(length($sequence) == length($firstSequenceAftergDNA));
				die if($justSawgDNA);
				for(my $i = 0; $i < length($firstSequenceAftergDNA); $i++)
				{
					my $c_seq = substr($sequence, $i, 1);
					if($c_seq eq '-')
					{
						substr($sequence, $i, 1) = substr($firstSequenceAftergDNA, $i, 1);
					}
				}
			}
			
			if(length($sequence) > 0)
			{
				if($preserveParts)
				{
					unless($sequence =~ /^[ACGTN_\*\.\|]+$/i)
					{
						my $sequence_rem = $sequence;
						$sequence_rem =~ s/[ACGTN_\*\.]+//gi;
						die "Weird characters in alignment $file for $locus $allele  -- $sequence_rem\n$sequence";
					}
				}
				else
				{
					unless($sequence =~ /^[ACGTN_\*\.]+$/i)
					{
						my $sequence_rem = $sequence;
						$sequence_rem =~ s/[ACGTN_\*\.]+//gi;
						die "Weird characters in alignment $file for $locus $allele  -- $sequence_rem\n$sequence";
					}				
				}
				
				$sequence =~ s/\./_/g;
				$sequence =~ s/\*/N/g;
				$sequences{$locus.'*'.$allele} .= $sequence;
			}
			
			if($justSawgDNA)
			{
				$justSawgDNA = 0;
				$firstSequenceAftergDNA = $sequence;
			}
		}
	}
	close(ALIGNMENT);
	die unless($inContent);
	
	my $expected_length;
	foreach my $seqID (keys %sequences)
	{
		my $sequence = $sequences{$seqID};
		if(not defined $expected_length)
		{
			$expected_length = length($sequence);
		}
		unless(length($sequence) == $expected_length)
		{
			if($seqID =~ /N/)
			{
				if(length($sequence) > $expected_length)
				{
					warn "Truncate allele sequence for $seqID";
					$sequences{$seqID} = substr($sequence, 0, $expected_length);
				}
				else
				{
					die;
				}
			}
			else
			{
				die "Length problem with $seqID: ".length($sequence). ' vs ' . $expected_length;
			}
		}
		
		if($preserveParts)
		{
			die unless($sequence =~ /^[ACGTN_\|]+$/);
		}
		else
		{
			die unless($sequence =~ /^[ACGTN_]+$/);
		}
	}
	
	if($verbose_reading)
	{
		print "\tread_HLA_alignment(..) for $file: ", scalar(keys %sequences), " sequences of length $expected_length.\n";
	}
	
	# print sequences for debugging
	if(1 == 0)
	{
		open(T, '>', '_a.txt') or die;
		foreach my $seqID (keys %sequences)
		{
			my $sequence = $sequences{$seqID};
			print T $sequence, "\n";
		}
		close(T);
	}
		
	return \%sequences;
}

sub kMerSharing
{
	my $seq1 = shift;
	my $seq2 = shift;
	my $k = shift;
	$k = 15 unless(defined $k);
	
	my %kmers_1 = map {$_ => 1} kmers($seq1, $k);
	my %kmers_2 = map {$_ => 1} kmers($seq2, $k);
	my %kmers_union = map {$_ => 1} ((keys %kmers_1), (keys %kmers_2));
	my $n_intersection = scalar(grep {exists $kmers_1{$_}} keys %kmers_2);
	my $n_union = scalar(keys %kmers_union);
	if($n_union == 0)
	{
		return 0;
	}
	else
	{
		return $n_intersection / $n_union;
	}
}

sub kMerSharingSubset
{
	my $seq1 = shift;
	my $seq2 = shift;
	my $k = shift;
	$k = 15 unless(defined $k);
	
	my %kmers_1 = map {$_ => 1} kmers($seq1, $k);
	my %kmers_2 = map {$_ => 1} kmers($seq2, $k);
	my $n_intersection = scalar(grep {exists $kmers_1{$_}} keys %kmers_2);
	if(scalar(keys %kmers_2) == 0)
	{
		return 0;
	}
	else
	{
		return $n_intersection / scalar(keys %kmers_2);
	}
}

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

sub globalAlignment
{
	my $S1 = shift;
	my $S2 = shift;
	my $temp_dir = shift;
	
	my $fn_ref = $temp_dir . '/gA_ref';
	my $fn_query = $temp_dir . '/gA_contig';
	
	writeFASTA($fn_ref, {ref => $S1});
	writeFASTA($fn_query, {contig => $S2});

	my $fn_output = $temp_dir . '/gA_alignment';
	my $edlib_cmd = qq(/data/projects/phillippy/software/edlib/build/bin/aligner -n 1 -m NW -p $fn_query $fn_ref > $fn_output);
	system($edlib_cmd) and die "Edlb error with command $edlib_cmd";

	my $TStart_0based;
	my $TStop_0based;
	my $QStart_0based;
	my $QStop_0based;
	my $T = '';
	my $Q = '';
	my $score;
	
	open(F, '<', $fn_output) or die "Cannot open $fn_output";
	while(<F>)
	{
		my $line = $_;
		if($line =~ /score = (\d+)/)
		{
			die "More than one alignment?" if(defined $score);
			$score = $1;
		}
		else
		{
			if($line =~ /^(T|Q)\: ([\w\-]+) \((\d+) - (\d+)\)/)
			{
				my $TQ = $1;
				my $alignment = $2;
				my $start = $3;
				my $stop = $4;
				
				if($TQ eq 'T')
				{
					$T .= $alignment;
					$TStart_0based = $start if(not defined $TStart_0based);
					$TStop_0based = $stop;
				}
				else
				{
					die unless($TQ eq 'Q');
					$Q .= $alignment;
					$QStart_0based = $start if(not defined $QStart_0based);
					$QStop_0based = $stop;							
				}
			}				
		}
	}
	close(F);
					
	$T =~ s/-/_/g;
	$Q =~ s/-/_/g;
	
	(my $T_noGaps = $T) =~ s/_//g;
	(my $Q_noGaps = $Q) =~ s/_//g;
	
	die unless($T_noGaps eq $S1);
	die unless($Q_noGaps eq $S2);
	
	return ($T, $Q);
}	

sub validate_assembly
{
	my $assembly_fn = shift;
	
	my $assembly_length = 0;
	my %contigIDs;
	open(F, '<', $assembly_fn) or die "Cannot open $assembly_fn";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		if(substr($line, 0, 1) eq '>')
		{
			my $contigID = substr($line, 1);
			$contigID =~ s/\s.+//;
			die "Contig ID $contigID used at least twice in $assembly_fn" if(exists $contigIDs{$contigID});
			$contigIDs{$contigID}++;
		}
		else
		{
			$assembly_length += length($line); 
		}
	}
	close(F);
	
	print "\nAssembly file $assembly_fn validated - total assembly length approximately ", sprintf("%.2f", $assembly_length / 1e9), "gb.\n\n";
}

sub produce_nucmer
{
	my $assembly_fn = shift;
	my $expected_nucmer_results = $assembly_fn . '.againstPrimary.delta.diff.mcoords';
	if(-e $expected_nucmer_results)
	{
		print "File $expected_nucmer_results found, will not re-run nucmer\n";
		return $expected_nucmer_results;
	}
	else
	{
		my $prefix1 = $assembly_fn . '.againstPrimary';
		my $cmd_nucmer = qq($nucmer_bin -c 1000 -l 20 -p $prefix1 $GRCh38_primary $assembly_fn);
		system($cmd_nucmer) and die "ncucmer command $cmd_nucmer failed";
		
		my $delta_file = $prefix1 . '.delta';
		die "Expected file $delta_file not existing" unless (-e $delta_file);
		
		my $prefix2 = $delta_file . '.diff';
		
		my $cmd_dnadiff = qq($dnadiff_bin -d $delta_file -p $prefix2);
		system($cmd_dnadiff) and die "dnadiff command $cmd_dnadiff failed";
	
		die "Expected file $expected_nucmer_results not existing" unless(-e $expected_nucmer_results);
		
		return $expected_nucmer_results;
	}
}

sub allele_list_to_G
{
	my $alleleList_aref = shift;
	my $locus = shift;
	
	die unless(keys %G_mapper_unambigious);
	my %g;
	foreach my $allele (@$alleleList_aref)
	{
		my $g = $G_mapper_unambigious{$allele};
		if(not defined $g)
		{
			if(defined $locus)
			{
				(my $locus_noHLA = $locus) =~ s/HLA-//;
				$g = $G_mapper_unambigious{$locus_noHLA . '*' . $allele};			
			}
			die Dumper("Can't translate allele $allele", ["Examples in the right format", [(keys %G_mapper_unambigious)[0 .. 9]]], $locus) if(not defined $g);
		}
			
		$g{$g}++;
			
	}
	return sort keys %g;
}

sub readG
{
	my $g_file = qq(reference_HLA_ASM/hla_nom_g.txt);
	open(GFILE, '<', $g_file) or die "Cannot open $g_file";
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

