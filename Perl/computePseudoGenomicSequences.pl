# Todo
# - Output PGF sequence with genes N-ed out
# - in output, replace "pgf" with "A*ref pgf:1--11111
# - include padding sequence in index() for coordinate determination
# - for non-pgf alleles, find exact reference sequence from qbl/mcf/mann etc.
# - do something about DRB5?
# - check that the new read-2-hla file is there and correc
# - filter reads using the new file
# 
# - sample some random reeads from PGF
# - increase padding
# - randomize read order
# - extra-fast read extraction

BEGIN {
	use FindBin;
	push(@INC, $FindBin::Bin);
	push(@INC, $FindBin::RealBin);
}

use warnings;
use strict;
use FindBin;
use File::Spec;
use Getopt::Long;
use Data::Dumper; 
use Sys::Hostname;
use Cwd qw/getcwd abs_path/;
use List::MoreUtils qw/mesh/;
use List::Util qw/all sum/;
use VCFFunctions;

my $graph;
my $n_padding = 1000;
my %padding_otherRefs = ('HLA-DRB3' => 'cox', 'HLA-DRB4' => 'ssto');
my %nonPGF_alleles = ('HLA-DRB3' => 'DRB3*01:01:02:01', 'HLA-DRB4' => 'DRB4*01:03:01:03');
GetOptions (
	'graph:s' => \$graph,
);

my $full_graph_dir = $FindBin::RealBin . '/../../graphs/' . $graph;
die "Directory $full_graph_dir not existing" unless(-e $full_graph_dir);

my $DRB5_PGF_allele = 'DRB5*01:01:01:01';
my $DRB5_firstPos_PGF;
my $DRB5_lastPos_PGF;
my $DRB5_alignment_href = read_HLA_alignment('static/DRB5_gen-1.txt', 0);
my $DRB5_alignment_referenceOrientation = { map {$_ => reverseComplement($DRB5_alignment_href->{$_})} keys %$DRB5_alignment_href };
my %ref_sequences;
{
	my $sequences_file = $full_graph_dir . '/sequences.txt';
	open(SEQUENCES, '<', $sequences_file) or die "Cannot open $sequences_file";
	my $header_sequences = <SEQUENCES>;
	$header_sequences =~ s/[\n\r]//g;	
	chomp($header_sequences);
	my @header_sequences_fields = split(/\t/, $header_sequences);
	die unless(($header_sequences_fields[0] eq 'SequenceID') and ($header_sequences_fields[1] eq 'Name'));
	my @additionalRefSequences;
	while(<SEQUENCES>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		next unless($line);
		my @line_fields = split(/\t/, $line, -1);
		die Dumper("Field mismatch") unless(scalar(@line_fields) == scalar(@header_sequences_fields));
		my %line = mesh @header_sequences_fields, @line_fields;

		sif($line{SequenceID} <= 9)
		{
			die if($line{Chr});
			push(@additionalRefSequences, [$line{Name}, $line{FASTAID}]);
		}
	}
	close(SEQUENCES);
	
	my $extendedReferenceGenome = $full_graph_dir . '/extendedReferenceGenome/extendedReferenceGenome.fa';
	my $referenceGenome_href = readFASTA($extendedReferenceGenome, 1);
	
	$ref_sequences{'pgf'} = VCFFunctions::getPGFSequence($full_graph_dir);

	
	foreach my $additionalRefSeqID (@additionalRefSequences)
	{
		die "Sequence $additionalRefSeqID->[1] missing" unless(exists $referenceGenome_href->{$additionalRefSeqID->[1]});
		$ref_sequences{$additionalRefSeqID->[0]} = $referenceGenome_href->{$additionalRefSeqID->[1]};
	}

	print "Acquired PGF sequence of length ", length($ref_sequences{'pgf'}), "\n";
	my $DRB5_ref_sequence = uc($DRB5_alignment_referenceOrientation->{$DRB5_PGF_allele});
	$DRB5_ref_sequence =~ s/_//g;
	
	$DRB5_firstPos_PGF = index(uc($ref_sequences{'pgf'}), $DRB5_ref_sequence);
	die if($DRB5_firstPos_PGF == -1);
	$DRB5_lastPos_PGF = $DRB5_firstPos_PGF + length($DRB5_ref_sequence) - 1;
}

my @files_in_order;

my %files_per_gene;
my %files_beforeAfter_gene;
my $fn_graph_segments = $full_graph_dir . '/PRG/segments.txt';
die "Missing file: $fn_graph_segments" unless(-e $fn_graph_segments);
open(SEGMENTS, '<', $fn_graph_segments) or die "Cannot open $fn_graph_segments";
my $lastLine_gene;
while(<SEGMENTS>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/[\n\r]//g;
	next unless($line);
	
	my $filename = $line;
	
	my $geneID;
	if($line =~ /\d+_gene_([\w\-]+)_\d+_((intron)|(exon))/)
	{
		$geneID = $1;
		push(@{$files_per_gene{$geneID}}, $filename);
		unless(defined $lastLine_gene)
		{	
			die if(defined $files_beforeAfter_gene{$geneID}[0]);
			$files_beforeAfter_gene{$geneID}[0] = $files_in_order[$#files_in_order];
		}
	}
	else
	{
		$geneID = undef;
	}
	
	push(@files_in_order, $filename);
	if((not defined $geneID) and (defined $lastLine_gene))
	{
		die if(defined $files_beforeAfter_gene{$lastLine_gene}[1]);
		$files_beforeAfter_gene{$lastLine_gene}[1] = $filename;
	}
	$lastLine_gene = $geneID;
	
}
close(SEGMENTS);

my $pgf_reconstruction_aref = VCFFunctions::find_reconstruction(\@files_in_order, $ref_sequences{'pgf'});
my $pgf_reconstruction_href = {mesh @files_in_order, @$pgf_reconstruction_aref};

my $pgf_sequence_with_N = '';
{
	my $pgf_sequence_noN_check = '';
	my $pgf_withGaps = '';
	my $lastLine_gene;
	my $padding_seq = 'N' x $n_padding;
	die unless(length($padding_seq) == $n_padding);
	
	foreach my $file (@files_in_order)
	{
		my $pgf_sequence_id = $pgf_reconstruction_href->{$file};
		die unless(defined $pgf_sequence_id);
		
		my $pgf_sequence_bit = acquireSequenceFromFile($full_graph_dir . '/PRG/' . $file, $pgf_sequence_id);
		$pgf_withGaps .= $pgf_sequence_bit;
		
		my $geneID;
		(my $pgf_sequence_bit_noGaps = $pgf_sequence_bit) =~ s/[\-\_]//g;
		my $pgf_sequence_bit_noGaps_withN = $pgf_sequence_bit_noGaps;
		if($file =~ /\d+_gene_([\w\-]+)_\d+_((intron)|(exon))/)
		{
			$geneID = $1;
			unless(defined $lastLine_gene)
			{	
				die unless(length($padding_seq) == $n_padding);
				die unless(length($pgf_sequence_with_N) >= $n_padding);				
				substr($pgf_sequence_with_N, length($pgf_sequence_with_N) - $n_padding, $n_padding) = $padding_seq;
			}
			$pgf_sequence_bit_noGaps_withN =~ s/./N/g;
		}
		else
		{
			$geneID = undef;
		}
		
		if((not defined $geneID) and (defined $lastLine_gene) and (length($pgf_sequence_bit_noGaps_withN)))
		{	
			die unless(length($padding_seq) == $n_padding);	
			die Dumper("Length padding problem", $file, length($pgf_sequence_bit_noGaps_withN)) unless(length($pgf_sequence_bit_noGaps_withN) >= $n_padding);
			substr($pgf_sequence_bit_noGaps_withN, 0, $n_padding) = $padding_seq;
		}
		
		$pgf_sequence_with_N .= $pgf_sequence_bit_noGaps_withN;
		$pgf_sequence_noN_check .= $pgf_sequence_bit_noGaps;
		die Dumper("Length difference (I)", $file, length($pgf_sequence_with_N), length($pgf_sequence_noN_check)) unless(length($pgf_sequence_with_N) == length($pgf_sequence_noN_check));
		$lastLine_gene = $geneID; # if (length($pgf_sequence_bit_noGaps_withN));
	}	
	die unless($pgf_sequence_noN_check eq $ref_sequences{'pgf'});
	die Dumper("Length difference (II)", length($pgf_sequence_with_N), length($pgf_sequence_noN_check), ($pgf_sequence_with_N =~ /^ACGTN+$/) ? 1 : 0, ($pgf_sequence_noN_check =~ /^ACGTN+$/) ? 1 : 0) unless(length($pgf_sequence_with_N) == length($pgf_sequence_noN_check));
	
	my $DRB_length = $DRB5_lastPos_PGF - $DRB5_firstPos_PGF + 1;
	my $DRB5_N = 'N' x ($DRB_length + 2 * $n_padding);
	die unless(length($DRB5_N) == ($DRB_length + 2 * $n_padding));
	substr($pgf_sequence_with_N, $DRB5_firstPos_PGF - $n_padding, $DRB_length + 2 * $n_padding) = $DRB5_N;
	die unless(length($pgf_sequence_with_N) == length($pgf_sequence_noN_check));
	
	my @pgf_rawCharacter_2_alignment;
	my $pgf_reachedCharacter = -1;
	for(my $i = 0; $i < length($pgf_withGaps); $i++)
	{
		my $alignment_character = substr($pgf_withGaps, $i, 1);
		my $isGap = (($alignment_character eq '-') or ($alignment_character eq '_'));
		unless($isGap)
		{
			$pgf_reachedCharacter++;
			$pgf_rawCharacter_2_alignment[$pgf_reachedCharacter] = $i;
		}
	}	
	die unless(scalar(@pgf_rawCharacter_2_alignment) == length($pgf_sequence_noN_check));
	die unless(all {defined $_} @pgf_rawCharacter_2_alignment);
	
	print "DRB5 in graph alignment coordinates: ", $pgf_rawCharacter_2_alignment[$DRB5_firstPos_PGF], " - ", $pgf_rawCharacter_2_alignment[$DRB5_lastPos_PGF], "\n";
}

my $output_dir = $full_graph_dir . '/pseudoGenomic_fullLengthMapping';
unless(-d $output_dir)
{
	mkdir($output_dir) or die "Cannot mkdir $output_dir";
}

my $fn_out_PGF_with_Ns = $output_dir . '/PGF_with_Ns.fa';
open(PGFWITHNS, '>', $fn_out_PGF_with_Ns) or die "Cannot open $fn_out_PGF_with_Ns";
print PGFWITHNS '>pgf_Ns', "\n", $pgf_sequence_with_N, "\n";

close(PGFWITHNS);

my @genes = sort keys %files_per_gene;
# @genes = ('HLA-DRA');

my %sequences_per_gene;
foreach my $gene (@genes)
{
	(my $gene_noHLA = $gene) =~ s/HLA\-//;
	
	my %all_alleles;
	my %sequences_per_gene_per_file;
	my $ref_allele;
	foreach my $file (@{$files_per_gene{$gene}})
	{
		# die unless($file =~ /\d+_gene_([\w\-]+)_\d+_((intron)|(exon))/);
		# my $geneID = $1;

		my $full_file_path = $full_graph_dir . '/PRG/' . $file;
		open(ONESEGMENT, '<', $full_file_path) or die "Cannot open $full_file_path";
		my $header_segment = <ONESEGMENT>;
		chomp($header_segment);
		$header_segment =~ s/[\n\r]//g;
		my @header_fields = split(/ /, $header_segment);
		die unless($header_fields[0] eq 'IndividualID');
		die unless(exists $pgf_reconstruction_href->{$file});
		my $ref_allele_sequence;
		while(<ONESEGMENT>)
		{
			my $line = $_;
			chomp($line);
			$line =~ s/[\n\r]//g;
			next unless($line); 
			my @line_fields = split(/ /, $line); 
			die unless(scalar(@line_fields) == scalar(@header_fields));
			my $hla_allele = shift(@line_fields);
			my $hla_allele_sequence = join('', @line_fields);
			die "Weird characters in HLA sequence: $hla_allele_sequence" unless($hla_allele_sequence =~ /^[ACGTN_]+$/);
			$all_alleles{$hla_allele}++;
			$sequences_per_gene_per_file{$file}{$hla_allele} = $hla_allele_sequence;	
			if($pgf_reconstruction_href->{$file} eq $hla_allele)
			{
				$ref_allele_sequence = $hla_allele_sequence;
			}
		}
		close(ONESEGMENT);
		
		die unless(defined $ref_allele_sequence);
		$ref_allele .= $ref_allele_sequence;
	}

	print "Gene $gene, know ", scalar(keys %all_alleles), " alleles.\n";
	foreach my $file (@{$files_per_gene{$gene}})
	{	
		#print "Gene $gene, file $file, have ", scalar(keys %{$sequences_per_gene_per_file{$file}}), " alleles.\n";
		
		my $allele_length = length((values %{$sequences_per_gene_per_file{$file}})[0]);
		die unless(all {length($_) == $allele_length} values %{$sequences_per_gene_per_file{$file}}); 
		
		my $augmented_sequences = 0;
		my $sequence_to_augment = 'N' x $allele_length;
		die unless(length($sequence_to_augment) == $allele_length);
		foreach my $allele (keys %all_alleles)
		{
			unless(exists $sequences_per_gene_per_file{$file}{$allele})
			{
				$sequences_per_gene_per_file{$file}{$allele} = $sequence_to_augment;
				$augmented_sequences++;
			}
		}
		
		#print "\tCarried out $augmented_sequences augmentations\n";
	}

	
	
	
	my %combined_allele_sequences;
	my $combined_allele_length;
	foreach my $allele (keys %all_alleles)
	{
		my $combined_allele_sequence = '';
		foreach my $file (@{$files_per_gene{$gene}})
		{
			die unless(exists $sequences_per_gene_per_file{$file}{$allele});
			$combined_allele_sequence .= $sequences_per_gene_per_file{$file}{$allele};
		}
		
		$combined_allele_length = length($combined_allele_sequence) unless(defined $combined_allele_length);
		die unless(length($combined_allele_sequence) == $combined_allele_length);
		
		$combined_allele_sequences{$allele} = $combined_allele_sequence;
	}
	die unless(length($ref_allele) == $combined_allele_length);
	
	my $whichRef;
	if($ref_allele =~ /^[\-_]+$/)
	{
		die unless(exists $nonPGF_alleles{$gene});
		die unless(exists $combined_allele_sequences{$nonPGF_alleles{$gene}});
		die unless(exists $padding_otherRefs{$gene});
		$ref_allele = $combined_allele_sequences{$nonPGF_alleles{$gene}};
		die unless($ref_allele !~ /N/i);
		$whichRef = $padding_otherRefs{$gene};
	}
	else
	{
		$whichRef = 'pgf';
	}
	
	my @full_genomic_allele_IDs = grep {$combined_allele_sequences{$_} !~ /N/i} keys %combined_allele_sequences;  
	my %alleleID_is_full_genomic = map {$_ => 1} @full_genomic_allele_IDs;
	print "\t ... of which ", scalar(@full_genomic_allele_IDs), " are full genomic sequences\n";
	die unless(scalar(@full_genomic_allele_IDs));
	
	my %allele_closest_full_genomic_allele;
	my %combined_allele_sequences_fullyGenomic;
	foreach my $allele (keys %all_alleles)
	{
		my $allele_sequence = $combined_allele_sequences{$allele};

		if($alleleID_is_full_genomic{$allele})
		{
			$combined_allele_sequences_fullyGenomic{$allele} = $allele_sequence;
		}
		else
		{
			my @allele_sequence_nonN_indices = grep {substr($allele_sequence, $_, 1) ne 'N'} (0 .. (length($allele_sequence) - 1)); 
			if(scalar(@allele_sequence_nonN_indices) == 0)
			{
				# warn "Allele $allele seems to be consistently N - $allele_sequence";
			}
			my %distances_to_full_genomic_alleles;
			foreach my $genomic_allele_ID (@full_genomic_allele_IDs)
			{
				my $genome_allele_sequence = $combined_allele_sequences{$genomic_allele_ID};
				my $distance = sum(map {(substr($allele_sequence, $_, 1) ne substr($genome_allele_sequence, $_, 1)) ? 1 : 0} @allele_sequence_nonN_indices);
				die if((not defined $distance) and scalar(@allele_sequence_nonN_indices));
				$distance = 0 unless(defined $distance);
				$distances_to_full_genomic_alleles{$genomic_allele_ID} = $distance;
			}
			die Dumper("Missing distances for some genomic allele IDs") unless(all {exists $distances_to_full_genomic_alleles{$_}} @full_genomic_allele_IDs);
			my @possible_other_alleles_sorted = sort { $distances_to_full_genomic_alleles{$a} <=> $distances_to_full_genomic_alleles{$b} } @full_genomic_allele_IDs;
			my $closest_allele = $possible_other_alleles_sorted[0];
			
			# print "\t Allele $allele, ", scalar(@allele_sequence_nonN_indices), " used positions, closest allele: ", $closest_allele, ", distance $distances_to_full_genomic_alleles{$closest_allele}\n";
			
			$allele_closest_full_genomic_allele{$allele} = $closest_allele;
						
			my $allele_fullyGenomicSequence = join('', 
				map {
					(substr($allele_sequence, $_, 1) ne 'N') ? (substr($allele_sequence, $_, 1)) : (substr($combined_allele_sequences{$closest_allele}, $_, 1))
				}
				(0 .. (length($allele_sequence) - 1))
			);
			
			die if($allele_fullyGenomicSequence =~ /N/i);

			$combined_allele_sequences_fullyGenomic{$allele} = $allele_fullyGenomicSequence;
		}
	}
	
	(my $ref_allele_noGaps = $ref_allele) =~ s/_//g;
	die unless(length($ref_allele_noGaps));	
	
	my ($paddingLeft, $paddingRight) = getPadding_pre_post($gene);	
	my $ref_allele_noGaps_withPadding = $paddingLeft . $ref_allele_noGaps . $paddingRight;
	
	my $whichRef_coordinate;	
	{
		my @positions_ref_allele_noGaps = find_all($ref_sequences{$whichRef}, $ref_allele_noGaps_withPadding);
		unless(scalar(@positions_ref_allele_noGaps) == 1)
		{
			die "Failed to find unambiguous position for allele sequence within $whichRef";
		}
		$whichRef_coordinate = $positions_ref_allele_noGaps[0];
	}	
	
	my $output_fn_alignments = $output_dir . '/alignments_' . $gene . '.fa';
	my $output_fn_raw = $output_dir . '/raw_' . $gene . '.fa';
	
	open(OUTPUT_ALIGNMENTS, '>', $output_fn_alignments) or die "Cannot open $output_fn_alignments";
	open(OUTPUT_RAW, '>', $output_fn_raw) or die "Cannot open $output_fn_raw";
	print OUTPUT_ALIGNMENTS '>', $gene_noHLA . '*ref', ' ', $whichRef . ':' . $whichRef_coordinate, "\n", $paddingLeft . $ref_allele . $paddingRight, "\n";
	print OUTPUT_RAW '>', $gene_noHLA . '*ref', ' ', $whichRef . ':' . $whichRef_coordinate, "\n", $ref_allele_noGaps_withPadding, "\n";
	
	foreach my $alleleID (keys %combined_allele_sequences_fullyGenomic)
	{
		print OUTPUT_ALIGNMENTS '>', $alleleID, "\n";
		print OUTPUT_RAW '>', $alleleID, "\n";
		
		print OUTPUT_ALIGNMENTS $paddingLeft . $combined_allele_sequences_fullyGenomic{$alleleID} . $paddingRight, "\n";
		(my $alleleSequence_noGaps = $combined_allele_sequences_fullyGenomic{$alleleID}) =~ s/_//g;
		print OUTPUT_RAW $paddingLeft . $alleleSequence_noGaps . $paddingRight, "\n";		
	}
	close(OUTPUT_ALIGNMENTS);
	close(OUTPUT_RAW);
}

# DRB5
{
	my $gene = 'HLA-DRB5';
	(my $gene_noHLA = $gene) =~ s/HLA\-//;

	my $ref_allele = uc($DRB5_alignment_referenceOrientation->{$DRB5_PGF_allele});
	(my $ref_allele_noGaps = $ref_allele) =~ s/_//g;
	die unless(length($ref_allele_noGaps));	
	
	my ($paddingLeft, $paddingRight);
	$paddingLeft = substr($ref_sequences{'pgf'}, $DRB5_firstPos_PGF - $n_padding, $n_padding);
	$paddingRight = substr($ref_sequences{'pgf'}, $DRB5_lastPos_PGF + 1, $n_padding);
	my $ref_allele_noGaps_withPadding = $paddingLeft . $ref_allele_noGaps . $paddingRight;
	
	my $whichRef = 'pgf';
	my $whichRef_coordinate;	
	{
		my @positions_ref_allele_noGaps = find_all($ref_sequences{$whichRef}, $ref_allele_noGaps_withPadding);
		unless(scalar(@positions_ref_allele_noGaps) == 1)
		{
			die "Failed to find unambiguous position for allele sequence within $whichRef";
		}
		$whichRef_coordinate = $positions_ref_allele_noGaps[0];
	}	
	
	my $output_fn_alignments = $output_dir . '/alignments_' . $gene . '.fa';
	my $output_fn_raw = $output_dir . '/raw_' . $gene . '.fa';
	
	open(OUTPUT_ALIGNMENTS, '>', $output_fn_alignments) or die "Cannot open $output_fn_alignments";
	open(OUTPUT_RAW, '>', $output_fn_raw) or die "Cannot open $output_fn_raw";
	print OUTPUT_ALIGNMENTS '>', $gene_noHLA . '*ref', ' ', $whichRef . ':' . $whichRef_coordinate, "\n", $paddingLeft . $ref_allele . $paddingRight, "\n";
	print OUTPUT_RAW '>', $gene_noHLA . '*ref', ' ', $whichRef . ':' . $whichRef_coordinate, "\n", $ref_allele_noGaps_withPadding, "\n";
	
	foreach my $alleleID (keys %$DRB5_alignment_href)
	{
		print OUTPUT_ALIGNMENTS '>', $alleleID, "\n";
		print OUTPUT_RAW '>', $alleleID, "\n";
		
		print OUTPUT_ALIGNMENTS $paddingLeft . $DRB5_alignment_href->{$alleleID} . $paddingRight, "\n";
		(my $alleleSequence_noGaps = $DRB5_alignment_href->{$alleleID}) =~ s/_//g;
		print OUTPUT_RAW $paddingLeft . $alleleSequence_noGaps . $paddingRight, "\n";		
	}
	close(OUTPUT_ALIGNMENTS);
	close(OUTPUT_RAW);
}

sub getPadding_pre_post
{
	my $geneID = shift;
	
	die unless(defined $files_beforeAfter_gene{$geneID});
	
	my $sequenceForPadding = (exists $padding_otherRefs{$geneID}) ? $padding_otherRefs{$geneID} : 'pgf';
	
	my $sequence_pre = acquireSequenceFromFile($full_graph_dir . '/PRG/' . $files_beforeAfter_gene{$geneID}[0], $sequenceForPadding);
	$sequence_pre =~ s/[\_\-]//g;
	die unless($sequence_pre);
	die unless(length($sequence_pre) >= $n_padding);
	
	my $sequence_post = acquireSequenceFromFile($full_graph_dir . '/PRG/' . $files_beforeAfter_gene{$geneID}[1], $sequenceForPadding);
	$sequence_post =~ s/[\_\-]//g;
	die unless($sequence_post);
	die unless(length($sequence_post) >= $n_padding);
	
	my $padding_left = substr($sequence_pre, length($sequence_pre) - $n_padding, $n_padding);
	die unless(length($padding_left) == $n_padding);
	
	my $padding_right = substr($sequence_post, 0, $n_padding);
	die unless(length($padding_right) == $n_padding);
	
	return ($padding_left, $padding_right);
}

sub acquireSequenceFromFile
{
	my $file = shift;
	my $whichSequence = shift;

	open(ONESEGMENT, '<', $file) or die "Cannot open $file";
	my $header_segment = <ONESEGMENT>;
	chomp($header_segment);
	$header_segment =~ s/[\n\r]//g;
	my @header_fields = split(/ /, $header_segment);
	die unless($header_fields[0] eq 'IndividualID');
	
	my $acquiredSequence;
	while(<ONESEGMENT>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		next unless($line); 
		my @line_fields = split(/ /, $line); 
		die unless(scalar(@line_fields) == scalar(@header_fields));
		my $sequenceID = shift(@line_fields);
		my $sequence = join('', @line_fields);
		if($sequenceID eq $whichSequence)
		{
			$acquiredSequence = $sequence;
		}
	}
	close(ONESEGMENT);
	
	unless(defined $acquiredSequence)
	{
		die "Could not find sequence $whichSequence in file $file";
	}
	
	return $acquiredSequence;
		
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

sub find_all
{
	my $findIn = shift;
	my $findWhat = shift;
	die unless(length($findIn) > length($findWhat));
	my @results;
	
	my $offset = 0;
	my $result = index($findIn, $findWhat, $offset);

	while ($result != -1)
	{
		if(scalar(@results) > 10)
		{
			die "Excessive number of matches in string of length " . length($findIn) . " for pattern: '$findWhat'";
		}
		push(@results, $result);
		$offset = $result + 1;
		$result = index($findIn, $findWhat, $offset);
	}
	
	return @results;
}

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	$kMer =~ tr/acgt/tgca/;
	return scalar(reverse($kMer));
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
	
	my $verbose_reading = 0;
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
