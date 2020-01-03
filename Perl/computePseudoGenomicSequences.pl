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

my $graph;
GetOptions (
	'graph:s' => \$graph,
);

my $full_graph_dir = $FindBin::RealBin . '/../../graphs/' . $graph;
die "Directory $full_graph_dir not existing" unless(-e $full_graph_dir);

my $pgf_sequence;
{
	my $sequences_file = $full_graph_dir . '/sequences.txt';
	my $pgf_id;
	my $pgf_chr;
	my $pgf_start;
	my $pgf_stop;
	open(SEQUENCES, '<', $sequences_file) or die "Cannot open $sequences_file";
	my $header_sequences = <SEQUENCES>;
	$header_sequences =~ s/[\n\r]//g;	
	chomp($header_sequences);
	my @header_sequences_fields = split(/\t/, $header_sequences);
	die unless(($header_sequences_fields[0] eq 'SequenceID') and ($header_sequences_fields[1] eq 'Name'));
	while(<SEQUENCES>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		next unless($line);
		my @line_fields = split(/\t/, $line, -1);
		die Dumper("Field mismatch") unless(scalar(@line_fields) == scalar(@header_sequences_fields));
		my %line = mesh @header_sequences_fields, @line_fields;
		if($line{Name} eq 'pgf')
		{
			$pgf_id = $line{FASTAID};
			$pgf_chr = $line{Chr};
			$pgf_start = $line{Start_1based};
			$pgf_stop = $line{Stop_1based};
		}
	}
	close(SEQUENCES);
	die unless(defined $pgf_id);
	
	my $extendedReferenceGenome = $full_graph_dir . '/extendedReferenceGenome/extendedReferenceGenome.fa';
	my $referenceGenome_href = readFASTA($extendedReferenceGenome, 1);
	die unless(exists $referenceGenome_href->{$pgf_chr});
	$pgf_sequence = substr($referenceGenome_href->{$pgf_chr}, $pgf_start - 1, $pgf_stop - $pgf_start + 1);
	
	print "Acquired PGF sequence of length ", length($pgf_sequence), "\n";
}

my @files_in_order;

my %files_per_gene;
my $fn_graph_segments = $full_graph_dir . '/PRG/segments.txt';
die "Missing file: $fn_graph_segments" unless(-e $fn_graph_segments);
open(SEGMENTS, '<', $fn_graph_segments) or die "Cannot open $fn_graph_segments";
while(<SEGMENTS>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/[\n\r]//g;
	next unless($line);
	
	my $filename = $line;
	
	if($line =~ /\d+_gene_([\w\-]+)_\d_((intron)|(exon))/)
	{
		my $geneID = $1;
		push(@{$files_per_gene{$geneID}}, $filename);
	}
	
	push(@files_in_order, $filename);
}
close(SEGMENTS);

my $pgf_reconstruction_aref = find_reconstruction(\@files_in_order, $pgf_sequence);
my $pgf_reconstruction_href = {mesh @files_in_order, @$pgf_reconstruction_aref};

my $output_dir = $full_graph_dir . '/pseudoGenomic_fullLengthMapping';
unless(-d $output_dir)
{
	mkdir($output_dir) or die "Cannot mkdir $output_dir";
}

my @genes = sort keys %files_per_gene;
# @genes = ('HLA-DRA');

my %sequences_per_gene;
foreach my $gene (@genes)
{
	(my $gene_noHLA = $gene) =~ s/HLA\-//;
	
	my %all_alleles;
	my %sequences_per_gene_per_file;
	my $pgf_allele;
	foreach my $file (@{$files_per_gene{$gene}})
	{
		my $full_file_path = $full_graph_dir . '/PRG/' . $file;
		open(ONESEGMENT, '<', $full_file_path) or die "Cannot open $full_file_path";
		my $header_segment = <ONESEGMENT>;
		chomp($header_segment);
		$header_segment =~ s/[\n\r]//g;
		my @header_fields = split(/ /, $header_segment);
		die unless($header_fields[0] eq 'IndividualID');
		die unless(exists $pgf_reconstruction_href->{$file});
		my $pgf_allele_sequence;
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
				$pgf_allele_sequence = $hla_allele_sequence;
			}
		}
		close(ONESEGMENT);
		
		die unless(defined $pgf_allele_sequence);
		$pgf_allele .= $pgf_allele_sequence;
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
	die unless(length($pgf_allele) == $combined_allele_length);
	
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
	
	(my $pgf_allele_noGaps = $pgf_allele) =~ s/_//g;
	my $pgf_allele_noGaps_startPos_0based_inPGF;
	if(length($pgf_allele_noGaps))
	{
		my @positions_pgf_allele_noGaps = find_all($pgf_sequence, $pgf_allele_noGaps);
		unless(scalar(@positions_pgf_allele_noGaps) == 1)
		{
			die "Failed to find unambiguous position for allele sequence within PGF";
		}
		$pgf_allele_noGaps_startPos_0based_inPGF = $positions_pgf_allele_noGaps[0];
	}
	else
	{
		$pgf_allele_noGaps_startPos_0based_inPGF = -1;
	}
	
	my $output_fn_alignments = $output_dir . '/alignments_' . $gene . '.fa';
	my $output_fn_raw = $output_dir . '/raw_' . $gene . '.fa';
	
	open(OUTPUT_ALIGNMENTS, '>', $output_fn_alignments) or die "Cannot open $output_fn_alignments";
	open(OUTPUT_RAW, '>', $output_fn_raw) or die "Cannot open $output_fn_raw";
	print OUTPUT_ALIGNMENTS '>', $gene_noHLA . '*pgf', ' ', $pgf_allele_noGaps_startPos_0based_inPGF, "\n", $pgf_allele, "\n";
	print OUTPUT_RAW '>', $gene_noHLA . '*pgf', ' ', $pgf_allele_noGaps_startPos_0based_inPGF, "\n", $pgf_allele_noGaps, "\n";
	
	foreach my $alleleID (keys %combined_allele_sequences_fullyGenomic)
	{
		print OUTPUT_ALIGNMENTS '>', $alleleID, "\n";
		print OUTPUT_RAW '>', $alleleID, "\n";
		
		print OUTPUT_ALIGNMENTS $combined_allele_sequences_fullyGenomic{$alleleID}, "\n";
		(my $alleleSequence_noGaps = $combined_allele_sequences_fullyGenomic{$alleleID}) =~ s/_//g;
		print OUTPUT_RAW $alleleSequence_noGaps, "\n";		
	}
	close(OUTPUT_ALIGNMENTS);
	close(OUTPUT_RAW);
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


sub find_reconstruction
{
	my $files_aref = shift;
	my $sequence_to_reconstruct = shift;
	
	print "Scanning for sequence of length ", length($sequence_to_reconstruct), "\n";
	my @running_reconstructions = ({sequence => [], reached_position => -1});
	foreach my $filename (@$files_aref)
	{
		print "\r\t $filename " . scalar(@running_reconstructions) . " alternatives.      ";
		die unless(@running_reconstructions);
		my @new_reconstructions;
		my $full_file_path = $full_graph_dir . '/PRG/' . $filename;
		open(ONESEGMENT, '<', $full_file_path) or die "Cannot open $full_file_path";
		my $header_segment = <ONESEGMENT>;
		chomp($header_segment);
		$header_segment =~ s/[\n\r]//g;
		my @header_fields = split(/ /, $header_segment);
		die unless($header_fields[0] eq 'IndividualID');
		my %processed_sequences;
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
			next if($processed_sequences{$hla_allele_sequence});
			$processed_sequences{$hla_allele_sequence}++;
			(my $hla_allele_sequence_noGaps = $hla_allele_sequence) =~ s/_//g;
			die unless($hla_allele_sequence_noGaps =~ /^[ACGTN]*$/);
			foreach my $running_reconstruction (@running_reconstructions)
			{
				# print $running_reconstruction->{reached_position}, "\n";
				if(length($hla_allele_sequence_noGaps))
				{
					my $supposed_sequence = uc(substr($sequence_to_reconstruct, $running_reconstruction->{reached_position}+1, length($hla_allele_sequence_noGaps)));
					if($hla_allele_sequence_noGaps eq $supposed_sequence)
					{
						my $new_reconstruction_href = {
							sequence => [@{$running_reconstruction->{sequence}}, $hla_allele],
							reached_position => $running_reconstruction->{reached_position} + length($hla_allele_sequence_noGaps),
						};
						push(@new_reconstructions, $new_reconstruction_href);
					}
				}
				else
				{
					my $new_reconstruction_href = {
						sequence => [@{$running_reconstruction->{sequence}}, $hla_allele],
						reached_position => $running_reconstruction->{reached_position},
					};
					push(@new_reconstructions, $new_reconstruction_href);				
				}
			}
		}
		
		@running_reconstructions = @new_reconstructions;
		
		close(ONESEGMENT);
	}
	
	print "\r\tDone. Found ", scalar(@running_reconstructions), " possible traversals.                                                 \n";
	die unless(scalar(@running_reconstructions) == 1);
	
	return $running_reconstructions[0]->{sequence};
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
