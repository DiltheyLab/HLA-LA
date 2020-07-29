package VCFFunctions;

use strict;
use List::MoreUtils qw/mesh/;

sub getFilesInOrder
{
	my $full_graph_dir = shift;
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
		push(@files_in_order, $filename);		
	}
	close(SEGMENTS);
	
	return \@files_in_order;
}

sub getPGFCoordinates
{
	my $full_graph_dir = shift;
	
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
	
	return ($pgf_chr, $pgf_start - 1, $pgf_stop - 1);
	
	my $extendedReferenceGenome = $full_graph_dir . '/extendedReferenceGenome/extendedReferenceGenome.fa';
	my $referenceGenome_href = readFASTA($extendedReferenceGenome, 1);
	die unless(exists $referenceGenome_href->{$pgf_chr});
	
	return substr($referenceGenome_href->{$pgf_chr}, $pgf_start - 1, $pgf_stop - $pgf_start + 1);;
}

sub getPGFSequence
{
	my $full_graph_dir = shift;

	(my $pgf_chr, my $pgf_start, my $pgf_stop) = getPGFCoordinates($full_graph_dir);

	my $extendedReferenceGenome = $full_graph_dir . '/extendedReferenceGenome/extendedReferenceGenome.fa';
	my $referenceGenome_href = readFASTA($extendedReferenceGenome, 1);
	die unless(exists $referenceGenome_href->{$pgf_chr});
	
	return substr($referenceGenome_href->{$pgf_chr}, $pgf_start, $pgf_stop - $pgf_start + 1);;
}

sub find_reconstruction
{
	my $full_graph_dir = shift;
	my $files_aref = shift;
	my $sequence_to_reconstruct = shift;
	
	print "Scanning for sequence of length ", length($sequence_to_reconstruct), "\n";
	my @running_reconstructions = ({sequence => [], reached_position => -1});
	foreach my $filename (@$files_aref)
	{
		# print "\r\t $filename " . scalar(@running_reconstructions) . " alternatives.     \n ";
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

1;
