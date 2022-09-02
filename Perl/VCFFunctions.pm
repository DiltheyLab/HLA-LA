package VCFFunctions;

use strict;
use List::MoreUtils qw/all mesh/;
use Data::Dumper;

my %classI = map {$_ => 1} qw/HLA-A HLA-B HLA-C HLA-E HLA-F HLA-G/;
my %classII = map {$_ => 1} qw/HLA-DQA1 HLA-DQB1 HLA-DRB1 HLA-DRB3 HLA-DRB4 HLA-DPB1 HLA-DPA1/;

sub readPseudoGenomicSequences
{
	my $directory = shift;
	my $forReturn_alignments = shift;
	my $forReturn_refPos = shift;
	
	%$forReturn_alignments = ();
	%$forReturn_refPos = ();
	
	my @files = glob($directory . '/*alignments*.fa');

	foreach my $f (@files)
	{
		next if($f =~ /controlUnmodified/);
		
		die "Can't parse file $f" unless($f =~ /((raw)|(alignments))_([\w\-]+).fa/);
		my $type = $1;
		my $gene = $4;
		die "Weird type $type" unless($type eq 'alignments');
		
		my $fasta_href = readFASTA($f, 0);
		my @ids_ref = grep {$_ =~ /\*ref/} keys %$fasta_href;
		die Dumper("Weird reference count", \@ids_ref) unless(scalar(@ids_ref) == 1);
		die "Can't parse ref - '$ids_ref[0]'" unless($ids_ref[0] =~ /^\S+ (\w+):(\d+)$/);
		
		$forReturn_refPos->{$gene} = [$1, $2];
	
		$fasta_href = readFASTA($f, 1);
		$forReturn_alignments->{$gene} = $fasta_href;
	}
}

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
	
my %_cache_getPGFCoordinates;

sub readPGFAlleles
{
	my $inputFile = shift;
	my $HLAtypes_PGF_href = shift;
	%$HLAtypes_PGF_href = ();
	
	open(PGF_LOCI_AND_ALLELES, '<', $inputFile) or die "Cannot open $inputFile";
	my $PGFalleles_header_line = <PGF_LOCI_AND_ALLELES>;
	$PGFalleles_header_line =~ s/[\n\r]//g;
	my @PGF_header_fields = split(/\t/, $PGFalleles_header_line);
	while(<PGF_LOCI_AND_ALLELES>)
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

		$HLAtypes_PGF_href->{$locus} = \%line;

		# unless((exists $alleles_to_fullGAmbiguity{$locus_with_HLA}) and (exists $alleles_to_fullGAmbiguity{$locus_with_HLA}{$allele_numeric}))
		# {
			# $allele_numeric .= 'G';
			# die Dumper("Undefined allele $line{PGFAllele}", [(keys %{$alleles_to_fullGAmbiguity{$locus_with_HLA}})[0 .. 10]]) unless((exists $alleles_to_fullGAmbiguity{$locus_with_HLA}) and (exists $alleles_to_fullGAmbiguity{$locus_with_HLA}{$allele_numeric}));
		# }
	}
	close(PGF_LOCI_AND_ALLELES);
	
}
sub getPGFCoordinates
{
	my $full_graph_dir = shift;
	if(exists $_cache_getPGFCoordinates{$full_graph_dir})
	{
		my @forReturn = @{$_cache_getPGFCoordinates{$full_graph_dir}};
		die unless(scalar(@forReturn) == 3);
		return $forReturn[0], $forReturn[1], $forReturn[2];
	}
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
	
	$_cache_getPGFCoordinates{$full_graph_dir} = [$pgf_chr, $pgf_start - 1, $pgf_stop - 1];

	return ($pgf_chr, $pgf_start - 1, $pgf_stop - 1);
}

my %_cache_getPGFSequence;
sub getPGFSequence
{
	my $full_graph_dir = shift;
	if(exists $_cache_getPGFSequence{$full_graph_dir})
	{
		my $forReturn = $_cache_getPGFSequence{$full_graph_dir};
		die unless(length($forReturn));
		return $forReturn;
	}
	
	(my $pgf_chr, my $pgf_start, my $pgf_stop) = getPGFCoordinates($full_graph_dir);

	my $extendedReferenceGenome = $full_graph_dir . '/extendedReferenceGenome/extendedReferenceGenome.fa';
	my $referenceGenome_href = readFASTA($extendedReferenceGenome, 1);
	die unless(exists $referenceGenome_href->{$pgf_chr});
	
	my $forReturn_PGF_sequence = substr($referenceGenome_href->{$pgf_chr}, $pgf_start, $pgf_stop - $pgf_start + 1);
	$_cache_getPGFSequence{$full_graph_dir} = $forReturn_PGF_sequence;
	return $forReturn_PGF_sequence;
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

sub outputToVCF
{
	my $VCF_fh = shift;
	my $start_sequence_0based = shift;
	my $PGF_allele_sequence_withGaps = shift;
	my $S1 = shift;
	my $S2 = shift;
	my $full_graph_dir = shift;
	my $geneName = shift;
	
	(my $PGF_chr, my $PGF_chr_start_0based, my $PGF_chr_stop_0based) = VCFFunctions::getPGFCoordinates($full_graph_dir);
	my $PGF_sequence = VCFFunctions::getPGFSequence($full_graph_dir);

	die unless(length($PGF_allele_sequence_withGaps) == length($S1));
	die unless(length($PGF_allele_sequence_withGaps) == length($S2));

	my $running_PGF_coordinate = $start_sequence_0based;		
	my $running_ref_start_0based = $start_sequence_0based;
	my $running_ref = '';
	my $running_S1 = '';
	my $running_S2 = '';
	
	my $flush_genotypes = sub {
		if(length($running_ref) > 1)
		{
			if((substr($running_ref, 0, 1) eq '_') or (substr($running_S1, 0, 1) eq '_') or (substr($running_S2, 0, 1) eq '_'))
			{
				warn "INDEL at the beginning - can't represent variant.";
			}
			else
			{
				if((substr($running_ref, 0, 1) eq substr($running_S1, 0, 1)) and (substr($running_ref, 0, 1) eq substr($running_S2, 0, 1)))
				{
					my $allele_S1 = substr($running_S1, 0, length($running_S1) - 1);
					my $allele_S2 = substr($running_S2, 0, length($running_S2) - 1);
					my $allele_ref = substr($running_ref, 0, length($running_ref) - 1);
					$allele_S1 =~ s/_//g;
					$allele_S2 =~ s/_//g;
					$allele_ref =~ s/_//g;
					my %alleles = map {$_ => 1} ($allele_S1, $allele_S2, $allele_ref);
					if(scalar(keys %alleles) > 1)
					{
						my @alleles_in_order = ($allele_ref);
						push(@alleles_in_order, grep {$_ ne $allele_ref} sort keys %alleles);
						my %allele_2_i = map {$alleles_in_order[$_] => $_} (0 .. $#alleles_in_order);
						my $GT = join('/', map {die unless(exists $allele_2_i{$_}); $allele_2_i{$_}} ($allele_S1, $allele_S2));
						
						my $adjustedOutputPos = $running_ref_start_0based + 1;
						if(all {length($_) == length($alleles_in_order[0])} @alleles_in_order)
						{
							$allele_ref = substr($allele_ref, 1);
							@alleles_in_order = map {substr($_, 1)} @alleles_in_order;
							$adjustedOutputPos++;
						}
						
						if(all {length($_) == length($alleles_in_order[0])} @alleles_in_order)
						{
							my $L = length($alleles_in_order[0]);
							for(my $i = 0; $i < $L; $i++)
							{
								my $thisPos = $adjustedOutputPos + $i;
								my $thisPos_allele_ref = substr($allele_ref, $i, 1),;
								my @thisPos_alleles_in_order = map {substr($_, $i, 1)} @alleles_in_order;
								print ${VCF_fh} join("\t",
									$PGF_chr,
									$thisPos + $PGF_chr_start_0based,
									'.',
									$thisPos_allele_ref,
									join(',', @thisPos_alleles_in_order[1 .. $#alleles_in_order]),
									'.',
									'PASS',
									'gene=' . $geneName,
									'GT',									
									$GT
								), "\n";
								die unless(substr($PGF_sequence, $thisPos - 1, length($thisPos_allele_ref)) eq $thisPos_allele_ref);	
							}
						}		
						else
						{
							print ${VCF_fh} join("\t",
								$PGF_chr,
								$adjustedOutputPos + $PGF_chr_start_0based,
								'.',
								$allele_ref,
								join(',', @alleles_in_order[1 .. $#alleles_in_order]),
								'.',
								'PASS',
								'gene=' . $geneName,
								'GT',
								$GT
							), "\n";
							die unless(substr($PGF_sequence, $adjustedOutputPos - 1, length($allele_ref)) eq $allele_ref);							
						}
					}
				}
				else
				{
					warn "First characters don't agree - weird";
				}
			}
		}
	};
	
	for(my $i = 0; $i < length($PGF_allele_sequence_withGaps); $i++)
	{
		my $c_ref = substr($PGF_allele_sequence_withGaps, $i, 1);
		my $c_S1 = substr($S1, $i, 1);
		my $c_S2 = substr($S2, $i, 1);
		
		$running_ref .= $c_ref;
		$running_S1 .= $c_S1;
		$running_S2 .= $c_S2;
		
		if(($c_ref eq $c_S1) and ($c_ref eq $c_S2))
		{
		
			$flush_genotypes->();
			
			$running_ref = $c_ref;
			$running_S1 = $c_S1;
			$running_S2 = $c_S2;
			
			$running_ref_start_0based = $running_PGF_coordinate;
		}
		
		$running_PGF_coordinate++ if($c_ref ne '_');			
	}
	
	$flush_genotypes->();
}

1;
