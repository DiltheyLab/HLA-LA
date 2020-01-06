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
my $inputSAM;
my $outputSAM;
my $bin_sam2alignment = '/gpfs/project/dilthey/software/NovoGraph_forJason/src/sam2alignment';
my $remap_reference = '/gpfs/project/dilthey/projects/HLA-LA-devel/working/NA12878_mini/remap/ref_for_remap.fa';
GetOptions (
	'graph:s' => \$graph,
	'inputSAM:s' => \$inputSAM,
	'outputSAM:s' => \$outputSAM,
);

die "Please specify --inputSAM" unless($inputSAM);
die "Please specify --outputSAM" unless($outputSAM);
die "bin_sam2alignment not existing: $bin_sam2alignment" unless(-e $bin_sam2alignment);

my $alignments_file = $inputSAM . '.alignments';
my $sam2alignment_cmd = qq($bin_sam2alignment $inputSAM $remap_reference > $alignments_file);
system($sam2alignment_cmd) and die "Could not execute: $sam2alignment_cmd";
print "Generated alignment file: $alignments_file\n";

my $full_graph_dir = $FindBin::RealBin . '/../../graphs/' . $graph;
die "Directory $full_graph_dir not existing" unless(-e $full_graph_dir);

my $raw_references_href = read_raw_references();
my $alignments_references_href = read_aligned_references();
my $alignments_translation_targets_href = compute_translation_targets($alignments_references_href);
my $raw_references_2AlignmentCoordinates_href = compute_alignment_coordinates($raw_references_href, $alignments_references_href);
my %translation_targets = map {$_ => 1} values %$alignments_translation_targets_href;

my $raw_references_doubleCheck_href = readFASTA($remap_reference);
foreach my $refID (keys %$raw_references_doubleCheck_href)
{
	die "Missing sequence: $refID" unless(exists $raw_references_href->{$refID});
	die "Sequence mismatch for sequence $refID" unless($raw_references_href->{$refID} eq $raw_references_doubleCheck_href->{$refID});
}

my $fn_output_ref = $outputSAM . '.ref.fa';
open(REFOUT, '>', $fn_output_ref) or die "Cannot open $fn_output_ref";

my $fh_output;
open($fh_output, '>', $outputSAM) or die "Cannot open $outputSAM";
print {$fh_output} "\@HD\tVN:1.6\n";
my @translation_targets_sorted = sort {length($raw_references_href->{$b}) <=> length($raw_references_href->{$a})} keys %translation_targets;
foreach my $translationTarget (@translation_targets_sorted)
{
	print {$fh_output} "\@SQ\tSN:$translationTarget\tLN:" . length($raw_references_href->{$translationTarget}) . "\n";
	print REFOUT '>', $translationTarget, "\n", $raw_references_href->{$translationTarget}, "\n";
}

close(REFOUT);

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
	
	
	if((defined $runningAlignmentInfo_readID) and ($readID_noPair ne $runningAlignmentInfo_readID))
	{
		processAlignments(\@runningAlignmentInfo, $fh_output);
		@runningAlignmentInfo = ();
	}
	$runningAlignmentInfo_readID = $readID_noPair;
	
	push(@runningAlignmentInfo, [\@headerLine_fields, $alignment_reference, $alignment_read, $qualities_read]);
	die Dumper(\@runningAlignmentInfo) if(scalar(@runningAlignmentInfo) > 2);

}
close(ALIGNMENTS);
processAlignments(\@runningAlignmentInfo, $fh_output);

print "\n";
print "n_secondary: $n_secondary\n";
print "n_supplementary: $n_supplementary\n";

sub processAlignments
{
	my $alignments_aref = shift;
	my $fh_output = shift;
	
	my @SAM_entries;
	foreach my $read_alignment (@$alignments_aref)
	{
		(my $readID_no12 = $read_alignment->[0][0]) =~ s/\/[12]//;
		
		(my $currentReferenceId, my $currentReferenceFrom_1based, my $currentReferenceTo_1based) = parseRefInfo($read_alignment->[0][1]);
		die unless($currentReferenceFrom_1based < $currentReferenceTo_1based);
		
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
		my $aligned_read_qual_noGaps = join('', map {substr($aligned_read, $_, 1)} @aligned_read_noGaps_indices);
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
		
		my $NM_new = 0;
		my $CIGAR = '';
		for(my $i = 0; $i < length($new_alignment_ref); $i++)
		{
			my $c_ref = substr($new_alignment_ref, $i, 1);
			my $c_read = substr($new_alignment_read, $i, 1);	
			
			my $ref_gap = (($c_ref eq '-') or ($c_ref eq '_'));
			my $read_gap = (($c_read eq '-') or ($c_read eq '_'));
			
			next if($ref_gap and $read_gap);
			
			if($ref_gap)
			{
				$CIGAR .= 'I';
				$NM_new++;
			}  
			elsif($read_gap)
			{
				$CIGAR .= 'D';  
				$NM_new++;
			}
			else
			{
				$CIGAR .= 'M';
				if($c_ref ne $c_read)
				{
					$NM_new++;
				}
			}
		}

		my $CIGAR_compressed = compressCIGAR($CIGAR);

		# print $CIGAR, "\n";		
		# print $CIGAR_compressed, "\n";
		# print $new_alignment_ref, "\n", $new_alignment_read, "\n";
		# print "nm $nm_before -> $NM_new \n";
		
		my $isPaired = ($FLAGs & 1) ? 1 : 0;
		my $isFirst = ($FLAGs & 64) ? 1 : 0;
		my $isSecond = ($FLAGs & 128) ? 1 : 0;
		if($isSecond)
		{
			die unless($read_alignment->[0][0] =~ /\/2/);
		}
		
		(my $sequence_prior_to_alignment_newRef = substr($alignments_references_href->{$projectOnto_sequenceId}, 0, $ref2ref_alignment_start_0based)) =~ s/[\-_]//g;
		die unless($sequence_prior_to_alignment_newRef =~ /^[ACGTN]*$/);
		
		my $newFlags = 0;
		$newFlags = ($newFlags | 1) if ($isPaired);
		$newFlags = ($newFlags | 16) if ($isReverseComplement);
		$newFlags = ($newFlags | 64) if ($isFirst);
		$newFlags = ($newFlags | 128) if ($isSecond);
		
		my @thisRead_CIGAR = (
			$readID_no12,
			$newFlags,
			$projectOnto_sequenceId,
			length($sequence_prior_to_alignment_newRef) + 1,
			60,
			$CIGAR_compressed,
			'*',
			0,
			length($aligned_read_noGaps),
			$aligned_read_noGaps,
			$aligned_read_qual_noGaps,
		);
		
		# print join("\t", @thisRead_CIGAR), "\n";
		push(@SAM_entries, \@thisRead_CIGAR);
	}
	
	if(scalar(@SAM_entries))
	{
		if(scalar(@SAM_entries) == 2)
		{
			my $read1_isReverseComplement = (($SAM_entries[0][1] & 16) ? 1 : 0);
			my $read2_isReverseComplement = (($SAM_entries[1][1] & 16) ? 1 : 0);
			
			# print $read1_isReverseComplement, " ", $read2_isReverseComplement, "\n";
			
			$SAM_entries[0][6] = $SAM_entries[1][2];
			$SAM_entries[1][6] = $SAM_entries[0][2];
			
			$SAM_entries[0][7] = $SAM_entries[1][3];
			$SAM_entries[1][7] = $SAM_entries[0][3];
			
			$SAM_entries[0][1] = ($SAM_entries[0][1] | 32) if ($read2_isReverseComplement);
			$SAM_entries[1][1] = ($SAM_entries[1][1] | 32) if ($read1_isReverseComplement);
			
			# print $SAM_entries[0][2], " ", $SAM_entries[1][2], "\n";
			if($SAM_entries[0][2] eq $SAM_entries[1][2])
			{
				if($read1_isReverseComplement == (! $read2_isReverseComplement))
				{
					my $fragment_left, my $fragment_right;
					
					if($SAM_entries[0][3] < $SAM_entries[1][3])
					{
						$fragment_left = $SAM_entries[0][3];
						$fragment_right = $SAM_entries[1][3] + $SAM_entries[1][8];
					}
					else
					{
						$fragment_left = $SAM_entries[1][3];
						$fragment_right = $SAM_entries[0][3] + $SAM_entries[0][8];
					}
					die unless($fragment_left < $fragment_right);
					my $fragment_size = $fragment_right - $fragment_left;
					
					# print "\t", $fragment_size, "\n";
					
					if($fragment_size <= 1000)
					{
						$SAM_entries[0][1] = ($SAM_entries[0][1] | 2);
						$SAM_entries[1][1] = ($SAM_entries[1][1] | 2);
					}
				}
			}
		}
		else
		{
			die unless(scalar(@SAM_entries) == 1);
		}
	}
	
	foreach my $SAMentry (@SAM_entries)
	{
		print {$fh_output} join("\t", @{$SAMentry}), "\n";
	}	
}

close($fh_output);

print "\n\nDone. Produced $outputSAM\n\n";

my $fn_output_BAM = $outputSAM . '.bam';
my $cmd_samtools_sort = qq(module load SamTools; samtools sort -o $fn_output_BAM $outputSAM; samtools index $fn_output_BAM);
system($cmd_samtools_sort) and die "Could not execute: $cmd_samtools_sort\n";

print "\n\nDone. Produced $fn_output_BAM\n\n";

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
	return $combined_raw_references_href;
}	

sub read_aligned_references
{
	my @filesToRead = glob($full_graph_dir . '/pseudoGenomic_fullLengthMapping/alignments_*');
	my $combined_aligned_references_href = {};
	foreach my $file (@filesToRead)
	{
		my $thisFile_href = readFASTA($file, 1);
		foreach my $seqID (keys %$thisFile_href)
		{
			next unless($seqID =~ /\*/);		
			die if(exists $combined_aligned_references_href->{$seqID});
			$combined_aligned_references_href->{$seqID} = $thisFile_href->{$seqID};
		}
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
		die unless($key =~ /^(\S+?)\*(\S+)$/);
		my $designated_target = $1 . '*pgf';
		if(exists $alignments_references_href->{$designated_target})
		{
			$translation_targets_href->{$key} = $designated_target;
		}
		else
		{
			warn "No target for $key";
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