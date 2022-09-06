package Util;

use strict;

my $verbose_reading = 0;
sub read_genomic_alignments
{
	my $directory = shift;
	my $preserveParts = shift;

	die unless(defined $directory);
	
	my %alignments;
	my @files =  glob ($directory . '/*_gen.txt'); 
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

1;