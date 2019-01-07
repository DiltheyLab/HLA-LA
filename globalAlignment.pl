#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::HTS;
use List::Util qw/max/;
use findPath;

my $S_match = 1;
my $S_mismatch = -1;
my $S_gap = -1;

my $bwa_bin;
my $samtools_bin;
my $reference;
my $query;
my $output;
my $endsFree_reference = 0;

GetOptions (
	'reference:s' => \$reference, 
	'query:s' => \$query,
	'output:s' => \$output,
	'bwa_bin:s' => \$bwa_bin,
	'samtools_bin:s' => \$samtools_bin,
);

die "Please specify valid --reference" unless(-e $reference);
die "Please specify valid --query" unless(-e $query);
die "Please specify valid --output" unless(defined $output);

# get correct paths of programs

$samtools_bin = findPath::find_path('samtools_bin', $samtools_bin, 'samtools');
$bwa_bin = findPath::find_path('bwa_bin', $bwa_bin, 'bwa');
findPath::check_samtools($samtools_bin);

my $reference_href = readFASTA($reference);
my $query_href = readFASTA($query);
die "Please specify a FASTA with exactly one entry for --reference" unless(scalar(values(%$reference_href)) == 1);
die "Please specify a FASTA with exactly one entry for --query" unless(scalar(values(%$query_href)) == 1);
my $querySequence_plus = (values %$query_href)[0];
my $querySequence_minus = reverseComplement($querySequence_plus);
my $referenceSequence = (values %$reference_href)[0];

my $prefix = $output;

my $forBWA_ref_fn = $prefix . '.ref';
writeFASTA($forBWA_ref_fn, {ref => $referenceSequence});

my $fromBWA_mapping_fn = $prefix . '.bam';
unless(-e $fromBWA_mapping_fn)
{
	my $index_cmd = qq($bwa_bin index $forBWA_ref_fn);
	system($index_cmd) and die "BWA index command $index_cmd failed";

	my $cmd_map = qq($bwa_bin mem -a -x pacbio $forBWA_ref_fn $query | $samtools_bin view -Sb - > $fromBWA_mapping_fn);
	system($cmd_map) and die "BWA mapping command $cmd_map failed";
}

my $sam = Bio::DB::HTS->new(-fasta => $forBWA_ref_fn, -bam => $fromBWA_mapping_fn);
my $iterator = $sam->features(-iterator=>1);

my $lenientOrder = 0;
my $warnedAboutOrder;
my %processedReadIDs;
my $runningReadID;
my %printed_complete_sequence;
my $queryID;
my $n_collected_alignments = 0;
my %alignments_perStrand;
while(my $alignment = $iterator->next_seq)
{
	next if($alignment->unmapped);
	
	die "Unknown reference ID '". $alignment->seq_id . "'"  unless($alignment->seq_id eq 'ref');
	$queryID = $alignment->query->name unless(defined $queryID);
	die unless($queryID eq $alignment->query->name);
	
	my $alignment_href = convertAlignmentToHash($alignment);	
	$n_collected_alignments++;
	push(@{$alignments_perStrand{$alignment_href->{strand}}}, $alignment_href) if($alignment_href);
	
	
}
print "Have collected $n_collected_alignments alignments for $queryID \n";

# find best global backtrace
if($n_collected_alignments)
{
	my @finalScores;
	my @finalScores_origins;
	my @finalScores_chromosome;
	my @finalScores_strand;
	
	my %chainAddr_2_i;

	foreach my $strand (keys %alignments_perStrand)
	{
		my @chains = enrichChains($alignments_perStrand{$strand});
		
		die unless(scalar(@chains));
		 
		foreach my $chain (@chains)
		{
			die unless((defined $chain->{n_matches}) and (defined $chain->{n_mismatches}) and (defined $chain->{n_gaps}));
			my $chain_score = $S_match * $chain->{n_matches} + $S_mismatch * $chain->{n_mismatches} + $S_gap * $chain->{n_gaps};
			$chain->{chainTraversalScore} = $chain_score;
		}
		
		@chains = sort {
			die unless(defined $a->{firstPos_reference});
			die unless(defined $b->{firstPos_reference});
			die unless(defined $a->{firstPos_read});
			die unless(defined $b->{firstPos_read});				
			
			if($a->{firstPos_reference} == $b->{firstPos_reference})
			{
				$a->{firstPos_read} <=> $b->{firstPos_read}
			}
			else
			{
				$a->{firstPos_reference} <=> $b->{firstPos_reference}				
			}
		} @chains;

		%chainAddr_2_i = ();
		if($strand eq '+')
		{
			print "Collected chains:\n";
			for(my $chainI = 0; $chainI <= $#chains; $chainI++)
			{
				my $chain = $chains[$chainI];
				print "Chain ${chainI}; reference ", $chain->{firstPos_reference}, " - ", $chain->{lastPos_reference}, ", read ", $chain->{firstPos_read}, " - ", $chain->{lastPos_read}, "; score ", $chain->{chainTraversalScore}, "; enrich $chain->{enrich}\n";
				$chainAddr_2_i{$chain} = $chainI;
			}
			print "";
			print "Sequence length: ", length($querySequence_plus), ", reference length: ", length($referenceSequence), "\n\n";
			
		}
		for(my $chainI = 0; $chainI <= $#chains; $chainI++)
		{
			my $chain = $chains[$chainI];
			
			my @inputScores;
			my @inputScoreOrigins;
			my @inputScores_preGaps;
			
			if($endsFree_reference)
			{
				push(@inputScores, $S_gap * $chain->{firstPos_read});
			}
			else
			{
				push(@inputScores, $S_gap * ($chain->{firstPos_read} + $chain->{firstPos_reference}));
			}
			
			push(@inputScoreOrigins, 0);
			push(@inputScores_preGaps, 0);
			
			for(my $chainII = 0; $chainII < $#chains; $chainII++)
			{
				if($chainII < $chainI)
				{
					my $chain2 = $chains[$chainII];
					if(($chain2->{lastPos_reference} < $chain->{firstPos_reference}) and ($chain2->{lastPos_read} < $chain->{firstPos_read}))
					{
						my $delta_reference = $chain->{firstPos_reference} - $chain2->{lastPos_reference} - 1;
						my $delta_read = $chain->{firstPos_read} - $chain2->{lastPos_read} - 1;
						die unless(defined $chain2->{chainOutputScore});
						push(@inputScores, $chain2->{chainOutputScore} + $S_gap * ($delta_read + $delta_reference));
						push(@inputScores_preGaps, $chain2->{chainOutputScore});
						push(@inputScoreOrigins, $chain2);
					}
				}
				else
				{
					my $chain2 = $chains[$chainII];
					die if(($chain2->{lastPos_reference} < $chain->{firstPos_reference}) and ($chain2->{lastPos_read} < $chain->{firstPos_read}));			
				}
			}
			die unless($#inputScores == $#inputScoreOrigins);
			
			if(0 and $chainI == 40)
			{
				print "Chain $chainI incoming scores:\n";	
				for(my $fromChainI = 0; $fromChainI <= $#inputScores; $fromChainI++)
				{
					my $origin = $inputScoreOrigins[$fromChainI];
					$origin = ($origin) ? $chainAddr_2_i{$origin} : "ORIGIN";
					print " - $origin $inputScores[$fromChainI] -- $inputScores_preGaps[$fromChainI] \n";
				}
			}
			my $max_input = which_max(\@inputScores);
			my $bestInputScore = $inputScores[$max_input];
			my $bestInputScore_origin = $inputScoreOrigins[$max_input];
			
			$chain->{chainOutputScore} = $bestInputScore + $chain->{chainTraversalScore};
			$chain->{chainOutputScoreOrigin} = $bestInputScore_origin;
		}
		
		my @finalOutputScores;
		my @finalOutputScoreOrigins;
		
		if($endsFree_reference)
		{
			push(@finalOutputScores, $S_gap * (length($querySequence_plus)));
		}
		else
		{
			push(@finalOutputScores, $S_gap * (length($querySequence_plus) + length($referenceSequence)));
		}
		push(@finalOutputScoreOrigins, 0);
		
		for(my $chainI = 0; $chainI <= $#chains; $chainI++)
		{
			my $chain = $chains[$chainI];
			die unless(defined $chain->{chainOutputScore});
			die unless(defined $chain->{chainOutputScoreOrigin});
			my $final_delta_read = length($querySequence_plus) - $chain->{lastPos_read} - 1;
			my $final_delta_ref = length($referenceSequence) - $chain->{lastPos_reference} - 1;
			if($endsFree_reference)
			{
				push(@finalOutputScores, $chain->{chainOutputScore} + $S_gap * $final_delta_read);
			}
			else
			{
				push(@finalOutputScores, $chain->{chainOutputScore} + $S_gap * ($final_delta_read + $final_delta_ref));
			}
			push(@finalOutputScoreOrigins, $chain);
		}			
		
		my $i_of_max_finalScore = which_max(\@finalOutputScores);
		my $maxFinalScore = $finalOutputScores[$i_of_max_finalScore];
		my $maxFinalScore_origin = $finalOutputScoreOrigins[$i_of_max_finalScore];

		push(@finalScores, $maxFinalScore);
		push(@finalScores_origins, $maxFinalScore_origin);
		push(@finalScores_chromosome, 'ref');
		push(@finalScores_strand, $strand);
	}
	
	die unless(scalar(@finalScores));
	my $i_of_max_finalfinalScore = which_max(\@finalScores);
	
	print "For read ${queryID}, found best alignment:\n";
	print "\t", "Score", ": ", $finalScores[$i_of_max_finalfinalScore], "\n";
	print "\t", "Chromosome", ": ", $finalScores_chromosome[$i_of_max_finalfinalScore], "\n";
	print "\t", "Strand", ": ", $finalScores_strand[$i_of_max_finalfinalScore], "\n";
	
	my $strand = $finalScores_strand[$i_of_max_finalfinalScore];
	my $chromosome = $finalScores_chromosome[$i_of_max_finalfinalScore];
	
	die unless(($strand eq '+') or ($strand eq '-'));
	my $useReadSequence = ($strand eq '+') ? $querySequence_plus : $querySequence_minus;
	
	my $bt_reference = '';
	my $bt_contig = '';
	
	my $next_bt_position = $finalScores_origins[$i_of_max_finalfinalScore];
	my $last_emitted_read_position = length($useReadSequence);
	my $last_emitted_reference_position = -1;
	my $min_emitted_reference_position;
	my $max_emitted_reference_position;
	
	my $leftOut_delta_reference = 0;
	my $n_thisBacktrace_chains = 0;
	while($last_emitted_read_position > 0)
	{
		die unless(defined $next_bt_position);
		if($next_bt_position == 0)
		{
			my $missing_contig_sequence = substr($useReadSequence, 0, $last_emitted_read_position);
			my $equivalent_reference_sequence_gap = ('-' x length($missing_contig_sequence));
			
			print "Left end, add delta ", $last_emitted_reference_position, "\n";
			$leftOut_delta_reference += ($last_emitted_reference_position);	
			
			$bt_contig .= reverse($missing_contig_sequence);
			$bt_reference .= $equivalent_reference_sequence_gap;
			
			$last_emitted_read_position = 0;
			$last_emitted_reference_position = -1;
		}
		else
		{
			my $chain_for_consumption = $next_bt_position;
			$n_thisBacktrace_chains++;
			
			if(exists $chainAddr_2_i{$chain_for_consumption})
			{
				print "Bt chain $chainAddr_2_i{$chain_for_consumption} \n";
			}
			my $delta_read = $last_emitted_read_position - $chain_for_consumption->{lastPos_read} - 1;
			if($n_thisBacktrace_chains == 1)
			{
				print "Right end, add delta ", (length($referenceSequence) - $chain_for_consumption->{lastPos_reference} - 1), "\n";
				
				$leftOut_delta_reference += (length($referenceSequence) - $chain_for_consumption->{lastPos_reference} - 1);
			}	
			die unless($delta_read >= 0);
			
			my $betweenChains_read = '';
			my $betweenChains_reference = '';
			
			my $jumpedOverRead_read = substr($useReadSequence, $chain_for_consumption->{lastPos_read} + 1, $last_emitted_read_position  - 1 - ($chain_for_consumption->{lastPos_read} + 1) + 1);
			die unless(length($jumpedOverRead_read) == $delta_read);			
			my $jumpedOverRead_reference = ('-' x $delta_read);
			die unless(length($jumpedOverRead_reference) == $delta_read);			
			
			$betweenChains_read .= $jumpedOverRead_read;
			$betweenChains_reference .= $jumpedOverRead_reference;
			
			if($last_emitted_reference_position != -1)
			{
				my $delta_reference = $last_emitted_reference_position - $chain_for_consumption->{lastPos_reference} - 1;
				
				my $jumpedOverReference_reference =  substr($reference_href->{$chromosome}, $chain_for_consumption->{lastPos_reference}+1, $last_emitted_reference_position - $chain_for_consumption->{lastPos_reference} - 1);
				die unless(length($jumpedOverReference_reference) == $delta_reference);				
				my $jumpedOverReference_read = ('-' x $delta_reference);
				die unless(length($jumpedOverReference_read) == $delta_reference);				
						
				$betweenChains_read .= $jumpedOverReference_read;
				$betweenChains_reference .= $jumpedOverReference_reference;				
			}

			die unless(length($betweenChains_read) == length($betweenChains_reference));
			
			$bt_contig .= reverse($chain_for_consumption->{alignment_read} . $betweenChains_read);
			$bt_reference .= reverse($chain_for_consumption->{alignment_reference} . $betweenChains_reference);
			
			die Dumper($chain_for_consumption->{firstPos_read}, $last_emitted_read_position) unless($chain_for_consumption->{firstPos_read} <= $last_emitted_read_position);
			$last_emitted_read_position = $chain_for_consumption->{firstPos_read};
			
			if($last_emitted_reference_position != -1)
			{
				die Dumper($chain_for_consumption->{firstPos_reference}, $last_emitted_reference_position)  unless($chain_for_consumption->{firstPos_reference} <= $last_emitted_reference_position);
			}
			$last_emitted_reference_position = $chain_for_consumption->{firstPos_reference};
			die Dumper($last_emitted_read_position)  unless($last_emitted_reference_position >= 0);
			
			$next_bt_position = $chain_for_consumption->{chainOutputScoreOrigin};
			
			if(not defined $max_emitted_reference_position)
			{
				$max_emitted_reference_position = $chain_for_consumption->{lastPos_reference};
			}
			
			if((not defined $min_emitted_reference_position) or ($chain_for_consumption->{firstPos_reference} < $min_emitted_reference_position))
			{
				$min_emitted_reference_position = $chain_for_consumption->{firstPos_reference};
			}
		}
	}
	
	if($last_emitted_reference_position != -1)
	{
		print "Final step with last_emitted_read_position = $last_emitted_read_position and last_emitted_reference_position = $last_emitted_reference_position :\n";
		print $last_emitted_reference_position, "\n";
		$leftOut_delta_reference += ($last_emitted_reference_position);	
	}
	
	$bt_contig = reverse($bt_contig);
	$bt_reference = reverse($bt_reference);
	
	my $bt_contig_noGaps = $bt_contig;
	$bt_contig_noGaps =~ s/\-//g;
	die Dumper($useReadSequence, $bt_contig_noGaps, length($bt_contig_noGaps), length($useReadSequence), 'II') unless($bt_contig_noGaps eq $useReadSequence);
	
	my $bt_reference_noGaps = $bt_reference;
	$bt_reference_noGaps =~ s/\-//g;
	my $reference_extract;
	if($bt_reference_noGaps)
	{
		die unless(defined $min_emitted_reference_position);
		die unless(defined $max_emitted_reference_position);
		die unless($min_emitted_reference_position <= $max_emitted_reference_position);
		$reference_extract = substr($reference_href->{$chromosome}, $min_emitted_reference_position, $max_emitted_reference_position - $min_emitted_reference_position + 1);
		my $reference_extract_expected_length = $max_emitted_reference_position - $min_emitted_reference_position + 1;
		 die unless(length($reference_extract) == $reference_extract_expected_length);
		if(index($reference_extract, "N") == -1)
		{
			if(($bt_reference_noGaps ne $reference_extract) and (length($bt_reference_noGaps) == length($reference_extract)))
			{
				for(my $i = 0; $i < length($bt_reference_noGaps); $i++)
				{
					my $c1 = substr($bt_reference_noGaps, $i, 1);
					my $c2 = substr($reference_extract, $i, 1);
					if($c1 ne $c2)
					{
						print "Mismatch position $i - $c1 vs $c2\n";
					}
				}
			}	
			die Dumper("Reference sequence mismatch", length($bt_reference_noGaps), length($reference_extract)) unless($bt_reference_noGaps eq $reference_extract);
			#die Dumper("Reference sequence mismatch",$bt_reference_noGaps, $reference_extract, length($bt_reference_noGaps), length($reference_extract)) unless($bt_reference_noGaps eq $reference_extract);
		}
		print "Check reference sequence correctness!\n";
	}
	
	die unless(length($bt_contig) == length($bt_reference));
	my $score_reconstructed = 0;
	my $n_mismatches = 0;
	for(my $i = 0; $i < length($bt_contig); $i++)
	{
		my $c1 = substr($bt_contig, $i, 1);
		my $c2 = substr($bt_reference, $i, 1);
		
		if(($c1 eq '-') or ($c2 eq '-'))
		{
			$score_reconstructed += $S_gap;
			$n_mismatches++;
		}
		else
		{
			if($c1 eq $c2)
			{
				$score_reconstructed += $S_match;
			}
			else
			{
				$score_reconstructed += $S_mismatch;
				$n_mismatches++;
			}
		}
	
	}
	
	my $score_before_endsFreeChange = $score_reconstructed;
	
	if(not $endsFree_reference)
	{
		$score_reconstructed += ($leftOut_delta_reference * $S_gap);
	}	
	unless($score_reconstructed == $finalScores[$i_of_max_finalfinalScore])
	{
		# die Dumper("Score mismatch", length($useReadSequence), $bt_reference, $bt_contig, $score_reconstructed, $finalScores[$i_of_max_finalfinalScore]);
		die Dumper("Score mismatch", [$score_reconstructed, $finalScores[$i_of_max_finalfinalScore]], ($score_reconstructed - $finalScores[$i_of_max_finalfinalScore]), $leftOut_delta_reference, $score_before_endsFreeChange, [$last_emitted_read_position, length($referenceSequence)],[$last_emitted_read_position,length($useReadSequence)]);
	}
	
	
	my $maxPos_contig = length($querySequence_plus) - 1;	
	open(OUTPUT, '>', $output) or die "Cannot open $output";
	print OUTPUT qq($n_mismatches ${min_emitted_reference_position}-${max_emitted_reference_position} ${strand}0-${maxPos_contig}), "\n";
	print OUTPUT $bt_reference, "\n";
	print OUTPUT $bt_contig, "\n";
	close(OUTPUT);	
}
else
{
	my $maxPos_contig = length($querySequence_plus) - 1;
	my $n_mismatches = length($querySequence_plus);
	my $min_emitted_reference_position = 0;
	my $max_emitted_reference_position = -1;
	my $strand = '+';
	my $bt_contig = $querySequence_plus;
	my $bt_reference = ('-' x $n_mismatches);
	open(OUTPUT, '>', $output) or die "Cannot open $output";
	print OUTPUT qq($n_mismatches ${min_emitted_reference_position}-${max_emitted_reference_position} ${strand}0-${maxPos_contig}), "\n";
	print OUTPUT $bt_reference, "\n";
	print OUTPUT $bt_contig, "\n";
	close(OUTPUT);	
}

print "\n\nDone.\n";

	
#unlink($forBWA_ref_fn);
#unlink($fromBWA_mapping_fn);



sub convertAlignmentToHash
{
	my $inputAlignment = shift;
	
	return undef if($inputAlignment->unmapped);
	
	my $readID = $inputAlignment->query->name;
	my $chromosome = $inputAlignment->seq_id;
	my $firstPos_reference = $inputAlignment->start - 1;
	die unless($firstPos_reference >= 0);
	
	my $lastPos_reference;
	my $firstPos_read = 0;
	my $lastPos_read;
	my $strand = $inputAlignment->strand;
	die unless(($strand == 1) or ($strand == -1));
	$strand = ($strand == 1) ? '+' : '-';
	
	my $n_matches = 0;
	my $n_mismatches = 0;
	my $n_insertions = 0;
	my $n_deletions = 0;
	my $n_gaps;
	my $alignment_reference;
	my $alignment_read;
	
	my $cigar  = $inputAlignment->cigar_str;
	my $remove_softclipping_front = 0;
	my $remove_softclipping_back = 0;	
	while($cigar =~ /^(\d+)([MIDHS])/)
	{
		my $number = $1;
		my $action = $2;
		$cigar =~ s/^(\d+)([MIDHS])//;		
		if(($action eq 'H') or ($action eq 'S'))
		{
			$firstPos_read += $number;
			if($action eq 'S')
			{
				$remove_softclipping_front += $number;
			}
		}
		else
		{
			last;
		}
	}
	
	{
		my $cigar  = $inputAlignment->cigar_str;
		while($cigar =~ /(\d+)([MIDHS])$/)
		{
			my $number = $1;
			my $action = $2;
			$cigar =~ s/(\d+)([MIDHS])$//;					
			if(($action eq 'H') or ($action eq 'S'))
			{
				if($action eq 'S')
				{
					$remove_softclipping_back += $number; 
				}
			}
			else
			{
				last;
			}
		}	
	}
	
	my ($ref, $matches, $query) = $inputAlignment->padded_alignment;
	if(length($ref) != length($query))
	{
		warn "Invalid output from padded_alignment() - " . length($ref) . " vs " . length($query) . ", ignore chain";
		return undef;
	}	
	
	unless(defined $ref)
	{
		warn "No reference sequence for $chromosome?";
		return undef;
	}
	
	if($remove_softclipping_front)
	{
		my $softclip_remove_ref = substr($ref, 0, $remove_softclipping_front);
		die "Weird softclipping for read" unless($softclip_remove_ref =~ /^\-+$/);
		substr($ref, 0, $remove_softclipping_front) = '';
		substr($query, 0, $remove_softclipping_front) = '';
	}

	if($remove_softclipping_back)
	{
		my $softclip_remove_ref = substr($ref, length($ref) - $remove_softclipping_back);
		die unless(length($softclip_remove_ref) == $remove_softclipping_back);
		
		die "Weird softclipping for read" unless($softclip_remove_ref =~ /^\-+$/);
		substr($ref, length($ref) - $remove_softclipping_back) = '';
		substr($query, length($query) - $remove_softclipping_back) = '';
	}
	
	$alignment_reference = $ref;
	$alignment_read = $query;
	
	my $alignment_length = length($ref);
	die Dumper(length($query), $alignment_length) unless(length($query) == $alignment_length);
	
	my $runningPos_query = 0;
	my $runningPos_reference = $firstPos_reference - 1;
	my $runningPos_read = $firstPos_read;
	for(my $alignmentPosI = 0; $alignmentPosI < $alignment_length; $alignmentPosI++)
	{
		my $c_ref = substr($ref, $alignmentPosI, 1);
		my $c_query = substr($query, $alignmentPosI, 1);
		die if(($c_ref eq '-') and ($c_query eq '-'));
		if($c_ref eq '-')
		{
			$n_insertions++;			
			$runningPos_read++;
		}
		elsif($c_query eq '-')
		{
			$n_deletions++;
			$runningPos_reference++;
			
		}
		else
		{
			$runningPos_reference++;
			$runningPos_read++;
			if($c_ref eq $c_query)
			{
				$n_matches++;
			}
			else
			{
				$n_mismatches++;
			}
		}
		
		$lastPos_reference = $runningPos_reference;
		$lastPos_read = $runningPos_read;
	}
	
	die unless(defined $lastPos_reference);
	die unless(defined $runningPos_read);
	$lastPos_read--;
	die unless($lastPos_read >= $firstPos_read);
	
	$n_gaps = $n_insertions + $n_deletions;
	
	die "Missing reference sequence for $chromosome" unless(exists $reference_href->{$chromosome});
	my $supposed_reference_sequence = substr($reference_href->{$chromosome}, $firstPos_reference, $lastPos_reference - $firstPos_reference + 1);
	
	my $read_sequence = $querySequence_plus;
	if($strand eq '-')
	{
		$read_sequence = reverseComplement($read_sequence);
	}
	my $supposed_read_sequence = substr($read_sequence, $firstPos_read, $lastPos_read - $firstPos_read + 1);
	
	my $alignment_reference_noGaps = $alignment_reference;
	$alignment_reference_noGaps =~ s/\-//g;
	
	my $alignment_read_noGaps = $alignment_read;
	$alignment_read_noGaps =~ s/\-//g;
	
	if(index($supposed_reference_sequence, "N") == -1)
	{
		unless($alignment_reference_noGaps eq $supposed_reference_sequence)
		{
			return undef;
		
			print "Reference mismatch for read $readID\n";
			print "\t", "CIGAR: ", $inputAlignment->cigar_str, "\n";
			print "\t", "Softclip remove: ", join("\t", $remove_softclipping_front, $remove_softclipping_back), "\n";			
			print "\t", $ref, "\n";
			print "\t", $query, "\n";
			print "\t", "\$inputAlignment->start: ", $inputAlignment->start, "\n";
			print "\t", "REF: ", $alignment_reference_noGaps, "\n";
			print "\t", "RE2: ", substr($reference_href->{$chromosome}, $firstPos_reference - 10, 20), "\n";
			print "\t", "ALG: ", $supposed_reference_sequence, "\n";
			print "\t", "strand: ", $strand, "\n";			
			print "\n";
			return undef;
		}

		unless($alignment_read_noGaps eq $supposed_read_sequence)
		{
			return undef;
			
			print "Sequence mismatch for read $readID\n";

			print "\t", "CIGAR: ", $inputAlignment->cigar_str, "\n";
			print "\t", "Softclip remove: ", join("\t", $remove_softclipping_front, $remove_softclipping_back), "\n";
			print "\t", "lastPos_read: ", $lastPos_read, "\n";
			print "\t", "firstPos_read: ", $firstPos_read, "\n";
			print "\t", "alignment_read_noGaps : ", $alignment_read_noGaps, "\n";
			print "\t", "supposed_read_sequence: ", $supposed_read_sequence, "\n\n";
			print "\t", "alignment_ref : ", $ref, "\n";
			print "\t", "alignment_query: ", $query, "\n";
			print "\t", "strand: ", $strand, "\n";
			print "\n";
			
			return undef;
		}
	}
	
	die Dumper("Mismatch query", $alignment_read_noGaps, $supposed_read_sequence, "Mismatch query") unless($alignment_read_noGaps eq $supposed_read_sequence);
	
	return {
		readID => $readID,
		chromosome => $chromosome,
		firstPos_reference => $firstPos_reference,
		lastPos_reference => $lastPos_reference,
		firstPos_read => $firstPos_read,
		lastPos_read => $lastPos_read,
		strand => $strand,
		n_matches => $n_matches,
		n_mismatches => $n_mismatches,
		n_gaps => $n_gaps,
		alignment_reference => $alignment_reference,
		alignment_read => $alignment_read,
		enrich => 0,		
	};
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

sub random_string
{
	my $l = 20;
	my $forReturn;
	for(my $i = 0; $i < $l; $i++)
	{
		$forReturn .= int(rand(9));
	}
	return $forReturn;
}

sub which_max
{
	my $input_aref = shift;
	die unless(@$input_aref);
	my $max = max(@$input_aref);
	for(my $i = 0; $i <= $#{$input_aref}; $i++)
	{
		return $i if($input_aref->[$i] == $max);
	}
	die;
}



sub enrichChains
{
	my $chainsIn_aref = shift;
	my @chainsIn = @$chainsIn_aref;
	
	my %beginCoordinates_reference;
	my %beginCoordinates_read;
	foreach my $chain (@chainsIn)
	{
		$beginCoordinates_reference{$chain->{firstPos_reference}}++;
		$beginCoordinates_read{$chain->{firstPos_read}}++;
	}
	
	
	my @chainsOut;	
	for(my $chainI = 0; $chainI <= $#chainsIn; $chainI++)
	{
		my %createdEndPointForReferenceCoordinate;
		my %createdEndPointForReadCoordinate;
		
		my $chain = $chainsIn[$chainI];
		
		my $useReadSequence = ($chain->{strand} eq '+') ? $querySequence_plus : $querySequence_minus;				
	
	
		{
			my $supposedReferenceSequence = substr($referenceSequence, $chain->{firstPos_reference}, $chain->{lastPos_reference} - $chain->{firstPos_reference} + 1);		
			if(index($supposedReferenceSequence, 'N') == -1)
			{
				my $chainSequence_reference_noGaps = $chain->{alignment_reference};
				$chainSequence_reference_noGaps =~ s/-//g;
				die Dumper(length($supposedReferenceSequence), length($chainSequence_reference_noGaps), $supposedReferenceSequence, $chainSequence_reference_noGaps) unless($supposedReferenceSequence eq $chainSequence_reference_noGaps);		
			}
		}
		push(@chainsOut, $chain);
		
		
		my $n_matches = 0;
		my $n_mismatches = 0;
		my $n_gaps = 0;
		my $alignment_reference = '';
		my $alignment_read = '';
		my $lastPos_reference = undef;
		my $lastPos_read = undef;
		for(my $alignmentPos = 0; $alignmentPos < length($chain->{alignment_reference}); $alignmentPos++)
		{
			my $character_reference = substr($chain->{alignment_reference}, $alignmentPos, 1);
			my $character_read = substr($chain->{alignment_read}, $alignmentPos, 1);
			
			$alignment_reference .= $character_reference;
			$alignment_read .= $character_read;
			
			if($character_reference eq $character_read)
			{
				$n_matches++;
			}
			else
			{
				if(($character_reference eq '-') or ($character_read eq '-'))
				{
					$n_gaps++;
				}
				else
				{
					$n_mismatches++;
				}
			}
			
			if($character_reference ne '-')
			{
				if(defined $lastPos_reference)
				{
					$lastPos_reference++;
				}
				else
				{
					$lastPos_reference = $chain->{firstPos_reference};
				}
				
			}
			if($character_read ne '-')
			{
				if(defined $lastPos_read)
				{
					$lastPos_read++;
				}
				else
				{
					$lastPos_read = $chain->{firstPos_read};
				}			
			}
			
			my $includeCurrentPositionAsEndpoint = 0;
			if((defined $lastPos_reference) and (defined $lastPos_read))
			{
				if(exists $beginCoordinates_reference{$lastPos_reference + 1})
				{
					if(not $createdEndPointForReferenceCoordinate{$lastPos_reference + 1})
					{
						$includeCurrentPositionAsEndpoint = 1;
					}
				}
				
				if(exists $beginCoordinates_read{$lastPos_read + 1})
				{
					if(not $createdEndPointForReadCoordinate{$lastPos_read + 1})
					{
						$includeCurrentPositionAsEndpoint = 1;
					}
				}
			}		

			if($includeCurrentPositionAsEndpoint)
			{
				$createdEndPointForReferenceCoordinate{$lastPos_reference + 1} = 1;
				$createdEndPointForReadCoordinate{$lastPos_read + 1} = 1;
				
				my $supposedReferenceSequence = substr($referenceSequence, $chain->{firstPos_reference}, $lastPos_reference - $chain->{firstPos_reference} + 1);
				if(index($supposedReferenceSequence, 'N') == -1)
				{								
					my $chainSequence_reference_noGaps = $alignment_reference;
					$chainSequence_reference_noGaps =~ s/-//g;
					die unless($supposedReferenceSequence eq $chainSequence_reference_noGaps);
				}
				
				my $supposedReadSequence = substr($useReadSequence, $chain->{firstPos_read}, $lastPos_read - $chain->{firstPos_read} + 1);
				my $chainSequence_read_noGaps = $alignment_read;
				$chainSequence_read_noGaps =~ s/-//g;
				die Dumper($character_reference, $character_read, $alignmentPos, $chain->{firstPos_read}, $chain->{lastPos_read}, $lastPos_read, length($supposedReadSequence), length($chainSequence_read_noGaps)) unless($supposedReadSequence eq $chainSequence_read_noGaps);				

				push(@chainsOut, {
					readID => $chain->{readID},
					chromosome => $chain->{chromosome},
					firstPos_reference => $chain->{firstPos_reference},
					lastPos_reference => $lastPos_reference,
					firstPos_read => $chain->{firstPos_read},
					lastPos_read => $lastPos_read,
					strand => $chain->{strand},
					n_matches => $n_matches,
					n_mismatches => $n_mismatches,
					n_gaps => $n_gaps,
					alignment_reference => $alignment_reference,
					alignment_read => $alignment_read,
					enrich => 1,
				});
			}			
		}
	}
	
	print "Chain enrichment: from ", scalar(@chainsIn), " to ", scalar(@chainsOut), " chains.\n";
	return @chainsOut;
}