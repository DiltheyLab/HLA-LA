#!/usr/bin/perl
use strict;
use Data::Dumper;
use List::MoreUtils qw/mesh/;

my $samtools_bin = findFileFromAlternatives('/home/dilthey/samtools-0.1.18/samtools', '/apps/well/samtools/0.1.19/bin/samtools');
my $KIR_info_file = qq(/gpfs1/well/gsk_hla/PRG-BWA/KIR_combined.txt);
my $BAM_in = qq(/gpfs1/well/gsk_hla/temp_mapping_2/NA19240/merged.bam);
my $BAM_out; qq(/gpfs1/well/gsk_hla/temp_mapping_2/NA19240/KIRonly.bam);
my @GRCH38_KIR_boundaries =  (54025634, 55084318);

# get KIR info

my %KIR_info;
my %known_KIR_lengths;
open(KIR, '<', $KIR_info_file) or die "Cannot open $KIR_info_file";
my $headerLine = <KIR>;
chomp($headerLine);
$headerLine =~ s/[\n\r]//g;
my @header_fields = split(/\t/, $headerLine);
while(<KIR>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/[\n\r]//g;
	my @line_fields = split(/\t/, $line);
	my %line = (mesh @header_fields, @line_fields);
	my $GI = $line{'HaplotypeID'};
	die unless($GI);
	die "Double ID $GI" if(exists $KIR_info{$GI});
	$KIR_info{$GI} = \%line;
	my $L = $line{'Length'};
	die "Double KIR length $L in $KIR_info_file" if(exists $known_KIR_lengths{$L});
	$known_KIR_lengths{$L} = $GI;
}
close(KIR);

my %contigs_by_length;

if(1 == 0)
{
	# get contigs from reference

	my $reference_FASTA = qq(/gpfs1/well/gsk_hla/GRCh38_with_alt/GRCh38_full_analysis_set_plus_decoy_hla.fa);
	my $ref_href = readFASTA($reference_FASTA);
	foreach my $contigID (keys %$ref_href)
	{
		my $L = length($ref_href->{$contigID});
		if(exists $contigs_by_length{$L})
		{
			die "Contig length $L observed more than once";
		}
		if(exists $known_KIR_lengths{$L})
		{
			$contigs_by_length{$L} = $contigID;
		}	
	}
}

# get BAM idx info
my $cmd_idx = qq($samtools_bin idxstats $BAM_in);
my $idx_output = `$cmd_idx`;
my @idx_output = split(/\n/, $idx_output);
my $chr19;
foreach my $l (@idx_output)
{
	next unless($l);
	my @f = split(/\t/, $l);
	die "Weird line $l" unless($#f == 3);
	my $L = $f[1];
	if(exists $contigs_by_length{$L})
	{
		die "Contig length $L observed more than once";
	}
	if(exists $known_KIR_lengths{$L})
	{
		$contigs_by_length{$L} = $f[0];
	}
	else
	{
		print $f[0], "\t-", $L, "-\n";
	}
	
	if(($f[0] eq '19') or ($f[0] eq 'chr19'))
	{
		die "Double definition for chr 19?" if (defined $chr19);
		$chr19 = $f[0];
	}
}
die "Can't determine chr19 contig ID from command $cmd_idx" unless(defined $chr19);

my @regions_to_extract = ($chr19.':'.$GRCH38_KIR_boundaries[0].'-'.$GRCH38_KIR_boundaries[1]);
foreach my $GI (keys %KIR_info)
{
	my $L = $KIR_info{$GI}{'Length'};
	die unless($L);
	if(exists $contigs_by_length{$L})
	{
		push(@regions_to_extract, $contigs_by_length{$L});
	}
	else
	{
		unless($GI eq 'ref')
		{
			die "No translation for $GI with length $L - is this GRCh38 + ALT input?";
		}
	}
}

my $BAM_out_unsorted = $BAM_out . '.unsorted';
my $regions_string = join(' ', @regions_to_extract);
my $command = qq($samtools_bin view -bh $BAM_in $regions_string > $BAM_out_unsorted);

print "Now executing:\n",$command,"\n\n";
my $ret = system($command);
unless($ret == 0)
{
	die "Command $command failed";
}
die unless(-e $BAM_out_unsorted);

my $command_sort = qq($samtools_bin sort $BAM_out_unsorted $BAM_out);
print "Now executing:\n",$command_sort,"\n\n";
my $ret_sort = system($command_sort);
unless($ret_sort == 0)
{
	die "Command $command_sort failed";
}
die unless(-e $BAM_out_unsorted);
my $BAM_out_sorted = $BAM_out . '.bam';
die "Sorted file $BAM_out_sorted not there" unless(-e $BAM_out_sorted);

my $cmd_mv = qq(mv $BAM_out_sorted $BAM_out);
print "Now executing:\n",$cmd_mv,"\n\n";
my $ret_mv = system($cmd_mv);
unless($ret_mv == 0)
{
	die "Command $cmd_mv failed";
}

my $command_idx = qq($samtools_bin index $BAM_out);

print "Now executing:\n",$command_idx,"\n\n";
my $ret_idx = system($command_idx);
unless($ret_idx == 0)
{
	die "Command $command_idx failed";
}

unlink($BAM_out_unsorted) or die "Cannot delete $BAM_out_unsorted";

die Dumper(\%contigs_by_length);




die Dumper(\%contigs_by_length);

sub findFileFromAlternatives
{
	foreach my $a (@_)
	{
		if(-e $a)
		{
			return $a;
		}
	}
	die "Couldn't find from alternatives: ".join(", ", @_);
}


sub readFASTA
{
	my $file = shift;	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{
		if(($. % 1000000) == 0)
		{
		# 	print "\r", $.;
		}
		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			$currentSequence = substr($line, 1);
		}
		else
		{
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}