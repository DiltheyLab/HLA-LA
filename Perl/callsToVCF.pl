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
use VCFFunctions;

my $callFile;
my $graph = 'PRG_MHC_GRCh38_withIMGT';
my $VCFoutput;
my $sampleID;
GetOptions (
	'callFile:s' => \$callFile,
	'graph:s' => \$graph,
	'VCFoutput:s' => \$VCFoutput,
	'sampleID:s' => \$sampleID,
);

die "Please specify parameter --callFile (specifying the path to the R1_bestguess.txt file you want to transform into a VCF)" unless($callFile);
die "Please specify parameter --graph (e.g. PRG_MHC_GRCh38_withIMGT)" unless($graph);
die "Please specify parameter --VCFoutput (output file)" unless($VCFoutput);
die "Please specify parameter --sampleID as \w+" unless($sampleID =~ /^\w+$/);

print "callsToVCF.pl\n";
print "\tcallFile: $callFile\n";
print "\tGraph: $graph\n";
print "\tOutput file: $VCFoutput\n";
print "\tSample ID (can set via --sampleID): $sampleID\n";
print "\n";

my $full_graph_dir = $FindBin::RealBin . '/../../graphs/' . $graph;
die "Directory $full_graph_dir not existing" unless(-e $full_graph_dir);

(my $PGF_chr, my $PGF_chr_start_0based, my $PGF_chr_stop_0based) = VCFFunctions::getPGFCoordinates($full_graph_dir);
print "PGF chromosomal coordinates: ${PGF_chr}:${PGF_chr_start_0based}-${PGF_chr_stop_0based}\n";

my $PGF_sequence = VCFFunctions::getPGFSequence($full_graph_dir);
print "Acquired PGF sequence of length " . length($PGF_sequence) . "\n";

my $files_in_order_aref = VCFFunctions::getFilesInOrder($full_graph_dir);
my $pgf_reconstruction_aref = VCFFunctions::find_reconstruction($full_graph_dir, $files_in_order_aref, $PGF_sequence);
my $pgf_reconstruction_href = {mesh @$files_in_order_aref, @$pgf_reconstruction_aref};
my $PGF_with_gaps = '';
foreach my $file (@$files_in_order_aref)
{
	my $allele = $pgf_reconstruction_href->{$file};
	(my $segment_href, undef, undef) = readGraphSegment($full_graph_dir . '/PRG/' . $file);
	die unless(exists $segment_href->{$allele});
	$PGF_with_gaps .= $segment_href->{$allele};
}

(my $PGF_noGaps_control = $PGF_with_gaps) =~ s/_//g;
die unless(uc($PGF_noGaps_control) eq uc($PGF_sequence));

my %classI = map {$_ => 1} qw/HLA-A HLA-B HLA-C HLA-E HLA-F HLA-G/;
my %classII = map {$_ => 1} qw/HLA-DQA1 HLA-DQB1 HLA-DRB1 HLA-DRB3 HLA-DRB4 HLA-DPB1 HLA-DPA1/;


my %HLAtypes_PGF;
VCFFunctions::readPGFAlleles('Perl/PGF_loci_and_alleles.txt', \%HLAtypes_PGF);

my %calledHLA;
{
	open(CALLS, '<', $callFile) or die "Cannot open $callFile";
	my $headerLine = <CALLS>;
	chomp($headerLine);
	my @header_fields = split(/\t/, $headerLine);
	while(<CALLS>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line, -1);
		die unless(scalar(@line_fields) == scalar(@header_fields));
		my %line = (mesh @header_fields, @line_fields);
		die unless(defined $line{'Locus'});
		die unless(defined $line{'Chromosome'});
		die unless(defined $line{'Allele'});
		my @alleles = split(/;/, $line{'Allele'});
		$calledHLA{$line{'Locus'}}{$line{'Chromosome'}} = \@alleles;
	}
	close(CALLS);
}

my $VCF_fh;
open($VCF_fh, '>', $VCFoutput) or die "Cannot open $VCFoutput for writing";
print ${VCF_fh} "##fileformat=VCFv4.2\n";
print ${VCF_fh} "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampleID", "\n";

my @loci_in_genomic_order = sort {
	die Dumper("Missing position data for locus", $a, \%HLAtypes_PGF) unless(defined $HLAtypes_PGF{$a}{FirstBase_B38_0based});
	die Dumper("Missing position data for locus", $b, \%HLAtypes_PGF) unless(defined $HLAtypes_PGF{$b}{FirstBase_B38_0based});
	$HLAtypes_PGF{$a}{FirstBase_B38_0based} <=> $HLAtypes_PGF{$b}{FirstBase_B38_0based}
} grep {exists $HLAtypes_PGF{$_}} keys %calledHLA;

foreach my $locus (@loci_in_genomic_order)
{
	die unless($HLAtypes_PGF{$locus});
	print "Processing locus $locus\n";
	print "\tCall file:\n";
	foreach my $chr (1, 2)
	{
		print "\t\t$chr " . $calledHLA{$locus}{$chr}[0], " (etc.)\n";
	}	
	print "\tPGF allele:\n\t\t" . $HLAtypes_PGF{$locus}{'PGFAllele'} . "\n";
	
	my @exons = (2);
	if($classI{'HLA-' . $locus})
	{
		push(@exons, 3);
	}	
	
	foreach my $exon (@exons)
	{
		print "\tProcessing exon $exon\n";
		my @potential_files = glob($full_graph_dir . '/PRG/*gene_HLA-' . $locus . '_*_exon_' . $exon . '.txt');
		
		die Dumper("Wrong number of retrieved exon files", $locus, $exon, $full_graph_dir) unless(scalar(@potential_files) == 1);
		
		(my $sequences_href, my $start_graph, my $stop_graph) = readGraphSegment($potential_files[0]);
		print "\t\tGraph coordinates: $start_graph - $stop_graph\n";
		
		my $PGFallele_as_in_graph = $HLAtypes_PGF{$locus}{'PGFAllele'};
		
		die "PGF allele $PGFallele_as_in_graph not found in graph" unless(exists $sequences_href->{$PGFallele_as_in_graph});
		
		my $start_graph_preceding_noGaps = substr($PGF_with_gaps, 0, $start_graph);
		$start_graph_preceding_noGaps =~ s/_//g;
		
		my $PGF_allele_sequence_withGaps = $sequences_href->{$PGFallele_as_in_graph};
		(my $PGF_allele_noGaps = $PGF_allele_sequence_withGaps) =~ s/_//g;
		
		my $start_exon_PGF_0based = length($start_graph_preceding_noGaps);
		
		unless(uc(substr($PGF_sequence, $start_exon_PGF_0based, length($PGF_allele_noGaps))) eq uc($PGF_allele_noGaps))
		{
			die Dumper("Sequence mismatch", uc(substr($PGF_sequence, $start_exon_PGF_0based, length($PGF_allele_noGaps))), uc($PGF_allele_noGaps), $start_exon_PGF_0based, length($PGF_allele_noGaps));
		}
		
		print "\t\tStart coordinate within PGF (0-based): $start_exon_PGF_0based\n";
		
		my $getSequence = sub {
			my $alleles_aref = shift;
			die unless($alleles_aref);
			my %S;
			foreach my $a (@$alleles_aref)
			{
				my $s = $sequences_href->{$a};
				die unless($s);
				$S{$s}++;
			}
			my @S = keys %S;
			if(scalar(@S) != 1)
			{
				die Dumper("Weird exon sequence count within G group", \@S);
			}
			return $S[0];
		};
		
		my $S1 = $getSequence->($calledHLA{$locus}{1});
		my $S2 = $getSequence->($calledHLA{$locus}{2});
		die unless(length($PGF_allele_sequence_withGaps) == length($S1));
		die unless(length($PGF_allele_sequence_withGaps) == length($S2));
		
		print "\t\tSequences for VCF output:\n";
		print "\t\t\tPGF ref: $PGF_allele_sequence_withGaps\n";
		print "\t\t\tChrom 1: $S1\n";
		print "\t\t\tChrom 2: $S2\n";
		
		VCFFunctions::outputToVCF($VCF_fh, $start_exon_PGF_0based, $PGF_allele_sequence_withGaps, $S1, $S2, $full_graph_dir, $locus . '-exon' . $exon);
	}
}
close($VCF_fh);
print "\n\nDone. Generated file $VCFoutput\n\n";

sub readGraphSegment
{
	my $segmentFile = shift;
	
	my %sequences_forReturn;
	
	open(ONESEGMENT, '<', $segmentFile) or die "Cannot open $segmentFile";
	my $header_segment = <ONESEGMENT>;
	chomp($header_segment);
	$header_segment =~ s/[\n\r]//g;
	my @header_fields = split(/ /, $header_segment);
	die unless($header_fields[0] eq 'IndividualID');
	my @column_coordinates = map {die unless($_ =~ /^(\d+)_/); $1} @header_fields[1 .. $#header_fields]; 
	
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
		
		(my $hla_allele_sequence_noGaps = $hla_allele_sequence) =~ s/_//g;
		die unless($hla_allele_sequence_noGaps =~ /^[ACGTN]*$/);

		$sequences_forReturn{$hla_allele} = $hla_allele_sequence; 
	}

	close(ONESEGMENT);	
	
	return (\%sequences_forReturn, min(@column_coordinates), max(@column_coordinates));
}

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	return reverse($kMer);
	return $kMer;
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


