#!/usr/bin/env perl

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
my $VCFin;
my $VCFout;
GetOptions (
	'graph:s' => \$graph,
	'VCFin:s' => \$VCFin,
	'VCFout:s' => \$VCFout,
);

my $full_graph_dir = $FindBin::RealBin . '/../../graphs/' . $graph;
die "Directory $full_graph_dir not existing" unless(-e $full_graph_dir);
die "Parameter --VCFin missing" unless($VCFin);
die "Parameter --VCFout missing" unless($VCFout);

my $PGF_sequence = VCFFunctions::getPGFSequence($full_graph_dir);

my %PGFcoordinates;
my @filesToRead = glob($full_graph_dir . '/pseudoGenomic_fullLengthMapping/alignments_*');
foreach my $file (@filesToRead)
{
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		if(substr($line, 0, 1) eq '>')
		{
			substr($line, 0, 1) = '';
			if($line =~ /^(\S+\*ref) pgf:(\d+)/)
			{
				my $contigID = $1;
				my $pgf_start = $2;
				
				my $seq_href = VCFFunctions::readFASTA($file, 1);
				die "Sequence $contigID missing from $file" unless(exists $seq_href->{$contigID});
				my $seq_contig = $seq_href->{$contigID};
				(my $seq_contig_noGaps = $seq_contig) =~ s/_//g;
				
				die Dumper("Sequence mismatch", substr($PGF_sequence, $pgf_start, length($seq_contig_noGaps)), uc($seq_contig_noGaps)) unless(uc(substr($PGF_sequence, $pgf_start, length($seq_contig_noGaps))) eq uc($seq_contig_noGaps));
				
				$PGFcoordinates{$contigID} = $pgf_start;				
			}	
		}
	}
	close(F);
}

(my $PGF_chr, my $PGF_chr_start_0based, my $PGF_chr_stop_0based) = VCFFunctions::getPGFCoordinates($full_graph_dir);
print "PGF chromosomal coordinates: ${PGF_chr}:${PGF_chr_start_0based}-${PGF_chr_stop_0based}\n";

open(VCFIN, '<', $VCFin) or die "Cannot open $VCFin";
open(VCFOUT, '>', $VCFout) or die "Cannot open $VCFout for writing";
while(<VCFIN>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	if(substr($line, 0, 2) eq '##')
	{
		next if($line =~ /contig=/);
		print VCFOUT $line, "\n";
	}
	elsif(substr($line, 0, 1) eq '#')
	{
		print VCFOUT $line, "\n";
	}
	else
	{
		my @line_fields = split(/\t/, $line);
		next if($line_fields[0] eq 'pgf_Ns');
		if(exists $PGFcoordinates{$line_fields[0]})
		{
			my $PGF_coordinate_1based = $line_fields[1] + $PGFcoordinates{$line_fields[0]};
			my $REF = $line_fields[3];
			
			die Dumper("Reference character mismatch", $line_fields[0], $PGF_coordinate_1based, uc(substr($PGF_sequence, $PGF_coordinate_1based, length($REF))), uc($REF)) unless(uc(substr($PGF_sequence, $PGF_coordinate_1based - 1, length($REF))) eq uc($REF));
			
			$line_fields[0] = 'chr6';			
			$line_fields[1] = ($PGF_coordinate_1based + $PGF_chr_start_0based);
			print VCFOUT join("\t", @line_fields), "\n";
		}
		else
		{
			print VCFOUT join("\t", @line_fields), "\n";			
		}
	}
}
close(VCFOUT);
close(VCFIN);

print "\n\nDone. Produced $VCFout\n\n";




