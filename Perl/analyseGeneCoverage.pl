#!/usr/bin/perl

my %coverager_per_gene;
my %positions_per_gene;

my $file = '/gpfs1/well/gsk_hla/MHC-PRG-2/working/NA12878_PLATINUM/reads_per_level.txt';
open(F, '<', $file) or die "Cannot open $file";
while(<F>)
{
	my $line = $_;
	chomp($line);
	next unless ($line);
	my @line_fields = split(/\t/, $line);
	die unless($#line_fields == 2);
	my $level = $line_fields[0];
	my $levelName = $line_fields[1];
	my $coverage = $line_fields[2];
	
	if($levelName =~ /(intron)|(exon)/)
	{
		die "Weird levelname $levelName" unless($levelName =~ /gene_(.+?)_/);
		my $geneName = $1;
		
		$coverager_per_gene{$geneName} += $coverage;
		$positions_per_gene{$geneName}++;
	}
}
close(F);

foreach my $geneName (sort keys %positions_per_gene)
{
	print join("\t", $geneName, $coverager_per_gene{$geneName}, $positions_per_gene{$geneName}, sprintf("%.2f", $coverager_per_gene{$geneName}/$positions_per_gene{$geneName})), "\n";
}