use warnings;
use strict;
use FindBin;
use Getopt::Long;
use Data::Dumper; 
use List::MoreUtils qw/all/;

my $assemblies_directory = '/home/dilthey/HLA-LA-devel/assemblies';
my $minimap2_bin = '/home/dilthey/minimap2-2.24_x64-linux/minimap2';
die unless(-e $minimap2_bin);

my $this_bin_dir = $FindBin::RealBin;
my $graphs_dir = $this_bin_dir . '/../../graphs/';

my $sampleID = 'NA12878';
my $prefix = '/home/dilthey/HLA-LA-devel/test/mock';
my $graph = 'PRG_MHC_GRCh38_withIMGT';
GetOptions (
	'sampleID:s' => \$sampleID,
	'prefix:s' => \$prefix,
	'graph:s' => \$graph,
);

my $outputFn = $prefix . '.evaluationSummary.txt';

die "Please specify --sampleID" unless($sampleID);
die "Please specify --prefix" unless($prefix);
die "Please specify --graph" unless($graph);

my $graph_dir = $graphs_dir . '/' . $graph;
die "Expected graph directory $graph_dir not found" unless(-d $graph_dir);

my $fastaFile = $prefix . '.fullLengthInference.fasta';
die "Expected file $fastaFile not existing" unless(-e $fastaFile);

my %assemblyFile_2_H;
my %H_2_assemblyFile;
my @assemblyFiles = glob($assemblies_directory . '/*' . $sampleID . '*');
unless(scalar(@assemblyFiles))
{
	die "Cannot identify exactly two assembly files for sample ID '$sampleID' in directory $assemblies_directory";
}
foreach my $asmFile (@assemblyFiles)
{
	if($asmFile =~ /\.h1-/)
	{
		$assemblyFile_2_H{$asmFile} = 1;
		$H_2_assemblyFile{1} = $asmFile;
	}
	
	if($asmFile =~ /\.h2-/)
	{
		die if($assemblyFile_2_H{$asmFile});
		$assemblyFile_2_H{$asmFile} = 2;
		$H_2_assemblyFile{2} = $asmFile;		
	}	
}

my $inferredSequences_href = readFASTA($fastaFile);

my %inferredGenes = map {my $seqID = $_; die unless($seqID =~ /^(.+)\-((H1)|(H2))$/); my $gene = $1; $gene => 1} keys %$inferredSequences_href;
warn Dumper("Some genes have only one haplotype in $fastaFile", [grep {!((exists $inferredSequences_href->{$_ . '-H1'}) and (exists $inferredSequences_href->{$_ . '-H2'}))} keys %inferredGenes]) unless(all {(exists $inferredSequences_href->{$_ . '-H1'}) and (exists $inferredSequences_href->{$_ . '-H2'})} keys %inferredGenes);

my %bestMappings;
foreach my $gene (keys %inferredGenes)
{
	my $gene_ref_fn = $graph_dir . '/pseudoGenomic_fullLengthMapping/raw_HLA-' . $gene . '.fa';
	die "Expected file $gene_ref_fn not found" unless(-e $gene_ref_fn);
	
	my $paf = qq(${fastaFile}.${gene}.paf);
	if(-e $paf)
	{
		unlink($paf) or die "Cannot unlink $paf";
	}
	my $cmd_map = qq($minimap2_bin -c $gene_ref_fn $fastaFile > $paf);
	print "Now mapping against $gene_ref_fn...\n";
	system($cmd_map) and die "Could not execute: $cmd_map";
	
	parsePAF($paf, \%bestMappings, 'pseudoGenomic_fullLengthMapping');
}

my $file0_h = $assemblyFile_2_H{$assemblyFiles[0]};
my $file1_h = $assemblyFile_2_H{$assemblyFiles[1]};
die if($file0_h eq $file1_h);

my $paf_asm_1 = qq(${fastaFile}.1.paf);
my $cmd_1 = qq($minimap2_bin -c $assemblyFiles[0] $fastaFile > $paf_asm_1);
print "Now mapping against $assemblyFiles[0]...\n";
system($cmd_1) and die "Could not execute: $cmd_1";
parsePAF($paf_asm_1, \%bestMappings, 'assembly-H' . $file0_h );

my $paf_asm_2 = qq(${fastaFile}.2.paf);
my $cmd_2 = qq($minimap2_bin -c $assemblyFiles[1] $fastaFile > $paf_asm_2);
print "Now mapping against $assemblyFiles[1]...\n";
system($cmd_2) and die "Could not execute: $cmd_2";
parsePAF($paf_asm_2, \%bestMappings, 'assembly-H' . $file1_h );

open(OUT, '>', $outputFn) or die "Cannot open $outputFn";
my @dataSources = ('pseudoGenomic_fullLengthMapping', 'assembly-H1', 'assembly-H2');
print OUT join("\t", "", @dataSources), "\n";
foreach my $gene (sort keys %inferredGenes)
{
	die "Results (H1) for $gene entirely missing" unless(exists $bestMappings{$gene . '-H1'});
	die "Results (H2) for $gene entirely missing" unless(exists $bestMappings{$gene . '-H2'});
	
	foreach my $H ('H1', 'H2')
	{
		my $resultsKey = $gene . '-' . $H;
		my @evaluationResults = ($resultsKey);
		
		foreach my $dS (@dataSources)
		{
			push(@evaluationResults, 
				join('  ', map {$_->[0] . " (" . sprintf("%.2f", $_->[1] * 100) . ' / ' . $_->[3] . ' NM)'} @{$bestMappings{$resultsKey}{$dS}})
			);
		}
		
		print OUT join("\t", @evaluationResults), "\n";
	}
}

close(OUT);

print "\n\nDone. Generated $outputFn\n\n";

sub parsePAF
{
	my $paf = shift;
	my $results_href = shift;
	my $results_prefix = shift;
	die unless(defined $results_prefix);
	
	my %bestResults;
	
	open(PAF, '<', $paf) or die "Cannot open $paf";
	while(<PAF>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		my $query_name = $line_fields[0];
		my $query_length = $line_fields[1];
		my $query_matched = $line_fields[9];
		my $matched_perc = $query_matched / $query_length;
		die unless($line =~ /NM:i:(\d+)/);
		my $nm = $1;
		
		die unless(($matched_perc >= 0) and ($matched_perc) <= 1);
		
		if((not exists $bestResults{$query_name}) or ($bestResults{$query_name}[0][1] < $matched_perc))
		{
			$bestResults{$query_name} = [[$line_fields[5], $matched_perc, $query_matched, $nm]];
		}
		elsif($bestResults{$query_name}[0][1] == $matched_perc)
		{
			push(@{$bestResults{$query_name}}, [$line_fields[5], $matched_perc, $query_matched, $nm]);  
		}
	}
	close(PAF);
	
	foreach my $queryID (keys %bestResults)
	{
		$results_href->{$queryID}{$results_prefix} = $bestResults{$queryID};
	}	
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


