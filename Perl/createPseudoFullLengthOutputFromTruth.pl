use warnings;
use strict;
use FindBin;
use Getopt::Long;
use Data::Dumper; 
use List::MoreUtils qw/all mesh/;

my $this_bin_dir = $FindBin::RealBin;
my $graphs_dir = $this_bin_dir . '/../../graphs/';

my $g_file = $this_bin_dir . '/../hla_nom_g.txt';
die "Missing file: $g_file" unless(-e $g_file);

my $truthFile = '/home/dilthey/HLAtypes.txt';
my $sampleID = 'NA12878';
my $prefixOutput = '/home/dilthey/HLA-LA-devel/test/mockFromTruth';
my $graph = 'PRG_MHC_GRCh38_withIMGT';
GetOptions (
	'truthFile:s' => \$truthFile,
	'sampleID:s' => \$sampleID,
	'prefixOutput:s' => \$prefixOutput,
	'graph:s' => \$graph,
);
	
my $graph_dir = $graphs_dir . '/' . $graph;

my %G_full_alleles;
my %G_mapper_unambigious;
my %G_mapper_multiples;
my %G_to_underlyingAlleles;
my %alleles_to_fullGAmbiguity;
readG();

my %hlaTypes_truth;
open(TRUTH, '<', $truthFile) or die "Cannot open $truthFile";
my $truth_header = <TRUTH>;
chomp($truth_header);
$truth_header =~ s/[\n\r]//g;
my @truth_header_fields = split(/\t/, $truth_header, -1);
my @hla_genes = grep {$_ =~ /^HLA/} @truth_header_fields; 
die unless(scalar(@hla_genes));

my %pseudoGenomicSequences;
foreach my $gene (@hla_genes)
{
	my $locus_for_seqLookup = $gene;
	$locus_for_seqLookup =~ s/HLA//;	
	my $gene_ref_fn = $graph_dir . '/pseudoGenomic_fullLengthMapping/raw_HLA-' . $locus_for_seqLookup . '.fa';
	die "Expected file $gene_ref_fn (from $gene / $locus_for_seqLookup) not found" unless(-e $gene_ref_fn);	 
	$pseudoGenomicSequences{$gene} = readFASTA($gene_ref_fn);
}

while(<TRUTH>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/[\n\r]//g;
	my @line_fields = split(/\t/, $line, -1);
	die unless(scalar(@truth_header_fields) == scalar(@line_fields));
	my %line_hash = (mesh @truth_header_fields, @line_fields);
	die "Missing IndividualID field in $truthFile" unless($line_hash{IndividualID});
	
	foreach my $gene (@hla_genes)
	{
		my @alleles = split(/\//, $line_hash{$gene});
		unless(scalar(@alleles) == 2)
		{
			die "Wrong format - expect 2 fields for gene $gene, separated by '/'";
		}
		my @alleles_tr;
		foreach my $allele (@alleles)
		{
			my $locus_for_gGroupLookup = $gene;
			$locus_for_gGroupLookup =~ s/HLA-/HLA/;
			my $allele_for_lookup = $allele;
			$allele_for_lookup =~ s/^\w+\*//;
			die "Unknown locus $line_fields[0]" unless(exists $G_to_underlyingAlleles{$locus_for_gGroupLookup});
			die Dumper("Can't resolve allele $allele / $allele_for_lookup at locus $gene at G-group resolution", $line_hash{$gene}, \@alleles) unless(defined $G_to_underlyingAlleles{$locus_for_gGroupLookup}{$allele_for_lookup});		
			my $alleles_tr = $G_to_underlyingAlleles{$locus_for_gGroupLookup}{$allele_for_lookup};
			#$alleles_tr =~ s/^\S+\*//;
			push(@alleles_tr, $alleles_tr);						
		}
		die unless(scalar(@alleles_tr) == 2);
		$line_hash{$gene} = \@alleles_tr;
	}
	$hlaTypes_truth{$line_hash{IndividualID}} = \%line_hash;
}
close(TRUTH);

unless(exists $hlaTypes_truth{$sampleID})
{
	die "Missing truth types for sample '$sampleID' in file $truthFile";
}

my %outputFasta;
foreach my $gene (@hla_genes)
{
	die unless(defined $hlaTypes_truth{$sampleID}{$gene});
	
	my @alleles = @{$hlaTypes_truth{$sampleID}{$gene}};
	die Dumper("Wrong number of alleles for $gene", \@alleles) unless(scalar(@alleles) == 2);
	
	my $translate_gene_alleleList = sub {
		my $alleleList_aref = shift;
		print "Gene $gene, translating list of " . scalar(@$alleleList_aref) . " alleles.\n";
		(my $gene_for_rawSequences = $gene) =~ s/HLA//;
		my @alleles_present = grep {exists $pseudoGenomicSequences{$gene}{$_}} map {$gene_for_rawSequences . '*' . $_} @$alleleList_aref;
		unless(scalar(@alleles_present))
		{
			warn "Cannot translate: $gene [" . join(', ', @$alleleList_aref) . "]";
			return undef;
		}
		print "\t" . scalar(@alleles_present) . " alternatives; choose $alleles_present[0]\n";
		return $alleles_present[0];
	};

	my $gene_for_fastaOutput = $gene;
	$gene_for_fastaOutput =~ s/HLA//;
			
	for(my $h = 0; $h < 2; $h++)
	{
		my $allele_list_aref = $alleles[$h];
		my $allele_use = $translate_gene_alleleList->($allele_list_aref);
		next unless(defined $allele_use);
		my $seq_id_for_output = $gene_for_fastaOutput . '-H' . ($h+1);
		$outputFasta{$seq_id_for_output} = $pseudoGenomicSequences{$gene}{$allele_use};
	}
}


my $outputFile = $prefixOutput . '.fullLengthInference.fasta';
writeFASTA($outputFile, \%outputFasta);

print "\n\nDone. Produced file $outputFile\n\n";


		
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



sub allele_list_to_G
{
	my $alleleList_aref = shift;
	my $locus = shift;
	
	die unless(keys %G_mapper_unambigious);
	my %g;
	foreach my $allele (@$alleleList_aref)
	{
		my $g = $G_mapper_unambigious{$allele};
		if(not defined $g)
		{
			if(defined $locus)
			{
				(my $locus_noHLA = $locus) =~ s/HLA-//;
				$g = $G_mapper_unambigious{$locus_noHLA . '*' . $allele};			
			}
			die Dumper("Can't translate allele $allele", ["Examples in the right format", [(keys %G_mapper_unambigious)[0 .. 9]]], $locus) if(not defined $g);
		}
			
		$g{$g}++;
			
	}
	return sort keys %g;
}

sub readG
{
	open(GFILE, '<', $g_file) or die "Cannot open $g_file";
	while(<GFILE>)
	{
		my $line = $_;
		chomp($line);
		next if(substr($line, 0, 1) eq '#');
		die unless($line =~ /^(.+?)\*;(.+)$/);
		my $locus = $1;
		my $alleles = $2;
		
		my $G_group;
		if($alleles !~ /;$/)
		{
			die "Weird alleles for $locus: '$alleles'" unless($alleles =~ /^(.+);(.+?G)/);
			$alleles = $1;
			$G_group = $2;
		}
		else
		{
			
			$alleles = substr($alleles, 0, length($alleles)-1);
			die if($alleles =~ /\//);
			die if($alleles =~ /G/);
			$G_group = $alleles;
		}

		my @alleles = split(/\//, $alleles);
		
		for(@alleles)
		{
			$G_full_alleles{'HLA'.$locus}{$_}++;
		}
		
		$G_full_alleles{'HLA'.$locus}{$G_group}++;	
				
		$G_to_underlyingAlleles{'HLA'.$locus}{$G_group} = \@alleles;

		$alleles_to_fullGAmbiguity{'HLA'.$locus}{$G_group} = \@alleles;
		for(@alleles)
		{
			$alleles_to_fullGAmbiguity{'HLA'.$locus}{$_} = \@alleles;		
		}
		
		my $alleleMaxIdx = $#alleles;	
		for(my $i = 0; $i <= $alleleMaxIdx; $i++)
		{
			my $allele = $alleles[$i];
			
			die if(exists $G_mapper_unambigious{$locus . '*' . $allele});
			$G_mapper_unambigious{$locus . '*' . $allele} = $locus . '*' .$G_group;		
			
			# die Dumper(\%G_mapper_unambigious);
			
			#$G_mapper_unambigious{'HLA'.$locus}{$allele} = $G_group;
			
			# $G_mapper_multiples{'HLA'.$locus}{$allele} = [$G_group];
			
			# my @allele_parts = split(/:/, $allele);
			# die "Weird alllele II: $allele $alleles $line" unless($#allele_parts >= 1);
			
			# my $fourDigit = join(':', @allele_parts[0, 1]);
					
			# if(exists $G_mapper_unambigious{'HLA'.$locus}{$fourDigit})
			# {
				# my $existing_G_group = $G_mapper_unambigious{'HLA'.$locus}{$fourDigit};
				# if($existing_G_group ne $G_group)
				# {
					# # warn "$locus $allele $fourDigit existing: $existing_G_group - now want to set $G_group";
					# $G_mapper_unambigious{'HLA'.$locus}{$fourDigit} = undef;
				# }
				
				# my @existing_multiples = @{$G_mapper_multiples{'HLA'.$locus}{$fourDigit}};
				# my %_existing_multiples = map {$_ => 1} @existing_multiples;
				# if(not $_existing_multiples{$G_group})
				# {
					# push(@{$G_mapper_multiples{'HLA'.$locus}{$fourDigit}}, $G_group);
				# }
				
			# }
			# else
			# {
				# $G_mapper_unambigious{'HLA'.$locus}{$fourDigit} = $G_group;
				# $G_mapper_multiples{'HLA'.$locus}{$fourDigit} = [$G_group];
			# }
		}
	}
	close(GFILE);
}


