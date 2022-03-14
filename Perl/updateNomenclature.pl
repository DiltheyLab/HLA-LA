use warnings;
use strict;
use FindBin;
use Getopt::Long;
use Data::Dumper; 
use List::MoreUtils qw/all mesh/;
use lib '.';

use Util;

my $this_bin_dir = $FindBin::RealBin;
my $graphs_dir = $this_bin_dir . '/../../graphs/';

my $g_file = $this_bin_dir . '/../hla_nom_g.txt';
die "Missing file: $g_file" unless(-e $g_file);

my $hla_nuc = '../hla_nuc.fasta';
my $graph = 'PRG_MHC_GRCh38_withIMGT';
GetOptions (
	'hla_nuc:s' => \$hla_nuc,
	'graph:s' => \$graph,
);
	
my $graph_dir = $graphs_dir . '/' . $graph;

unless($hla_nuc and (-e $hla_nuc))
{
	die "No IMGT/HLA nucleotide file specified (--hla_nuc), or specified file '$hla_nuc' not present";
}

my $fn_sequences = $graph_dir . '/sequences.txt';
die "Missing file $fn_sequences" unless(-e $fn_sequences); 

my $fn_sequences_old = $fn_sequences . '.old';
my $fn_sequences_new = $fn_sequences . '.new';
if(-e $fn_sequences_old)
{
	die "File $fn_sequences_old existing, abort";
}

my %knownAlleles;
my $nuc_fasta_href = Util::readFASTA($hla_nuc, 0);
foreach my $alleleID (keys %$nuc_fasta_href)
{
	my @alleleID_parts = split(/ /, $alleleID);
	die unless(scalar(@alleleID_parts) == 4);
	$knownAlleles{$alleleID_parts[1]}++;
}

system("cp $fn_sequences $fn_sequences_old") and die "sequences.txt backup command filed";

my $updates_sequence_txt = 0;
my $unknownAlleles_in_sequence_txt = 0;
my %allele_translation;
open(SEQUENCES_NEW, '>', $fn_sequences_new) or die "Cannot open $fn_sequences_new";
open(SEQUENCES_OLD, '<', $fn_sequences_old) or die "Cannot open $fn_sequences_old";
my $sequences_headerLine = <SEQUENCES_OLD>;
print SEQUENCES_NEW $sequences_headerLine;
chomp($sequences_headerLine);
my @sequences_headerFields = split(/\t/, $sequences_headerLine);
die unless($sequences_headerFields[1] eq 'Name');
while(<SEQUENCES_OLD>)
{
	my $line = $_;
	chomp($line);
	my @line_fields = split(/\t/, $line, -1);
	my $allele = $line_fields[1];
	if($allele =~ /^(\w+?)\*(.+)$/)
	{
		my $locus = $1;
		my $allele_numeric = $2;
		unless(exists $knownAlleles{$allele})
		{
			# print "Unknwon allele $allele";
			my $extendedAllele = $allele . ':01';
			my $extendedAllele_II = $allele . ':01:01';
			if(exists $knownAlleles{$extendedAllele})
			{
				# print " -> $extendedAllele"; 
				$allele_translation{$allele} = $extendedAllele;
				$line_fields[1] = $extendedAllele;
				$updates_sequence_txt++;
			}
			elsif(exists $knownAlleles{$extendedAllele_II})
			{
				# print " -> $extendedAllele_II"; 
				$allele_translation{$allele} = $extendedAllele_II;
				$line_fields[1] = $extendedAllele;
				$updates_sequence_txt++;
			}
			else
			{
				$unknownAlleles_in_sequence_txt++;
			}
			
			# print "\n";
		}
	#	die "Unknown locus $locus" unless(exists $G_full_alleles{'HLA' . $locus});
	}
	print SEQUENCES_NEW join("\t", @line_fields), "\n";
}
close(SEQUENCES_NEW);
close(SEQUENCES_OLD);

print "\n\nGenerated file: $fn_sequences_new - $updates_sequence_txt updates - $unknownAlleles_in_sequence_txt unknown in the IMGT input data\n\n";
my @copy_tuples = ([$fn_sequences_new, $fn_sequences]);
my @prg_files = grep {$_ =~ /(exon_\d+\.txt)|(intron\d+\.txt)/} glob("${graph_dir}/PRG/*.txt");
foreach my $prg_file_old (@prg_files)
{
	my $n_alleles_unchanged = 0;
	my $n_updates_prgFile = 0;
	my $n_alleles_unknown = 0;
	print "\tUpdating $prg_file_old ";
	my $prg_file_new = $prg_file_old . '.new';
	open(PRG_NEW, '>', $prg_file_new) or die "Cannot open $fn_sequences_new";
	open(PRG_OLD, '<', $prg_file_old) or die "Cannot open $prg_file_old";
	my $prg_headerLine = <PRG_OLD>;
	print PRG_NEW $prg_headerLine;
	while(<PRG_OLD>)
	{
		my $line = $_;
		if($line =~ /^(\S+) /)
		{
			my $allele = $1;
			
			my $extendedAllele = $allele . ':01';
			my $extendedAllele_II = $allele . ':01:01';
			my $updateAllele;
			if(exists $knownAlleles{$allele})
			{
				$n_alleles_unchanged++;
			}
			elsif(exists $knownAlleles{$extendedAllele})
			{
				$updateAllele = $extendedAllele;
			}
			elsif(exists $knownAlleles{$extendedAllele_II})
			{
				$updateAllele = $extendedAllele_II;
			}
			else
			{
				$n_alleles_unknown++;
			}
			
			if($updateAllele)
			{
				$n_updates_prgFile++;
				$line =~ s/^(\S+) /${updateAllele} /;
			}
		}
		print PRG_NEW $line;
	}
	close(PRG_OLD);
	close(PRG_NEW);
	print " -> $n_alleles_unchanged, $n_updates_prgFile updates, $n_alleles_unknown unknown.\n";
	push(@copy_tuples, [$prg_file_new, $prg_file_old]);
}

foreach my $cT (@copy_tuples)
{
	my $copy_cmd = qq(cp $cT->[0] $cT->[1]);
	system($copy_cmd) and die "Could not execute: '$copy_cmd'";
}

print "\nDone. Updated " . scalar(@copy_tuples) . " files.\n\n";



