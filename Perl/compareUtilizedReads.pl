#!/usr/bin/perl

use strict;
use Data::Dumper;

# my $sample_new = "NA19240_1000G";
# my $sample_old = "I3_NA19240";
my $sample_new = "NA18939_1000G_RED";
my $sample_old = "I3_NA18939";
my @loci = qw/A/;

my $dir_new = qq(/gpfs1/well/gsk_hla/MHC-PRG-2/working/$sample_new/hla/);
my $dir_old = qq(/Net/birch/data/dilthey/MHC-PRG/tmp/hla/$sample_old/);

print "\n", join(' ', qw/New Shared Old/), "\n\n";

my %all_new;
my %all_old;

foreach my $locus (@loci)
{
	my %readIDs_new;
	my %readIDs_old;
	
	my $f_new = $dir_new . '/R1_readIDs_' . $locus . '.txt';
	my $f_old = $dir_old . '/R1_readIDs_' . $locus . '.txt';
	
	open(NEW, '<', $f_new) or die "Cannot open $f_new";
	while(<NEW>)
	{
		my $line = $_;
		chomp($line);
		$readIDs_new{$line}++;
		$all_new{$line}++;
	}
	close(NEW);

	open(OLD, '<', $f_old) or die "Cannot open $f_old";
	while(<OLD>)
	{
		my $line = $_;
		chomp($line);
		die unless($line =~ /\@\@(.+?):normal/);
		my $readID = $1;
		$readIDs_old{$readID}++;
		$all_old{$readID}++;
	}
	close(OLD);
	
	print $locus, "\n";
	print "\t", join(' ', intersectionStats([keys %readIDs_new], [keys %readIDs_old])), "\n";
}

print "\nAll\n";
print "\t", join(' ', intersectionStats([keys %all_new], [keys %all_old])), "\n";

sub intersectionStats
{
	my $aref_1 = shift;
	my $aref_2 = shift;
	my %h_1 = map {$_ => 1} @$aref_1;
	my %h_2 = map {$_ => 1} @$aref_2;
	
	my @excl_1 = grep {not exists $h_2{$_}} keys %h_1;
	my @shared = grep {exists $h_2{$_}} keys %h_1;
	my @excl_2 = grep {not exists $h_1{$_}} keys %h_2;
	
	my @forReturn = (scalar(@excl_1), scalar(@shared), scalar(@excl_2));

	return @forReturn;
}