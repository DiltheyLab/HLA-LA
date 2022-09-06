package GGroups;

use strict;

sub readG
{
	my $inputFile = shift;
	my $G_full_alleles_href = shift;
	my $G_to_underlyingAlleles_href = shift;
	my $alleles_to_fullGAmbiguity_href = shift;
	my $G_mapper_unambigious_href = shift;
	
	die unless(-e $inputFile);
	die unless($G_full_alleles_href);
	die unless($G_to_underlyingAlleles_href);
	die unless($G_mapper_unambigious_href);
	
	%$G_full_alleles_href = ();
	%$G_to_underlyingAlleles_href = ();
	%$alleles_to_fullGAmbiguity_href = ();
	%$G_mapper_unambigious_href = ();
	
	open(GFILE, '<', $inputFile) or die "Cannot open $inputFile";
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
			$G_full_alleles_href->{'HLA'.$locus}{$_}++;
		}
		
		$G_full_alleles_href->{'HLA'.$locus}{$G_group}++;	
				
		$G_to_underlyingAlleles_href->{'HLA'.$locus}{$G_group} = \@alleles;

		$alleles_to_fullGAmbiguity_href->{'HLA'.$locus}{$G_group} = \@alleles;
		for(@alleles)
		{
			$alleles_to_fullGAmbiguity_href->{'HLA'.$locus}{$_} = \@alleles;		
		}
		
		my $alleleMaxIdx = $#alleles;	
		for(my $i = 0; $i <= $alleleMaxIdx; $i++)
		{
			my $allele = $alleles[$i];
			
			die if(exists $G_mapper_unambigious_href->{$locus . '*' . $allele});
			$G_mapper_unambigious_href->{$locus . '*' . $allele} = $locus . '*' .$G_group;		
			
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


1;
