package simpleHLA;

sub modernHLA_is_missing
{
	my $hla = shift;
	return($hla =~ /\?/);
}

sub is_missing
{
	my $hla = shift;
	
	if(ref($hla) eq 'ARRAY')
	{
		my $v1 = $hla->[0];
		my $v2 = $hla->[1];
		
		my $m1 = 0, my $m2 = 0;
		
		
		if ($v1 =~ /\?/)
		{
			$m1 = 1;
		}
		
		if ($v1 eq '9999')
		{
			$m1 = 1;
		}
		
		if (! $v1)
		{
			$m1 = 1;
		}
		
		
		
		if ($v2 =~ /\?/)
		{
			$m2 = 1;
		}
		
		if ($v2 eq '9999')
		{
			$m2 = 1;
		}
		
		if (! $v2)
		{
			$m2 = 1;
		}	
		
		
		unless($m1 == $m2)
		{
			die "Inconsistency - both alleles should either be known or unknown";
		}	
		
		return $m1;
				
	}
	else
	{
		if ($hla =~ /\?/)
		{
			return 1;
		}
		
		if ($hla eq '9999')
		{
			return 1;
		}
		
		if (! $hla)
		{
			return 1;
		}
		
		if ($hla =~ /^(0+)$/)
		{
			return 1;
		}
	}
	
	return 0;
}

sub is_24compatible
{
	my $old = shift;
	my $new = shift;
	
	$old = &HLA_4digit($old);
	$new = &HLA_4digit($new);
	
	if(&HLA_is4digit($old))
	{
		if(&HLA_is4digit($new))
		{
			return &is_compatible($old, $new);
		}
		elsif(&HLA_is2digit($new))
		{
			return 0;
			#die "Old $old, new $new: resolution should improve, at least!";
		}
		else
		{
			die;
		}
	}
	elsif(&HLA_is2digit($old))
	{
		if(&HLA_is2digit($new))
		{
			return &is_compatible($old, $new);
		}
		elsif(&HLA_is4digit($new))
		{
			my $d2_new = &HLA_reduce_to_2($new);
			return &is_compatible($old, $d2_new);
		}
		else
		{
			die;
		}
	}
	else
	{
		die;
	}
}

sub is_compatible
{
	my $one = shift;
	my $two = shift;
	
	if(ref($one) eq 'ARRAY')
	{
		(ref($two) eq 'ARRAY') || die "Both arrays, please!".Dumper($one, $two);
		
		if(&is_missing($one) and &is_missing($two))
		{
			return 1;
		}
		
		my $v1_1 = &HLA_4digit($one->[0]);
		my $v1_2 = &HLA_4digit($one->[1]);
		
		my $v2_1 = &HLA_4digit($two->[0]);
		my $v2_2 = &HLA_4digit($two->[1]);	
		
		return 
			(
			 	(($v1_1 eq $v2_1) and ($v1_2 eq $v2_2))
				or
				(($v1_1 eq $v2_2) and ($v1_2 eq $v2_1))
			);	
				
	}
	else
	{
		my $v1 = &HLA_4digit($one);
		my $v2 = &HLA_4digit($two);
		
		if(&is_missing($v1) and &is_missing($v2))
		{
			return 1;
		}
		
		return ($v1 eq $v2);
	}
}


sub autoHLA_2digit
{
	my $hla = shift;
	if($hla =~ /:/)
	{
		return modernHLA_2digit($hla);
	}
	else
	{
		return HLA_2digit($hla);
	}
}

sub modernHLA_4digit
{
	my $hla = shift;
	die unless($hla);
	my @elements = split(/:/, $hla);
	foreach my $e (@elements)
	{
		die if(length($e) < 2);
		die "Problem with HLA $hla" unless($e =~ /^\d\d/);
		die "Problem with HLA $hla $e" if(length($e) > 3);
	}
	
	die unless(scalar(@elements) >= 2);
	
	return join(':', @elements[0, 1]);
}


sub modernHLA_2digit
{
	my $hla = shift;
	die unless($hla);
	$hla =~ s/N//;
	$hla =~ s/Q//;
	my @elements = split(/:/, $hla);
	foreach my $e (@elements)
	{
		die if(length($e) < 2);
		if($e eq $elements[0])
		{  
			die "Problem with HLA $hla $e" unless($e =~ /^(\w+\*)?(\d\d)/);
			die "Problem with HLA $hla $e" if(length($2) > 3);

		}
		else
		{
			die "Problem with HLA $hla $e" unless($e =~ /^\d\d/);
			die "Problem with HLA $hla $e" if(length($e) > 3);			
		}
		
	}
	
	die unless(scalar(@elements) >= 2);
	
	if(scalar(@elements) == 1)
	{
		push(@elements, '00');
	}
	else
	{
		$elements[1] = '00';
		@elements = @elements[0, 1];
	}
	
	return join(':', @elements[0, 1]);
}

sub HLA_2digit
{
	my $hla = shift;
	
	my $four_dig = &HLA_4digit($hla);
	if(is_missing($four_dig))
	{
		return '????';
	}
	
	return HLA_reduce_to_2($hla);
	
	# return substr($four_dig, 0, 2).'00';
}

sub HLA_is_g
{
	my $hla = shift;
	
	if((length($hla) == 5) and ((substr($hla, 4, 1) eq 'g') or (substr($hla, 4, 1) eq 'G')))
	{	
		die unless(HLA_4digit(substr($hla, 0, 4)));
		return 1;
	}
	else
	{
		return 0;
	}
}	

sub HLA_strip_g
{
	my $hla = shift;
	die unless(HLA_is_g($hla));
	return substr($hla, 0, 4);
}	

sub HLA_4digit
{
	my $hla = shift;
	
	if($hla =~ /\:/)
	{
		my @components = split(/\:/, $hla);
		die unless($#components >= 1);
		if($#components == 0)
		{
			return $components[0].':00';
		}
		else
		{
			return $components[0].':'.$components[1];
		}
	}
	if(is_missing($hla))
	{
		return '????';
	}
	
	if(HLA_is_g($hla))
	{
		return $hla;
	}
	elsif(length($hla) == 1)
	{
		return '0'.$hla.'00';
	}
	elsif(length($hla) == 2)
	{
		return $hla.'00';
	}
	elsif(length($hla) == 3)
	{
		return '0'.$hla;
	}
	elsif(length($hla) == 4)
	{
		return $hla;
	}
	elsif($hla =~ /(\w*)\*(\d+?)\:(\d+)/)
	{
		$hla = $1 . '*' . $2 . ':' . $3;
	}
	else
	{
		print "This " . length($hla) . " shall be an HLA code? ---- $hla ----\n\n\n\n\n";
		die "This  " . length($hla) . " shall be an HLA code? $hla\n\n\n\n";
	}
}

sub HLA_reduce_to_2
{
	my $hla = shift;
	my $four_digit = &HLA_4digit($hla);
	
	if(HLA_is_g($four_digit))
	{
		return HLA_reduce_to_2(HLA_strip_g($four_digit));
	}
	elsif($four_digit =~ /(\w*)\*(\d+?)\:(\d+)/)
	{
		return $1 . '*' . $2 . ':' . '00';	
	}
	else
	{
		return substr($four_digit, 0, 2).'00';
	}
}


sub autoHLA_is2digit
{
	my $hla = shift;
	if($hla =~ /:/)
	{
		return modernHLA_is2digit($hla);
	}
	else
	{
		return HLA_is2digit($hla);
	}
}


sub HLA_is2digit
{
	my $hla = shift;
	$hla = &HLA_4digit($hla);
	return (not HLA_is_g($hla)) && (substr($hla, 2, 2) eq '00');	
}

sub modernHLA_is2digit
{
	my $hla = shift;
	my @components = split(/\:/, $hla);
	die unless($#components >= 1);
	return($components[1] =~ /^0+$/);		
}

sub HLA_is4digit
{
	my $hla = shift;
	if($hla =~ /\:/)
	{
		my @components = split(/\:/, $hla);
		if($#components >= 1)
		{
			return ($components[1] !~ /^0+$/);
		}
		else
		{
			return 0;
		}
	}
	else
	{
		$hla = &HLA_4digit($hla);
		if(HLA_is_g($hla))
		{
			return HLA_is4digit(HLA_strip_g($hla));
		}
		else
		{
			return ((!is_missing($hla)) and (substr($hla, 2, 2) ne '00'));	
		}
	}
}

sub HLA_determineInformationLevel
{
	my $hla = shift;
	$hla = &HLA_4digit($hla);
	
	if(is_missing($hla))
	{
		return 0;
	}
	
	if(&HLA_is2digit($hla))
	{
		return 2;
	}
	
	if(&HLA_is4digit($hla))
	{
		return 4;
	}
	
	die "Cannot determine HLA information level of $hla!\n";
	
	
}

1;