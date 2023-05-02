use strict;
use warnings;

while(my $pileup = <STDIN>) # look at file line by line
	{
	my @entry = split("\t", $pileup);
	my $nt_ref = $entry[2];
	if($nt_ref eq "t" or $nt_ref eq "T")
		{	
		my $nt_read = $entry[4];
		my $c_count = $nt_read =~ tr/Cc//;
		if( $c_count > 0)
			{
			my $cov = $entry[3];
			if($c_count >= $cov/5 and $c_count >= 5)
				{ 
				my $ref = $entry[0];
				my $pos = $entry[1];
				print STDOUT $ref . ":" . $pos . "-" . $pos;
				print STDOUT ";";
				print STDOUT $c_count . "; " . $cov . "\n";
				}
			}
		}
	}


