use strict;
use warnings;

while(my $sam = <STDIN>) # look at file line by line
	{
	if(substr($sam, 0, 1) eq "@")
		{
		print STDOUT $sam;
		next;
		}
	my @entry = split("\t", $sam);
	my $CIGAR = $entry[5];

	if($CIGAR =~ /N/) # if there is an N in the CIGAR string
		{
		print STDOUT $sam;
		}
	}
