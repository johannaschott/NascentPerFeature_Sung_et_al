use strict;
use warnings;

my $MD_field = $ARGV[0] - 1; # field in bam file that contains MD tag
my $featureCounts = $ARGV[1]; # name of the corresponding  featureCounts file from command-line
my $trans_from = $ARGV[2]; # nucleotide that is converted (T for T>C transitions)
my $trans_to = $ARGV[3]; # result of conversion (C for T>C transition)
my $strandedness = $ARGV[4]; # strandedness of library (forward or reverse)
my $max_ref = $ARGV[5]; # maximum number of covered Ts in reference to consider
my $max_trans = $ARGV[6]; # maximum number of T>C transitions to consider


open(FILE, $featureCounts); # open featureCounts file 

my %result; # generate array for the result of gene-wise transition counts

while(my $line = <FILE> and my $bam = <STDIN>) # bam file content is read from STDIN; parse through featureCounts output and bam file in parallel
	{
	
	my @line = split("\t", $line);
	my $status = $line[1]; # extracts status of the read featureCounts file (e.g. "Assigned")
	if($status ne "Assigned") # if read was not assigned to any feature...
		{
		next; # ... skip read.
		}
	
	my $feature = $line[2]; # extracts name of the feature (e.g. gene ID) from featureCounts file


########## If a T position in the read is among the mismatching positions of the MD tag, don't count it for the matching Ts.
########## If a C position in the read is among the mismatching T positions of the MD tag, count it for the T>C transitions.
########## Total Ts are the sum of matches and transitions.

	my @bam = split("\t", $bam); 
	my $MD_tag = $bam[$MD_field]; # MD tag is extracted from the field given as command line argument. 
	my $CIGAR = $bam[5]; # CIGAR string is 5th field in bam file.
	my $read = $bam[9]; # Read sequence is 10th field in bam file.
	my $flag = $bam[1]; # FLAG is second field in the bam file.

	# If the library is reverse-stranded, transitions in reads aligning to the plus-strand (FLAG = 0) need to be complemented 
	# (e.g. a T>C transition becomes an A>G transition). If the library is forward-stranded, reads aligning to the 
	# minus-strand (FLAG = 16) need to be complemented.
	my $trans_from_plus;
	my $trans_to_plus;
	my $trans_from_minus;
	my $trans_to_minus;	

	if($strandedness eq "reverse")
		{
		$trans_from_minus = $trans_from;
		$trans_to_minus = $trans_to;	
		$trans_from_plus = ($trans_from =~ tr/GATC/CTAG/);
		$trans_to_plus = ($trans_to =~ tr/GATC/CTAG/);
		}  

	if($strandedness eq "forward")
		{
		$trans_from_plus = $trans_from;
		$trans_to_plus = $trans_to;
		$trans_from_minus = ($trans_from =~ tr/GATC/CTAG/);
		$trans_to_minus = ($trans_to =~ tr/GATC/CTAG/);	
		}  

	my $from;
	my $to;

	if($flag == 16)
		{
		$from = $trans_from_minus;
		$to = $trans_to_minus;
		}

	if($flag == 0)
		{
		$from = $trans_from_plus;
		$to = $trans_to_plus;
		}

	my $trans_number;
	my $ref_number;
	($trans_number, $ref_number) = count_transitions($trans_from, $trans_to, $read, $CIGAR, $MD_tag);

	if(exists $result{$feature}{$ref_number}{$trans_number})
		{
		$result{$feature}{$ref_number}{$trans_number}++;
		}else{
		$result{$feature}{$ref_number}{$trans_number} = 1;
		}
	}



close(FILE);
		
# print header line:
print STDOUT "feature;";
foreach my $ref_number (0..$max_ref)
	{
		
	# The number of transitions cannot be higher than the number of target bases in the reference.
	my $max_trans_printed;
	if($max_trans > $ref_number)
	{
		$max_trans_printed = $ref_number;
	}else{
		$max_trans_printed = $max_trans;
	}
		
	foreach my $trans_number (0..$max_trans_printed)
		{
		print STDOUT "of_" . $ref_number . "_" . $trans_number;

		if( $ref_number == $max_ref and $trans_number == $max_trans_printed)	
			{
			print STDOUT "\n";
			}else{
			print STDOUT ";";
			}

		}
		
	}

# print results of transition count (one line per gene):
foreach ( keys %result ) 
	{
	print STDOUT "$_;";
	foreach my $ref_number (0..$max_ref)
		{	
		# The number of transitions cannot be higher than the number of target bases in the reference.
		my $max_trans_printed;
		if($max_trans > $ref_number)
		{
			$max_trans_printed = $ref_number;
		}else{
			$max_trans_printed = $max_trans;
		}
		
		foreach my $trans_number (0..$max_trans_printed)  
			{
			if(exists $result{$_}{$ref_number}{$trans_number})
				{
				print STDOUT "$result{$_}{$ref_number}{$trans_number}";
				}else{
				print STDOUT "0";
				}
			if( $ref_number == $max_ref and $trans_number == $max_trans_printed)	
				{
				print STDOUT "\n";
				}else{
				print STDOUT ";";
				}		
			}
		}
	}		

##############################SUBROUTINES##################################

######### Purpose: 
######### This subroutine counts transitions from a specific nucleotide to another
######### between a read sequence and its alignment to a reference sequence. In addition,
######### it returns the number of aligning positions that could potentially display the 
######### transition of interest (e.g. for T>C transitions: How many T>C transitions are present
######### in the read sequence among the aligning nucleotides, and how many Ts are present in 
######### the reference sequence among the aligning nucleotides?).

######### Arguments: 
######### As arguments, it requires the nucleotide that is converted (e.g. "T" for T>C transitions),
######### the result of the transition ("C" for T>C transitions), the read sequence, the CIGAR string and the MD tag. 

######### Background:
######### The CIGAR string provides information about the position of aligning bases (M), but not about the 
######### nature of the alignment (match or mismatch). This information is contained in the MD tag.
######### The MD tag has to be read from the perspective of the the reference. 
######### Nucleotides that are missing from the reference are not considered 
######### (hard-clipped, soft-clipped, padded or inserted nucleotides, 
######### represented by H, S, P and I in the CIGAR string). 
######### Deletions with respect to the reference(D in the CIGAR string) are contained in the MD tag (as \^[A-Z]+),
######### while larger omissions (such as introns, N in the CIGAR string) are not part of the MD tag.
######### As a consequence, only nucleotides that align (M in the CIGAR string) are present in the read sequence,
######### the CIGAR string and the MD tag. 

######### Strategy:
######### Step 1: Positions that do not align to the reference (insertions, soft-clipped bases), 
######### are removed from the read sequence.
######### Step 2: The positions of nucleotides that can potentially transition or result from a 
######### transition (e.g. Ts and Cs for T>C transitions) within the aligning nucleotides from Step 1 
######### are determined.
######### Step 3: If there are no transitions in the MD tag, 
######### zero transitions and the number of aligning positions that could potentially undergo the 
######### transition of interest within the reference are returned as result.
######### Step 4: If there are transitions, positions that correspond to the converted or unconverted 
######### nucleotide are checked for (mis-)matches and counted.


sub count_transitions {
	my $from = $_[0]; # nucleotide that is converted ("T" for T>C transitions)
	my $to = $_[1]; # result of conversion ("C" for T>C transition)
	my $read = $_[2]; # Read sequence 
	my $CIGAR = $_[3]; # CIGAR string 
	my $MD_tag = $_[4]; # MD tag 
	
	my $trans_number;
	my $ref_number;


	# Step 1: 
	
	$CIGAR =~ s/[0-9]+[HPND]//gi; # remove hard-clipped, padded, omitted and deleted bases from CIGAR
	my @CIGAR_numbers = split(/[MSI]/, $CIGAR); # split CIGAR along characters to extract numbers
	my @CIGAR_chars = split(/[0-9]+/, $CIGAR); # split CIGAR along numbers to extract characters
	my @matching_stretch;
	my $read_pos = 0;
	my $match_pos = 0; # position among matching bases (for comparison to MD tag)

	for(my $i = 0; $i < scalar(@CIGAR_numbers); $i++)
		{
		if($CIGAR_chars[$i + 1] eq "M") # If positions are matches, extract nucleotides from read
			{
			my $new_bases = substr($read, $read_pos, $CIGAR_numbers[$i]); # extract bases of matching stretch from read sequence
			push(@matching_stretch, $new_bases); 
			$read_pos = $read_pos + $CIGAR_numbers[$i];
			}else{ # If positions are insertions or soft-clipped bases, only add number of bases to positions on read
			$read_pos = $read_pos + $CIGAR_numbers[$i];
			}
		}

	my $matching_bases = join("", @matching_stretch);
	
	# Step 2: 
	
	my @matching_bases = split(//, $matching_bases);
	my @from_pos; # position of bases that correspond to unconverted base
	my @to_pos; # position of bases that correspond to potentially converted base

	for(my $i = 0; $i < scalar(@matching_bases); $i++)
		{
		if($matching_bases[$i] eq $from)
			{
			push(@from_pos, $i + 1);
			}
		if($matching_bases[$i] eq $to)
			{
			push(@to_pos, $i + 1);
			}
		}
	
	# Step 3: 
		
	if($MD_tag !~ /[ATGC]/)
		{
		$ref_number = scalar(@from_pos);
		return 0, $ref_number;
		}
		
	# Step 4:
	
	if($MD_tag =~ /[ATGC]/)
		{
		my $MD = substr($MD_tag, 5); # removes the MD:Z: part of the tag
		$MD =~ s/\^[ACTG]+/D/gi; # replaces deletions given in MD tag by "D"

		$trans_number = 0; # number of matching bases that correspond to converted base
		$ref_number = scalar(@from_pos); # number of matching bases that correspond to unconverted base
	
		my @MD_array = split(//, $MD);
		my $pos = 0;
		my @to_add;

		foreach(@MD_array)
			{
			if($_ =~ "D") # If there is a deletion, ignore it.
				{
				if(scalar(@to_add) > 0)
					{
					my $to_add = join("", @to_add);
					$pos = $pos + $to_add;
					@to_add = (); # empty array of digits again
					}
				next;
				}
			if($_ =~ /[ATCG]/)
				{
				if(scalar(@to_add) > 0)
					{
					my $to_add = join("", @to_add);
					$pos = $pos + $to_add;
					@to_add = (); # empty array of digits again
					}
				$pos++;
			
				foreach my $from_pos (@from_pos)
					{
					if($pos == $from_pos) # If an unconverted base in the read is a mismatch...
						{
						$ref_number--; #... subtract 1 from the number of unconverted bases.
						}
					}

				if($_ eq $from) # If the mismatch corresponds to the unconverted base in the reference...
					{
					foreach my $to_pos (@to_pos) 
					{
					if($pos == $to_pos) # ... and the position corresponds to a converted base in the read ...
						{
						$trans_number++; # ... add 1 to the number of transitions.
						}
					}	
				
				}			 
			next;
			}
			if($_ =~ /[0-9]/)
				{
				push(@to_add, $_);
				}
			}

		if(scalar(@to_add) > 0)
			{
			my $to_add = join("", @to_add);
			$pos = $pos + $to_add;
			}
	
		my $T_total = $trans_number + $ref_number; # total number of Ts is matching Ts plus T>C transitions 
		
		return $trans_number, $T_total;

		}
	
}
