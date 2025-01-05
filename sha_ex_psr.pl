# File: sha_ex_psr.pl
# Author: Ryan Jahnige
#
# Description: Parse the SHA-3 examples contained in the
# 		Test_Examples directory. Example were obtained 
# 		from the NIST website and copied into a text 
# 		document (see FIPS-202:A-B).

#!/usr/bin/perl

use strict;
use warnings;

my $fileName = 'Test_Examples/SHA3-224_Msg0.txt';

open(FH_TEST, '<', $fileName) or die "ERROR :: $fileName not found\n";
#open(FH_MSG, '>', 'bin_msg.txt') or die "ERROR :: could not open bin_msg.txt\n";
open(FH_IN, '>', 'hex_input.txt') or die "ERROR :: could not open hex_input.txt\n";
open(FH_OUT, '>', 'hex_output.txt') or die "ERROR :: could not open hex_output.txt\n";

while (<FH_TEST>) 
{

=if ($_ eq "Msg as bit string\n") 
	{
		my $msg = <FH_TEST>;			# Read next line in file
		$msg =~ s/\s//g;			# Remove white space
		if ($msg =~ m/[01]+/)			# If binary string is provided
		{
			print FH_MSG "$msg";
		}
	}
=cut	

	if ($_ eq "Data to be absorbed\n") 
	{
		my $inData = <FH_TEST>;
		$inData =~ s/\s//g;
		print FH_IN "$inData";
	}	

	if ($_ eq "After Theta\n") 
	{
		my $outData = <FH_TEST>;
		$outData =~ s/\s//g;
		print FH_OUT "$outData";
	}
}

close(FH_OUT);
close(FH_IN);
#close(FH_MSG);
close(FH_TEST);
