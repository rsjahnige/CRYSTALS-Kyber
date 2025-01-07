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

my $numArgs = $#ARGV + 1;
my $fileName = "";

if ($numArgs == 1)
{
	if ($ARGV[0] =~ m/(\.txt)$/) 
	{
		$fileName = $ARGV[0];
	}
	else
	{
		$fileName = "ERROR";
	}
}
else
{
	print "ERROR :: No command line arguments provided\n";
	print "Usage: perl sha_ex_psr.pl [flags] fileName.txt\n";
	print "\t-t	: hex_output.txt is populated with Round #0 - After Theta hex data\n";
	exit;
}

if ($fileName eq "ERROR")
{
	print "ERROR :: No input file path found :: File format must be '.txt'\n";
	exit;
}

open(FH_TEST, '<', $fileName) or die "ERROR :: $fileName not found\n";
#open(FH_MSG, '>', 'bin_msg.txt') or die "ERROR :: could not open bin_msg.txt\n";
open(FH_IN, '>', 'hex_input.txt') or die "ERROR :: could not open hex_input.txt\n";
open(FH_OUT, '>', 'hex_output.txt') or die "ERROR :: could not open hex_output.txt\n";


sub print_output 
{
	my $outData = <FH_TEST>;
	$outData =~ s/\s//g;
	print FH_OUT "$outData";
}

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

=if ($_ eq "After Theta\n") 
	{
		print_output();
	}
=cut

	if ($_ eq "After Rho\n")
	{
		print_output();
	}
}

close(FH_OUT);
close(FH_IN);
#close(FH_MSG);
close(FH_TEST);

exit;
