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
	print "Usage: perl sha_ex_psr.pl [fileName].txt\n";
	exit;
}

if ($fileName eq "ERROR")
{
	print "ERROR :: No input file path found :: File format must be '.txt'\n";
	exit;
}

open(FH_TEST, '<', $fileName) or die "ERROR :: $fileName not found\n";
open(FH_MSG, '>', 'bin_msg.txt') or die "ERROR :: could not open bin_msg.txt\n";
open(FH_OUT, '>', 'exp_output.txt') or die "ERROR :: could not open exp_output.txt\n";

while (<FH_TEST>) 
{

	if ($_ eq "Msg as bit string\n") 
	{
		my $msg = <FH_TEST>;			# Read next line in file
		while (<FH_TEST>)			# Read in all lines of input data 
		{
			if ($_ eq "\n") 		# Break from the loop when an empty line is reached
			{
				last;
			}
			else
			{
				$msg = $msg.$_;
			}
		}	
		$msg =~ s/\s//g;			# Remove white space
		if ($msg =~ m/[01]+/)			# If binary string is provided
		{
			print FH_MSG "$msg";		# Print binary string to bin_msg.txt
		}
		else
		{
			print FH_MSG "";		# Print empty string to bin_msg.txt
		}
	}	

	if (($_ eq "Hash val is\n") or ($_ eq "Output val is\n"))
	{
		my $hash = <FH_TEST>;			# Read next line in file
		while (<FH_TEST>) 
		{
			$hash = $hash.$_;		# Concatenate lines until end of file
		}
		$hash =~ s/\s//g;			# Remove white space
		print FH_OUT "$hash"			# Print hash value to exp_output.txt
	}
}

close(FH_OUT);
close(FH_MSG);
close(FH_TEST);

exit;
