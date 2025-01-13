#!/bin/bash

declare -a inputFiles

inputFiles=(Test_Examples/SHA/*)

for file in "${inputFiles[@]}"
do
	perl sha_ex_psr.pl $file
	expResult=$(cat exp_output.txt)

	input=$(echo "$file" | cut -d '/' -f 3)
	input=$(echo "$input" | cut -d '_' -f 1)
	
	fnct=$(echo "$input" | cut -d '-' -f 1) 
	dlen=$(echo "$input" | cut -d '-' -f 2) 

	if [ "${fnct}" = "XOF" ] 
	then 
		./test05 "$fnct" "$dlen" 4096
	else
		./test05 "$fnct" "$dlen"
	fi

	actResult=$(cat act_output.txt)

	if [ "${expResult,,}" = "${actResult,,}" ];
	then
		echo "Test successful - $file"
	else
		echo "Test failed - $file"
		exit
	fi
	
	rm -f *.txt
done

exit
