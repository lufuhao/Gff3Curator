#!/bin/bash

if [ -z "$@" ]; then
	echo "Error: no scaffold id detected" >&2
	exit 1
fi

scaffoldid=$1
if [ -z "$scaffoldid" ]; then
	echo "Error: no scaffold id detected" >&2
	exit 1
fi

wheat_annotation_extractor20150914.sh $scaffoldid

fstart=0
fend=0
test=0


while [ $test -eq 0 ]; do
	echo "Please input feature start: "
	read fstart
	echo "Detected: $fstart"
	if [[ "$fstart" =~ ^[0-9]+$ ]] && [ $fstart -gt 0 ]; then
		echo "Accepted: $fstart"
		test=1
	fi
done

test=0
while [ $test -eq 0 ]; do
	echo "Please input feature end: "
	read fend
	echo "Detected: $fend"
	if [[ "$fend" =~ ^[0-9]+$ ]] && [ $fend -gt $fstart ] ; then
		echo "Accepted: $fend"
		test=1
	fi
done
echo "gff_sum_exons_for_annotation20150914.pl $scaffoldid:$fstart-$fend"
gff_sum_exons_for_annotation20150914.pl $scaffoldid:$fstart-$fend
