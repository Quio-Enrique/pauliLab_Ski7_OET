#!/bin/perl -w

use strict;

my $abundance = $ARGV[0];


#Script to sum read_counts and tpm from a sorted by transcript kallisto abundance output

#First open file and declare variables
open(ABUNDANCE, "$abundance");

my $i = 0;
my (@transcript, @length, @eff_length, @counts, @tpm) = ("", "", "", "", "");

while(<ABUNDANCE>){
	chomp;
	my @line = split(/\t/);
	my @isoform = split /_/, $line[0];
	($transcript[$i], $length[$i], $eff_length[$i], $counts[$i], $tpm[$i]) = ($isoform[0], $line[1], $line[2], $line[3], $line[4]);
	if($i != 0){
		if($transcript[$i] eq $transcript[$i-1]){
			if($length[$i] > $length[$i-1]){
				$length[$i-1] = $length[$i];
			}
			if($eff_length[$i] > $eff_length[$i-1]){
				$eff_length[$i-1] = $eff_length[$i];
			}
			$counts[$i-1] += $counts[$i];
			$tpm[$i-1] += $tpm[$i];
		}else{
			$i++;
		}
	}else{
		$i++;
	}
}

for($i=0; $i<scalar @transcript; $i++){
	print "$transcript[$i]\t$length[$i]\t$eff_length[$i]\t$counts[$i]\t$tpm[$i]\n";
}

close(ABUNDANCE);
