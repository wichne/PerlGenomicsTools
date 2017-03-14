#!/usr/bin/perl

$decision = "Y";
while (my $l = <STDIN>) {
    if ($l =~ /^>(\S+)/) {
	$id = $1;
    } elsif ($l =~ /max. C\s+(\d+)\s+([\d\.]+).*(YES|NO)/) {
	$Cpos = $1;
	$Cmax = $2;
	$decision = $3 eq "NO" ? "N" : $decision;
    } elsif ($l =~ /max. Y\s+(\d+)\s+([\d\.]+).*(YES|NO)/) {
	$Ypos = $1;
	$Ymax = $2;
	$decision = $3 eq "NO" ? "N" : $decision;
    } elsif ($l =~ /max. S\s+(\d+)\s+([\d\.]+).*(YES|NO)/) {
	$Spos = $1;
	$Smax = $2;
	$decision = $3 eq "NO" ? "N" : $decision;
    } elsif ($l =~ /mean S\s+\d+\-\d+\s+([\d\.]+).*(YES|NO)/) {
	$Smean = $1;
	$decision = $2 eq "NO" ? "N" : $decision;
    } elsif ($l =~ /max. D\s+\d+\-\d+\s+([\d\.]+).*(YES|NO)/) {
	$Dmax = $1;
	$decision = $2 eq "NO" ? "N" : $decision;
    } elsif ($l =~ /^\/\//) {
	print "$id\t$Cmax\t$Cpos\t$Ymax\t$Ypos\t$Smax\t$Spos\t$Smean\t$Dmax\t$decision\n";
	$id=$Cmax=$Cpos=$Ymax=$Ypos=$Smax=$Spos=$Smean=$Dmax="";
	$decision = "Y";
    }
}
    
