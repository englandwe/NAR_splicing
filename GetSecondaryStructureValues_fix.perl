#!/usr/bin/perl

######################
# (c) Michael Hiller #
######################

# call: perl -w GetSecondaryStructureValues.perl -f inputfile [-o outputfile -l motiflen -flankUp sequence -flankDown sequence -method [EF,PU] ]\n";
# inputfile has to have the following format:
#		>name
#		sequence
#		>name2
#		sequence2
#	   ...
# if outputfile is not specified, the output is written to 'inputputfile + .sec'
# motiflen gives the length of the motif (you have to use this parameter later for MEME -w)
# flankUp specifies a sequence as a constant upstream flank for all sequences in inputfile that should be included for the folding but not for the motif search 
#   (this can be for example a constant primer annealing site); 
# flankDown specifies a sequence as a constant downstream flank
# method gives the probabilites to be computed; either the expected fraction of a motif that is paired (EF) or the probability that the complete motif is unpaired (PU)
#
# NOTE: all sequences (inputfile, flankUp, flankDown) are converted to upper case and U is replaced by T
#
# NOTE: you may have to edit the $RNAfold variable below


use strict;
use warnings;
use Cwd;

my $RNAfold = "/data/apps/ViennaRNA/2.2.8/bin";			# directory of RNAfold without '/' at the end
my $KT = 0.616320775;											# used for probability calculation (constant value)
my @matrix;															# used to store the dot plot matrix

#my $globalcount = 0;

# compute pair-prob matrix with RNAfold and read the matrix
sub getProbMatrix {
	# read parameter
	my $seq = $_[0];
	
	# write fasta file and call RNAfold
	my $name = "./~tmpSeq.fa";		
	open(file_local, ">$name");
	print file_local ">$name\n";
	print file_local "$seq\n";
	close(file_local);	
	
	my $name1 = "./~tmpSeq.out";		
	# call RNAfold with input and output filename
	system "$RNAfold/RNAfold -p < $name > $name1";
	
	# create one line vector for the matrix
	my (@zeile);
	my $len = length($seq);
	for (my $j=0; $j<$len; $j++) {
		push(@zeile,0);
	}

	# read dot plot
	$name = "./~tmpSeq.fa_dp.ps";		
	open(file_local, "$name");
	while (my $line1 = <file_local>) {

		# search for sequence line
		if ($line1 =~/\/sequence \{ \(\\/) {
			# read next line --> sequence
			$line1 = <file_local>;
			chomp($line1);
			my $seq1;
			while (!($line1 =~ /def/)){
				$line1 =~/([ACGU]+)\\/;		# remove \ at line end
				$seq1 .= $1;
				$line1 = <file_local>;
				chomp($line1);
			}

			my $seqlen = length($seq1);
			$seq1 =~ tr/U/T/;		# our input seq has T instead of U
			# security check
			if ((! ($seq eq $seq1)) || (! ($len == $seqlen))) {
				print "ERROR: sequence in ps file does not match input sequence\n$seq\n$seq1\n";
			}

			# create matrix
			for (my $i=0; $i<$seqlen; $i++) {
				$matrix[$i] = [@zeile];		# push reference of a copy-array to matrix
			}

		# search for lines like 1 8 0.00971 ubox
		}elsif ($line1 =~/(\d+) (\d+) ([\d\.]+) ubox\n/) {
			my $value = $3*$3;						# need to square the values	
			$matrix[$1-1][$2-1] = $value;		# indices begin with 1 in dot.ps
		}
	}
	close(file_local);

	# for easier evaluation copy upper left part to lower left part
	for (my $i=0; $i<$len; $i++) {
		for (my $j=$i; $j<$len; $j++) {
			$matrix[$j][$i] = $matrix[$i][$j];	
		}	
	}		

}

# compute pi (prob that position i pairs) and compute the expected fraction of the given region that is paired
sub getExpectedFraction {
	# read parameter
	my ($start,$end,$seqlen) = ($_[0],$_[1],$_[2]);

	# we sum up lines (does'nt matter since the matrix is symmetrical)
	my $EF = 0;
	for (my $i=$start; $i<=$end; $i++) {
		my $sum = 0;
		for (my $j=0; $j<$seqlen; $j++) {
			$sum += $matrix[$i][$j];
		}
		$EF += $sum;
	}
	return $EF /= ($end - $start + 1);
}


# compute energy of the complete ensembl
sub getEnergyAll {
	# read parameter
	my ($seq) = ($_[0]);
	my $len = length($seq);
	
	# compute energy of whole ensemble
	# write fasta file and call RNAfold
	my $name = "./~tmpSeq.fa";		
	open(file_local, ">$name");
	print file_local ">$name\n";
	print file_local "$seq\n";
	close(file_local);
	
	my $name1 = "./~tmpSeq.out";		
	system "$RNAfold/RNAfold -p0 < $name > $name1";
	
	# read energy of whole ensemble (all structures)
	my $energy_all = 10000000;		# just init with a big number
	open(file_local, "$name1");
	while (my $line1 = <file_local>) {
		# search for line with ensemble energy
		if ($line1 =~/ free energy of ensemble = +([-\d\.]+) kcal\/mol.*/) {
			$energy_all = $1;	
		}
	}
	close(file_local);

	#print "energy of all-ensemble: $energy_all\n";

	# sometimes RNAfold with constraints breaks --> catch those cases and return '-'
	if ($energy_all == 10000000) {
		my $return_value = "-";
		return $return_value;
	}
	
	return $energy_all;
}


# compute prob that the given region is unpaired (using RNAfold with constraints)
sub getUnpairedProb {
	# read parameter
	my ($seq,$motif_start,$motif_end,$energy_all) = ($_[0],$_[1],$_[2],$_[3]);
	my $len = length($seq);
	
	# compute energy of whole ensemble with unpair-constraints
	# write fasta file and call RNAfold
	my $name = "./~tmpSeq.fa";		
	open(file_local, ">$name");
	print file_local ">$name\n";
	print file_local "$seq\n";
	# create constraints line
	for (my $i=0; $i<$motif_start; $i++) {			# '.' means no constraint (for the positions upstream)
		#print ".";
		print file_local ".";
	}
	for (my $i=$motif_start; $i<=$motif_end; $i++) {			# 'x' means unpaired (for the positions in the motif)
		#print "x";
		print file_local "x";
	}
	for (my $i=$motif_end+1; $i<$len; $i++) {			# '.' means no constraint (for the positions downstream)
		#print ".";
		print file_local ".";
	}
	#print "\n";
	print file_local "\n";
	close(file_local);
	
	my $name1 = "./~tmpSeq.out";
	system "$RNAfold/RNAfold -p0 -C < $name > $name1";
	# read energy of ensembl with unpaired constraints
	my $energy_constraint = 10000000;		# just init with a big number
	open(file_local, "$name1");
	while (my $line1 = <file_local>) {
		# search for line with ensemble energy
		if ($line1 =~/ free energy of ensemble = +([-\d\.]+) kcal\/mol.*/) {
			$energy_constraint = $1;	
		}
	}
	close(file_local);

	print "energy of constraint-ensemble: $energy_constraint\n";

	# sometimes RNAfold with constraints breaks --> catch those cases and return '-'
	if ($energy_constraint == 10000000) {
        #$globalcount+=1;
        #my $outname = $name1 . $globalcount . '_bad';		
     	#system "cp $name1 $outname";
		my $return_value = "-";
		return $return_value;
	}
	
	# compute prob of set of all structures with unpaired constraints (PU = Z_unpaired / Z where Z is the partition function)
	my $prob_unpaired = exp( ((-1*$energy_constraint) - (-1*$energy_all)) / $KT);
    #$globalcount+=1;
    #my $outname = $name1 . $globalcount . '_good';		
   	#system "cp $name1 $outname";

	return $prob_unpaired;
}




###################################################################################################
# 														main
###################################################################################################

my $inputfile = "";			# where to read the sequence from
my $outputfile = "";			# where to write results to 
my $method = "EF";			# default method to be computed is EF
my $flankUpstream = "";		# upstream flank (default none)
my $flankDownstream = "";	# downstream flank (default none)
my $motiflen = 6; 			# default value is motif length 6

if (@ARGV == 0) {
	print "usage: perl -w GetSecondaryStructureValues.perl -f inputfile [-o outputfile -l motiflen -flankUp sequence -flankDown sequence -method [EF,PU] ]\n";
	print "this perl file uses RNAfold to compute either the expected fraction of base pairing for a motif or the probability that the complete motif is unpaired\n";
	print "\t-l          gives the length of the motif (default is 6)\n";
	print "\t-method     gives the method to be used (either expected fraction = 'EF' or probability unpaired = 'PU', default EF)\n";
	print "\t-flankUp    gives the sequence of the upstream flank\n";
	print "\t-flankDown  gives the sequence of the downstream flank\n";
	exit(0);
}


# determine input parameters
for (my $i = 0; $i < @ARGV; $i++){
	if ($ARGV[$i] eq "-f") {
		if (defined($ARGV[$i+1])){
			$inputfile = $ARGV[$i+1];
			$i++;
		}
	}
	elsif ($ARGV[$i] eq "-o") {
		if (defined($ARGV[$i+1])){
			$outputfile = $ARGV[$i+1];
			$i++;
		}
	}
	elsif ($ARGV[$i] eq "-flankUp") {
		if (defined($ARGV[$i+1])){
			$flankUpstream = uc $ARGV[$i+1];		# upper case
			$i++;
			$flankUpstream =~ tr/U/T/;				# T instead of U
		}
	}
	elsif ($ARGV[$i] eq "-flankDown") {
		if (defined($ARGV[$i+1])){
			$flankDownstream = uc $ARGV[$i+1];	# upper case
			$i++;
			$flankDownstream =~ tr/U/T/;			# T instead of U
		}
	}
	elsif ($ARGV[$i] eq "-l") {
		if (defined($ARGV[$i+1])){
			$motiflen = $ARGV[$i+1];
			$i++;
		}
	}
	elsif ($ARGV[$i] eq "-method") {
		if (defined($ARGV[$i+1])){
			if ($ARGV[$i+1] eq "EF") {
				$method = "EF";
			}elsif ($ARGV[$i+1] eq "PU") {
				$method = "PU";
			}else{	
				print "ERROR: expect 'EF' or 'PU' after -method\n";
				exit(-1);
			}
			$i++;
		}else{
			print "ERROR: expect 'EF' or 'PU' after -method\n";
			exit(-1);
		}
	}else{
		print "ERROR: unknown parameter: $ARGV[$i]\n";
		exit(-1);
	}
}
if ($inputfile eq ""){
	print "You have to specify an input file: -f file\n";
	exit(-1);
}
if ($outputfile eq ""){
	$outputfile = $inputfile . ".sec";
}

# test if RNAfold binary exist
my $RNAfoldBin = "$RNAfold/RNAfold";
$RNAfoldBin =~ s/~/$ENV{HOME}/;		# replace ~ 
if (! (-e $RNAfoldBin)) {
	print "ERROR: RNAfold binary does not exist ($RNAfoldBin)\n";
	exit(-1);
}


my $flankLen = length($flankUpstream);
my $energy_all = 0;					# free energy of complete ensembl

# print parameter
print "reading:          $inputfile\n";
print "write results to: $outputfile\n";
print "upstream flank:   \"$flankUpstream\"\n";
print "downstream flank: \"$flankDownstream\"\n";
print "motif length:     $motiflen\n";
print "method:           $method\n\n";


#####################################################################
# open input and output file
open(INPUT, $inputfile);
open(OUTPUT, ">$outputfile");

# read all sequences
while (my $line1 = <INPUT>) {
	chomp($line1); 					# fasta header
	my $line2 = <INPUT>; 			# sequence
	chomp($line2); 

	my $seq = uc $line2;				# sequence
	$seq =~ tr/U/T/;					# T instead of U
	my $len = length($seq);
	$line1 =~ />(.*)/;				# get sequence name
	my $seqName = $1;


	print OUTPUT "$seqName;";

	my $completeSeq = $flankUpstream . $seq . $flankDownstream;		# sequence with the up/downstream flank
	print ">$seqName: $completeSeq\n";

	# in case of EF compute matrix only once
	if ($method eq "EF") {
		# compute matrix with pair probabilities
		getProbMatrix($completeSeq);
	}else {		# in case of PU compute energy of the complete ensembl only once
		$energy_all = getEnergyAll($completeSeq);
	}

	# compute value for all motif starts (0 ... seqlen - motiflen + 1)
	for (my $motif_start = 0; $motif_start <= length($seq)-$motiflen; $motif_start++) {
		# start and end of motif in cur_seq
		my $motif_end = $motif_start + $motiflen - 1;

		my $ExpectedFraction = 0;
		my $ProbUnpaired = 0;
		if ($method eq "EF") {					# compute EF
			$ExpectedFraction = getExpectedFraction($motif_start+$flankLen, $motif_end+$flankLen, length($completeSeq));
		}elsif ($method eq "PU") {				# compute prob unpaired  
			$ProbUnpaired = getUnpairedProb($completeSeq, $motif_start+$flankLen, $motif_end+$flankLen, $energy_all);
		}

		# print into outfile
		if ($method eq "EF") {
			printf OUTPUT "%1.4f;", 1-$ExpectedFraction;					# $ExpectedFraction = 1 means complete motif is paired --> for prior probability use 1 - $ExpectedFraction
			printf "motif_start: %03d  motif_end: %03d  EF: %1.4f\n",$motif_start+1,$motif_end+1,1-$ExpectedFraction;
		}elsif ($method eq "PU") {
            if ($ProbUnpaired == "-") {
    			printf OUTPUT "%s;", $ProbUnpaired;
    			printf "motif_start: %03d  motif_end: %03d  PU: %s\n",$motif_start+1,$motif_end+1,$ProbUnpaired;
            }else {    
			    printf OUTPUT "%1.4f;", $ProbUnpaired;
    			printf "motif_start: %03d  motif_end: %03d  PU: %1.4f\n",$motif_start+1,$motif_end+1,$ProbUnpaired;
            }
		}
	}
	print OUTPUT "\n";

}	# all sequences in input file

close(INPUT);	
close(OUTPUT);

# delete temporary files
system "rm -f ./~tmpSeq.fa ./~tmpSeq.out ./~tmpSeq.fa_dp.ps ./~tmpSeq.fa_ss.ps";



