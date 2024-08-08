#!/usr/bin/perl -w

# PARSER: A Model for Word Segmentation
#
# proposed in:
# Perruchet, P. & Vinter, A. (1998). PARSER: A Model for Word Segmentation.
# Journal of Memory and Language, 39, 246-263.
#
# Implemented by Ansgar D. Endress (ansgar.endress@polytechnique.org)

# Since P&V paper was unclear in certain respects, there are several incertainties in our simulations:
# 
# - How can "primitives" (i.e., syllables) be prevented from disappearing from the lexicon, 
# which is an requirement to parse a text? According to P&V, primitives are initially assigned
# a weight of 1, and a then subject to forgetting. However, since units can only be activated if 
# their weight is at least 1, words with syllables that have not been parsed at the first time 
# step could not possibly be recognized. Moreover, P&V report primitives with unchanged weights 
# (p.253). We presume therefore that P&V forgot to specify that the weight of primitives couldn't 
# drop under 1.0.

#############################################################################################################
# VERSION JAN 12th, 2023										    #
# This version forces the model to read in a fixed number of units in the initial iterations.   	    #
# We force the network to read a total of 3 chunks on the first iteration to bias it to consider part-words #

our @N_UNITS_INITIAL = (3);

#############################################################################################################

######################################################
#                                                    #
#                  GLOBAL AND CONSTANTS              #
#                                                    #
######################################################

$INCREMENT_ON_MATCH = 0.5;
$DECREMENT_THROUGH_FORGETTING = 0.05;
$DECREMENT_THROUGH_INTERFERENCE = 0.005;
$THRESHOLD = 1;
$INITIAL_WEIGHTS = 1;

# Variables that will be global
# @lexicon;
# @allSylls
# %weights;
# $nPrimitives;

sub randint
{
    # Return random number from 1 to $_[0]
    return int (rand ($_[0])) + 1;
}

sub min
{
    # Return minimum of the arguments
    if ($_[0] < $_[1]){
	return $_[0];
    }
    return $_[1];
}

sub initLexicon
{
    # Find all syllables in the stream, place them 
    # (the 'primitives') in the lexicon
    # and initialize their weights
    my @mystream = @{$_[0]};
    my %alreadySeen = ();
    foreach (@mystream) {
    	$_ =~ s/\s+//g;
    	$alreadySeen{$_}++;
    }
    @lexicon = keys %alreadySeen;
    @allSylls = keys %alreadySeen;
    
    $nPrimitives = @lexicon;
    foreach (@lexicon){
	$weights{$_} = $INITIAL_WEIGHTS;
    }
}

sub isWord
{
    # Test whether $_[0] has a lexical entry IRRESPECTIVELY of its weights
    if (exists $weights{$_[0]}){
	return 1;
    }
    return -1;
}

sub getCandidates
{
    # Return list of lexical entries that start with $_[0]
    my $currentChunk = $_[0];
    my @candidates = ();

    foreach (@lexicon){
	if (length ($_) >= length ($currentChunk)){
	    if ($_ =~ /^($currentChunk)/) {
		if ($weights{$_} >= $THRESHOLD){
		    push (@candidates, $_);
		}
	    }
	}
    }
    return @candidates;
}

sub verifyCandidates 
{
    # Check whether candidates agree only with the first part of the 
    # forthcoming chunks or whether they can really be extracted up to their end
    # For example, kipulikipuli would be reject for stream part like kipuliki#
    #
    # candidateCandidates are candidates to be lexical candidates
    # candidates are words that can really be extracted
    # 
    # The resulting list is sorted by decreasing length

    my @candidateCandidates = @{$_[0]}; 
    my $currentChunk = $_[1];
    my @stream = @{$_[2]};

    my ($l, $finalChunk, $cl, @stream2);
    my @candidates = ();

    foreach (@candidateCandidates) {
	$l = length ($_);
	if (length ($currentChunk) == $l){
	    if ($currentChunk eq $_){
	    	push (@candidates, $_);
	    }
	} else {
	    $finalChunk = $currentChunk;

	    @stream2 = @stream;
	    while ((@stream2) &&
		   (length ($finalChunk) < $l)){
		$cl = shift (@stream2);
		$finalChunk .= (split (/\s+/, $cl))[0];
	    }
	    if ($finalChunk eq $_){
	    	push (@candidates, $_);
	    }
	}
    }

    @candidates = sort {length($b) <=> length ($a)} @candidates;
    
    return @candidates;
}


sub strengthenUnitWeights
{
# Strengthen weights of the words in the argument array
    foreach (@{$_[0]}) {
	$weights{$_} += $INCREMENT_ON_MATCH;	    
    }
}

sub getSyllableList
{
# Return list of syllables that are present in the actual units 
    my @units = @{$_[0]};
    my ($as, $u, @sl);
        
   ALL_SYLLS_LOOP: foreach $as (@allSylls){    
	foreach $u (@units){
	    if ($u =~ /$as/){
	    	push (@sl, $as);
		next ALL_SYLLS_LOOP;
	    }
	}
    }        
    return @sl;
}

sub interference
{
    # First interfere with all chunks 
    # Then UNDO interference for current units
    my @units = @{$_[0]};
    my ($l, $syll, %units2, $v);
    my @sl = &getSyllableList (\@units);

    LEXICON_LOOP: foreach $l (@lexicon)  {
	foreach $syll (@sl){	
	    if ($l =~ /$syll/){
	    	$weights{$l} -= $DECREMENT_THROUGH_INTERFERENCE;
	    	next LEXICON_LOOP;
	    }
	}
    }
    #Undo decrements for units - but make sure to increase each entry only ONCE

    # Make table that contains each entry only once
    foreach (@units){
	$units2{$_} = $_;
    }
    foreach $v (values %units2) { #     foreach $i (keys %units2) would also be OK
	$weights{$v} += $DECREMENT_THROUGH_INTERFERENCE;
    }
}

sub forgetting
{
    foreach (@lexicon){
	$weights{$_} -= $DECREMENT_THROUGH_FORGETTING;
    }
}

sub addNewOrBelowThresholdEntry
{
    # add weight to the concatentation of the extracted units, arguments is @units
    my $tmpStr = join ("", @{$_[0]});

    if (exists $weights{$tmpStr}){
	$weights{$tmpStr} += $INCREMENT_ON_MATCH;
	$weights{$tmpStr} += $DECREMENT_THROUGH_INTERFERENCE; # undo interference that has already been applied
    } else {
	my $nWords = @lexicon;
	$lexicon[$nWords] = $tmpStr;
	$weights{$lexicon[$nWords]} = $INITIAL_WEIGHTS;
    }
}

sub keepPrimitives
{
    # Make sure that the weight of the initial syllables (= 'primitves')
    # doesn't fall below 1

    foreach (@allSylls) {
	if ($weights{$_} < 1){
	    $weights{$_} = 1;
	}
    }
}

sub clearLexicon
{
    # Remove lexical entries with weights <= 0
    my ($i);

    for ($i=(@lexicon-1); $i>($nPrimitives-1); $i--){ # Primitives cannot be removed
	if ($weights{$lexicon[$i]} <= 0){
	    delete $weights{$lexicon[$i]};
	    splice (@lexicon, $i, 1);
	}
    }
}

######################################################
#                                                    #
#                        MAIN                        #
#                                                    #
######################################################

my ($i, $currentLine, @stream, $nUnits, $currentChunk, $nWords, $complainStr, @candidates, @units);
srand;

if (@ARGV < 1){
    $complainStr = "Usage: parser.pl <file name> [INCREMENT_ON_MATCH DECREMENT_THROUGH_FORGETTING DECREMENT_THROUGH_INTERFERENCE THRESHOLD]\n\n";
    $complainStr .= "Note: The program assumes that syllables ";
    $complainStr .= "contain exactly two letters and that each syllable is in a separate line.\n\n";
    die $complainStr;
}

$complainStr = "File doesn't exist - exiting.\n";
open (FD, $ARGV[0]) or die "$complainStr: $ARGV[0]";
@stream = <FD>;
close (FD);
chop (@stream);
&initLexicon (\@stream);

if (@ARGV == 5){
    # Modify default values of constants
    $INCREMENT_ON_MATCH = $ARGV[3];
    $DECREMENT_THROUGH_FORGETTING = $ARGV[1];
    $DECREMENT_THROUGH_INTERFERENCE = $ARGV[2];
    $THRESHOLD = $ARGV[4];
}

######################################################
#                                                    #
#                   PROCESS LINES                    #
#                                                    #
######################################################

while (@stream){
  $nWords = @lexicon;

  # Changed Jan 12th, 2023 - Start
  #$nUnits = randint (3);
  $nUnits = (@N_UNITS_INITIAL) ? shift @N_UNITS_INITIAL : randint (3);
  # Changed Jan 12th, 2023 - End

    for ($i=0; (($i<$nUnits) && (@stream)); $i++){
	$currentLine = shift (@stream);
	($currentChunk) = split (/\s+/, $currentLine);
	@candidates = &getCandidates ($currentChunk);
	@candidates = &verifyCandidates (\@candidates, 
						 $currentChunk, \@stream);
	if (@candidates > 0){
	    $units[$i] = $candidates[0];
	    if ($candidates[0] ne $currentChunk){
		my $removeString = $currentChunk;
		while (length ($removeString) < length ($candidates[0])){
		    $removeString .= shift (@stream);
		}
	    } 
	}
    }
    $nUnits = @units;

    &strengthenUnitWeights(\@units);
    &interference (\@units); 
    &forgetting; 
    &keepPrimitives;
    if ($nUnits > 1){
	&addNewOrBelowThresholdEntry (\@units);
    }
    &clearLexicon;
    @units = ();
}

######################################################
#                                                    #
#                  PRINT RESULTS                     #
#                                                    #
######################################################

@lexicon = sort {$weights{$b} <=> $weights{$a}} @lexicon;

foreach (@lexicon) {
    if ($weights{$_}	> $THRESHOLD) {
	print $_ . "\t" . $weights{$_} . "\n";
    }
}







