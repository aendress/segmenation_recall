#!/usr/bin/perl -wl

use strict;

$, = "\t";

my ($result_file, $test_list) = @ARGV;

my ($result_file_base, $test_list_base) = map {s/^.*\///r} ($result_file, $test_list);
my $outputFile = "$test_list_base.$result_file_base.results.txt";

#my @result_files = qw/EM1lang8.syll EM1lang8.syll.faded.2syll/;
#my @test_lists = qw/testList.forParser testList.PhWPW.forParser/;

my (%lexicon, $nSubj, @trials);

# Hash of hash references:
# stream -> subject -> item -> weight
%lexicon = read_results ($result_file . ".results.txt");

# Array of arrays. Outer array reflects lines, inner array columns
@trials = read_test_list ($test_list . ".txt");


$nSubj = (sort {$a <=> $b} (keys %{$lexicon{$result_file_base}}))[-1];




open (OUT, ">output/$outputFile")
    or die "Cannot open file output/$outputFile: $!";

foreach my $subj (1..$nSubj){

  my @current_results = analyze_subject ($subj, $result_file_base, $test_list_base, 
			       \%lexicon, \@trials);
  
    foreach (@current_results){
      	print OUT join ("\t", @$_);      
    }
}

close (OUT);



sub read_results 
{
  my ($file) = @_;

  my @line;
  my %lexicon = ();
  
  open (FD, $file)
    or die "Cannot open file $file: $!";
  while (<FD>){
    chomp;
    
    @line = split (/\t/);

    $lexicon{$line[0]} = {}
      unless (exists ($lexicon{$line[0]}));
    
    $lexicon{$line[0]}->{$line[1]} = {}
      unless (exists ($lexicon{$line[0]}->{$line[1]}));

    # stream -> subject -> item -> weight
    $lexicon{$line[0]}->{$line[1]}->{$line[2]} = $line[3];
    
  }
  close (FD);

  return %lexicon;
}

sub read_test_list
{
  my ($file) = @_;

  my @trials = ();
  open (FD, "$file")
    or die "Cannot open file $file: $!";

  while (<FD>){
    chomp;
    push (@trials, [split (/\t/)]);
  }

  close (FD);

  return @trials;

}

sub analyze_subject 
{
  my ($subj, 
      $result_file, 
      $test_list, 
      $lexicon, 
      $trials) = @_;

  die "There is no lexicon for result file $result_file, exiting."
    unless (exists ($lexicon{$result_file}));

  # Each entry is an array reference including the following information
  # * stream name
  # * subject number
  # * trial info from as line in @trial (i.e, the two items, the test type and the test sub type
  # * accuracy variable indicating whether the answer on the trial was correct.
  my @results = ();

  # weights of the two items in a test trial
  my (@scores);

  
  if (exists ($lexicon{$result_file}->{$subj})) {
    # There are items with 6 phonemes for $subj. We thus can determine their weights.
    
    # This hash has the structure
    # item -> weight
    # We will test these item weights for each trial in the test list below
    my %sub_lexicon = %{$lexicon{$result_file}->{$subj}};
    
    
    foreach my $t (@trials){
      push (@results, [$result_file, $subj]);
      push (@{$results[-1]}, @$t);
      
      foreach my $i (0..1){
	$scores[$i] = (exists ($sub_lexicon{$t->[$i]})) ? $sub_lexicon{$t->[$i]} : 0;
	
      }

      # Determine weight of the two test items
      if ($scores[0] > $scores[1]){
	push (@{$results[-1]}, 1);
      } elsif ($scores[0] < $scores[1]){
	push (@{$results[-1]}, 0);
      } else {
	push (@{$results[-1]}, .5);      
      }    
    }
    
  } else {
    # There are NO items with 6 phonemes for $subj. We thus can determine their weights.
    # We will set all comparisons to .5 as the model has to guess

    foreach my $t (@trials){
      push (@results, [$result_file, $subj]);
      push (@{$results[-1]}, @$t);
      # Set accuracy to .5
      push (@{$results[-1]}, .5);            
    }
    

  }
  
  return @results;
}
