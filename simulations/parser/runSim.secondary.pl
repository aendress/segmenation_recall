#!/usr/bin/perl -wl


$, = "\t";

my ($file, $nSubj, $forgetting, $interference) = @ARGV;

my @output;
my $streamId;

my $streamFile;
my $fileBase = $file;
$fileBase =~ s/^.*\///;



open (OUT, ">output/$fileBase.results.txt")
    or die "Cannot open file output/$file.results.txt: $!";

foreach my $s (1..$nSubj){

  # Parameters:
  # forgetting (original rate: .05)
  # interference (original rate: .005)
  # increment on match
  # threshold

  # Randomly choose a familiarization stream
  $streamId = 1 + int rand(40); 
  $streamFile = "$file.$streamId.txt";
  
  @output = `./parser.cleaned.fixed.initial.window.pl $streamFile $forgetting $interference .5 1`;
  chomp @output;
  @output = grep (/^\S{6}\s/, @output);
  
  foreach (@output){
    
    print OUT $fileBase, $s, $_;
    
  }
}


