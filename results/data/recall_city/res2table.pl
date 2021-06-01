#!/Users/endress/bin/perl -0015 -l012

use strict;
use warnings;

# Set working directory and output file
my $wd = (@ARGV) ? $ARGV[0] : ".";
$wd =~ s/\/$//;

my $cond_str = $wd;
$cond_str =~ s/^res\.//;


my $out_file = (@ARGV > 1) ? $ARGV[1] : $cond_str . ".tab";

my $print_header = (-e "$out_file") ? 0 : 1;

# Get results
opendir(DIR, "$wd")
    or die;
my @resFiles = grep(/\.res$/, readdir(DIR));
closedir(DIR);

open (OUT, ">>$out_file")
    or die "Cannot open $wd/$out_file: $!";


# Get header from first file in directory
if ($print_header) {
  my @header = qw/condDir subj gender age/;
  push (@header, get_header_row ("$wd/$resFiles[0]"));
  print  OUT join ("\t", @header);
}


 

map {analyze_file ("$wd/$_", $cond_str)} @resFiles;

close (OUT);


sub get_header_row{
  
  my $resFile = shift;

  open (IN, $resFile)
    or die "Cannot open file $resFile: $!";

  
  while (<IN>){

    if (/^Trial/){
      my @header = split (/\s/);
      
      return @header;
    }
    
  }

  return ();
}
sub analyze_file 
{  

  my ($in_file, $cond_str) = @_;
  
  # Read in current file
  open (IN, $in_file)
    or die "Cannot open $in_file: $!";
  
  my @lines = <IN>;  
  close (IN);
  
  chomp @lines;
  
  # Go through lines to find subject info; stop at trials
  my ($current_line, $subjnumber, $initials, $age, $gender, $subj, @out_line);
  
  while (($current_line = shift @lines) !~ /^Trial/){
    $subjnumber = $1
      if ($current_line =~ /^SubjectNumber\:\s+(\d+)/i);
    
    $initials = $1
      if ($current_line =~ /^Initials\:\s+(\S+)/i);

    $age = $1
      if ($current_line =~ /^age\:\s+(\d+)/i);

    $gender = $1
      if ($current_line =~ /^sex\:\s+(\S+)/i);
    
  }

  $subj = "$subjnumber.$initials";
  
  # now we arrive at the actual trials
  @lines = grep (/^\d/, @lines);
  print STDERR scalar @lines . " trials found in file $in_file."; 
  
  
  foreach $current_line (@lines){
    $current_line =~ s/[\[\]]//g;
    $current_line =~ s/\_/\t/g;
    

    @out_line = ($cond_str, $subj, $gender, $age);
    push (@out_line, split (/\t/, $current_line));
    print OUT join ("\t", @out_line)
      
  }
}






