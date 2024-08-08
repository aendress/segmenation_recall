#!/usr/bin/perl

use strict;
use warnings;

use File::Slurper qw (read_text);
use Digest::SHA qw (sha1_hex);
use File::Basename;

=pod

=head1 NAME



=head1 SYNOPSIS

anonymize_testable [path]

=head1 DESCRIPTION
 
This script reads in all files in a directory, replaces testable minds IDs with SHA1 hashes, and renames the files. 

It also prints a file ids.txt to recover the correspondance between ids and hashes
 
=cut


# Find current directory
my $scriptDir = dirname(__FILE__);

# Read arguments
my $dir = (@ARGV < 1) ? '.' : $ARGV[0];
$dir = "$scriptDir/$dir";

# Read files
opendir (DIR, $dir);
my @files = grep(/\_M\d+\.csv$/,readdir(DIR));
closedir (DIR);

# Create hashes 
my %ids = map {
  my $k = $_;
  $k =~ s/^.*\_(M\d+)\.csv$/$1/;
  $k => sha1_hex($k);
} @files;

# Modify and rename files
foreach my $oldFile (@files){

  my $content = read_text("$dir/$oldFile")
    or die "Cannot read file: $!";

  $content =~ /\,\"?(M\d+)\"?\,/;

  my $oldId = $1;
  my $newId = $ids{$oldId};

  $content =~ s/$oldId/$newId/;
  my $newFile = $oldFile;
  $newFile =~ s/$oldId/$newId/;    

  open (OUT, ">$dir/$newFile")
    or die "Cannnot open file: $!";
  print OUT $content;
  close (OUT);  

}


# Print correspondance table
open (OUT, ">$dir/ids.csv")
  or die "Cannot open file: $!";


print OUT "id,hash";

foreach my $k (keys %ids){
  
  print  OUT "$k,$ids{$k}";
  
}

close (OUT);


