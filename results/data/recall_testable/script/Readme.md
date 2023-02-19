# anonymize_testable

## Description 
This script reads in all files in a directory, replaces testable minds IDs with SHA1 hashes, and renames the files. 

It also prints a file ids.txt to recover the correspondance between ids and hashes.

## Usage 
### All systems
anonymize_testable.pl [path]

### Mac OS X
Change the extension to .command instead of .pl, place the file in the directory you want to anonymize and double click on it after installation of perl and its modules (see below).

### Windows 
Place file in the directory you want to anonymize and double click on it after installation of perl and its modules (see below).



## Installation 
### Mac OS X/Linux 
Use the built-in perl version or install ActivePerl from https://www.activestate.com/products/perl/downloads/

Install the following modules
* File::Slurper
* Digest::SHA
* File::Basename

If you use the build-in version of perl, follow the instructions here: 
http://triopter.com/archive/how-to-install-perl-modules-on-mac-os-x-in-4-easy-steps/


### Windows
Install ActivePerl from https://www.activestate.com/products/perl/downloads/

Install the following modules
* File::Slurper
* Digest::SHA
* File::Basename













 