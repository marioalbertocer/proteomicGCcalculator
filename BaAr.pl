#usr/bin/perl;

use File::Basename; # Alway read the paths since Users

opendir (mono_euks, "/Users/katzlab32/Dropbox/for David-Miguel/mono_euk-euk_only"); # open the folder in dropbox
@files = readdir mono_euks;
open (REP, ">>/Users/katzlab32/Desktop/report.txt") or die "can't open report\n\n"; #open an empty file for writing the report


@files = grep (/\.fas$/, @files); # filtering temporary files


my $count = 0;  #a counting... in order to know how many files have been processed
foreach $file (@files) { #looping in the list of files that compose the folder in dropbox
    
    $count++; 
    $OG = $file; 
#    $OG =~ s/\_alignment.*$//; 
    open (F, "/Users/katzlab32/Dropbox/for David-Miguel/mono_euk-euk_only/$file\n"); #open each file
    @file = <F>; #putting the information of each file in an array
    $file = join ('', @file); # transforming the array in a line for easy handeling
    
    $BaAr = "-"; # As a default there is absence of bacteria or archaea until ...
    $BaAr = "+" if ($file =~ />Ba_/); # it finds a bacteria or ...
    $BaAr = "+" if ($file =~ />Za_/); # it finds an Archaea
    
    print "$count\t$OG\t$BaAr\n"; #print results
    print REP "$count\t$OG\t$BaAr\n";
}

