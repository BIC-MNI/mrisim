#!/usr/bin/perl

die "Usage $0 <sequence> <output_base>\n" if @ARGV != 2;
$seq         = shift(@ARGV);
$output_base = shift(@ARGV);

open(SEQ_HANDLE, $seq) || die "Couldn't read $seq: $!\n";
$SEQ_FILE = join('', <SEQ_HANDLE>);
close(SEQ_HANDLE);

for ($count=0; $count<70; $count++){
   $TR = 300 + $count * 50;

   # Get filename for temporary sequence file
   $temp = "temp_${count}.seq";
 
   # Set Repetition time to TR
   $SEQ_FILE =~ s/(Repetition.*:\s*)\d+\n/${1}${TR}\n/;

   # Save temporary sequence file
   open(SEQ_HANDLE, ">$temp") || die "Couldn't open $temp: $!\n";
   print SEQ_HANDLE $SEQ_FILE;
   close(SEQ_HANDLE);

   $command = "se.pl $temp ${output_base}_${count}";
   system ( "$command" );

   #unlink($temp);
}


