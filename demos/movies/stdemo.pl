#!/usr/bin/perl

die "Usage $0 <sequence> <output_base>\n" if @ARGV != 2;
$seq         = shift(@ARGV);
$output_base = shift(@ARGV);

open(SEQ_HANDLE, $seq) || die "Couldn't read $seq: $!\n";
$SEQ_FILE = join('', <SEQ_HANDLE>);
close(SEQ_HANDLE);

for ($count=0; $count<64; $count++){

   $ST = sprintf("%2.1f",(1.0 + 0.25*$count));

   # Get filename for temporary sequence file
   $temp = "temp_${count}.seq";
 
   # Set slice thickness to ST
   $SEQ_FILE =~ s/(Slice thickness.*:\s*)(\d+\.?\d*|\.\d+)\n/${1}${ST}\n/;

   # Save temporary sequence file
   open(SEQ_HANDLE, ">$temp") || die "Couldn't open $temp: $!\n";
   print SEQ_HANDLE $SEQ_FILE;
   close(SEQ_HANDLE);

   $command = "se.pl $temp ${output_base}_${count}";
   system ( "$command" );

   #unlink($temp);
   
}


