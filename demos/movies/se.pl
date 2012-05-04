#!/usr/bin/env perl -w

# sim_rf.pl
#
# Simulates a volume with receive, transmit, and both receive and transmit
# inhomogeneities.

require 5.001;

use JobControl;
use Batch;

$Verbose     = 1;
$Execute     = 1;
$Batch       = 1;
$MergeErrors = 1;

&JobControl::SetOptions("Verbose", $Verbose, "Execute", $Execute,
                        "Batch", $Batch, "MergeErrors", $MergeErrors);

die "Usage $0 <sequence> <output_base>\n" if @ARGV != 2;
$seq_file    = shift(@ARGV);
$output_base = shift(@ARGV);

$mrisim      = "mrisim";

$phantom  = "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_bck.mnc:0," .
            "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_csf.mnc:1," .
            "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_gry.mnc:2," .
            "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_wht.mnc:3," .
            "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_fat.mnc:4," .
            "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_m+s.mnc:5," .
            "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_skn.mnc:6," .
            "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_skl.mnc:7," .
            "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_gli.mnc:8," .
            "/data/icbm/sbd/phantoms/1.0mm/phantom_1.0mm_normal_mit.mnc:9";
$tissue   = "/data/icbm/sbd/sequences/tissue_normal.prm";
$sequence = "$seq_file";
$coil     = "coil.rf";
#$coil     = "noiseless.rf";

$mrisim_command = "$mrisim -fuzzy $phantom -tissue $tissue -coil $coil " .
            "-sequence $sequence " .
            "-clob ${output_base}.mnc";

&Spawn ( $mrisim_command, "${output_base}.log" );
#&Spawn ( "$reshape ${output_base}_tmp.mnc ${output_base}.mnc " .
#         "-start 0,0,0 -count 181,217,181", "${output_base}.log" );

