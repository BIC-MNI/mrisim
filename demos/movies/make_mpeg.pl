#!/usr/local/bin/perl

    $output = shift;

    if( $output )
    {
        $out_file = " -s $output ";
    }
    else
    {
        die( "Usage: $0 output.mpeg basename nframes ...\n" );
    }

    $base = shift;
    $nframes = shift;

    @files = ();
    for ($count=1; $count<=$nframes; $count++){
       @files = (@files,"${base}${count}.tiff");
    }

    $n_frames = 0;
 
    while( $file = shift(@files) )
    {
        ++$n_frames;

        if( $n_frames == 1 )
        {
            $out = `imginfo $file | grep Dimensions`;
            $out =~ /\s+(\d+),\s+(\d+)/;
            $x_size = $1;
            $y_size = $2;
        }

        $dir = $file;
        $dir =~ /(.*\/)/;
        $dir = $1;
        $ppm_file = "${dir}frame_${n_frames}.ppm";
        $yuv_file_prefix = "${dir}frame_${n_frames}";

        @yuvs = ( @yuvs, $yuv_file_prefix );
 
        print( "File: $file -> $ppm_file\n" );
        system( "tifftopnm $file > $ppm_file" ); 
        system( "ppmtoyuvsplit $yuv_file_prefix $ppm_file" ); 
        unlink( $ppm_file );
    }

    system( "mpeg $out_file -h $x_size -v $y_size -q 1 -a 1 -b $n_frames ${dir}frame_" );

    foreach $file ( @yuvs )
    {
        unlink( "${file}.Y" );
        unlink( "${file}.U" );
        unlink( "${file}.V" );
    }
