#!/usr/bin/perl

### Updated script to perturb atomic coordinates

### Takes as arguments mol2 file and maximum displacement


if( @ARGV != 2 )
	{	
		die "USAGE: perturbCoords.pl <mol2 file> <max disp> \n" ;
	}
	
$mol2File = $ARGV[0] ;
$maxDisp = $ARGV[1] ;


open( INPUT, $mol2File ) || die "COULD NOT OPEN $mol2FIle - Exit!\n" ;

$STATE = LookForAtom ;

while( $line = <INPUT> )
	{
		if( $STATE eq LookForAtom )
			{
				print $line ;
				
				if( $line =~  /@\<TRIPOS\>ATOM/ )
					{
						$STATE = CollectAtoms ;
					}
			}
		elsif( $STATE eq CollectAtoms )
			{	
				if( $line =~ /^@/ )
					{
						$STATE = Finish ;
						print $line ;
						
						next ;
					}
				
					
				### Compute scaled coordinates
				
				$first = substr( $line, 0, 16 ) ;
				
				$x = substr( $line, 16, 12 ) ;
				$y = substr( $line, 28, 12 ) ;
				$z = substr( $line, 40, 12 ) ;
				
				$last = substr( $line, 52 ) ;
				
				$disp = 2.*$maxDisp*rand(1.) - $maxDisp ;
				$x += $disp ;
				$disp = 2.*$maxDisp*rand(1.) - $maxDisp ;
				$y += $disp ;
				$disp = 2.*$maxDisp*rand(1.) - $maxDisp ;
				$z += $disp ;
				
				$coord = sprintf( "%10.4f%10.4f%10.4f", $x, $y, $z ) ;
				
				print $first, $coord, $last ;
			}
		else
			{
				print $line ;
			}
	}
	
close INPUT ;


