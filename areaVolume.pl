#!/usr/bin/perl

### Compute area and volume for a flats surface
### I think that either new or old format should work with this script

if( @ARGV != 1 )
	{
		die "USAGE: areaVolume.pl <flats file>\n" ;
	}
	
open INPUT, $ARGV[0] || die "COULD NOT OPEN INPUT FILE - Exit!\n" ;

### Get number of vertices

$line = <INPUT> ;

$line =~ s/^\s+// ;

$nVertices = ( split( /\s+/, $line ) )[0] ;

@vertices = () ;

for( $i = 0 ; $i < $nVertices ; ++$i )
	{
		$line = <INPUT> ;
		
		$line =~ s/^\s+// ;
		
		@words = split( /\s+/, $line ) ;
		
		$x = $words[0] ;
		$y = $words[1] ;
		$z = $words[2] ;
		
		push @vertices, [ $x, $y, $z ] ;
	}
	
### Get number of elements

$line = <INPUT> ;

$line =~ s/^\s+// ;

$nElements = ( split( /\s+/, $line ) )[0] ;

$area - 0. ;
$volume = 0. ;

for( $i = 0 ; $i < $nElements ; ++$i )
	{
		$line = <INPUT> ;
		
		$line =~ s/^\s+// ;
		
		@words = split( /\s+/, $line ) ;
		
		$v1 = $words[0] ;
		$v2 = $words[1] ;
		$v3 = $words[2] ;
		
		$midPointX = 0. ;
		$midPointY = 0. ;
		$midPointZ = 0. ;
		
		foreach $v( $v1, $v2, $v3 )
			{
				$midPointX += $vertices[$v][0] ;
				$midPointY += $vertices[$v][1] ;
				$midPointZ += $vertices[$v][2] ;
			}
			
		$midPointX /= 3. ;
		$midPointY /= 3. ;
		$midPointZ /= 3. ;
		
		$nX = $words[3] ;
		$nY = $words[4] ;
		$nZ = $words[5] ;
		
		### area from size of cross produce
		
		$delta12X = $vertices[$v2][0] - $vertices[$v1][0] ;
		$delta12Y = $vertices[$v2][1] - $vertices[$v1][1] ;
		$delta12Z = $vertices[$v2][2] - $vertices[$v1][2] ;
		
		$delta13X = $vertices[$v3][0] - $vertices[$v1][0] ;
		$delta13Y = $vertices[$v3][1] - $vertices[$v1][1] ;
		$delta13Z = $vertices[$v3][2] - $vertices[$v1][2] ;
		
		$crossX = $delta12Y*$delta13Z - $delta13Y*$delta12Z ;
		$crossY = $delta12Z*$delta13X - $delta13Z*$delta12X ;
		$crossZ = $delta12X*$delta13Y - $delta13X*$delta12Y ;
		
		$A = (1./2.) * sqrt( $crossX*$crossX + $crossY*$crossY + $crossZ*$crossZ ) ;
		
		$area += $A ;
		
		$dot = $midPointX*$nX + $midPointY*$nY + $midPointZ*$nZ ;
		
		$volume += $A * $dot ;
	}
	
$volume /= 3. ;
		
		
		
print "\nAREA = $area ANG^2\n" ;
print "\nVOLUME = $volume ANG^3\n" ;


		
		
		
		