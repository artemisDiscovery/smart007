#!/usr/bin/perl

## This script tests my arc-arc intersection algorithm. 

## Two arcs are specified by 1) Axis vector (x, y, z), 2) polar angle from axis

if( @ARGV != 8 )
	{
		die "USAGE: testIntAlgo.pl ax1 ay1 az1 angle1 ax2 ay2 az2 angle2\n" ;
	}
	
	
$a1x = $ARGV[0] ;
$a1y = $ARGV[1] ;
$a1z = $ARGV[2] ;

$angle1 = $ARGV[3] ;

$angle1 *= (180./acos(-1.) ) ; #/

### Normalize

$s = sqrt( $a1x*$a1x + $a1y*$a1y + $a1z*$a1z ) ;

$a1x /= $s ;
$a1y /= $s ;
$a1z /= $s ;

### Assume radius of 1.

$d = cos($angle1) ;

$c1x = $d * $a1x ;
$c1y = $d * $a1y ;
$c1z = $d * $a1z ;

$r1 = sin($angle1) ;



$a2x = $ARGV[4] ;
$a2y = $ARGV[5] ;
$a2z = $ARGV[6] ;

$angle2 = $ARGV[7] ;

$angle2 *= (180./acos(-1.) ) ; #/

### Normalize

$s = sqrt( $a2x*$a2x + $a2y*$a2y + $a2z*$a2z ) ;

$a2x /= $s ;
$a2y /= $s ;
$a2z /= $s ;

### Assume radius of 1.

$d = cos($angle2) ;

$c2x = $d * $a2x ;
$c2y = $d * $a2y ;
$c2z = $d * $a2z ;

$r2 = sin($angle2) ;



