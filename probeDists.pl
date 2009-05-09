

open( INPUT, $ARGV[0] ) || die "COULD NOT OPEN INPUT FILE!\n" ;

$line = <INPUT> ;

while( !( $line =~ /ATOM/ ) )
	{
		$line = <INPUT> ;
	}


@coords = () ;

$line = <INPUT> ;

while( !( $line =~ /BOND/ ) )
	{
		$line =~ s/^\s+// ;

		@words = split( /\s+/, $line ) ;

		push @coords, [ $words[2], $words[3], $words[4] ] ;

		$line = <INPUT> ;
	}

@dists = () ;

for( $i = 0 ; $i < @coords-1 ; ++$i )
	{
		for( $j = $i + 1 ; $j < @coords ; ++$j )
			{
				$dx = $coords[$i][0] - $coords[$j][0] ;
				$dy = $coords[$i][1] - $coords[$j][1] ;
				$dz = $coords[$i][2] - $coords[$j][2] ;

				$d = sqrt( $dx*$dx + $dy*$dy + $dz*$dz ) ;

				push @dists, [ $d, $i, $j ] ;
			}
	}

@distsSort = sort { $$a[0] <=> $$b[0] } @dists ;

print "DIST\t\tPROBES\n" ;

for( $i = 0 ; $i < @distsSort ; ++$i )
	{
	
		print "$distsSort[$i][0]\t\t$distsSort[$i][1]\t$distsSort[$i][2]\n" ;
	}

 
