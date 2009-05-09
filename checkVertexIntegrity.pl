### This script will check the integrity of the vertices of a flats file 
### It will check if any vertices have same coordinate

open( INPUT, $ARGV[0] ) || die "USAGE: checkVertexIntegrity.pl <flats file> \n" ;

@vertices = () ;

%coordsToVertices = () ;


$line = <INPUT> ;

@words = split( /\s+/, $line ) ;

$nVert = $words[0] ;

$count = 0 ;

for( $i = 0 ; $i < $nVert ; ++$i )
	{
		$line = <INPUT> ;
		
		$line =~ s/^\s+// ;
		
		@words = split( /\s+/, $line ) ;
		
		push( @vertices, [ $words[0], $words[1], $words[2] ] ) ;
		
		$key = $words[0]."_".$words[1]."_".$words[2] ;
		
		if( $coordsToVertices{$key} )
			{
				push @{$coordsToVertices{$key}}, $count ;
			}
		else
			{
				$coordsToVertices{$key} = [ $count ] ;
			}
		
		++$count ;
		
		
	}
	
	

foreach $key( keys %coordsToVertices )
	{
		@verts = @{$coordsToVertices{$key}} ;
		
		if( @verts > 1 )
			{
				print "Vertices ", "@verts", " share coord: ", $key, "\n" ;
			}
	}
		
		