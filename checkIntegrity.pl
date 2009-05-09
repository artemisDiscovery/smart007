### This script will check the integrity of a flat surface. 
### It will generate a list of edges, and check for any that are not associated with 
### exactly two elements

open( INPUT, $ARGV[0] ) || die "USAGE: checkIntegrity.pl <flats file> \n" ;

@vertices = () ;

%edgesToElems = () ;

@elems = () ;

$line = <INPUT> ;

@words = split( /\s+/, $line ) ;

$nVert = $words[0] ;

for( $i = 0 ; $i < $nVert ; ++$i )
	{
		$line = <INPUT> ;
		
		$line =~ s/^\s+// ;
		
		@words = split( /\s+/, $line ) ;
		
		push( @vertices, [ $words[0], $words[1], $words[2] ] ) ;
	}
	
$line = <INPUT> ;

@words = split( /\s+/, $line ) ;

$nElem = $words[0] ;

for( $i = 0 ; $i < $nElem ; ++$i )
	{
		$line = <INPUT> ;
		
		$line =~ s/^\s+// ;
		
		@words = split( /\s+/, $line ) ;
		
		push( @elems, [ $words[0], $words[1], $words[2] ] ) ;
		
		if( $words[0] < $words[1] )
			{
				$edge = $words[0]."_".$words[1] ;
			}
		else
			{
				$edge = $words[1]."_".$words[0] ;
			}
			
		if( $edgesToElems{$edge} )
			{
				push( @{$edgesToElems{$edge}}, $i ) ;
			}
		else
			{
				$edgesToElems{$edge} = [ $i ] ;
			}
			
		if( $words[1] < $words[2] )
			{
				$edge = $words[1]."_".$words[2] ;
			}
		else
			{
				$edge = $words[2]."_".$words[1] ;
			}
			
		if( $edgesToElems{$edge} )
			{
				push( @{$edgesToElems{$edge}}, $i ) ;
			}
		else
			{
				$edgesToElems{$edge} = [ $i ] ;
			}
			
		if( $words[0] < $words[2] )
			{
				$edge = $words[0]."_".$words[2] ;
			}
		else
			{
				$edge = $words[2]."_".$words[0] ;
			}
			
		if( $edgesToElems{$edge} )
			{
				push( @{$edgesToElems{$edge}}, $i ) ;
			}
		else
			{
				$edgesToElems{$edge} = [ $i ] ;
			}
			
			
	}
	
### Check for any edges with # elems not 2!

foreach $edge( keys %edgesToElems )
	{
		if( @{$edgesToElems{$edge}} != 2 )
			{
				$num = @{$edgesToElems{$edge}} ;
				
				print "Edge $edge assoc with $num elements:\n" ;
				
				for( $i = 0 ; $i < $num ; ++$i )
					{
						print "ELEM: $edgesToElems{$edge}[$i] VERTICES: $elems[$edgesToElems{$edge}[$i]][0]  $elems[$edgesToElems{$edge}[$i]][1] $elems[$edgesToElems{$edge}[$i]][2]\n" ;
					}
					
			}
	}
	
	

			
		
		