//
//  SMVertex.m
//  4SMART
//
//  Created by zauhar on 18/5/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SMVertex.h"
#import "SMArc.h"
#import "SMMol.h"
#import "SMCycle.h"
#import "SMArcEnd.h"

@implementation SMVertex

- (id) init
	{
		self = [ super init ] ;
		
		arcEnds = [ [ NSMutableSet alloc ] initWithCapacity:5 ] ;
		
		position = [ [ MMVector3 alloc ] initX:0. Y:0. Z:0. ] ; 
		normal = [ [ MMVector3 alloc ] initX:0. Y:0. Z:1. ] ;
		
		// DEBUGGING
		
		mergeVertices =  [ [ NSMutableArray alloc ] initWithCapacity:100 ] ;
		
		oldIndex = -1 ;
		
		index = -1 ;
		
		subsurface = -1 ;
		subsurfaceIndex = -1 ;
		
		return self ;
	}
	
- (void) dealloc
	{
		[ arcEnds release ] ;
		
		[ super dealloc ] ;
		
		return ;
	}
	
- (NSSet *) arcEnds
	{
		return arcEnds ;
	}
	
- (void) addArc:(SMArc *)a start:(BOOL)s 
	{
		SMArcEnd *endToAdd = [ [ SMArcEnd alloc ] initWithArc:a atStart:s ] ;
		
		[ arcEnds addObject:endToAdd ] ;
		
		if( s == YES )
			{
				[ a setStartVertex:self ] ;
			}
		else
			{
				[ a setEndVertex:self ] ;
			}
	
		return ;
	}
	
- (void) mergeVertex:(SMVertex *) v
	{
		//printf( "START MERGE of %p and %p \n", self, v ) ;
		//int nHere, nToMerge, nAfterMerge  ;
		
		//nHere = [ arcEnds count ] ;
		//nToMerge = [ [ v arcEnds ] count ] ;
		
		//printf( "\tInitial %d arcEnds:\n", nHere ) ;
		//printf( "%s\n", [ [ arcEnds description ] cString ] ) ; 
		
		//printf( "\t%d Arc Ends to Merge:\n", nToMerge ) ;
		//printf( "%s\n",  [ [ [ v arcEnds ] description ] cString ] ) ; 
		
		if( [ v->mergeVertices count ] == 0 )
			{
				// Merged vertex has no history
				
				[ mergeVertices addObject:v ]  ;
			}
		else
			{
				[ mergeVertices addObject:(v->mergeVertices) ] ;
			}
		
		[ arcEnds unionSet:[ v arcEnds ] ] ;
		
		//nAfterMerge = [ arcEnds count ] ;
		
		//printf( "\t%d AFTER MERGE:\n", nAfterMerge ) ;
		//printf( "%s\n",  [ [ arcEnds description ] cString ] )  ;
		
		//if( nAfterMerge < nHere )
		//	{
		//		printf( "WARNING: MERGE LEADS TO REDUCTION IN CONNECTIONS!\n" ) ;
		//	}
		
		
		
		
		// Update all arcs to point at new vertex
		
		NSEnumerator *arcEndEnumerator ;
		SMArcEnd *nextArcEnd ;
		
		arcEndEnumerator = [ arcEnds objectEnumerator ] ;
		
		while( ( nextArcEnd = [ arcEndEnumerator nextObject ] ) )
			{
				if( [ nextArcEnd atStart ] == YES )
					{
						[ [ nextArcEnd arc ] setStartVertex:self ] ;
					}
				else
					{
						[ [ nextArcEnd arc ] setEndVertex:self ] ;
					}
					
			}
		
		// Average position
		
		[ position setX:( [ position X ] + [ [ v vertexPosition ] X ] )/2. ] ;
		[ position setY:( [ position Y ] + [ [ v vertexPosition ] Y ] )/2. ] ;
		[ position setZ:( [ position Z ] + [ [ v vertexPosition ] Z ] )/2. ] ;
		
		// Average normals
		
		[ normal setX:( [ normal X ] + [ [ v normal ] X ] ) ] ;
		[ normal setY:( [ normal Y ] + [ [ v normal ] Y ] ) ] ;
		[ normal setZ:( [ normal Z ] + [ [ v normal ] Z ] ) ] ;
		
		[ normal normalize ] ;
		
		return ;
	}
		
- (int) printMergeHistory
	{
		NSEnumerator *mergeObjectEnumerator ;
		
		mergeObjectEnumerator = [ mergeVertices objectEnumerator ] ;
		
		id nextMergeObject ;
		
		printf( "(" ) ;
		
		while( ( nextMergeObject = [ mergeObjectEnumerator nextObject ] ) )
			{
				if( [ nextMergeObject isKindOfClass:[ SMVertex class ] ] == YES )
					{
						// Just print the old index
						
						printf( "%d ", ((SMVertex *)nextMergeObject)->oldIndex ) ;
					}
				else
					{						
						[ nextMergeObject printMergeHistory ] ;
					}
			}
			
		printf( ")" ) ;
		
		return 1 ;
	}
		
- (void) setOldIndex:(int)i
	{
		oldIndex = i ;
	
		return ;
	}
	
- (void) computePositionForSaddleUsingMolecule:(SMMol *)m
	{
		// Get a corner - associated arc will provide phi, theta, and torus
		
		SMArcEnd *theArcEnd = [ arcEnds anyObject ] ;
		
		SMArc *theArc = [ theArcEnd arc ] ;
		
		double phi, theta ;
		
		if( [ theArcEnd atStart ] == YES )
			{
				phi = [ theArc phiStart ] ;
				theta = [ theArc thetaStart ] ;
			}
		else
			{
				phi = [ theArc phiEnd ] ;
				theta = [ theArc thetaEnd ] ;
			}
			
		MMVector3 *computePosition = [ theArc computePositionForTheta:theta andPhi:phi usingMolecule:m allowSelfIntersection:m->allowSelfIntersection normal:nil ] ;
		
		//if( fabs( [ computePosition X ] ) > 1000. )
			//{
			//	printf( "WARNING: Generated big coordinate!\n" ) ;
			//}
		
		[ position setX:[ computePosition X ] ] ;
		[ position setY:[ computePosition Y ] ] ;
		[ position setZ:[ computePosition Z ] ] ;
		
		[ computePosition release ] ;
		
		MMVector3 *probePos = [ theArc probePositionForPhi:phi ] ;
		
		[ normal setX:([ probePos X ] - [ position X ]) ] ;
		[ normal setY:([ probePos Y ] - [ position Y ]) ] ;
		[ normal setZ:([ probePos Z ] - [ position Z ]) ] ;
		
		[ normal normalize ] ;
		
		[ probePos release ] ;
		
		return ;
	}
		
		
		
- (void) computePositionForReentrantCycle:(SMCycle *)cyc usingMolecule:(SMMol *)m
	{
		if( m->allowSelfIntersection == NO )
			{
				// Need to check if the vertex is assocaciated with a saddle arc
				
				NSEnumerator *arcEndEnumerator ;
				SMArcEnd *nextArcEnd ;
				SMArc *saddleArc ;
				
				arcEndEnumerator = [ arcEnds objectEnumerator ] ;
				
				double phi, theta ;
				
				phi = theta = -1. ;
				saddleArc = nil ;
				
				while( ( nextArcEnd = [ arcEndEnumerator nextObject ] ) )
					{
						SMTorus *torus ;
						
						if( ( torus = [ [ nextArcEnd arc ] torusSection ] ) != nil )
							{
								saddleArc = [ nextArcEnd arc ] ;
								
								if( [ nextArcEnd atStart ] == YES )
									{
										phi = [ [ nextArcEnd arc ] phiStart ] ;
										theta = [ [ nextArcEnd arc ] thetaStart ] ;
									}
								else
									{
										phi = [ [ nextArcEnd arc ] phiEnd ] ;
										theta = [ [ nextArcEnd arc ] thetaEnd ] ;
									}
									
								break ;
							}
					}
							
				if( saddleArc )
					{
						position = [ saddleArc computePositionForTheta:theta andPhi:phi usingMolecule:m allowSelfIntersection:NO normal:nil ] ;
					}
				else
					{
						// Need to use limit plane, if available, and if reentrant cycle is marked as self-intersecting
						
						SMArcEnd *theArcEnd = [ arcEnds anyObject ] ;
		
						SMArc *theArc = [ theArcEnd arc ] ;
						
						if( cyc->theLimitPlane && cyc->selfIntersection == YES )
							{
								
								if( [ theArcEnd atStart ] == YES )
									{
										[ position release ] ;
										
										position = [ [ MMVector3 alloc ] initUsingVector:[ theArc startPosition ] ] ;
										
										[ cyc->theLimitPlane adjustPosition:position ] ;
									}
								else
									{
										[ position release ] ;
										
										position = [ [ MMVector3 alloc ] initUsingVector:[ theArc endPosition ] ] ;
										
										[ cyc->theLimitPlane adjustPosition:position ] ;
									}
								}
							else
								{
									[ position release ] ;
		
									if( [ theArcEnd atStart ] == YES )
										{
											position = [ [ MMVector3 alloc ] initUsingVector:[ theArc startPosition ] ] ;
										}
									else
										{
											position = [ [ MMVector3 alloc ] initUsingVector:[ theArc endPosition ] ] ;
										}
								}
					}
					
				// Should have something special here!
				
		
				MMVector3 *probePos ;
				
				SMArcEnd *theArcEnd = [ arcEnds anyObject ] ;
		
				SMArc *theArc = [ theArcEnd arc ] ;
									
				probePos = [ theArc hostCenter ] ;
					
				[ normal setX:( [ probePos X ] - [ position X ] ) ] ;
				[ normal setY:( [ probePos Y ] - [ position Y ] ) ] ;
				[ normal setZ:( [ probePos  Z ] - [ position X ] ) ] ;
				
				[ normal normalize ] ;
				
			}
		else
			{
						
						
				SMArcEnd *theArcEnd = [ arcEnds anyObject ] ;
				
				SMArc *theArc = [ theArcEnd arc ] ;
				
				MMVector3 *probePos ;
				
				[ position release ] ;
				
				if( [ theArcEnd atStart ] == YES )
					{
						position = [ [ MMVector3 alloc ] initUsingVector:[ theArc startPosition ] ] ;
					}
				else
					{
						position = [ [ MMVector3 alloc ] initUsingVector:[ theArc endPosition ] ] ;
					}
					
				probePos = [ theArc hostCenter ] ;
					
				[ normal setX:( [ probePos X ] - [ position X ] ) ] ;
				[ normal setY:( [ probePos Y ] - [ position Y ] ) ] ;
				[ normal setZ:( [ probePos  Z ] - [ position X ] ) ] ;
				
				[ normal normalize ] ;
			}
			
		if( fabs( [ position X ] ) > 1000. )
			{
				printf( "WARNING: GENERATED BIG COORDINATE!\n" ) ;
			}
		
		return ;
	}
		
		
- (void) computePositionForContact 
	{
		SMArcEnd *theArcEnd = [ arcEnds anyObject ] ;
		
		SMArc *theArc = [ theArcEnd arc ] ;
		
		MMVector3 *atomPos ;
		
		[ position release ] ;
		
		if( [ theArcEnd atStart ] == YES )
			{
				position = [ [ MMVector3 alloc ] initUsingVector:[ theArc startPosition ] ] ;
			}
		else
			{
				position = [ [ MMVector3 alloc ] initUsingVector:[ theArc endPosition ] ] ;
			}
			
		atomPos = [ theArc hostCenter ] ;
			
		[ normal setX:([ position X ] - [ atomPos X ]) ] ;
		[ normal setY:([ position Y ] - [ atomPos Y ]) ] ;
		[ normal setZ:([ position Z ] - [ atomPos Z ]) ] ;
		
		[ normal normalize ] ;
		
		return ;
	}
		
		
- (MMVector3 *)vertexPosition 
	{
		return position ;
	}
	
- (MMVector3 *)normal 
	{
		return normal ;
	}
		
- (int) index 
	{
		return index ;
	}
	
- (void) setIndex:(int)i 
	{
		index = i ;
	
		return ;
	}
	
- (int) subsurface 
	{
		return subsurface ;
	}
	
- (void) setSubsurface:(int)i 
	{
		subsurface = i ;
		
		return ;
	}
			
		
- (int) subsurfaceIndex 
	{
		return subsurfaceIndex ;
	}
	
- (void) setSubsurfaceIndex:(int)i 
	{
		subsurfaceIndex = i ;
		
		return ;
	}
		
- (NSComparisonResult) compareSubsurfaceVertices:(SMVertex *)v
	{
		if( subsurfaceIndex < [ v subsurfaceIndex ] )
			{
				return NSOrderedAscending ;
			}
		else if( subsurfaceIndex > [ v subsurfaceIndex ] )
			{
				return NSOrderedDescending ;
			}
		else
			{
				return NSOrderedSame ;
			}
			
		return NSOrderedSame ;
	}
		
- (NSString *) description
	{
		NSString *returnString = [ NSString stringWithFormat:@"Vertex:%p RetainCount:%d\n\tArc Ends:\n",self, [ self retainCount ] ] ;
		
		returnString = [ returnString stringByAppendingString:[ arcEnds description ] ]  ;
		
		return returnString ;
	}
		
@end
