//
//  SMArc.m
//  4SMART
//
//  Created by Randy Zauhar on 1/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SMArc.h"
#import "SMMol.h"
#import "SMTorus.h"
#include <math.h>
#import "SMCycle.h"


@implementation SMArc

- (id) initWithHostCenter:(MMVector3 *)hc radius:(double)hr torusSection:(SMTorus *)ts arcType:(arcClass)at
	{
		self = [ super init ] ;
		
		hostCenter = [ [ MMVector3 alloc ] initUsingVector:hc ] ;
		
		hostRadius = hr ;
		
		arcCenter = arcAxis = nil ;
		
		arcRadius = 0. ;
		
		startU = endU = nil ;
		
		lineSegment = NO ;
		
		length = angle = 0.  ;
		
		selfIntersectionType = UNKNOWN ; 
	
		hostProbe = nil ;
		
		// Weak reference, torus section is never dealloced
		torusSection = ts ;
		
		arcType = at ;
		
		parentArc = nil ;
		parentCycles = nil ;
		
		startPosition = nil ;
		endPosition = nil ;
		
		startTangent = nil ;
		endTangent = nil ;
		
		startVertex = nil ;
		endVertex = nil ;
		
		thetaStart = -1. ;
		thetaEnd = -1. ; 
		phiStart = -1. ; 
		phiEnd 	= -1. ;
				
		twin = nil ;
		
		skip = NO ;
		
		//collapseArc = nil ;
		
		startConnections = [ [ NSMutableArray alloc ] initWithCapacity:5 ]  ;
		endConnections  = [ [ NSMutableArray alloc ] initWithCapacity:5 ] ;
	
		startConnectStart = [ [ NSMutableArray alloc ] initWithCapacity:5 ]  ;
		endConnectStart = [ [ NSMutableArray alloc ] initWithCapacity:5 ]  ;

		
		return self ;
	}

		
- (id) copyArc
	{
		SMArc *arcCopy = [ [ SMArc alloc ] init ] ;
		
		arcCopy->hostCenter = [ [ MMVector3 alloc ] initUsingVector:hostCenter ] ;
				
		arcCopy->hostRadius = hostRadius ;
		
		arcCopy->arcCenter = [ [ MMVector3 alloc ] initUsingVector:arcCenter ] ;
		
		arcCopy->arcAxis = [ [ MMVector3 alloc ] initUsingVector:arcAxis ] ;
		
		arcCopy->arcRadius = arcRadius ;
		
		arcCopy->startU = [ [ MMVector3 alloc ] initUsingVector:startU ] ;
		
		arcCopy->endU = [ [ MMVector3 alloc ] initUsingVector:endU ] ;
		
		arcCopy->uPerp = [ [ MMVector3 alloc ] initUsingVector:uPerp ] ;
		
		arcCopy->lineSegment = lineSegment ;
		
		arcCopy->length = length ;
		
		arcCopy->angle = angle  ;
				
		
		arcCopy->thetaStart = thetaStart ;
		arcCopy->thetaEnd = thetaEnd ; 
		arcCopy->phiStart = phiStart ; 
		arcCopy->phiEnd 	= phiEnd ;
		
		
		arcCopy->skip = skip ;
	
		arcCopy->selfIntersectionType = selfIntersectionType ;
							
		arcCopy->hostProbe = hostProbe ;
		[ hostProbe retain ] ;
		
		arcCopy->torusSection = torusSection ;
		[ torusSection retain ] ;
		
		arcCopy->arcType = arcType ;
		
		arcCopy->parentArc = parentArc ;
		[ parentArc retain ] ;
		
		if( parentCycles )
			{
				arcCopy->parentCycles = [ [ NSMutableArray alloc ] initWithArray:parentCycles ] ;
			}
		else
			{
				arcCopy->parentCycles = nil ;
			}
		
		arcCopy->startPosition = [ [ MMVector3 alloc ] initUsingVector:startPosition ] ;
		arcCopy->endPosition = [ [ MMVector3 alloc ] initUsingVector:endPosition ] ;
		
		arcCopy->startTangent = [ [ MMVector3 alloc ] initUsingVector:startTangent ] ;
		arcCopy->endTangent = [ [ MMVector3 alloc ] initUsingVector:endTangent ] ;
		
		arcCopy->startVertex = startVertex ;
		arcCopy->endVertex = endVertex ;
		
		if( startConnections )
			{
				arcCopy->startConnections = [ [ NSMutableArray alloc ] initWithArray:startConnections ] ;
			}
	
		if( endConnections )
			{
				arcCopy->endConnections = [ [ NSMutableArray alloc ] initWithArray:endConnections ] ;
			}
			
		if( startConnectStart )
			{
				arcCopy->startConnectStart = [ [ NSMutableArray alloc ] initWithArray:startConnectStart ] ;
			}

		if( endConnectStart )
			{
				arcCopy->endConnectStart = [ [ NSMutableArray alloc ] initWithArray:endConnectStart ] ;
			}

		
		return arcCopy ;
	}
	
	
- (void) reverse
	{
		// This method reverses the physical orientation of an arc. It currently has one application - when combining contact cycles on the same atom,
		// we need to link them with a pair of edges with opposite physical orientation. These are "twinned" so that vertices can be correctly generated
		// at the end.
		
		MMVector3 *swap ;
		double angleSwap ;
		SMVertex *vertexSwap ;
		
		[ arcAxis reverse ] ;
		
		swap = startU ;
		startU = endU ;
		endU = swap ;
		
		if( uPerp ) [ uPerp release ] ;
		
		uPerp = [ [ MMVector3 alloc ] initByCrossing:arcAxis and:startU ] ;
		
		[ uPerp normalize ] ;
		
		swap = startPosition ;
		startPosition = endPosition ;
		endPosition = swap ;
		
		swap = startTangent ;
		startTangent = endTangent ;
		endTangent = swap ;
		
		// Need to reverse tangents!
		
		[ startTangent reverse ] ;
		[ endTangent reverse ] ;
		
		angleSwap = thetaStart ;
		thetaStart = thetaEnd ;
		thetaEnd = angleSwap ;
		
		angleSwap = phiStart ;
		phiStart = phiEnd ;
		phiEnd = angleSwap ;
		
		vertexSwap = startVertex ;
		startVertex = endVertex ;
		endVertex = vertexSwap ;
		
		return ;
	}
		
		

- (void) initializeWithArcCenter:(MMVector3 *)ac arcRadius:(double)ar axis:(MMVector3 *)ax start:(MMVector3 *)s end:(MMVector3 *)e hostProbe:(SMProbe *)hp
	{
		arcRadius = ar ;
		
		hostProbe = hp ;	// Should be no need to retain
		
		[ ac retain ] ;
		
		if( arcCenter ) [ arcCenter release ] ;
		
		arcCenter = ac ;
	
		[ ax retain ] ;
		
		if( arcAxis ) [ arcAxis release ] ;
		
		arcAxis = ax ;
		
		[ s retain ] ;
		
		if( startU ) [ startU release ] ;
		
		startU = s ;

		[ e retain ] ;
		
		if( endU ) [ endU release ] ;
		
		endU = e ;
		
		// Normalize just in case
		
		[ startU normalize ] ;
		[ endU normalize ] ;
		
		// Compute arc length (assume limit plane not in effect)
		
		MMVector3 *crossProduct = [ [ MMVector3 alloc ] initByCrossing:startU and:endU ] ;
		
		double dot ;
		
		dot = [ startU dotWith:endU ] ;
		
		if( fabs(dot) > 1. )
			{
				dot = dot/fabs(dot) ;
			}
		
		if( [ crossProduct dotWith:arcAxis ] > 0. )
			{
				angle = acos( dot ) ;
			}
		else
			{
				angle = 2.*acos(-1.) - acos( dot ) ;
			}
			
		length = arcRadius * angle ;
		
		[ crossProduct release ] ;
		
		uPerp = [ [ MMVector3 alloc ] initByCrossing:arcAxis and:startU ] ;
		
		[ uPerp normalize ] ;
		
		// Tangents
		
		startTangent = [ [ MMVector3 alloc ] initByCrossing:arcAxis and:startU ] ;
		[ startTangent normalize ] ;

		endTangent = [ [ MMVector3 alloc ] initByCrossing:arcAxis and:endU ] ;
		[ endTangent normalize ] ;
		
		// Start & end positions
		
		startPosition = [ [ MMVector3 alloc ] initX:([ arcCenter X ] + arcRadius*[ startU X ])
			Y:([ arcCenter Y ] + arcRadius*[ startU Y ])
			Z:([ arcCenter Z ] + arcRadius*[ startU Z ]) ] ;
			
		endPosition = [ [ MMVector3 alloc ] initX:([ arcCenter X ] + arcRadius*[ endU X ])
			Y:([ arcCenter Y ] + arcRadius*[ endU Y ])
			Z:([ arcCenter Z ] + arcRadius*[ endU Z ]) ] ;
			
		
		
		return ;
	}
		
		
- (id) initWithTorusSection:(SMTorus *)ts molecule:(SMMol *)mol phiStart:(double)phiS thetaStart:(double)thetaS phiEnd:(double)phiE thetaEnd:(double)thetaE
	{
		// This method creates a new saddle arc
		
		self = [ super init ] ;
		
		hostCenter = nil ;
		
		hostRadius = 0. ;
		
		arcCenter = arcAxis = nil ;
		
		arcRadius = 0. ;
		
		startU = endU = nil ;
		
		lineSegment = NO ;
		
		length = angle = 0.  ;
		
		selfIntersectionType = UNKNOWN  ;
	
		hostProbe = nil ;
		
		torusSection = ts ;
		
		arcType = INTERIORSADDLE ;
		
		parentArc = nil ;
		parentCycles = nil ;
		
		startPosition = nil ;
		endPosition = nil ;
		
		startTangent = nil ;
		endTangent = nil ;
		
		startVertex = nil ;
		endVertex = nil ;
		
		thetaStart = thetaS ;
		thetaEnd = thetaE ; 
		phiStart = phiS ; 
		phiEnd 	= phiE ;
		
		twin = nil ;
		
		startConnections = [ [ NSMutableArray alloc ] initWithCapacity:5 ]  ;
		endConnections  = [ [ NSMutableArray alloc ] initWithCapacity:5 ] ;
	
		startConnectStart = [ [ NSMutableArray alloc ] initWithCapacity:5 ]  ;
		endConnectStart = [ [ NSMutableArray alloc ] initWithCapacity:5 ]  ;
	
		
		
		// Compute length by numerical approximation. I will not take limit plane into account 
		
		startPosition = [ self computePositionForTheta:thetaStart andPhi:phiStart usingMolecule:mol allowSelfIntersection:mol->allowSelfIntersection normal:nil ] ;
		
		// Use five divisions of the arc - use first and last as a rough approximation of start and end tangents
		
		MMVector3 *prePos, *pos, *diff ;
		double theta, phi, deltaTheta, deltaPhi ;
		int i ;
		
		diff = [ [ MMVector3 alloc ] init ] ;
		
		deltaTheta = (thetaEnd - thetaStart)/5. ;
		deltaPhi = (phiEnd - phiStart)/5. ;
		
		prePos = startPosition ;
		
		for( i = 0 ; i < 5 ; ++i )
			{
				phi = phiStart + (i + 1)*deltaPhi ;
				theta = thetaStart + (i + 1)*deltaTheta ;
				
				pos = [ self computePositionForTheta:theta andPhi:phi usingMolecule:mol allowSelfIntersection:mol->allowSelfIntersection  normal:nil ] ;
				
				[ diff setX:([ pos X ] - [ prePos X ]) ] ;
				[ diff setY:([ pos Y ] - [ prePos Y ]) ] ;
				[ diff setZ:([ pos Z ] - [ prePos Z ]) ] ;
				
				length += [ diff length ] ;
				
				if( i == 0 )
					{
						startTangent = diff ;
						[ startTangent normalize ] ;
						
						diff = [ [ MMVector3 alloc ] init ] ;
						
						prePos = pos ;
					}
				else if( i == 4)
					{
						endTangent = diff ;
						[ endTangent normalize ] ;
						
						endPosition = pos ;
						
						[ prePos release ] ;
					}
				else
					{
						[ prePos release ] ;
						prePos = pos ;
					}
						
			}
			
		// For debugging purposed, I am going to reset the length to simply the distance between start and end points
		/*
		diff = [ [ MMVector3 alloc ] init ] ;
		
		[ diff setX:([ endPosition X ] - [ startPosition X ]) ] ;
		[ diff setY:([ endPosition Y ] - [ startPosition Y ]) ] ;
		[ diff setZ:([ endPosition Z ] - [ startPosition Z ]) ] ;
		
		length = [ diff length ] ;
		
		[ diff release ] ;
		*/
		
		
		return self ;
	}
		
		

- (void) dealloc
	{
	
	
		if( hostCenter ) [ hostCenter release ] ;
		
		if( arcCenter ) [ arcCenter  release ] ;
		if( arcAxis ) [ arcAxis release ] ;
		
		if( startU ) [ startU release ] ;
		if( endU ) [ endU release ] ;
		if( uPerp ) [ uPerp release ] ;
		
		if( startPosition ) [ startPosition release ] ; 
		if( endPosition ) [ endPosition release ]  ;
		
		if( startTangent ) [ startTangent release ] ; 
		if( endTangent ) [ endTangent release ]  ;
		
		[ super dealloc ] ;
		
		return ;
	}
	
- (NSArray *) interiorArcsForGeodesicWithCenter:(MMVector3 *)cntr limitAxis:(MMVector3 *)lmaxis minimumDisplacement:(double) minDisp arcSize:(double) size
	{
		// This method will work for geodesic arcs, whether on atom or probe surface. It can handle, if needed, a limit plane or torus
		// buffer zone. 
	
		double dx, dy, dz ;
	
		if ( selfIntersectionType == BOTHTOUCH ) {
			// Treat this as a straight line
			
			// This class of "geodesic" arc is treated as a simple straight line between 'startPosition' and 'endPosition' vectors
			
			MMVector3 *line = [ [ MMVector3 alloc ] initBySubtracting:endPosition minus:startPosition ] ;
			
			double len = [ line length ] ;
			
			[ line normalize ] ;
			
			int nDiv = len / size ;
			int iDiv ;
			
			if (nDiv <= 1 ) {
				return nil ;
			}
			
			NSMutableArray *returnArray = [ NSMutableArray arrayWithCapacity:5 ] ;
			
			double div = [ a length ] / nDiv ;
			
			for( iDiv = 0 ; iDiv < nDiv ; ++iDiv  ) {
				MMVector3 *s, *e ;
				
				if (iDiv == 0 ) {
					s = startU ;
					
					e = [ [ MMVector3 alloc ] initUsingVector:line ] ;
					[ e scaleBy:div ] ;
					[ e addVector:startPosition ] ;
					
					SMArc *newArc = [ self copyArc ] ;
					newArc->
					
					
					
					
				}
			}
			
			
			
			
			
		}
	
		
	}
- (NSArray *) subdivisionDataWithDivisionSize:(double) div 
	{
		// This critical function finds the division points to be used when subdividing an arc. It works with either geodesic arcs in reentrant 
		// regions, or arcs in saddle sections. 
	
		// For simple circular arcs/geodesics, the return is an array of normalized vectors directed at the division points. 
	
		// For reentrant geodesics, the return array includes, for each division point, an array containing (1) a vector directed from arc center to 
		// division point, 2) a flag indicating whether the point touches either a saddle buffer zone or limit plane.
	
		// For interior saddle arcs, the retuen array contains, for each division point, an array of three NSNumbers - theta, phi, and 
		// a flag indicating touching the buffer zone of the saddle.
	
		NSMutableArray *returnArray = [ NSMutableArray arrayWithCapacity:10 ] ;
	
		switch ( a->arcType ) {
			case CONTACTI:
			case CONTACTJ:
				
				int nDiv ;
				
				if ( div <= 0. ) {
					 nDiv = 2 ;
				}
				else {
					nDiv = a->length / div ;
				}
				
				if (nDiv <= 1 ) {
					return nil ;
				}

				double angleIncrement ;
				
				angleIncrement = ((double)[ a angle ]) / nDiv ;
				
				double dx, dy, dz ;
				
				int iDiv ;
				
				for( iDiv = 0 ; iDiv < nDiv ; ++iDiv ) {
					
					double ang = iDiv * angleIncrement ;
					
					dx = cos(ang)*[ [ a startU ] X ] + sin(ang)*[ [ a uPerp ] X ] ;
					dy = cos(ang)*[ [ a startU ] Y ] + sin(ang)*[ [ a uPerp ] Y ] ;
					dz = cos(ang)*[ [ a startU ] Z ] + sin(ang)*[ [ a uPerp ] Z ] ;
					
					MMVector3 *divideVector ;
					
					divideVector = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
					
					[ divideVector normalize ] ;
					
					[ returnArray addObject:divideVector ] ;
					[ divideVector release ] ;
				}
			break;

				case REENTRANTL:
				case REENTRANTR:
					
					// Need to check for intersection with buffer zone - collect subarcs to process individually
				
					// Consider that we may have theta start > theta end, in which case we reverse order of thetaBufferLo and thetaBufferHi
					// To make things a little easier, work in "delta thetas" from one end of the other
					double theta1, theta2, thetaBuff1, thetaBuff2, sign ;
					BOOL span ; // Does arc enclose the self-intersection zone?
				
					if (a->torusSection->thetaStart < a->torusSection->thetaEnd ) {
						theta1 = 0. ;
						theta2 = a->torusSection->thetaEnd - a->torusSection->thetaStart ;
						thetaBuff1 = a->torusSection->thetaBufferLo - a->torusSection->thetaStart ;
						thetaBuff2 = a->torusSection->thetaBufferHi - a->torusSection->thetaStart ;
						sign = 1. ;
					}
					else {
						theta1 = 0. ;
						theta2 = a->torusSection->thetaStart - a->torusSection->thetaEnd ;
						thetaBuff1 = a->torusSection->thetaBufferHi - a->torusSection->thetaEnd ;
						thetaBuff2 = a->torusSection->thetaBufferLo - a->torusSection->thetaEnd ;
						sign = -1. ;
					}

					BOOL inside1, inside2 ;
				
					if (thetaBuff1 <= theta1 && theta1 <= thetaBuff2 ) {
						inside1 = YES ;
					}
					else {
						inside1 = NO ;
					}
				
					if (thetaBuff1 <= theta2 && theta2 <= thetaBuff2 ) {
						inside2 = YES ;
					}
					else {
						inside2 = NO ;
					}
				
					if (theta1 <= thetaBuff1 && theta2 >= thetaBuff2 ) {
						span = YES ;
					}
					else {
						span = NO ;
					}


					NSMutableArray *subArcsByTheta = [ NSMutableArray arrayWithCapacity:3 ] ;
				
					if (inside1 == YES && inside2 == YES ) {
						[ subArcsByTheta 
						 addObject:[ NSArray arrayWithObjects:[ NSNumber numberWithDouble:theta1 ], 
																[ NSNumber numberWithDouble:theta2 ], [ NSNumber numberWithInt:(int)BOTHTOUCH ], nil ] ] ;
					}
					else if (inside1 == YES && inside2 == NO ) {
						[ subArcsByTheta 
						 addObject:[ NSArray arrayWithObjects:[ NSNumber numberWithDouble:theta1 ], 
									[ NSNumber numberWithDouble:thetaBuff2 ], [ NSNumber numberWithInt:(int)BOTHTOUCH ], nil ] ] ;
						[ subArcsByTheta 
						 addObject:[ NSArray arrayWithObjects:[ NSNumber numberWithDouble:thetaBuff2 ], 
									[ NSNumber numberWithDouble:theta2 ], [ NSNumber numberWithInt:(int)STARTTOUCH ], nil ] ] ;
					}
					else if (inside1 == NO && inside2 == YES ) {
						[ subArcsByTheta 
						 addObject:[ NSArray arrayWithObjects:[ NSNumber numberWithDouble:theta1 ], 
									[ NSNumber numberWithDouble:thetaBuff1 ], [ NSNumber numberWithInt:(int)ENDTOUCH ], nil ] ] ;
						[ subArcsByTheta 
						 addObject:[ NSArray arrayWithObjects:[ NSNumber numberWithDouble:thetaBuff1 ], 
									[ NSNumber numberWithDouble:theta2 ], [ NSNumber numberWithInt:(int)BOTHTOUCH ], nil ] ] ;
					}
					else {
						// Both flags are NO
						
						if (span == YES ) {
							[ subArcsByTheta 
							 addObject:[ NSArray arrayWithObjects:[ NSNumber numberWithDouble:theta1 ], 
										[ NSNumber numberWithDouble:thetaBuff1 ], [ NSNumber numberWithInt:(int)ENDTOUCH ], nil ] ] ;
							[ subArcsByTheta 
							 addObject:[ NSArray arrayWithObjects:[ NSNumber numberWithDouble:thetaBuff1 ], 
										[ NSNumber numberWithDouble:thetaBuff2 ], [ NSNumber numberWithInt:(int)BOTHTOUCH ], nil ] ] ;
							[ subArcsByTheta 
							 addObject:[ NSArray arrayWithObjects:[ NSNumber numberWithDouble:thetaBuff2 ], 
										[ NSNumber numberWithDouble:theta2 ], [ NSNumber numberWithInt:(int)STARTTOUCH ], nil ] ] ;
						}
					}

				// Now that subarcs are defined, process individually to set interior theta values
				
				// Need theta reference - this actually only works for a reentrant boundary arc (probe fixed)
				
				MMVector3 *thetaRef, *thetaPerp ;
				
				double dx = [ [ a hostProbe ] X ] - xAtom[ a->torusSection->atomI ] ;
				double dy = [ [ a hostProbe ] Y ] - yAtom[ a->torusSection->atomI ] ;
				double dz = [ [ a hostProbe ] Z ] - zAtom[ a->torusSection->atomI ] ;
				
				thetaRef = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
				[ thetaRef normalize ] ;
				
				thetaPerp = [ [ MMVector3 alloc ] initAlong:a->torusSection->axis perpTo:thetaRef ] ;
				
				
				
				
				
					
				
					
				
			break;

			
				break;
			default:
				break;
		}
	
		
	}


- (double) angle
	{
		return angle ;
	}
	
- (BOOL) skip 
	{
		return skip ;
	}

- (void) setSkip:(BOOL)s 
	{
		skip = s ;
	
		return ;
	}
	
	
- (MMVector3 *) hostCenter 
	{
		return hostCenter ;
	}

- (double) hostRadius 
	{
		return hostRadius ;
	}
	
- (double) length
	{	
		return length ;
	}
	
- (SMProbe *) hostProbe
	{
		return hostProbe ;
	}
	
- (double) arcRadius
	{
		return arcRadius ;
	}
	
- (SMTorus *) torusSection 
	{
		return torusSection ;
	}
	
	
	
- (arcClass) arcType 
	{
		return arcType ;
	}
	
	

- (void) setStartVertex: (SMVertex *)s 
	{
		startVertex = s ;
		
		return ;
	}

- (void) setEndVertex: (SMVertex *)e 
	{
		endVertex = e ;
				
		return ;
	}
	
- (void) setStartU:(MMVector3 *)sU 
	{
		[ sU retain ] ;
		
		if( startU ) [ startU release ] ;
		
		startU = sU ;
		
		return ;
	}
	
- (void) setEndU:(MMVector3 *)eU 
	{
		[ eU retain ] ;
		
		if( endU ) [ endU release ] ;
		
		endU = eU ;
		
		return ;
	}

- (MMVector3 *) startTangent 
	{
		return startTangent ;
	}
	
- (MMVector3 *) endTangent 
	{
		return endTangent ;
	}

	

- (SMVertex *) startVertex 
	{
		return startVertex ;
	}
	
- (SMVertex *) endVertex 
	{
		return endVertex ;
	}

- (NSMutableArray *) parentCycles 
	{
		return parentCycles ;
	}
	
- (void) addParentCycle:(SMCycle *)c 
	{
		if( ! parentCycles )
			{
				parentCycles = [ [ NSMutableArray alloc ] initWithCapacity:2 ] ;
			}
			
		if( [ parentCycles containsObject:c ] == NO )
			{
				[ parentCycles addObject:c ] ;
			}
		
		if( [ parentCycles count ] > 2 )
			{
				printf( "WARNING - ARC PARTICIPATES IN %d CYCLES\n", (int)[ parentCycles count ] ) ;
			}
		
		return ;
	}
	
- (void) mergeParentCyclesFromArc:(SMArc *)a
	{
		NSMutableSet *cycleSet = [ [ NSMutableSet alloc ] initWithArray:[ a parentCycles ] ] ;
		
		[ cycleSet addObjectsFromArray:parentCycles ] ;
		
		[ parentCycles setArray:[ cycleSet allObjects ] ] ;
		
		[ cycleSet release ] ;
		
		return ;
	}
		
		
	
- (void) removeParentCycle:(SMCycle *)c 
	{
		if( ! parentCycles )
			{
				return ;
			}
			
		[ parentCycles removeObject:c ] ;
				
		return ;
	}
	

- (MMVector3 *) startPosition 
	{
		return startPosition ;
	}
	
- (MMVector3 *) endPosition ;
	{
		return endPosition ;
	}
	
- (MMVector3 *) startU
	{
		return startU ;
	}
	
- (MMVector3 *) endU 
	{
		return endU ;
	}
	
- (MMVector3 *) uPerp 
	{
		return uPerp ;
	}
	
- (MMVector3 *) arcAxis
	{
		return arcAxis ;
	}
	
-(MMVector3 *) arcCenter 
	{
		return arcCenter ;
	}
	
	
- (void) setStartVertexPosition:(MMVector3 *)s endVertexPosition:(MMVector3 *)e
	{
		[ s retain ] ;
		
		if( startPosition ) [ startPosition release ] ;
		
		startPosition = s ;
		
		[ e retain ] ;
		
		if( endPosition ) [ endPosition release ] ;
		
		endPosition = e ;
		
		return ;
	}

- (void) computeVertexPositions 
	{
		double dx, dy, dz ;
		
		dx = [ arcCenter X ] + (arcRadius * [ startU X ]) ; 
		dy = [ arcCenter Y ] + (arcRadius * [ startU Y ]) ; 
		dz = [ arcCenter Z ] + (arcRadius * [ startU Z ]) ; 
		
		startPosition = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
		
		dx = [ arcCenter X ] + (arcRadius * [ endU X ]) ; 
		dy = [ arcCenter Y ] + (arcRadius * [ endU Y ]) ; 
		dz = [ arcCenter Z ] + (arcRadius * [ endU Z ]) ; 
		
		endPosition = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
		
		return ;
	}
		
- (void) restoreVectorsAndAngle
	{
		// Compute arc length (assume limit plane not in effect), vertex positions, and uPerp direction
		// This is needed after arc subdivision!
		
		MMVector3 *crossProduct = [ [ MMVector3 alloc ] initByCrossing:startU and:endU ] ;
		
		double dot ;
		
		dot = [ startU dotWith:endU ] ;
		
		if( fabs(dot) > 1. ) dot = dot/fabs(dot) ;
		
		if( [ crossProduct dotWith:arcAxis ] > 0. )
			{
				angle = acos( dot ) ;
			}
		else
			{
				angle = 2.*acos(-1.) - acos( dot ) ;
			}
			
		length = arcRadius * angle ;
		
		[ crossProduct release ] ;
		
		[ self computeVertexPositions ] ;
		
		if( uPerp ) [ uPerp release ] ;
		
		uPerp = [ [ MMVector3 alloc ] initByCrossing:arcAxis and:startU ] ;
		
		[ uPerp normalize ] ;
		
		// Tangents
		
		if( startTangent ) [ startTangent release ] ;
		if( endTangent ) [ endTangent release ] ;
		
		startTangent = [ [ MMVector3 alloc ] initByCrossing:arcAxis and:startU ] ;
		[ startTangent normalize ] ;

		endTangent = [ [ MMVector3 alloc ] initByCrossing:arcAxis and:endU ] ;
		[ endTangent normalize ] ;
		

		return ;
	}
	
- (int) computeThetaAndPhiAnglesUsingMolecule:(SMMol *)mol
	{
		// This function will compute phi and theta start and end angles for a saddle-associated arc
		
		// Theta measures the angle along the I->J axis, using the probe center as reference. The vector reference for the angle is a 
		// unit vector directed from probe center toward atom I.
		
		// Phi measures the current position of the center of the probe. The moving vector is from the torus base point to the current probe, the 
		// reference vector is directed from the torus bp to the right reference probe position. moving X ref is directed along the torus axis, then
		// the angle is < 180 degrees, else the angle lies between 180 and 360. 
		
		// Don't try if no associated torus section
		
		// NOTE that in the current implementation this method is ONLY called to compute theta and phi angles for initial saddle edges. 
		// This information will be used in detecting a vector reversal that may be required in presence of self-intersecting surface
		
		if( torusSection == nil ) return 0 ;
		
		// Compute current phi angle
		
		// NOTICE: Phi and theta angles only computed if not already set!!
		
		double dx, dy, dz, dot, dotRefPerp, saddleRadius ;
		int failure ;
		
		
		// Dot in the direction of axis
		
		MMVector3 *diff, *axis, *base, *perp, *refR, *refL, *refPerp ;
		
		failure = 0 ;
		
		axis = [ torusSection axis ] ;
		
		base = [ torusSection base ] ;
		
		refR = [ torusSection refR ] ;
		refL = [ torusSection refL ] ;
		refPerp = [ torusSection refPerp ] ;
		
		saddleRadius = (double) [ torusSection saddleRadius ] ;
		
		diff = [ [ MMVector3 alloc ] initX:0. Y:0. Z:0. ] ;
		
		dx = [ startPosition X ] - [ [ torusSection base ] X ] ;
		dy = [ startPosition Y ] - [ [ torusSection base ] Y ] ;
		dz = [ startPosition Z ] - [ [ torusSection base ] Z ] ;
		
		[ diff setX:dx ] ;
		[ diff setY:dy ] ;
		[ diff setZ:dz ] ;
		
		[ diff normalize ] ;
		
		
		dot = [ axis dotWith:diff ] ;
		
		// Perpendicular vector
		
		dx = [ diff X ] - dot * [ axis X ] ;
		dy = [ diff Y ] - dot * [ axis Y ] ;
		dz = [ diff Z ] - dot * [ axis Z ] ;
		
		perp = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
		
		if( [ perp length ] < 1e-6 ) 
			{
				failure = 1 ;
				goto SKIP_START ;
			}
			
		
		// There is a danger here that perp could be too short (this because of a probe that almost falls through the I-J gap), but
		// I will asssume for now that this is precluded by narrow-torus skipping
		
		// NOTE - the perp vector can also go to zero owing to presence of self-intersecting surface (see test above)
		
		[ perp normalize ] ;
		
		// See if perp needs to be reversed
		
		if( arcType == REENTRANTR )
			{
				// Check if vector in direction of host probe
				
				if( [ refR dotWith:perp ] < 0. )
					{
						[ perp reverse ] ;
					}
			}
		else if( arcType == REENTRANTL )
			{
				// Check if vector in direction of host probe
				
				if( [ refL dotWith:perp ] < 0. )
					{
						[ perp reverse ] ;
					}
			}
		else if( phiStart != 0. && phiStart != [ torusSection phiAngleSubtended ] && [ torusSection withinTorusPhiLimits:perp ] == NO )
			{
				[ perp reverse ] ;
			}
		
		dot = [ perp dotWith:refR ] ;
		dotRefPerp = [ perp dotWith:refPerp ] ;
		
		if( fabs(dot) > 1. ) dot = dot/fabs(dot) ;
		if( fabs(dotRefPerp) > 1. ) dotRefPerp = dotRefPerp/fabs(dotRefPerp) ;
		
		//testCross = [ [ MMVector3 alloc ] initByCrossing:perp and:refR ] ;
		
		/*
		if( [ testCross dotWith:axis ] >= 0. )
			{
				phiStart = acos( dot ) ;
			}
		else
			{
				phiStart = 2.*acos(-1.) - acos(dot) ;
			}
			
		// TEMPORARY HACK - don't let angle get too close to 2 PI
		
		if( 2.*acos(-1.) - phiStart < 0.0001 ) phiStart = 0. ;
		*/
		
		if( phiStart < 0. )
			{
				phiStart = acos( dot ) ;
				
				if( dotRefPerp < 0. )
					{
						phiStart = 2.*acos(-1.) - phiStart ;
					}
			}
	
		
		double probeX, probeY, probeZ ;
		MMVector3 *thetaRef, *pointRef ;
		
		thetaRef = [ [ MMVector3 alloc ] initX:0. Y:0. Z:0. ] ;
		pointRef = [ [ MMVector3 alloc ] initX:0. Y:0. Z:0. ] ;
		
		if( thetaStart < 0. )
			{
				probeX = [ base X ] + saddleRadius*[ perp X ] ;
				probeY = [ base Y ] + saddleRadius*[ perp Y ] ;
				probeZ = [ base Z ] + saddleRadius*[ perp Z ] ;
				
				dx = ((double *)mol->xAtom)[ [ torusSection atomI ] ] - probeX ;
				dy = ((double *)mol->yAtom)[ [ torusSection atomI ] ] - probeY ;
				dz = ((double *)mol->zAtom)[ [ torusSection atomI ] ] - probeZ ;
				
				[ thetaRef setX:dx ] ;
				[ thetaRef setY:dy ] ;
				[ thetaRef setZ:dz ] ;
				
				[ thetaRef normalize ] ;
				
				
				dx = [ startPosition X ] - probeX ;
				dy = [ startPosition Y ] - probeY ;
				dz = [ startPosition Z ] - probeZ ;
				
				[ pointRef setX:dx ] ;
				[ pointRef setY:dy ] ;
				[ pointRef setZ:dz ] ;
				
				[ pointRef normalize ] ;
				
				dot = [ pointRef dotWith:thetaRef ] ;
				
				if( fabs(dot) > 1. ) dot = dot/fabs(dot) ;
				
				thetaStart = acos( dot ) ;
		}
		
SKIP_START:
		
		// Now do the whole thing over for the end point - should have done this in the torus object?
		
		dx = [ endPosition X ] - [ [ torusSection base ] X ] ;
		dy = [ endPosition Y ] - [ [ torusSection base ] Y ] ;
		dz = [ endPosition Z ] - [ [ torusSection base ] Z ] ;
		
		[ diff setX:dx ] ;
		[ diff setY:dy ] ;
		[ diff setZ:dz ] ;
		
		[ diff normalize ] ;
		
		
		dot = [ axis dotWith:diff ] ;
		
		// Perpendicular vector
		
		dx = [ diff X ] - dot * [ axis X ] ;
		dy = [ diff Y ] - dot * [ axis Y ] ;
		dz = [ diff Z ] - dot * [ axis Z ] ;
		
		[ perp setX:dx ] ;
		[ perp setY:dy ] ;
		[ perp setZ:dz ] ;
		
		if( [ perp length ] < 1e-6 ) 
			{
				failure = 2 ;
				goto SKIP_END ;
			}

		
		[ perp normalize ] ;
		
		if( arcType == REENTRANTR )
			{
				// Check if vector in direction of host probe
				
				if( [ refR dotWith:perp ] < 0. )
					{
						[ perp reverse ] ;
					}
			}
		else if( arcType == REENTRANTL )
			{
				// Check if vector in direction of host probe
				
				if( [ refL dotWith:perp ] < 0. )
					{
						[ perp reverse ] ;
					}
			}
		else if( phiEnd != 0. && phiEnd != [ torusSection phiAngleSubtended ] && [ torusSection withinTorusPhiLimits:perp ] == NO )
			{
				[ perp reverse ] ;
			}
		
		dot = [ perp dotWith:refR ] ;
		dotRefPerp = [ perp dotWith:refPerp ] ;
		
		if( fabs(dot) > 1. ) dot = dot/fabs(dot) ;
		if( fabs(dotRefPerp) > 1. ) dotRefPerp = dotRefPerp/fabs(dotRefPerp) ;
		
		/*	
		[ testCross release ] ;
		
		testCross = [ [ MMVector3 alloc ] initByCrossing:perp and:refR ] ;
		
		if( [ testCross dotWith:axis ] >= 0. )
			{
				phiEnd = acos( dot ) ;
			}
		else
			{
				phiEnd = 2.*acos(-1.) - acos(dot) ;
			}
			
		// TEMPORARY HACK - don't let angle get too close to 2 PI
		
		if( 2.*acos(-1.) - phiEnd < 0.0001 ) phiEnd = 0. ;
		
		*/
		
		if( phiEnd < 0. )
			{
				phiEnd = acos( dot ) ;
				
				if( dotRefPerp < 0. )
					{
						phiEnd = 2.*acos(-1.) - phiEnd ;
					}
			}
		
		if( thetaEnd < 0. )
			{
				probeX = [ base X ] + saddleRadius*[ perp X ] ;
				probeY = [ base Y ] + saddleRadius*[ perp Y ] ;
				probeZ = [ base Z ] + saddleRadius*[ perp Z ] ;
				
				dx = mol->xAtom[ [ torusSection atomI ] ] - probeX ;
				dy = mol->yAtom[ [ torusSection atomI ] ] - probeY ;
				dz = mol->zAtom[ [ torusSection atomI ] ] - probeZ ;
				
				[ thetaRef setX:dx ] ;
				[ thetaRef setY:dy ] ;
				[ thetaRef setZ:dz ] ;
				
				[ thetaRef normalize ] ;
				
				
				dx = [ endPosition X ] - probeX ;
				dy = [ endPosition Y ] - probeY ;
				dz = [ endPosition Z ] - probeZ ;
				
				[ pointRef setX:dx ] ;
				[ pointRef setY:dy ] ;
				[ pointRef setZ:dz ] ;
				
				[ pointRef normalize ] ;
				
				dot = [ pointRef dotWith:thetaRef ] ;
				
				if( fabs(dot) > 1. ) dot = dot/fabs(dot) ;
				
				thetaEnd = acos( dot ) ;
			}
		
		// Sanity check - start and stop positions for arc should be regenerated from phi and theta
		
		
		// Start
		/*
		pos = [ self computePositionForTheta:thetaStart andPhi:phiStart usingMolecule:mol allowSelfIntersection:mol->allowSelfIntersection normal:nil ] ;
		
		[ diff setX:([ pos X] - [ startPosition X ]) ] ;
		[ diff setY:([ pos Y] - [ startPosition Y ]) ] ;
		[ diff setZ:([ pos Z] - [ startPosition Z ]) ] ;
		
		
		// The following always occurs when self-int surface is present - I will remove for now
		//if( [ diff length ] > 0.0001 )
			//{
			//	printf( "WARNING: START POSITION OF ARC NOT REGENERATED FROM THETA,PHI! DIFF = %f\n", [ diff length ] ) ;
			//}
			
		[ pos release ] ;
		
		pos = [ self computePositionForTheta:thetaEnd andPhi:phiEnd usingMolecule:mol allowSelfIntersection:mol->allowSelfIntersection  normal:nil ] ;
		
		[ diff setX:([ pos X] - [ endPosition X ]) ] ;
		[ diff setY:([ pos Y] - [ endPosition Y ]) ] ;
		[ diff setZ:([ pos Z] - [ endPosition Z ]) ] ;
		
		//if( [ diff length ] > 0.0001 )
			//{
			//	printf( "WARNING: END POSITION OF ARC NOT REGENERATED FROM THETA,PHI! DIFF = %f\n", [ diff length ] ) ;
			//}
			
		[ pos release ] ;
		*/
		
		
		[ pointRef release ] ;
		[ thetaRef release ] ;
		[ diff release ] ;
		[ perp release ] ;
		//[ testCross release ] ;
		
SKIP_END:
		return failure ;
	
	}
		
- (void) arcPoint:(MMVector3 *)p atFraction:(double)f usingMolecule:(SMMol *)m
	{
		// This function computes the position of an interior point of an arc, fraction f (0 <= f <= 1) from the start
		// It will take into account twinning
		
		static MMVector3 *twinPoint = nil ;
		
		if( ! twinPoint )
			{
				twinPoint = [ [ MMVector3 alloc ] initX:0. Y:0. Z:0. ] ;
			}
		
		if( twin )
			{
				[ self theArcPoint:p atFraction:f usingMolecule:m ] ;
				[ twin theArcPoint:twinPoint atFraction:(1. - f) usingMolecule:m ] ;
				
				[ p setX:( [ p X ] + [ twinPoint X ] )/2. ] ;
				[ p setY:( [ p Y ] + [ twinPoint Y ] )/2. ] ;
				[ p setZ:( [ p Z ] + [ twinPoint Z ] )/2. ] ;
				
			}
		else
			{
				[ self theArcPoint:p atFraction:f usingMolecule:m ] ;
			}
			
		return ;
	}
		
- (void) theArcPoint:(MMVector3 *)p atFraction:(double) f usingMolecule:(SMMol *)m
	{
		
		double deltaAngle ;
		
		double cX, cY, cZ, dx, dy, dz ;
		
		
		
		MMVector3 *nextPos ;
		
		
		double phi1, theta1 ;
		
		// Note - use saddle-type subdivision on reentrant arcs (reentrantR and reentrantL) to handle self-intersecting surface
		
		switch( arcType )
			{
				case CONTACTI: // Contact I
				case CONTACTJ: // Contact J
				case INTERIORGEODESIC: // Interior geodesic
		
					deltaAngle = angle * f ;
									
					dx = cos(deltaAngle)*[ startU X ] + sin(deltaAngle)*[ uPerp X ] ;
					dy = cos(deltaAngle)*[ startU Y ] + sin(deltaAngle)*[ uPerp Y ] ;
					dz = cos(deltaAngle)*[ startU Z ] + sin(deltaAngle)*[ uPerp Z ] ;
					
					cX = [ arcCenter X ] ;
					cY = [ arcCenter Y ] ;
					cZ = [ arcCenter Z ] ;
					
					
					
					[ p setX:( cX + arcRadius * dx ) ] ;
					[ p setY:( cY + arcRadius * dy ) ] ;
					[ p setZ:( cZ + arcRadius * dz ) ] ;
					
					// Check if we have reentrant arc, limit plane
					
					SMCycle *parentCycle = [ parentCycles lastObject ] ;
					SMLimitPlane *limitPlane ;
					
					if( (limitPlane = parentCycle->theLimitPlane ) && parentCycle->selfIntersection == YES )
						{
							[ limitPlane adjustPosition:p andNormal:nil ] ;
						}
					
					break ;
					
				case REENTRANTR: // Reentrant R
				case REENTRANTL: // Reentrant L
				case INTERIORSADDLE: // Interior saddle
				
			
					phi1 = phiStart + (phiEnd - phiStart)*f ;
					theta1 = thetaStart + (thetaEnd - thetaStart)*f ;
					
					nextPos = [ self computePositionForTheta:theta1 andPhi:phi1 usingMolecule:m 
						allowSelfIntersection:m->allowSelfIntersection normal:nil ] ;
						
					// this is clunky
					
					[ p setX:[ nextPos X ] ] ;
					[ p setY:[ nextPos Y ] ] ;
					[ p setZ:[ nextPos Z ] ] ;
					
					[ nextPos release ] ;
					
					break ;
				
				default:
					printf( "UNRECOGNIZED ARC TYPE IN METHOD arcPoint - Exit!\n" ) ;
					exit(1) ;
			
				
						
			}
			
		return ;
	}
						
		
		
- (MMVector3 *) computePositionForTheta:(double)t andPhi:(double)p usingMolecule:(SMMol *)mol allowSelfIntersection:(BOOL)SIFlag normal:(MMVector3 *)norm
	{
		// Inverse function - find position given phi and theta
		
		if( torusSection == nil )
			{
				return nil ;
			}
			
		// First place probe
		
		double dx, dy, dz, dispX, dispY, dispZ ;
		
		double saddleRadius ;
		double probePosX, probePosY, probePosZ ;
		double perpRefX, perpRefY, perpRefZ ;
		double thetaRefX, thetaRefY, thetaRefZ ;
		double dispSaveX, dispSaveY, dispSaveZ ;
		
		saddleRadius = [ torusSection saddleRadius ] ;
		
		dispSaveX = saddleRadius * ( cos(p) * [ [ torusSection refR ]  X ] + sin(p) * [ [ torusSection refPerp ] X ] ) ;
		dispSaveY = saddleRadius * ( cos(p) * [ [ torusSection refR ]  Y ] + sin(p) * [ [ torusSection refPerp ] Y ] ) ;
		dispSaveZ = saddleRadius * ( cos(p) * [ [ torusSection refR ]  Z ] + sin(p) * [ [ torusSection refPerp ] Z ] ) ;
		
		//dispVector = [ [ MMVector3 alloc ] initX:dispX Y:dispY Z:dispZ ] ;
	
		probePosX = [ [ torusSection base ] X ] + dispSaveX ;
		probePosY = [ [ torusSection base ] Y ] + dispSaveY ;
		probePosZ = [ [ torusSection base ] Z ] + dispSaveZ ;
			
		//probePos = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
		
		// Theta reference vector
		
			
		thetaRefX = mol->xAtom[ [ torusSection atomI ] ] - probePosX ;
		thetaRefY = mol->yAtom[ [ torusSection atomI ] ] - probePosY ;
		thetaRefZ = mol->zAtom[ [ torusSection atomI ] ] - probePosZ ;
		
		//thetaRef = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
		
		//[ thetaRef normalize ] ;
		
		double size = sqrt( thetaRefX*thetaRefX + thetaRefY*thetaRefY + thetaRefZ*thetaRefZ ) ;
		
		thetaRefX /= size ;
		thetaRefY /= size ;
		thetaRefZ /= size ;
			
		dx = mol->xAtom[ [ torusSection atomJ ] ] - probePosX ;
		dy = mol->yAtom[ [ torusSection atomJ ] ] - probePosY  ;
		dz = mol->zAtom[ [ torusSection atomJ ] ] - probePosZ  ;
		
		//tempVec = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
		
		double dot ;
		
		dot = dx*thetaRefX + dy*thetaRefY + dz*thetaRefZ ;
		
		perpRefX = dx - dot*thetaRefX ;
		perpRefY = dy - dot*thetaRefY ;
		perpRefZ = dz - dot*thetaRefZ ;
		
		size = sqrt( perpRefX*perpRefX + perpRefY*perpRefY + perpRefZ*perpRefZ ) ;
		
		perpRefX  /= size ;
		perpRefY  /= size ;
		perpRefZ  /= size ;
		
		//perpRef = [ [ MMVector3 alloc ] initAlong:tempVec perpTo:thetaRef ] ;
		
		//[ perpRef normalize ] ;
		
		//[ tempVec release ] ;
		
		double probeRadius = mol->probeRadius ;
		
		BOOL deformed = NO ;
		
		if( mol->allowSelfIntersection == YES || torusSection->selfIntersection == NO )
			{
				dx = probePosX + (  probeRadius*( cos(t)*thetaRefX + sin(t)*perpRefX ) ) ;
				dy = probePosY + (  probeRadius*( cos(t)*thetaRefY + sin(t)*perpRefY ) ) ;
				dz = probePosZ + (  probeRadius*( cos(t)*thetaRefZ + sin(t)*perpRefZ ) ) ;
			}
		else
			{
				dispSaveX = -dispSaveX ;
				dispSaveY = -dispSaveY ;
				dispSaveZ = -dispSaveZ ;
				
				//[ dispVector reverse ] ;
				
				//[ dispVector normalize ] ;
				
				size = sqrt( dispSaveX*dispSaveX + dispSaveY*dispSaveY + dispSaveZ*dispSaveZ ) ;
				
				dispSaveX /= size ;
				dispSaveY /= size ;
				dispSaveZ /= size ;
				
				dispX = probeRadius*( cos(t)*thetaRefX + sin(t)*perpRefX ) ;
				dispY = probeRadius*( cos(t)*thetaRefY + sin(t)*perpRefY ) ;
				dispZ = probeRadius*( cos(t)*thetaRefZ + sin(t)*perpRefZ ) ;
				
				dot = dispX * dispSaveX + dispY * dispSaveY + dispZ * dispSaveZ ;
				
				double h = torusSection->saddleRadius - dot ;
				
				double factor ;
				
				if( h < torusSection->bufferIJ )
					{
						factor = (saddleRadius - torusSection->bufferIJ)/dot ;
						deformed = YES ;
					}
				else
					{
						factor = 1. ;
					}
					
				dx = probePosX + factor*dispX ;
				dy = probePosY + factor*dispY ;
				dz = probePosZ + factor*dispZ ;
			}
				
				
		
		MMVector3 *returnVec = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
		
		if( norm )
			{
				// This is wasteful 
				
				dx = dx - probePosX ;
				dy = dy - probePosY ;
				dz = dz - probePosZ ;
				
				if( deformed == YES )
					{
						// Make orthogonal to torus axis
						
						MMVector3 *a = [ torusSection axis ] ;
						double ax, ay, az ;
						
						ax = [ a X ] ; ay = [ a Y ] ; az = [ a Z ] ;
						
						dot = dx*ax + dy*ay + dz*az ;
						
						dx = dx - dot*ax ;
						dy = dy - dot*ay ;
						dz = dz - dot*az ;
						
					}
					
				[ norm setX:(-dx) ] ;
				[ norm setY:(-dy) ] ;
				[ norm setZ:(-dz) ] ;
				
				[ norm normalize ] ;
			}
			
		
		//[ probePos release ] ;
		//[ perpRef release ] ;
		//[ thetaRef release ] ;
		//[ dispVector release ] ;
		
		return returnVec ;
	}
		
			
- (MMVector3 *) probePositionForPhi:(double) p
	{
		// Compute normal using probe position
		
		if( torusSection == nil )
			{
				return nil ;
			}
			
		// First place probe
		
		double dx, dy, dz, dispX, dispY, dispZ ;
		
		MMVector3 *probePos ;
		
		dispX = [ torusSection saddleRadius ] * ( cos(p) * [ [ torusSection refR ]  X ] + sin(p) * [ [ torusSection refPerp ] X ] ) ;
		dispY = [ torusSection saddleRadius ] * ( cos(p) * [ [ torusSection refR ]  Y ] + sin(p) * [ [ torusSection refPerp ] Y ] ) ;
		dispZ = [ torusSection saddleRadius ] * ( cos(p) * [ [ torusSection refR ]  Z ] + sin(p) * [ [ torusSection refPerp ] Z ] ) ;
	
		dx = [ [ torusSection base ] X ] + dispX ;
		dy = [ [ torusSection base ] Y ] + dispY ;
		dz = [ [ torusSection base ] Z ] + dispZ ;
			
		probePos = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
		
		return probePos ;
	}
		
		
- (BOOL) intersectWith:(SMArc *)t
	{
		// Test for intersection with other general circular arc
		
		// Strategy: Use startU and uPerp of the home arc as a frame; intersect points are 
		// expressed as alpha*uStart + beta*uPerp,  with alpha^2 + beta^2 = arcRadius^2
		//
		// Let a2 be the axis of the test ("with") arc, c1 the center of the home arc and c2 the center of the 
		// test arc. Then for intersection points
		//
		//	(alpha*uStart + beta*uPerp + c1 - c2).a2 = 0 
		//
		// alpha*(uStart.a2) + beta*(uPerp.a2)  = (c2 - c1).a2 OR
		//
		// alpha*A + beta*B = C
		//
		// Next,  beta = +/- sqrt( arcRadius^2 - alpha^2), so
		//
		//	alpha*A +/- sqrt(arcRadius^2 - alpha^2)*B = C
		//
		// Rearranging and squaring, 
		//
		//	(arcRadius^2 - alpha^2)*B^2 = (C - alpha*A)^2 = C^2 - 2*A*C*alpha + A^2*alpha^2
		//
		//	alpha^2*(A^2 + B^2) - 2*A*C*alpha + (C^2 - B^2*arcRadius^2) = 0
		// 
		// Solve for alpha - 
		//						alpha = (2*A*C +/- Sqrt( 4*A^2*C^2 - 4*(A^2 + B^2)*(C^2 - B^2*arcRadius^2) ) ) / (2*(A^2 + B^2))
		//
		// For each alpha, have +/- beta. Need to compute location of intersection point and confirm it is valid (lies on the home arc), and then test if it
		// lies inside the target arc
		
		MMVector3 *axis2, *centDiff ;
		
		axis2 = [ t arcAxis ] ;
		
		centDiff = [ [ MMVector3 alloc ] initX:( [ [ t arcCenter ] X ] - [ arcCenter X ] )
			Y:( [ [ t arcCenter ] Y ] - [ arcCenter Y ] )
			Z:( [ [ t arcCenter ] Z ] - [ arcCenter Z ] ) ] ;
		
		double A, B, C ;
		
		A = [ startU dotWith:axis2 ] ;
		B = [ uPerp dotWith:axis2 ] ;
		C = [ centDiff dotWith:axis2 ] ;
		
		[ centDiff release ] ;
		
		double alpha[2], beta[2], disc ;
		
		disc = 4.*( A*A*C*C - (A*A + B*B)*(C*C - B*B*arcRadius*arcRadius) ) ;
		
		if( disc < 0. ) return NO ;
		
		disc = sqrt( disc ) ;
		
		alpha[0] = (2.*A*C + disc)/(2.*(A*A + B*B)) ; 
		alpha[1] = (2.*A*C - disc)/(2.*(A*A + B*B)) ; 
		
		int i ;
		MMVector3 *testPos ;
		
		testPos =  [ [ MMVector3 alloc ] initX:0. Y:0. Z:0. ] ;
		
		[ testPos autorelease ] ;
		
		for( i = 0 ; i < 2 ; ++i )
			{
				double argu ;
				
				argu = arcRadius*arcRadius - alpha[i]*alpha[i] ;
				
				if( argu < 0. ) argu = 0. ;
				
				beta[0] = sqrt( argu ) ;
				beta[1] = -beta[0] ;
				
				int j ;
				
				for( j = 0 ; j < 2 ; ++j )
					{
						// Check for legal alpha and beta
						
						
						double testRad ;
						
						testRad = alpha[i]*alpha[i] + beta[j]*beta[j] ;
						
						if( fabs( testRad - arcRadius*arcRadius ) > 1.e-06 )
							{
								continue ;
							}
							
						[ testPos setX:( [ arcCenter X ] + alpha[i]*[ startU X ] + beta[j]*[ uPerp X ] ) ] ;
						[ testPos setY:( [ arcCenter Y ] + alpha[i]*[ startU Y ] + beta[j]*[ uPerp Y ] ) ] ;
						[ testPos setZ:( [ arcCenter Z ] + alpha[i]*[ startU Z ] + beta[j]*[ uPerp Z ] ) ] ;
						
						MMVector3 *diff ;
						double testDot ;
						
						diff = [ [ MMVector3 alloc ] initX:([ testPos X ] - [ arcCenter X ])
							Y:([ testPos Y ] - [ arcCenter Y ])
							Z:([ testPos Z ] - [ arcCenter Z ]) ]  ;
							
						testDot = [ diff dotWith:arcAxis ] ;
						
						// I will comment this error out for now... It can occur easily in very tight spots, but is generally innocuous
						if( fabs(testDot) > 1e-06 )
							{
								printf( "WARNING: Intersection test error!\n" ) ;
							}
							
						[ diff setX:([ testPos X ] - [ [ t arcCenter ] X ]) ] ;
						[ diff setY:([ testPos Y ] - [ [ t arcCenter ] Y ]) ] ;
						[ diff setZ:([ testPos Z ] - [ [ t arcCenter ] Z ]) ] ;
						
						testDot = [ diff dotWith:[ t arcAxis ] ] ;
						
						[ diff release ] ;
						
						if( fabs(testDot) > 1e-06 )
							{
								// Invalid alpha/beta pair
								//printf( "WARNING: Intersection test error!\n" ) ;
								continue ;
							}
							
						if( [ self withinArc:testPos ] == YES && [ t withinArc:testPos ] == YES ) return YES ;
					}
			}
			
		return NO ;
	}
	
- (BOOL) intersectWith2:(SMArc *)arc2
	{
		// This is modelled after the test in 3smart
		
		// In what follows, arc1 is the "self" arc, arc2 is the target
		
         double dotab, dotsa, dotva ;
         double alpha1, beta1, alpha2, beta2, alphamax, betamax ;
         double a, b, c, disc, arg, ang, a1, a2, a3, a4 ;
         MMVector3 *v ;
         short in11, in12, in21, in22 ;
		 static MMVector3 *w = nil ;
		 double px, py, pz, tx, ty, tz, sx, sy, sz, wx, wy, wz, vx, vy, vz ;
		 double axis2x, axis2y, axis2z ;
		 double sUx, sUy, sUz, eUx, eUy, eUz ;

		 //v = [ [ MMVector3 alloc ] initByCrossing:arcAxis and:startU ] ;
         // cross( arc1->axis, arc1->ustart, & v ) ;
		 
		 //[ v normalize ] ;
         //normalize( & v ) ;
		 
		 v = uPerp ;
		 
		 vx = [ uPerp X ] ;
		 vy = [ uPerp Y ] ;
		 vz = [ uPerp Z ] ;
		 
		 sUx = [ startU X ] ;
		 sUy = [ startU Y ] ;
		 sUz = [ startU Z ] ;
		 
		 eUx = [ endU X ] ;
		 eUy = [ endU Y ] ;
		 eUz = [ endU Z ] ;

    /*     printf( "v = %f %f %f\n", v.x, v.y, v.z ) ;  */
	
		/*
		if( ! w )
			{
				w = [ [ MMVector3 alloc ] initX:([ [ arc2 arcCenter ] X ] - [ arcCenter X ] )
					Y:([ [ arc2 arcCenter ] Y ] - [ arcCenter Y ] )
					Z:([ [ arc2 arcCenter ] Z ] - [ arcCenter Z ] ) ] ;
			}
		else
			{
				[ w setX:([ [ arc2 arcCenter ] X ] - [ arcCenter X ] ) ] ;
				[ w setY:([ [ arc2 arcCenter ] Y ] - [ arcCenter Y ] ) ] ;
				[ w setZ:([ [ arc2 arcCenter ] Z ] - [ arcCenter Z ] ) ] ;
			}
		*/
			
		wx = [ [ arc2 arcCenter ] X ] - [ arcCenter X ] ;
		wy = [ [ arc2 arcCenter ] Y ] - [ arcCenter Y ] ;
		wz = [ [ arc2 arcCenter ] Z ] - [ arcCenter Z ] ;
		
		axis2x = [ [ arc2 arcAxis ] X ] ;
		axis2y = [ [ arc2 arcAxis ] Y ] ;
		axis2z = [ [ arc2 arcAxis ] Z ] ;

         //w.x = arc2->base.x - arc1->base.x ;
         //w.y = arc2->base.y - arc1->base.y ;
         //w.z = arc2->base.z - arc1->base.z ;

    /*     printf( "w = %f %f %f\n", w.x, w.y, w.z ) ;  */
		 
		 dotab = axis2x*wx + axis2y*wy + axis2z*wz ;
		 //dotab = [ [ arc2 arcAxis ] dotWith:w ] ;
         //dotab = dot(arc2->axis, w ) ;
		 
		 dotsa =  ( axis2x*sUx + axis2y*sUy + axis2z*sUz ) * arcRadius ;
		 //dotsa = [ [ arc2 arcAxis ] dotWith:startU ] * arcRadius ;
         //dotsa = dot(arc2->axis, arc1->ustart) * 
         //            arc1->radius ;
		 
		 dotva =  ( axis2x*vx + axis2y*vy + axis2z*vz ) * arcRadius ;
		 
	     //dotva = [ [ arc2 arcAxis ] dotWith:v ] * arcRadius ;
         //dotva = dot(arc2->axis, v) *
         //            arc1->radius ;

     /*    printf( "dotab = %f, dotsa = %f, dotva = %f\n",
           dotab, dotsa, dotva ) ;  */
		   
		 alphamax = eUx*sUx + eUy*sUy + eUz*sUz ;
		 //alphamax = [ endU dotWith:startU ] ;
		 //betamax = [ endU dotWith:v ] ;
		 betamax = eUx*vx + eUy*vy + eUz*vz ;
		 
         //alphamax = dot( arc1->uend, arc1->ustart ) ;
         //betamax = dot( arc1->uend, v ) ;

      /*
         printf( "alphamax = %f   betamax = %f\n", alphamax, betamax ) ;

         printf( "dot(w,w) = %f\n", dot(w,w) ) ; 
      */

         if( (dotsa == 0.) && (dotva == 0.) )
            {
               /* Arcs are parallel */

               //if( dot(w,w) < 1.e-12 )
			   if( [ w length ] < 1.e-6 )
                  {
                     alpha1 = 1. ;
                     beta1 = 0. ;
 
                     alpha2 = alphamax ;
                     beta2 = betamax ;
                  }
               else
                  {
                     return NO ;
                  }
            }
         else
            {
               if( fabs(dotva) > fabs(dotsa) )
                  {
                     /* Solve in terms of alpha */

                     a = 1. + ((dotsa/dotva)*(dotsa/dotva)) ;

                     b = -2.*(dotab/dotva)*(dotsa/dotva) ;

                     c = ((dotab/dotva)*(dotab/dotva)) - 1. ;
                  }
               else 
                  {
                     /* Solve in terms of beta */

                     a = 1. + ((dotva/dotsa)*(dotva/dotsa)) ;

                     b = -2.*(dotab/dotsa)*(dotva/dotsa) ;

                     c = ((dotab/dotsa)*(dotab/dotsa)) - 1. ;
                  }

               disc = (b*b) - (4.*a*c) ;

               if( disc < 0. ) return NO ;

               if( fabs(dotva) > fabs(dotsa) )
                  {
                     alpha1 = (-b + sqrt(disc))/(2.*a) ;
                     alpha2 = (-b - sqrt(disc))/(2.*a) ;
 
                     beta1 = (dotab/dotva) - ((dotsa/dotva)*alpha1) ;
                     beta2 = (dotab/dotva) - ((dotsa/dotva)*alpha2) ;
                  }
               else
                  {
                     beta1 = (-b + sqrt(disc))/(2.*a) ;
                     beta2 = (-b - sqrt(disc))/(2.*a) ;

                     alpha1 = (dotab/dotsa) - ((dotva/dotsa)*beta1) ;
                     alpha2 = (dotab/dotsa) - ((dotva/dotsa)*beta2) ;
                  } 
            }

   /*      printf( "alpha,beta1 = %f %f   alpha,beta2 = %f %f\n",
            alpha1, beta1, alpha2, beta2 ) ; */

         /* Have angles for two points - see if either falls inside 
            both arcs */

         in11 = FALSE ;
         in21 = FALSE ;
         in12 = FALSE ;
         in22 = FALSE ;

         a1 = 0. ;
		 
		 /*
		 d = alphamax ;
		// d = [ startU dotWith:endU ] ;
         //d = dot( arc1->ustart, arc1->uend ) ;

         //if(fabs(d) > 1. ) d = d/fabs(d) ;
		 
		 if( fabs(d) > 1. )
			{
				if( d < 0. )
					{
						d = -1. ;
					}
				else
					{
						d = +1. ;
					}
			}

         a2 = acos( d ) ;
	
		
		 t = [ [ MMVector3 alloc ] initByCrossing:startU and:endU ] ;
         //cross( arc1->ustart, arc1->uend, & t ) ;

		 if( [ t dotWith:arcAxis ] < 0. ) a2 = 2.*acos(-1.) - a2 ;
         //if( dot( t, arc1->axis ) < 0. ) a2 = 2.*PI - a2 ;
		 */
		 
		 a2 = angle ;

         a3 = 0. ;
		 
		 /*
		 d = [ [ arc2 startU ] dotWith:[ arc2 endU ] ] ;
         //d = dot( arc2->ustart, arc2->uend ) ;
		 
         //if(fabs(d) > 1. ) d = d/fabs(d) ;
		 
		 if( fabs(d) > 1. )
			{
				if( d < 0. )
					{
						d = -1. ;
					}
				else
					{
						d = +1. ;
					}
			}
		 

         a4 = acos( d ) ;
		 
		 [ t release ] ;
		 t = [ [ MMVector3 alloc ] initByCrossing:[ arc2 startU ] and:[ arc2 endU ] ] ;
         //cross( arc2->ustart, arc2->uend, & t ) ;
		 
		 if( [ t dotWith:[ arc2 arcAxis ] ] < 0. ) a4 = 2.*acos(-1.) - a4 ;
         //if( dot( t, arc2->axis ) < 0. ) a4 = 2.*PI - a4 ;
		 */
		 
		 a4 = [ arc2 angle ] ;

    /*     printf( "a1 = %f  a2 = %f  a3 = %f  a4 = %f  \n",
          a1, a2, a3, a4 ) ;  */

         /* Test 1 & 2 inside 1  */ 
               
         if( (alpha1 == 1.) && (beta1 == 0.) ) 
            {
               in11 = TRUE ;
               //f11 = 0. ;
            }
         else
            {
               arg = alpha1/sqrt(alpha1*alpha1 + beta1*beta1) ;
      
               //if( fabs(arg) > 1. ) arg = arg/fabs(arg) ;
			   
				 if( fabs(arg) > 1. )
					{
						if( arg < 0. )
							{
								arg = -1. ;
							}
						else
							{
								arg = +1. ;
							}
					}

               if( beta1 >= 0. ) 
                  {
                     ang = acos( arg ) ;
                  }
               else
                  {
                     ang = 2.*acos(-1.) - acos( arg ) ;
                  }

           /*    printf( "Angle 1 in 1 = %f\n", ang ) ;  */

               if( (a1 <= ang) && (ang <= a2) )
                  {
                     in11 = TRUE ; 
					/*
                     if( (ang - a1) < (a2 - ang) )
                        {
                           fmin = ang - a1 ;
                        }
                     else
                        {
                           fmin = a2 - ang ;
                        }
                     
                     f11 = fmin/a2 ;
					*/
                  }
            }

         if( (alpha2 == alphamax) && (beta2 == betamax) )
            {
               in21 = TRUE ;
               //f21 = 0. ;
            }
         else
            {
               arg = alpha2/sqrt(alpha2*alpha2 + beta2*beta2) ;

               //if( fabs(arg) > 1. ) arg = arg/fabs(arg) ;
			   
				 if( fabs(arg) > 1. )
					{
						if( arg < 0. )
							{
								arg = -1. ;
							}
						else
							{
								arg = +1. ;
							}
					}
			   

               if( beta2 >= 0. )
                  {
                     ang = acos( arg ) ;
                  }
               else
                  {
                     ang = 2.*acos(-1.) - acos( arg ) ;
                  }

          /*     printf( "Angle 2 in 1 = %f\n", ang ) ; */

               if( (a1 <= ang) && (ang <= a2) )
                  {
                     in21 = TRUE ; 
					 
					/*
                     if( (ang - a1) < (a2 - ang) )
                        {
                           fmin = ang - a1 ;
                        }
                     else
                        {
                           fmin = a2 - ang ;
                        }
                     
                     f21 = fmin/a2 ;
				*/
                  }
            }      

         /* Test 1 & 2 inside 2 */ 
		 
		px = [ arcCenter X ] + 
					arcRadius*( alpha1*sUx + beta1*vx ) ;
		py = [ arcCenter Y ] + 
					arcRadius*( alpha1*sUy + beta1*vy ) ;
		pz = [ arcCenter Z ] + 
					arcRadius*( alpha1*sUz + beta1*vz ) ;
			
		/*
		if( ! p )
			{
				 p = [ [ MMVector3 alloc ] initX:( [ arcCenter X ] + 
					arcRadius*( alpha1*[ startU X ] + beta1*[ v X ] ) )
					Y:( [ arcCenter Y ] + 
					arcRadius*( alpha1*[ startU Y ] + beta1*[ v Y ] ) )
					Z:( [ arcCenter Z ] + 
					arcRadius*( alpha1*[ startU Z ] + beta1*[ v Z ] ) ) ] ;
			}
		else
			{
				[ p setX:( [ arcCenter X ] + 
					arcRadius*( alpha1*[ startU X ] + beta1*[ v X ] ) ) ] ;
					
				[ p setY:( [ arcCenter Y ] + 
					arcRadius*( alpha1*[ startU Y ] + beta1*[ v Y ] ) ) ] ;
					
				[ p setZ:( [ arcCenter Z ] + 
					arcRadius*( alpha1*[ startU Z ] + beta1*[ v Z ] ) ) ] ;
			}
		*/
				
		/*
         p.x = arc1->base.x + arc1->radius * (alpha1*arc1->ustart.x
                                                  +  beta1*v.x ) ; 
         p.y = arc1->base.y + arc1->radius * (alpha1*arc1->ustart.y
                                                  +  beta1*v.y ) ; 
         p.z = arc1->base.z + arc1->radius * (alpha1*arc1->ustart.z
                                                  +  beta1*v.z ) ; 
		*/

   /*      printf( "p1 = %f %f %f\n", p.x, p.y, p.z ) ; */
   
		tx = px - [ [ arc2 arcCenter ] X ] ;
		ty = py - [ [ arc2 arcCenter ] Y ] ;
		tz = pz - [ [ arc2 arcCenter ] Z ] ;
   
		//[ t setX:( [ p X ] - [ [ arc2 arcCenter ] X ] ) ] ;
		//[ t setY:( [ p Y ] - [ [ arc2 arcCenter ] Y ] ) ] ;
		//[ t setZ:( [ p Z ] - [ [ arc2 arcCenter ] Z ] ) ] ;

         //t.x = p.x - arc2->base.x ;
         //t.y = p.y - arc2->base.y ;
         //t.z = p.z - arc2->base.z ;

    /*     printf( "t = %f %f %f\n", t.x, t.y, t.z ) ; */
	
		 double size = sqrt( tx*tx + ty*ty + tz*tz ) ;
		 
		 tx /= size ;
		 ty /= size ;
		 tz /= size ;

		 
		 //[ t normalize ] ;
         //normalize( & t ) ; 
		 
		 double sU2x, sU2y, sU2z ;
		 
		 sU2x = [ [ arc2 startU ] X ] ;
		 sU2y = [ [ arc2 startU ] Y ] ;
		 sU2z = [ [ arc2 startU ] Z ] ;
		 
		 arg = tx*sU2x + ty*sU2y + tz*sU2z ;
		 //arg = [ t dotWith:[ arc2 startU ] ] ; 
         //arg = dot( t, arc2->ustart ) ;

        // if( fabs(arg) > 1. ) arg = arg/fabs(arg) ;
		 
		 if( fabs(arg) > 1. )
			{
				if( arg < 0. )
					{
						arg = -1. ;
					}
				else
					{
						arg = +1. ;
					}
			}
		 
		 
		 //s = [ [ MMVector3 alloc ] initByCrossing:[ arc2 startU ] and:t ] ;
         //cross( arc2->ustart, t, & s ) ;
		 
		 sx = sU2y*tz - ty*sU2z ;
		 sy = sU2z*tx - tz*sU2x ;
		 sz = sU2x*ty - tx*sU2y ;
		 
		 if( ( sx*axis2x + sy*axis2y + sz*axis2z ) > 0. )
		 //if( [ s dotWith:[ arc2 arcAxis ] ] > 0. )
         //if( dot(s, arc2->axis) > 0. )
            {
               ang = acos(arg) ;
            }
         else
            {
               ang = 2.*acos(-1.) - acos(arg) ;
            }

    /*     printf( "Angle 1 in 2 = %f\n", ang ) ; */

         if( (a3 <= ang) && (ang <= a4) )
            {
               in12 = TRUE ;
			
			/*
               if( (ang - a3) < (a4 - ang) )
                  {
                     fmin = ang - a3 ;
                  }
               else
                  {
                     fmin = a4 - ang ;
                  }
                   
               f12 = fmin/a4 ;
			*/
            }
			
		//[ p release ] ;
		
		// p is alloced at this point
		/*
		[ p setX:( [ arcCenter X ] + 
			arcRadius*( alpha2*[ startU X ] + beta2*[ v X ] ) ) ] ;
		[ p setY:( [ arcCenter Y ] + 
			arcRadius*( alpha2*[ startU Y ] + beta2*[ v Y ] ) ) ] ;
		[ p setZ:( [ arcCenter Z ] + 
			arcRadius*( alpha2*[ startU Z ] + beta2*[ v Z ] ) ) ] ;
		*/
		
		px = [ arcCenter X ] + 
					arcRadius*( alpha2*sUx + beta2*vx ) ;
		py = [ arcCenter Y ] + 
					arcRadius*( alpha2*sUy + beta2*vy ) ;
		pz = [ arcCenter Z ] + 
					arcRadius*( alpha2*sUz + beta2*vz ) ;
		
		 //p = [ [ MMVector3 alloc ] initX:( [ arcCenter X ] + 
		//	arcRadius*( alpha2*[ startU X ] + beta2*[ v X ] ) )
		//	Y:( [ arcCenter Y ] + 
		//	arcRadius*( alpha2*[ startU Y ] + beta2*[ v Y ] ) )
		//	Z:( [ arcCenter Z ] + 
		//	arcRadius*( alpha2*[ startU Z ] + beta2*[ v Z ] ) ) ] ;

		/*
         p.x = arc1->base.x + arc1->radius * (alpha2*arc1->ustart.x
                                            +  beta2*v.x ) ;
         p.y = arc1->base.y + arc1->radius * (alpha2*arc1->ustart.y
                                            +  beta2*v.y ) ;
         p.z = arc1->base.z + arc1->radius * (alpha2*arc1->ustart.z
                                            +  beta2*v.z ) ;
		*/

      /*   printf( "p2 = %f %f %f\n", p.x, p.y, p.z ) ;  */
	  
		tx = px - [ [ arc2 arcCenter ] X ] ;
		ty = py - [ [ arc2 arcCenter ] Y ] ;
		tz = pz - [ [ arc2 arcCenter ] Z ] ;
	  
		//[ t setX:( [ p X ] - [ [ arc2 arcCenter ] X ] ) ] ;
		//[ t setY:( [ p Y ] - [ [ arc2 arcCenter ] Y ] ) ] ;
		//[ t setZ:( [ p Z ] - [ [ arc2 arcCenter ] Z ] ) ] ;

         //t.x = p.x - arc2->base.x ;
         //t.y = p.y - arc2->base.y ;
         //t.z = p.z - arc2->base.z ;

      /*   printf( "t = %f %f %f\n", t.x, t.y, t.z ) ; */
	  
		 //[ t normalize ] ;
		 
		 size = sqrt( tx*tx + ty*ty + tz*tz ) ;
		 
		 tx /= size ;
		 ty /= size ;
		 tz /= size ;
		 
         //normalize( & t ) ;
		
		 arg = tx*sU2x + ty*sU2y + tz*sU2z ;
		 //arg = [ t dotWith:[ arc2 startU ] ] ;
         //arg = dot( t, arc2->ustart ) ;

         //if( fabs(arg) > 1. ) arg = arg/fabs(arg) ;
		 
		 if( fabs(arg) > 1. )
			{
				if( arg < 0. )
					{
						arg = -1. ;
					}
				else
					{
						arg = +1. ;
					}
			}
		 
		 
		 //[ s release ] ;
		 //s = [ [ MMVector3 alloc ] initByCrossing:[ arc2 startU ] and:t ] ;
         //cross( arc2->ustart, t, & s ) ;
		 
		 sx = sU2y*tz - ty*sU2z ;
		 sy = sU2z*tx - tz*sU2x ;
		 sz = sU2x*ty - tx*sU2y ;
		 
		 if( ( sx*axis2x + sy*axis2y + sz*axis2z ) > 0. )
		 //if( [ s dotWith:[ arc2 arcAxis ] ] > 0. )
         //if( dot(s, arc2->axis) > 0. )
            {
               ang = acos(arg) ;
            }
         else
            {
               ang = 2.*acos(-1.) - acos(arg) ;
            }

       /*  printf( "Angle 2 in 2 = %f\n", ang ) ; */

         if( (a3 <= ang) && (ang <= a4) )
            {
               in22 = TRUE ;
			
			/*
               if( (ang - a3) < (a4 - ang) )
                  {
                     fmin = ang - a3 ;
                  }
               else
                  {
                     fmin = a4 - ang ;
                  }
                   
               f22 = fmin/a4 ;
			*/
            }

     /*    printf( "in 11, 12, 21, 22 = %hd %hd %hd %hd\n",
           in11, in12, in21, in22 ) ; */
		   
		//[ p release ] ;
		//[ s release ] ;
		//[ t release ] ;

         if( in11 && in12 )
            {
               //if( mode && ((f11 < 0.01) || (f12 < 0.01)) )
               //   {
               //      return FALSE ;
               //   }
               //else
              //    {
                     return YES ;
               //   }
            }

         if( in21 && in22 )
            {
               //if( mode && ((f21 < 0.01) || (f22 < 0.01)) )
               //   {
               //      return FALSE ;
               //   }
              // else
               //   {
                     return YES ;
               //   }
            }

         return NO ;

      }
		
		
		
- (BOOL) withinArc:(MMVector3 *) pos
	{
		// Check if position (which is assumed to fall on the extended arc circle) is within the limits of the arc
		
		MMVector3 *disp, *crossTestS, *crossTestE ;
		
		disp = [ [ MMVector3 alloc ] initX:( [ pos X ] - [ arcCenter X ] )
					Y:( [ pos Y ] - [ arcCenter Y ] )
					Z:( [ pos Z ] - [ arcCenter Z ] ) ] ;
					
		[ disp normalize ] ;
					
		crossTestS = [ [ MMVector3 alloc ] initByCrossing:startU and:disp ] ;
		crossTestE = [ [ MMVector3 alloc ] initByCrossing:disp and:endU ] ;
		
		// TEST - if we are too close to an arc edge, assume intersection
		/*
		if( [ disp dotWith:startU ] > 0.999999 || [ disp dotWith:endU ] > 0.999999 )
			{
				[ disp release ] ;
				[ crossTestS release ] ;
				[ crossTestE release ] ;
				
				return YES ;
			}
		*/
		 
		double dotTestS, dotTestE ;
		
		dotTestS = [ crossTestS dotWith:arcAxis ] ;
		dotTestE = [ crossTestE dotWith:arcAxis ] ;
		
		[ disp release ] ;
		[ crossTestS release ] ;
		[ crossTestE release ] ;
			
		if( angle <=  acos(-1.) )
			{
				if( dotTestS >= 0. &&  dotTestE >= 0. )
					{
						return YES ;
					}
				else
					{
						return NO ;
					}
			}
		else
			{
				if( dotTestS <= 0. &&  dotTestE <= 0. )
					{
						return NO ;
					}
				else
					{
						return YES;
					}
			}
			
		return NO ;
	}
			
		
	
- (SMArc *) twin 
	{
		return twin ;
	}
	
- (void) setTwin:(SMArc *)t 
	{
		twin = t ;
		return ;
	}

- (double) phiStart 
	{
		return phiStart ;
	}
	
- (double) phiEnd 
	{
		return phiEnd ;
	}

- (double) thetaStart
	{
		return thetaStart ;
	}
	
- (double) thetaEnd 
	{
		return thetaEnd ;
	}

- (void) setPhiStart:(double)ps
	{
		phiStart = ps ;
		
		return ; 
	}
	
- (void) setPhiEnd:(double)pe 
	{
		phiEnd = pe ;
		
		return ;
	}

- (void) setThetaStart:(double)ts 
	{
		thetaStart = ts ; 
		
		return ;
	}
	
- (void) setThetaEnd:(double)te 
	{
		thetaEnd = te ;
		
		return ;
	}

- (NSMutableArray *) startConnections 
	{
		return startConnections ;
	}
	
- (NSMutableArray *) endConnections 
	{
		return endConnections ;
	}
	
	
- (NSMutableArray *) startConnectStart 
	{
		return startConnectStart ;
	}
	
- (NSMutableArray *) endConnectStart 
	{
		return endConnectStart ;
	}

- (void) addConnectionToArc:(SMArc *)a descend:(BOOL)desc
	{
		// This method connects the current arc to the argument. It also connects the argument to all arcs that the 
		// host is already connected to. The method is also called on the argument arc, using the host as argument
		
		// First, try to add connection to arcs already connected
		
		double distHostStartToArgumentStart, distHostEndToArgumentStart, distHostStartToArgumentEnd, distHostEndToArgumentEnd ;
		
		// Don't try to connect if already connected
		
		if( [ startConnections indexOfObject:a ] != NSNotFound ) return ;
		if( [ endConnections indexOfObject:a ] != NSNotFound ) return ;
		
		// Note that I am using a simple metric, not Euclidean distance
		
		distHostStartToArgumentStart =	fabs( [ startPosition X ] - [ [ a startPosition ] X ] ) +
										fabs( [ startPosition Y ] - [ [ a startPosition ] Y ] ) +
										fabs( [ startPosition Z ] - [ [ a startPosition ] Z ] ) ;
		
		distHostEndToArgumentStart =	fabs( [ endPosition X ] - [ [ a startPosition ] X ] ) +
										fabs( [ endPosition Y ] - [ [ a startPosition ] Y ] ) +
										fabs( [ endPosition Z ] - [ [ a startPosition ] Z ] ) ;

		distHostStartToArgumentEnd =	fabs( [ startPosition X ] - [ [ a endPosition ] X ] ) +
										fabs( [ startPosition Y ] - [ [ a endPosition ] Y ] ) +
										fabs( [ startPosition Z ] - [ [ a endPosition ] Z ] ) ;
										
		distHostEndToArgumentEnd =		fabs( [ endPosition X ] - [ [ a endPosition ] X ] ) +
										fabs( [ endPosition Y ] - [ [ a endPosition ] Y ] ) +
										fabs( [ endPosition Z ] - [ [ a endPosition ] Z ] ) ;
										
		enum { S_TO_S, E_TO_S, S_TO_E, E_TO_E } connectType ;
		double minDist ;
		
		minDist = distHostStartToArgumentStart ;
		connectType = S_TO_S ;
		
		if( distHostEndToArgumentStart < minDist )
			{
				minDist = distHostEndToArgumentStart ;
				connectType = E_TO_S ;
			}
			
		if( distHostStartToArgumentEnd < minDist )
			{
				minDist = distHostStartToArgumentEnd ;
				connectType = S_TO_E ;
			}
		
		if( distHostEndToArgumentEnd < minDist )
			{
				minDist = distHostEndToArgumentEnd ;
				connectType = E_TO_E ;
			}
								
										
		switch( connectType )
			{
				case S_TO_S:
					[ startConnections addObject:a ] ;
					[ startConnectStart addObject:[ NSNumber numberWithBool:YES ] ] ;
				
				break ;
				
				case E_TO_S:
				
					[ endConnections addObject:a ] ;
					[ endConnectStart addObject:[ NSNumber numberWithBool:YES ] ] ;
				
				break ;
				
				case S_TO_E:
				
					[ startConnections addObject:a ] ;
					[ startConnectStart addObject:[ NSNumber numberWithBool:NO ] ] ;
				
				break ;
				
				case E_TO_E:
				
					[ endConnections addObject:a ] ;
					[ endConnectStart addObject:[ NSNumber numberWithBool:NO ] ] ;
				
				break ;
			}
			
		if( minDist > 0.001 )
			{
				printf( "WARNING: ADDING SUSPICIOUS CONNECTION BETWEEN ARCS;  DIST = %f\n", minDist ) ;
			}
			
		NSEnumerator *arcEnumerator ;
		SMArc *nextArc ;
		
		if( desc == YES )
			{
							
				switch( connectType )
					{
						case S_TO_S:
						case S_TO_E:
						
						arcEnumerator = [ startConnections objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								if( nextArc == a ) continue ;
								
								[ nextArc addConnectionToArc:a descend:NO ] ;
							}
							
						break ;
						
						case E_TO_S:
						case E_TO_E:
						
						arcEnumerator = [ endConnections objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								if( nextArc == a ) continue ;
								
								[ nextArc addConnectionToArc:a descend:NO ] ;
							}
							
						break ;
					}
			}
						
		return ;
	}
						

- (NSString *) description
	{
		NSString *returnString ;
		
		int startVertexIndex, endVertexIndex ;
		
		if( startVertex )
			{
				startVertexIndex = [ startVertex index ] ;
			}
		else
			{
				startVertexIndex = -1 ;
			}
			
		if( endVertex )
			{
				endVertexIndex = [ endVertex index ] ;
			}
		else
			{
				endVertexIndex = -1 ;
			}
		
		returnString = [ NSString 
			stringWithFormat:@"ID:%p Start:%f %f %f End:%f %f %f  PhiS:%f ThetaS:%f PhiE:%f ThetaE:%f startV:%p endV:%p startVIdx:%d endVIdx:%d Length:%f\n",self,[startPosition X],[startPosition Y],[startPosition Z],
																			[endPosition X],[endPosition Y],[endPosition Z], phiStart, thetaStart, phiEnd, thetaEnd,
																			startVertex, endVertex, startVertexIndex, endVertexIndex, length  ] ;
																			
		return returnString ;
	}
		

@end
