//
//  SMTorus.m
//  4SMART
//
//  Created by Randy Zauhar on 1/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SMTorus.h"


@implementation SMTorus

- (id) initWithProbeR:(SMProbe *)pi probeL:(SMProbe *)pj atomI:(int)i atomJ:(int)j free:(BOOL)f usingMolecule:(SMMol *)m usingBase:(MMVector3 *)b
				usingTorusAxis:(MMVector3 *)a usingSaddleRadius:(double)sr usingProbeRadius:(double)probeRad
	{
		self = [ super init ] ;
		
		freeTorus = f ;
		
		atomI = i ;
		atomJ = j ;
		
		probeR = pi ;
		probeL = pj ;
		
		// This is shady, but I'm doing it anyway. 
		
		double *xAtom = m->xAtom ;
		double *yAtom = m->yAtom ;
		double *zAtom = m->zAtom ;
		
		double *radii = m->radii ;
		
		contactArcI = reentrantArcR = contactArcJ = reentrantArcL = nil ;
		
		double dij, rI, rJ ;
		
		rI = radii[i] ;
		rJ = radii[j] ;
		
		dij = sqrt( pow( (xAtom[j] - xAtom[i]), 2) + pow( (yAtom[j] - yAtom[i]), 2) + pow( (zAtom[j] - zAtom[i]), 2) ) ;
		
		//axis = [ [ MMVector3 alloc ] initX:(xAtom[j] - xAtom[i]) Y:(yAtom[j] - yAtom[i]) Z:(zAtom[j] - zAtom[i]) ] ;
		
		// Maximum theta angle
		
		thetaMax = acos( (pow( (probeRad + rI), 2) +  pow( (probeRad + rJ), 2) - dij*dij)/(2.*(probeRad + rI)*(probeRad + rJ)) ) ;
		
		[ a retain ] ;
		
		axis = a ;
		
		[ b retain ] ;
		
		base = b ;
		
		saddleRadius = sr ;
		
		double dx, dy, dz, testDot, phiDot ;
		
		if( freeTorus == NO )
			{
				dx = [ probeR X ] - [ base X ] ;
				dy = [ probeR Y ] - [ base Y ] ;
				dz = [ probeR Z ] - [ base Z ] ;
				
				refR = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
				
				[ refR normalize ] ;
				
				dx = [ probeL X ] - [ base X ] ;
				dy = [ probeL Y ] - [ base Y ] ;
				dz = [ probeL Z ] - [ base Z ] ;
				
				refL = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
				
				[ refL normalize ] ;
				
				// Compute total phi angle subtended
				
				phiDot = [ refR dotWith:refL ] ;
				
				MMVector3 *testCross = [ [ MMVector3 alloc ] initByCrossing:refL and:refR ] ;
				
				testDot = [ testCross dotWith:axis ] ;
				
				if( testDot >= 0. )
					{
						phiAngleSubtended = acos( phiDot ) ;
					}
				else
					{
						phiAngleSubtended = 2.*acos(-1.) - acos( phiDot ) ;
					}
					
				[ testCross release ] ;
				
				
			}
		else
			{
				MMVector3 *testX, *testY ;
				
				testX = [ [ MMVector3 alloc ] initByCrossing:axis and:[ MMVector3 xAxis ] ] ;
				testY = [ [ MMVector3 alloc ] initByCrossing:axis and:[ MMVector3 yAxis ] ] ;
				
				if( [ testX length ] > [ testY length ] )	
					{
						refR = testX ;
						[ testY release ] ;
					}
				else
					{
						refR = testY ;
						[ testX release ] ;
					}
					
				[ refR normalize ] ;
				
				refL = refR ;
				
				phiAngleSubtended = 2.*acos(-1.) ;
			}
			
		refPerp = [ [ MMVector3 alloc ] initByCrossing:refR and:axis ] ;
		
		[ refPerp normalize ] ;
				
		
		//[ axis normalize ] ;
		
		skipTorus = NO ;
		
		// Determine if this torus presents self-intersecting surface
		
		if( saddleRadius > probeRad )
			{
				// No chance!
				
				selfIntersection = NO ;
				
				thetaSelfLo = -1. ;
				thetaSelfHi = -1. ;
				
				bufferI = -1. ;
				bufferJ = -1. ;
			}
		else
			{
				// Need another test - compute thetaSelfLo, thetaSelfHi and check against thetaMax for this torus section
				
				// Omega is the angle beween 1) the ray from the probe center to atom I and 2) the perpendicular from probe center to base point
				// Tau is the angle between 1) the ray from probe center to first intersection between probe and axis
				
				double omega, tau ;
				
				omega = acos( saddleRadius / ( probeRad + radii[atomI] ) ) ;
				tau = acos( saddleRadius / probeRad ) ;
				
				thetaSelfLo = omega - tau ;
				thetaSelfHi = omega + tau ;
				
				// NOTE that we must have the range (thetaSelfLo,thetaSelfHi) totally inside the range (0, thetaMax) or totally outside!
				
				if( thetaSelfLo < thetaMax && thetaSelfHi < thetaMax )
					{
						selfIntersection = YES ;
						
						bufferI = ( saddleRadius - probeRad*cos(omega) ) / 2. ;
						bufferJ = ( saddleRadius - probeRad*cos(thetaMax - omega) ) / 2. ;
						
						// Compute theta angles corresponding to entry into the buffer zone
						
						// FOR SIMPLICITY, I am going to use one buffer value, the minimum of bufferI and bufferJ
						
						if( bufferI < bufferJ )
							{
								bufferIJ = bufferI ;
							}
						else
							{
								bufferIJ = bufferJ ;
							}
						
						thetaBufferI = omega - acos( ( saddleRadius - bufferIJ ) / probeRad ) ;
						thetaBufferJ = omega + acos( ( saddleRadius - bufferIJ ) / probeRad ) ;
					}
				else if( thetaSelfLo > thetaMax && thetaSelfHi > thetaMax )
					{
						selfIntersection = NO ;
						
						bufferI = -1. ;
						bufferJ = -1. ;
						bufferIJ = -1. ;
						
						thetaBufferI = -1. ;
						thetaBufferJ = -1. ;
					}
				else
					{
						printf( "WARNING: SELF-INTERSECTION ANGLE RANGE COMPUTATION BROKEN IN SADDLE INITIALIZATION!\n" ) ;
						selfIntersection = NO ;
						
						bufferI = -1. ;
						bufferJ = -1. ;
						bufferIJ = -1. ;
						
						thetaBufferI = -1. ;
						thetaBufferJ = -1. ;
					}
			}
				
				
		
		contactIArcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ; 
		contactJArcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ; 
		reentrantRArcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
		reentrantLArcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
		
		return self ;
	}
		
- (int) atomI 
	{
		return atomI ;
	}
	
- (int) atomJ 
	{
		return atomJ ;
	}

- (SMProbe *) probeR 
	{
		return probeR ;
	}
	
- (SMProbe *) probeL 
	{
		return probeL ;
	}

- (MMVector3 *) base
	{
		return base ;
	}
	
- (MMVector3 *) axis
	{
		return axis ;
	}

- (MMVector3 *)refR 
	{
		return refR ;
	}
	
- (MMVector3 *)refL 
	{
		return refL ;
	}
	
- (MMVector3 *)refPerp 
	{
		return refPerp ;
	}

- (double) saddleRadius
	{
		return saddleRadius ;
	}
	
- (double) phiAngleSubtended
	{
		return phiAngleSubtended ;
	}
	
- (void) computeContactAndReentrantArcsUsingMolecule:(SMMol *)m andSkipWidth:(double)skipW
	{
		// Compute arcs using geometry in hand (right and left probe, atoms, molecule info)
		
		double rI, rJ ;
		
		rI = m->radii[atomI] ;
		rJ = m->radii[atomJ] ;
		
		double dx, dy, dz, dIToPR, dIToPL, dJToPL, dJToPR ;
		
		// Contact Arc I
		
		MMVector3 *uIToPR, *uJToPL, *uIToPL, *uJToPR, *startVertex, *endVertex ;
		
		dx = [ probeR X ] - m->xAtom[atomI] ;
		dy = [ probeR Y ] - m->yAtom[atomI] ;
		dz = [ probeR Z ] - m->zAtom[atomI] ;
		
		dIToPR = sqrt( dx*dx + dy*dy + dz*dz ) ;
		
		uIToPR = [ [ [ MMVector3 alloc ] initX:(dx/dIToPR) Y:(dy/dIToPR) Z:(dz/dIToPR) ] autorelease ] ;
		
		double dotWithAxis = [ uIToPR dotWith:axis ] ;
		
		dx =  m->xAtom[atomI] + rI*dotWithAxis*[axis X] ;
		dy =  m->yAtom[atomI] + rI*dotWithAxis*[axis Y] ;
		dz =  m->zAtom[atomI] + rI*dotWithAxis*[axis Z] ;
		
		MMVector3 *arcCenter = [ [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] autorelease ] ;
		
		double arcR = rI * sqrt( 1. - dotWithAxis*dotWithAxis ) ;
		
		dx = [ probeL X ] - m->xAtom[atomI] ;
		dy = [ probeL Y ] - m->yAtom[atomI] ;
		dz = [ probeL Z ] - m->zAtom[atomI] ;
		
		dIToPL = sqrt( dx*dx + dy*dy + dz*dz ) ;
		
		startVertex = [ [ [ MMVector3 alloc ] 
			initX:(m->xAtom[atomI] + dx*(rI/dIToPL))
				Y:(m->yAtom[atomI] + dy*(rI/dIToPL))
				Z:(m->zAtom[atomI] + dz*(rI/dIToPL)) ] autorelease ] ;
		
		
		MMVector3 *sU = [ [ [ MMVector3 alloc ] initX:( [ startVertex X ] - [ arcCenter X ] )
			Y:( [ startVertex Y ] - [ arcCenter Y ] ) 
			Z:( [ startVertex Z ] - [ arcCenter Z ] ) ] autorelease ] ;
		
		[ sU normalize ] ;
		
		dx = [ probeR X ] - m->xAtom[atomI] ;
		dy = [ probeR Y ] - m->yAtom[atomI] ;
		dz = [ probeR Z ] - m->zAtom[atomI] ;
		
		dIToPR = sqrt( dx*dx + dy*dy + dz*dz ) ;
		
		endVertex = [ [ [ MMVector3 alloc ] 
			initX:(m->xAtom[atomI] + dx*(rI/dIToPR))
				Y:(m->yAtom[atomI] + dy*(rI/dIToPR))
				Z:(m->zAtom[atomI] + dz*(rI/dIToPR)) ] autorelease ] ;
		
		
		MMVector3 *eU = [ [ [ MMVector3 alloc ] initX:( [ endVertex X ] - [ arcCenter X ] )
			Y:( [ endVertex Y ] - [ arcCenter Y ] ) 
			Z:( [ endVertex Z ] - [ arcCenter Z ] ) ] autorelease ] ;
		
		[ eU normalize ] ;
		
		
		MMVector3 *hc = [ [ [ MMVector3 alloc ] initX:m->xAtom[atomI] Y:m->yAtom[atomI] Z:m->zAtom[atomI] ] autorelease ] ;
		
		MMVector3 *arcAxis = [ [ [ MMVector3 alloc ] initX:[ axis X] Y:[ axis Y] Z:[ axis Z] ] autorelease ] ;
		
		contactArcI = [ [ SMArc alloc ] initWithHostCenter:hc radius:rI torusSection:self arcType:0 ] ;
		
		[ contactArcI initializeWithArcCenter:arcCenter arcRadius:arcR axis:arcAxis start:sU end:eU hostProbe:nil  ] ;
		
		[ contactArcI setStartVertexPosition:startVertex endVertexPosition:endVertex ] ;
		
		[ contactArcI setPhiStart:phiAngleSubtended ] ;
		[ contactArcI setPhiEnd:0. ] ;
		[ contactArcI setThetaStart:0. ] ;
		[ contactArcI setThetaEnd:0. ] ;
		
		
		
		// Contact Arc J
				
		
		dx = [ probeL X ] - m->xAtom[atomJ] ;
		dy = [ probeL Y ] - m->yAtom[atomJ] ;
		dz = [ probeL Z ] - m->zAtom[atomJ] ;
		
		dJToPL = sqrt( dx*dx + dy*dy + dz*dz ) ;
		
		uJToPL = [ [ [ MMVector3 alloc ] initX:(dx/dJToPL) Y:(dy/dJToPL) Z:(dz/dJToPL) ] autorelease ] ;
		
		dotWithAxis = [ uJToPL dotWith:axis ] ;
		
		dx =  m->xAtom[atomJ] + rJ*dotWithAxis*[axis X] ;
		dy =  m->yAtom[atomJ] + rJ*dotWithAxis*[axis Y] ;
		dz =  m->zAtom[atomJ] + rJ*dotWithAxis*[axis Z] ;
		
		arcCenter = [ [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] autorelease ] ;
		
		arcR = rJ * sqrt( 1. - dotWithAxis*dotWithAxis ) ;

		dx = [ probeR X ] - m->xAtom[atomJ] ;
		dy = [ probeR Y ] - m->yAtom[atomJ] ;
		dz = [ probeR Z ] - m->zAtom[atomJ] ;
		
		dJToPR = sqrt( dx*dx + dy*dy + dz*dz ) ;
		
		startVertex = [ [ [ MMVector3 alloc ] 
			initX:(m->xAtom[atomJ] + dx*(rJ/dJToPR))
				Y:(m->yAtom[atomJ] + dy*(rJ/dJToPR))
				Z:(m->zAtom[atomJ] + dz*(rJ/dJToPR)) ] autorelease ] ;
		
		
		sU = [ [ [ MMVector3 alloc ] initX:( [ startVertex X ] - [ arcCenter X ] )
			Y:( [ startVertex Y ] - [ arcCenter Y ] ) 
			Z:( [ startVertex Z ] - [ arcCenter Z ] ) ] autorelease ] ;
		
		[ sU normalize ] ;
		
		dx = [ probeL X ] - m->xAtom[atomJ] ;
		dy = [ probeL Y ] - m->yAtom[atomJ] ;
		dz = [ probeL Z ] - m->zAtom[atomJ] ;
		
		dJToPL = sqrt( dx*dx + dy*dy + dz*dz ) ;
		
		endVertex = [ [ [ MMVector3 alloc ] 
			initX:(m->xAtom[atomJ] + dx*(rJ/dJToPL))
				Y:(m->yAtom[atomJ] + dy*(rJ/dJToPL))
				Z:(m->zAtom[atomJ] + dz*(rJ/dJToPL)) ] autorelease ] ;
		
		
		eU = [ [ [ MMVector3 alloc ] initX:( [ endVertex X ] - [ arcCenter X ] )
			Y:( [ endVertex Y ] - [ arcCenter Y ] ) 
			Z:( [ endVertex Z ] - [ arcCenter Z ] ) ] autorelease ] ;
		
		[ eU normalize ] ;
				
		
		hc = [ [ [ MMVector3 alloc ] initX:m->xAtom[atomJ] Y:m->yAtom[atomJ] Z:m->zAtom[atomJ] ] autorelease ] ;
		
		arcAxis = [ [ [ MMVector3 alloc ] initX:(-[ axis X]) Y:(-[ axis Y]) Z:(-[ axis Z]) ] autorelease ] ;
		
		contactArcJ = [ [ SMArc alloc ] initWithHostCenter:hc radius:rJ torusSection:self arcType:2 ] ;
		
		[ contactArcJ initializeWithArcCenter:arcCenter arcRadius:arcR  axis:arcAxis start:sU end:eU hostProbe:nil  ] ;
		
		[ contactArcJ setStartVertexPosition:startVertex endVertexPosition:endVertex ] ;
		
		[ contactArcJ setPhiStart:0. ] ;
		[ contactArcJ setPhiEnd:phiAngleSubtended ] ;
		[ contactArcJ setThetaStart:thetaMax ] ;
		[ contactArcJ setThetaEnd:thetaMax ] ;
		
		
		// I will put in the "limit plabe" later - for now just want to get this working ...
				
		// Reentrant Arc R
		
		dx =  m->xAtom[atomI] - [ probeR X ] ;
		dy =  m->yAtom[atomI] - [ probeR Y ] ;
		dz =  m->zAtom[atomI] - [ probeR Z ] ;
		
		sU = [ [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] autorelease ] ;
		
		[ sU normalize ] ;
		
		dx =  m->xAtom[atomJ] - [ probeR X ] ;
		dy =  m->yAtom[atomJ] - [ probeR Y ] ;
		dz =  m->zAtom[atomJ] - [ probeR Z ] ;
		
		eU = [ [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] autorelease ] ;
		
		[ eU normalize ] ;		
		
		arcAxis = [ [ [ MMVector3 alloc ] initByCrossing:sU and:eU ] autorelease ] ;
		
		[ arcAxis normalize ] ;
		
		arcR = m->probeRadius ;
		
		
		hc = [ [ [ MMVector3 alloc ] initX:[ probeR X ] Y:[ probeR Y ] Z:[ probeR Z ] ] autorelease ] ;
		
		reentrantArcR = [ [ SMArc alloc ] initWithHostCenter:hc radius:m->probeRadius torusSection:self arcType:1  ] ;
		
		[ reentrantArcR initializeWithArcCenter:hc arcRadius:arcR axis:arcAxis start:sU end:eU hostProbe:probeR  ] ;

		[ reentrantArcR setPhiStart:0. ] ;
		[ reentrantArcR setPhiEnd:0. ] ;
		[ reentrantArcR setThetaStart:0. ] ;
		[ reentrantArcR setThetaEnd:thetaMax ] ;
		
		
		// Reentrant Arc L
		
		dx =  m->xAtom[atomJ] - [ probeL X ] ;
		dy =  m->yAtom[atomJ] - [ probeL Y ] ;
		dz =  m->zAtom[atomJ] - [ probeL Z ] ;
		
		sU = [ [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] autorelease ] ;
		
		[ sU normalize ] ;
		
		dx =  m->xAtom[atomI] - [ probeL X ] ;
		dy =  m->yAtom[atomI] - [ probeL Y ] ;
		dz =  m->zAtom[atomI] - [ probeL Z ] ;
		
		eU = [ [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] autorelease ] ;
		
		[ eU normalize ] ;		
		
		arcAxis = [ [ [ MMVector3 alloc ] initByCrossing:sU and:eU ] autorelease ] ;
		
		[ arcAxis normalize ] ;
		
		arcR = m->probeRadius ;
		
		
		hc = [ [ [ MMVector3 alloc ] initX:[ probeL X ] Y:[ probeL Y ] Z:[ probeL Z ] ] autorelease ] ;
		
		reentrantArcL = [ [ SMArc alloc ] initWithHostCenter:hc radius:m->probeRadius torusSection:self arcType:3 ] ;
		
		[ reentrantArcL initializeWithArcCenter:hc arcRadius:arcR axis:arcAxis start:sU end:eU hostProbe:probeL  ] ;
		
		[ reentrantArcL setPhiStart:phiAngleSubtended ] ;
		[ reentrantArcL setPhiEnd:phiAngleSubtended ] ;
		[ reentrantArcL setThetaStart:thetaMax ] ;
		[ reentrantArcL setThetaEnd:0. ] ;
		
		
		// Assign skip flags to arcs
		
		if( freeTorus == NO )
			{
				if( [ contactArcI length ] <= skipW )
					{
						[ contactArcI setSkip:YES ] ;
						[ contactArcJ setSkip:YES ] ;
						[ reentrantArcR setTwin:reentrantArcL ] ;
						[ reentrantArcL setTwin:reentrantArcR ] ;
						skipTorus = YES ;
					}
				else if( [ reentrantArcR length ] <= skipW )
					{
						[ reentrantArcR setSkip:YES ] ;
						[ reentrantArcL setSkip:YES ] ;
						[ contactArcI setTwin:contactArcJ ] ;
						[ contactArcJ setTwin:contactArcI ] ;				
						skipTorus = YES ;
					}
			}
		else
			{
				[ reentrantArcR setTwin:reentrantArcL ] ;
				[ reentrantArcL setTwin:reentrantArcR ] ;
			}
			
			
		
		// Compute vertex positions
		
		//[ contactArcI computeVertexPositions ] ;
		//[ contactArcJ computeVertexPositions ] ;
		
		//[ reentrantArcR computeVertexPositions ] ;
		//[ reentrantArcL computeVertexPositions ] ;
		
		// FOR DEBUGGING - 
		
		// PRINT OUT ARC START AND END POINTS, CCW ORDER
		/*
		printf( "CI %d %f %f %f %f %f %f \n", 
			atomI, 
			[ [ contactArcI startPosition ] X ] ,   
			[ [ contactArcI startPosition ] Y ] , 
			[ [ contactArcI startPosition ] Z ] , 
			[ [ contactArcI endPosition ] X ] ,   
			[ [ contactArcI endPosition ] Y ] , 
			[ [ contactArcI endPosition ] Z ]  ) ;
			
		printf( "RR %f %f %f %f %f %f \n",
			[ [ reentrantArcR startPosition ] X ] ,
			[ [ reentrantArcR startPosition ] Y ] ,
			[ [ reentrantArcR startPosition ] Z ] ,
			[ [ reentrantArcR endPosition ] X ] ,
			[ [ reentrantArcR endPosition ] Y ] ,
			[ [ reentrantArcR endPosition ] Z ]  ) ;
			
			
		printf( "CJ %d %f %f %f %f %f %f \n", 
			atomJ,
			[ [ contactArcJ startPosition ] X ] ,   
			[ [ contactArcJ startPosition ] Y ] , 
			[ [ contactArcJ startPosition ] Z ] , 
			[ [ contactArcJ endPosition ] X ] ,   
			[ [ contactArcJ endPosition ] Y ] , 
			[ [ contactArcJ endPosition ] Z ]  ) ;
			
		printf( "RL %f %f %f %f %f %f \n",
			[ [ reentrantArcL startPosition ] X ] ,
			[ [ reentrantArcL startPosition ] Y ] ,
			[ [ reentrantArcL startPosition ] Z ] ,
			[ [ reentrantArcL endPosition ] X ] ,
			[ [ reentrantArcL endPosition ] Y ] ,
			[ [ reentrantArcL endPosition ] Z ]  ) ;
			*/
			

		// Register arcs
		
		NSArray *cIArray, *cJArray, *rRArray, *rLArray ;
		
		cIArray = [ NSArray arrayWithObjects:contactArcI,nil ] ;
		cJArray = [ NSArray arrayWithObjects:contactArcJ,nil ] ;
		
		rRArray = [ NSArray arrayWithObjects:reentrantArcR,nil ] ;
		rLArray = [ NSArray arrayWithObjects:reentrantArcL,nil ] ;
		
		[ self registerContactIArcs:cIArray parentArc:nil ] ;
		[ self registerContactJArcs:cJArray parentArc:nil ] ;
		[ self registerReentrantRArcs:rRArray parentArc:nil ] ;
		[ self registerReentrantLArcs:rLArray parentArc:nil ] ;
		
		
		// That's it
		
		return ;
	}

		
		
		
- (SMArc *) contactArcI 
	{
		return contactArcI ;
	}
	
- (SMArc *) contactArcJ 
	{
		return contactArcJ ;
	}
	
- (SMArc *) reentrantArcR 
	{
		return reentrantArcR ;
	}
			
- (SMArc *) reentrantArcL 
	{
		return reentrantArcL ;
	}

/*
- (void) setSkipTorus:(BOOL)s 
	{
		skipTorus = s ;
		
		return ;
	}
*/

- (BOOL) skipTorus 
	{
		return skipTorus ;
	}
	
- (BOOL) freeTorus
	{
		return freeTorus ;
	}

	
- (void) registerArcs:(NSArray *)arcs parentArc:(SMArc *)p 
	{
		// See which array parent arc is in
		
		int idx ;
		
		if( ( idx = [ contactIArcs indexOfObject:p ] ) != NSNotFound )
			{
				[ self registerContactIArcs:arcs parentArc:p ] ;
			}
		else if( ( idx = [ contactJArcs indexOfObject:p ] ) != NSNotFound )
			{
				[ self registerContactJArcs:arcs parentArc:p ] ;
			}
		else if( ( idx = [ reentrantRArcs indexOfObject:p ] ) != NSNotFound )
			{
				[ self registerReentrantRArcs:arcs parentArc:p ] ;
			}
		else if( ( idx = [ reentrantLArcs indexOfObject:p ] ) != NSNotFound )
			{
				[ self registerReentrantLArcs:arcs parentArc:p ] ;
			}
		else
			{
				printf( "COULD NOT REGISTER ARC WITH TORUS - Exit!\n" ) ;
				exit(1) ;
			}
	}

- (void) registerContactIArcs:(NSArray *)cI parentArc:(SMArc *)p 
	{
		// Add to contactIArcs ; replace if parent arc specified
		
		if( p )
			{
				int idx ;
				
				idx = [ contactIArcs indexOfObject:p ] ;
				
				if( idx == NSNotFound )
					{
						printf( "BROKEN POINTER TO CONTACTI ARC - Exit!\n" ) ;
						exit(1) ;
					}
					
				NSEnumerator *arcEnumerator ;
				arcEnumerator = [ cI objectEnumerator ] ;
				
				SMArc *nextArc ;
				
				[ contactIArcs removeObjectAtIndex:idx ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						[ contactIArcs insertObject:nextArc atIndex:idx ] ;
						++idx ;
					}
				
			}
		else
			{
				// Sanity check - array should be empty
				
				if( [ contactIArcs count ] != 0 )
					{
						printf( "WARNING: ARC ADDED TO NONEMPTY CONTACTI ARRAY, AND WITH NULL PARENT!\n" ) ;
					}
				
				[ contactIArcs addObjectsFromArray:cI ] ;
			}
			
		return ;
	}
				
			
- (void) registerContactJArcs:(NSArray *)cJ parentArc:(SMArc *)p 
	{
		// Add to contactJArcs ; replace if parent arc specified
		
		if( p )
			{
				int idx ;
				
				idx = [ contactJArcs indexOfObject:p ] ;
				
				if( idx == NSNotFound )
					{
						printf( "BROKEN POINTER TO CONTACTJ ARC - Exit!\n" ) ;
						exit(1) ;
					}
					
				NSEnumerator *arcEnumerator ;
				arcEnumerator = [ cJ objectEnumerator ] ;
				
				SMArc *nextArc ;
				
				[ contactJArcs removeObjectAtIndex:idx ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						[ contactJArcs insertObject:nextArc atIndex:idx ] ;
						++idx ;
					}
				
			}
		else
			{
				// Sanity check - array should be empty
				
				if( [ contactJArcs count ] != 0 )
					{
						printf( "WARNING: ARC ADDED TO NONEMPTY CONTACTJ ARRAY, AND WITH NULL PARENT!\n" ) ;
					}
				
				[ contactJArcs addObjectsFromArray:cJ ] ;
			}
			
		return ;
	}

- (void) registerReentrantRArcs:(NSArray *)rR parentArc:(SMArc *)p
	{
		// Add to reentrantRArcs ; replace if parent arc specified
		
		if( p )
			{
				int idx ;
				
				idx = [ reentrantRArcs indexOfObject:p ] ;
				
				if( idx == NSNotFound )
					{
						printf( "BROKEN POINTER TO REENTRANTR ARC - Exit!\n" ) ;
						exit(1) ;
					}
					
				NSEnumerator *arcEnumerator ;
				arcEnumerator = [ rR objectEnumerator ] ;
				
				SMArc *nextArc ;
				
				[ reentrantRArcs removeObjectAtIndex:idx ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						[ reentrantRArcs insertObject:nextArc atIndex:idx ] ;
						++idx ;
					}
				
			}
		else
			{
				// Sanity check - array should be empty
				
				if( [ reentrantRArcs count ] != 0 )
					{
						printf( "WARNING: ARC ADDED TO NONEMPTY REENTRANTR ARRAY, AND WITH NULL PARENT!\n" ) ;
					}
				
				[ reentrantRArcs addObjectsFromArray:rR ] ;
			}
			
		return ;
	}

- (void) registerReentrantLArcs:(NSArray *)rL parentArc:(SMArc *)p 
	{
		// Add to reentrantRArcs ; replace if parent arc specified
		
		if( p )
			{
				int idx ;
				
				idx = [ reentrantLArcs indexOfObject:p ] ;
				
				if( idx == NSNotFound )
					{
						printf( "BROKEN POINTER TO REENTRANTL ARC - Exit!\n" ) ;
						exit(1) ;
					}
					
				NSEnumerator *arcEnumerator ;
				arcEnumerator = [ rL objectEnumerator ] ;
				
				SMArc *nextArc ;
				
				[ reentrantLArcs removeObjectAtIndex:idx ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						[ reentrantLArcs insertObject:nextArc atIndex:idx ] ;
						++idx ;
					}
				
			}
		else
			{
				// Sanity check - array should be empty
				
				if( [ reentrantLArcs count ] != 0 )
					{
						printf( "WARNING: ARC ADDED TO NONEMPTY REENTRANTL ARRAY, AND WITH NULL PARENT!\n" ) ;
					}
				
				[ reentrantLArcs addObjectsFromArray:rL ] ;
			}
			
		return ;
	}

- (NSArray *)contactIArcs 
	{
		return contactIArcs ;
	}
	
- (NSArray *)contactJArcs 
	{
		return contactJArcs ;
	}
	
- (NSArray *)reentrantRArcs 
	{
		return reentrantRArcs ;
	}
	
- (NSArray *)reentrantLArcs 
	{
		return reentrantLArcs ;
	}

- (BOOL) withinTorusPhiLimits:(MMVector3 *) disp
	{
		// Check if vector (which is assumed to be a unit vector directed from the torus base) is within the phi limits of the torus
		
		MMVector3  *crossTestS, *crossTestE ;
							
		crossTestS = [ [ MMVector3 alloc ] initByCrossing:refL and:disp ] ;
		crossTestE = [ [ MMVector3 alloc ] initByCrossing:disp and:refR ] ;
		
		double dotTestS, dotTestE ;
		
		dotTestS = [ crossTestS dotWith:axis ] ;
		dotTestE = [ crossTestE dotWith:axis ] ;
		
		[ crossTestS release ] ;
		[ crossTestE release ] ;
			
		if( phiAngleSubtended <=  acos(-1.) )
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

@end
