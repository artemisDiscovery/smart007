//
//  SMLimitPlane.m
//  4SMART
//
//  Created by zauhar on 8/6/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SMLimitPlane.h"


@implementation SMLimitPlane

- (id) initWithReentrantArcs:(NSArray *)arcs usingMolecule:(SMMol *)mol 
	{
		self = [ super init ] ;
		
		selfIntersection = NO ;
		
		bordersSelfintersectingTorus = NO ;
		
		// Cycle probe position
		
		double dx, dy, dz ;
		
		dx = dy = dz = 0. ;
		
		cycleAtoms = [ [ NSMutableSet alloc ] initWithCapacity:3 ] ;
		
		NSEnumerator *arcEnumerator ;
		SMArc *nextArc ;
		
		arcEnumerator = [ arcs objectEnumerator ] ;
		
		// NOTE: Our assumption here is that this is an initial, undivided cycle - all arcs have associated tori and probes
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				SMProbe *arcProbe = [ nextArc hostProbe ] ;
				
				SMTorus *torus = [ nextArc torusSection ] ;			
				
				if( torus->selfIntersection == YES )
					{
						bordersSelfintersectingTorus = YES ;
					}
				
				if( ! arcProbe || ! torus )
					{
						printf( "WARNING: COULD NOT CONSTRUCT LIMIT PLANE FOR REENTRANT CYCLE!\n" ) ;
						return nil ;
					}
					
				dx += [ arcProbe X ] ;
				dy += [ arcProbe Y ] ;
				dz += [ arcProbe Z ] ;
				
				[ cycleAtoms addObject:[ NSNumber numberWithInt:[ torus atomI ] ] ] ;
				[ cycleAtoms addObject:[ NSNumber numberWithInt:[ torus atomJ ] ] ] ;
				
			}
			
		dx /= [ arcs count ] ;
		dy /= [ arcs count ] ;
		dz /= [ arcs count ] ;
		
		cycleProbe = [ [ SMProbe alloc ] initWithPosition:[ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ] ;
		
		// Plane through atoms
		
		dx = dy = dz = 0. ;
		
		NSEnumerator *atomEnumerator ;
		NSNumber *nextAtom ;
		double *xAtom, *yAtom, *zAtom ;
		int a ;
		
		xAtom = mol->xAtom ;
		yAtom = mol->yAtom ;
		zAtom = mol->zAtom ;
		
		// Compute normal to plane
		
		MMVector3 *v12, *v23, *v31, *vCross, *vCrossMax ;
		
		// Use first three atoms
		
		if( [ cycleAtoms count ] < 3 )
			{
				printf( "WARNING: ONLY %d ATOMS IN REENTRANT CYCLE - CANNOT BUILD LIMIT PLANE!\n" ) ;
				return nil ;
			}
			
		NSArray *atoms = [ cycleAtoms allObjects ] ;
			
		int atom1 = [ [ atoms objectAtIndex:0 ] intValue ] ;
		int atom2 = [ [ atoms objectAtIndex:1 ] intValue ] ;
		int atom3 = [ [ atoms objectAtIndex:2 ] intValue ] ;
		
		v12 = [ [ MMVector3 alloc ] initX:(xAtom[atom2] - xAtom[atom1])
			Y:(yAtom[atom2] - yAtom[atom1])
			Z:(zAtom[atom2] - zAtom[atom1]) ] ;
			
		v23 = [ [ MMVector3 alloc ] initX:(xAtom[atom3] - xAtom[atom2])
			Y:(yAtom[atom3] - yAtom[atom2])
			Z:(zAtom[atom3] - zAtom[atom2]) ] ;
			
		v31 = [ [ MMVector3 alloc ] initX:(xAtom[atom1] - xAtom[atom3])
			Y:(yAtom[atom1] - yAtom[atom3])
			Z:(zAtom[atom1] - zAtom[atom3]) ] ;
			
		vCrossMax = [ [ MMVector3 alloc ] initByCrossing:v12 and:v23 ] ;
		
		vCross = [ [ MMVector3 alloc ] initByCrossing:v23 and:v31 ] ;
		
		if( [ vCross length ] > [ vCrossMax length ] )
			{
				[ vCrossMax release ] ;
				vCrossMax = vCross ;
			}
		else
			{
				[ vCross release ] ;
			}
			
		
		vCross = [ [ MMVector3 alloc ] initByCrossing:v31 and:v12 ] ;
		
		if( [ vCross length ] > [ vCrossMax length ] )
			{
				[ vCrossMax release ] ;
				vCrossMax = vCross ;
			}
		else
			{
				[ vCross release ] ;
			}
			
		[ vCrossMax normalize ] ;
		
		// May need to reverse ...
		
		MMVector3 *disp ;
		
		disp = [ [ MMVector3 alloc ] initX:( [ cycleProbe X ] - xAtom[atom1] )
			Y:( [ cycleProbe Y ] - yAtom[atom1] )
			Z:( [ cycleProbe Z ] - zAtom[atom1] ) ] ;
			
		if( [ vCrossMax dotWith:disp ] < 0. )
			{
				[ vCrossMax reverse ] ;
			}
			
		planeNormal = vCrossMax ;
		
		// Find plane center as average using all three atoms
		
		double dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, dot ;
		
		dot = [ planeNormal dotWith:disp ] ;
		
		dx1 = [ cycleProbe X ] - dot*[ planeNormal X ] ;
		dy1 = [ cycleProbe Y ] - dot*[ planeNormal Y ] ;
		dz1 = [ cycleProbe Z ] - dot*[ planeNormal Z ] ;
		
		[ disp setX:( [ cycleProbe X ] - xAtom[atom2] ) ] ;
		[ disp setY:( [ cycleProbe Y ] - yAtom[atom2] ) ] ;
		[ disp setZ:( [ cycleProbe Z ] - zAtom[atom2] ) ] ;
		
		dot = [ planeNormal dotWith:disp ] ;
		
		dx2 = [ cycleProbe X ] - dot*[ planeNormal X ] ;
		dy2 = [ cycleProbe Y ] - dot*[ planeNormal Y ] ;
		dz2 = [ cycleProbe Z ] - dot*[ planeNormal Z ] ;
		
		[ disp setX:( [ cycleProbe X ] - xAtom[atom3] ) ] ;
		[ disp setY:( [ cycleProbe Y ] - yAtom[atom3] ) ] ;
		[ disp setZ:( [ cycleProbe Z ] - zAtom[atom3] ) ] ;
		
		dot = [ planeNormal dotWith:disp ] ;
		
		dx3 = [ cycleProbe X ] - dot*[ planeNormal X ] ;
		dy3 = [ cycleProbe Y ] - dot*[ planeNormal Y ] ;
		dz3 = [ cycleProbe Z ] - dot*[ planeNormal Z ] ;
		
		dx = (dx1 + dx2 + dx3)/3. ;
		dy = (dy1 + dy2 + dy3)/3. ;
		dz = (dz1 + dz2 + dz3)/3. ;
				
		planeCenter = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
		
		[ disp setX:([ cycleProbe X ] - [ planeCenter X ]) ] ;
		[ disp setY:([ cycleProbe Y ] - [ planeCenter Y ]) ] ;
		[ disp setZ:([ cycleProbe Z ] - [ planeCenter Z ]) ] ;
		
		probeHeight = [ disp length ] ;
		
		[ disp release ] ;
		[ v12 release ] ;
		[ v23 release ] ;
		[ v31 release ] ;
							
		// Cycle buffer height
		
		arcEnumerator = [ arcs objectEnumerator ] ;
		
		double projBufferMin ;
		
		projBufferMin = 1.e6 ; 
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				SMTorus *torus = [ nextArc torusSection ] ;
				
				if( torus->bufferIJ <= 0. ) continue ;
				
				dx = [ planeCenter X ] - [ [ torus base ] X ] ;
				dy = [ planeCenter Y ] - [ [ torus base ] Y ] ;
				dz = [ planeCenter Z ] - [ [ torus base ] Z ] ;
				
				double b = sqrt( dx*dx + dy*dy + dz*dz ) ;
				
				double r = sqrt( b*b + probeHeight*probeHeight ) ;
				
				double projBuffer = (probeHeight/r)*(torus->bufferIJ ) ;
				
				if( projBuffer < projBufferMin )
					{
						projBufferMin = projBuffer ;
					}
			}
			
		if( projBufferMin == 1.e6)
			{
				// No selfintersecting tori!
				
				//bufferHeight = probeHeight - mol->probeRadius ;
				
				// Find the height of each probe-atom contact above the plane, and take 1/2 the
				// minimum as the buffer height
				
				NSEnumerator *atomEnumerator = [ cycleAtoms objectEnumerator ] ;
				
				int nextAtomIndex ; NSNumber *nextAtomNumber ;
				
				double minHeight = 1.e6 ;
				double height ;
				
				disp = [ [ MMVector3 alloc ] initX:0. Y:0. Z:0. ] ;
				
				while( ( nextAtomNumber = [ atomEnumerator nextObject ] ) ) 
					{
						nextAtomIndex = [ nextAtomNumber intValue ] ;
						
						[ disp setX:( ( mol->xAtom[nextAtomIndex] - [ cycleProbe X ])*(mol->probeRadius/(mol->probeRadius + mol->radii[nextAtomIndex])) ) ] ;
						[ disp setY:( ( mol->yAtom[nextAtomIndex] - [ cycleProbe Y ])*(mol->probeRadius/(mol->probeRadius + mol->radii[nextAtomIndex])) ) ] ;
						[ disp setZ:( ( mol->zAtom[nextAtomIndex] - [ cycleProbe Z ])*(mol->probeRadius/(mol->probeRadius + mol->radii[nextAtomIndex])) ) ] ;
						
						height = fabs( [ disp dotWith:planeNormal ] ) ;
						
						if( height < minHeight ) minHeight = height ;
						
						
					}
					
				[ disp release ] ;
				
				// bufferHeight = probeHeight / 2 ;  // This is not optimal
				
				bufferHeight = minHeight / 2 ;
			}
		else
			{
				bufferHeight = projBufferMin ;
				
				selfIntersection = YES ;
			}
		
		return self ;
	}
	
- (BOOL) adjustPosition:(MMVector3 *)p
	{
		// Adjust position of supplied point to obey the limit plane
		
		double dx, dy, dz,  dot ;
		BOOL adjust = NO ;
		
		dx = [ p X ] - [ cycleProbe X ] ;
		dy = [ p Y ] - [ cycleProbe Y ] ;
		dz = [ p Z ] - [ cycleProbe Z ] ;
		
		dot = ( dx*[ planeNormal X ] + dy*[ planeNormal Y ] + dz*[ planeNormal Z ] ) ;
		
		if( dot < -(probeHeight - bufferHeight) )
		//if( ( probeHeight - dot ) < bufferHeight )
			{
				double factor = -(probeHeight - bufferHeight) / dot ;
				
				[ p setX:( [ cycleProbe X ] + factor*dx ) ] ;
				[ p setY:( [ cycleProbe Y ] + factor*dy ) ] ;
				[ p setZ:( [ cycleProbe Z ] + factor*dz ) ] ;
				
				adjust = YES ;
			}
			
		return adjust ;
	}
	
	
		
		

				
				
		
		
		
		
				

@end
