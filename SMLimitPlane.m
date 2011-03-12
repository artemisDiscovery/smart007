//
//  SMLimitPlane.m
//  4SMART
//
//  Created by zauhar on 8/6/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SMLimitPlane.h"
#include <math.h>


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
		
		double *xAtom, *yAtom, *zAtom ;
		
		xAtom = mol->xAtom ;
		yAtom = mol->yAtom ;
		zAtom = mol->zAtom ;
		
		// Compute normal to plane
		
		MMVector3 *v12, *v23, *v31, *vCross, *vCrossMax ;
		
		// Use first three atoms
		
		if( [ cycleAtoms count ] < 3 )
			{
				printf( "WARNING: ONLY %d ATOMS IN REENTRANT CYCLE - CANNOT BUILD LIMIT PLANE!\n", [ cycleAtoms count ] ) ;
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
		
		double projHeightMin ;
		
		projHeightMin = 1.e6 ; 
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				SMTorus *torus = [ nextArc torusSection ] ;
			
				
				dx = [ planeCenter X ] - [ [ torus base ] X ] ;
				dy = [ planeCenter Y ] - [ [ torus base ] Y ] ;
				dz = [ planeCenter Z ] - [ [ torus base ] Z ] ;
				
				double b = sqrt( dx*dx + dy*dy + dz*dz ) ;
				
				double r = sqrt( b*b + probeHeight*probeHeight ) ;
				
				double projHeight ;
			
				if (torus->bufferIJ > 0. ) {
					projHeight = (probeHeight/r)*(torus->bufferIJ ) ;
				}
				else {
					projHeight = (probeHeight/r)*(torus->heightIJ ) ;
				}

				if( projHeight < projHeightMin )
					{
						projHeightMin = projHeight ;
					}
			}
			
				
		bufferHeight = projHeightMin / 2 ;
				
		// BUT is this self-intersecting or not?
		
		if( probeHeight < mol->probeRadius )
			{
				selfIntersection = YES ;
			
				// Compute theta values for limit-plane intersection
			
				// This is based on a rather complicated diagram that I will "verbally" describe - 
				// We are computing the line of intesection between the plane described by the arc normal and 
				// the limit plane (which lies distance 'h' below the probe center 'C' ). As the arc traces out its path, there 
				// are three points of interest: p1 = first intersection with limit plane, p2 = second, m = midpoint between p1 and p2,
				// all on the line of intersection between the planes. 
				// 'a' is the arc normal ; if the "official" arc normal has negative inner product with limit plane normal, a is the reverse ;
				// 'omega' is the angle between limit plane normal and 'a' ( < 90 deg)
				// 'phi' is angle between arc ray that passes through 'm' and the limit plane perpendicular ; phi = 180 - 90 - omega = 90 - omega
				// Accessory variables :
				// t = displacement from limit plane center to midpoint m : t = h*tan(phi)
				// w = displacement from probe center to 'm' along arc ray : w = sqrt( h^2 + t^2 )
				// Delta = displacement along line of intersection from m to p1 or p2 : Delta = sqrt( probeRad^2 - w^2 )
				//
				// Make local coord system : Nlp (limit plane normal), e = a (schmit ortho to ) Nlp , f = Nlp X e
				// Then in global coords, 
				// p1 = C - h*Nlp + t*e - Delta*f
				// p2 = C - h*Nlp + t*e + Delta*f
				//
				// NOTE that if Rp cos(phi) <= h there is no intersection between the arc plane and the limit plane
			
				MMVector3 *A, *E, *F, *p1, *p2 ;
			
				for( SMArc *nextArc in arcs ) {
					
					A = [ [ MMVector3 alloc ] initUsingVector:nextArc->arcAxis ] ;
					
					if ( [ planeNormal dotWith:nextArc->arcAxis ] < 0. ) {
						[ A reverse ] ; 
					}
					
					double omega = acos([ A dotWith:planeNormal ]) ;
					double phi = acos(-1.)/2. - omega ;
					
					if ( mol->probeRadius*cos(phi) <= probeHeight ) {
						[ A release ] ;
						continue ;
					}
					
					double t = probeHeight * tan(phi) ;
					double w = sqrt(probeHeight*probeHeight + t*t ) ;
					
					double delta = sqrt(mol->probeRadius*mol->probeRadius - w*w ) ;
					
					E = [ [ MMVector3 alloc ] initAlong:A perpTo:planeNormal ] ;
					F = [ [ MMVector3 alloc ] initByCrossing:planeNormal and:E ] ;
					
					double dx = [ cycleProbe X ] - probeHeight * [ planeNormal X ] + t * [ E X ] - delta * [ F X ] ;
					double dy = [ cycleProbe Y ] - probeHeight * [ planeNormal Y ] + t * [ E Y ] - delta * [ F Y ] ;
					double dz = [ cycleProbe Z ] - probeHeight * [ planeNormal Z ] + t * [ E Z ] - delta * [ F Z ] ;
					
					
					p1 = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
					
					// Other side 
					dx = dx + 2.*delta*[ F X ] ;
					dy = dy + 2.*delta*[ F Y ] ;
					dz = dz + 2.*delta*[ F Z ] ;
					
					p2 = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
					
					// Convert to displacement from probe center
					
					[ p1 setX:( [ p1 X ] - [ nextArc->arcCenter X ] ) ] ;
					[ p1 setY:( [ p1 Y ] - [ nextArc->arcCenter Y ] ) ] ;
					[ p1 setZ:( [ p1 Z ] - [ nextArc->arcCenter Z ] ) ] ;
					
					[ p1 normalize ] ;
					
					[ p2 setX:( [ p2 X ] - [ nextArc->arcCenter X ] ) ] ;
					[ p2 setY:( [ p2 Y ] - [ nextArc->arcCenter Y ] ) ] ;
					[ p2 setZ:( [ p2 Z ] - [ nextArc->arcCenter Z ] ) ] ;
					
					[ p2 normalize ] ;
					
					double sgn = nextArc->thetaEnd > nextArc->thetaStart ? 1. : -1. ;
					
					double t1 = acos([ nextArc->startU dotWith:p1 ]) ;
					double t2 = acos([ nextArc->startU dotWith:p2 ]) ;
					
					if (t1 < t2 ) {
						nextArc->thetaLPStart = thetaStart + sgn * t1 ;
						nextArc->thetaLPEnd = thetaStart + sgn * t2 ;
					}
					else {
						nextArc->thetaLPStart = thetaStart + sgn * t2 ;
						nextArc->thetaLPEnd = thetaStart + sgn * t1 ;
					}
					
					// Sanity check
					
					if (nextArc->thetaLPStart < 0. || nextArc->thetaLPEnd < 0. ||
						nextArc->thetaLPStart > nextArc->torusSection->thetaMax || nextArc->thetaLPEnd > nextArc->torusSection->thetaMax) {
						printf("LIMIT PLANE CONSTRUCTION FAILED - Could not assign arc-limit plane interctions - Exit!\n" );
						exit(1) ;
					}
					
					
					
				[ A release ] ;
				[ E release ] ;
				[ F release ] ;
				[ p1 release ] ;
				[ p2 release ] ;
					
					
				}
			
			}
						
		return self ;
	}
	
- (BOOL) adjustPosition:(MMVector3 *)p andNormal:(MMVector3 *)n
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
				
				if( n )
					{
						[ n setX:[ planeNormal X ] ] ;
						[ n setY:[ planeNormal Y ] ] ;
						[ n setZ:[ planeNormal Z ] ] ;
					}
			}
			
		return adjust ;
	}
	
	
		
		

				
				
		
		
		
		
				

@end
