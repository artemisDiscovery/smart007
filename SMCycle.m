//
//  SMCycle.m
//  4SMART
//
//  Created by zauhar on 6/2/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SMCycle.h"


@implementation SMCycle

- (id) init
	{
		self = [ super init ] ;
		
		active = YES ;
		
		
		arcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
		
		forward = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
		
		atoms = [ [ NSMutableSet alloc ] initWithCapacity:10 ] ;
		probe = nil ;
		
		parentCycle = nil ;
		
		theLimitPlane = nil ;
		
		selfIntersection = NO ;
		
		reentrantRegionContact = NO ;
		
		subsurface = -1 ;
		
		return self ;
	}
	
- (void) dealloc
	{
		[ arcs release ] ;
		[ forward release ] ;
		[ atoms release ] ;
		
		[ super dealloc ] ;
		
		return ;
	}

- (void) killCycle 
	{
		// Kill this cycle - remove from parentCycle array of each arc 
		
		NSEnumerator *arcEnumerator = [ arcs objectEnumerator ] ;
		
		SMArc *nextArc ;
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				[ [ nextArc parentCycles ] removeObject:self ] ;
			}
			
		active = NO ;
		
		return ;
	
		
	}
	
- (BOOL) active 
	{
		return active ;
	}
	
	

- (void) addArc:(SMArc *)a forward:(BOOL)f 
	{
		[ arcs addObject:a ] ;
		
		[ forward addObject:[ NSNumber numberWithBool:f ] ] ;
		
		[ a addParentCycle:self ] ;
		
		return ;
	}
	
- (BOOL) replaceArc:(SMArc *)oldArc with:(SMArc *)newArc
	{
		// Start position of newArc should co-incide with start position of oldArc (in which case it will
		// adopt the same orientation) or the end (in which case orientation is the reverse)
		
		double dStartToStart, dStartToEnd ;
		
		dStartToStart = [ [ newArc startPosition ] distWith:[ oldArc startPosition ] ] ;
		dStartToEnd = [ [ newArc startPosition ] distWith:[ oldArc endPosition ] ] ;
		
		int replaceIndex = [ arcs indexOfObject:oldArc ] ;
		
		if( replaceIndex == NSNotFound ) 
			{
				return NO ;
			}
		
		BOOL oldOrientation = [ [ forward objectAtIndex:replaceIndex ] boolValue ] ;
		
		BOOL newOrientation ;
		
		[ arcs replaceObjectAtIndex:replaceIndex withObject:newArc  ] ;
		
		if( dStartToStart < dStartToEnd )
			{
				newOrientation = oldOrientation ;
			}
		else
			{
				if( oldOrientation == YES )
					{
						newOrientation = NO ;
					}
				else
					{
						newOrientation = YES ;
					}
			}
			
		[ forward replaceObjectAtIndex:replaceIndex withObject:[ NSNumber numberWithBool:newOrientation ] ] ;
		
		[ newArc addParentCycle:self ] ;
		
		return YES ;
	}

- (NSMutableArray *) arcs
	{
		return arcs ;
	}
	
- (NSMutableArray *) forward
	{
		return forward ;
	}
	
- (NSMutableSet *) atoms
	{
		return atoms ;
	}
	
- (SMProbe *) probe 
	{
		return probe ;
	}
	
	
- (void) setParentCycle:(SMCycle *)p 
	{
		parentCycle = p ;
		
		return ;
	}
	
- (SMCycle *) parentCycle 
	{
		return parentCycle ;
	}
	
- (void) setThetaMax:(double)tM phiMax:(double)pM 
	{
		thetaMax = tM ;
		phiMax = pM ;
		
		return ;
	}

- (double) thetaMax
	{
		return thetaMax ;
	}
	
- (double) phiMax 
	{
		return phiMax ;
	}

- (void) addAtomWithIndex:(int)a
	{
		[ atoms addObject:[ NSNumber numberWithInt:a ] ] ;
		
		return ;
	}

- (void) setAtoms:(NSSet *)s 
	{
		[ atoms release ] ;
		
		atoms = [ [ NSMutableSet alloc ] initWithCapacity:[ s count ] ] ;
		
		[ atoms addObjectsFromArray:[ s allObjects ] ] ;
		
		return ;
	}
	
- (void) setProbe:(SMProbe *)p 
	{
		probe = p ;
		
		[ p retain ] ;
		
		return ;
	}
	
	
- (NSString *) description
	{
		// Find top parent
		
		SMCycle *theParent ;
		
		theParent = self ;
		
		while( [ theParent parentCycle ] )
			{
				theParent = [ theParent parentCycle ] ;
			}
			
		NSString *returnString = [ NSString stringWithFormat:@"CYCLE ID: %p TOP ANCESTOR: %p \n",self, theParent ] ;
		
		returnString = [ returnString stringByAppendingString:[ arcs description ] ] ;
		
		return returnString ;
	}
	
	
- (int) subsurface 
	{
		return subsurface ;
	}
	
- (void) setSubsurface:(int) i 
	{
		subsurface = i ;
		
		return ;
	}
	
- (void) contactElementMidPoint:(MMVector3 *)p andNormal:(MMVector3 *)norm forMolecule:(SMMol *)m 
	{
		// Do this simply using associated atom
		
		double xA, yA, zA, xMid, yMid, zMid  ;
		int iAtom ;
		
		iAtom = [ [ atoms anyObject ] intValue ] ;
		
		xA =  m->xAtom[iAtom] ;
		yA =  m->yAtom[iAtom] ;
		zA =  m->zAtom[iAtom] ;
		
		// Collect vertices
		
		static NSMutableSet *vertexSet ;
		
		if( ! vertexSet )
			{
				vertexSet = [ [ NSMutableSet alloc ] initWithCapacity:10 ] ;
			}
			
		[ vertexSet removeAllObjects ] ;
		
		NSEnumerator *arcEnumerator = [ arcs objectEnumerator ] ;
		
		SMArc *nextArc ;
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				[ vertexSet addObject:[ nextArc startVertex ] ] ;
				[ vertexSet addObject:[ nextArc endVertex ] ] ;
			}
			
		NSEnumerator *vertexEnumerator = [ vertexSet objectEnumerator ] ;
		
		SMVertex *nextVertex ;
		
		xMid = yMid = zMid = 0. ;
		
		while( ( nextVertex = [ vertexEnumerator nextObject ] ) )
			{
				xMid += [ [ nextVertex vertexPosition ] X ] - xA ;
				yMid += [ [ nextVertex vertexPosition ] Y ] - yA ;
				zMid += [ [ nextVertex vertexPosition ] Z ] - zA ;
			}
			
		// Midpoint position
		
		[ norm setX:xMid ] ;
		[ norm setY:yMid ] ;
		[ norm setZ:zMid ] ;
		
		[ norm normalize ] ;
		
		double rad = m->radii[iAtom] ;
		
		[ p setX:( xA + rad*[ norm X ] ) ] ;
		[ p setY:( yA + rad*[ norm Y ] ) ] ;
		[ p setZ:( zA + rad*[ norm Z ] ) ] ;
	
		return ;
	}
		
- (void) reentrantElementMidPoint:(MMVector3 *)p andNormal:(MMVector3 *)norm forMolecule:(SMMol *)m 
	{
		// Do this simply using associated atom
		
		double xP, yP, zP, xMid, yMid, zMid  ;
		int iAtom ;
		
		iAtom = [ [ atoms anyObject ] intValue ] ;
		
		xP =  [ probe X ] ;
		yP =  [ probe Y ] ;
		zP =  [ probe Z ] ;
		
		// Collect vertices
		
		static NSMutableSet *vertexSet ;
		
		if( ! vertexSet )
			{
				vertexSet = [ [ NSMutableSet alloc ] initWithCapacity:10 ] ;
			}
			
		[ vertexSet removeAllObjects ] ;
		
		NSEnumerator *arcEnumerator = [ arcs objectEnumerator ] ;
		
		SMArc *nextArc ;
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				[ vertexSet addObject:[ nextArc startVertex ] ] ;
				[ vertexSet addObject:[ nextArc endVertex ] ] ;
			}
			
		NSEnumerator *vertexEnumerator = [ vertexSet objectEnumerator ] ;
		
		SMVertex *nextVertex ;
		
		xMid = yMid = zMid = 0. ;
		
		while( ( nextVertex = [ vertexEnumerator nextObject ] ) )
			{
				xMid += [ [ nextVertex vertexPosition ] X ] - xP ;
				yMid += [ [ nextVertex vertexPosition ] Y ] - yP ;
				zMid += [ [ nextVertex vertexPosition ] Z ] - zP ;
			}
			
		// Midpoint position
		
		[ norm setX:xMid ] ;
		[ norm setY:yMid ] ;
		[ norm setZ:zMid ] ;
		
		[ norm normalize ] ;
		
		double rad = m->probeRadius ;
		
		// Initial position 
		
		[ p setX:( xP + rad*[ norm X ] ) ] ;
		[ p setY:( yP + rad*[ norm Y ] ) ] ;
		[ p setZ:( zP + rad*[ norm Z ] ) ] ;
		
		// What if limit plane?
		
		BOOL didAdjust = NO ;
		
		if( theLimitPlane && selfIntersection == YES )
			{
				didAdjust = [ theLimitPlane adjustPosition:p andNormal:norm ] ;
			}
		
		if( didAdjust == YES )
			{
				// Use normal from limitplane
				
				//[ norm setX:( [ theLimitPlane->planeNormal X ] ) ] ;
				//[ norm setY:( [ theLimitPlane->planeNormal Y ] ) ] ;
				//[ norm setZ:( [ theLimitPlane->planeNormal Z ] ) ] ;
			}
		else
			{
				
				// Normal was not modified by limit plan method - but we do need to reverse normal!
				
				[ norm reverse ] ;
			}
	
		return ;
	}

- (void) saddleElementMidPoint:(MMVector3 *)p andNormal:(MMVector3 *)norm forMolecule:(SMMol *)m
	{
		// Collect vertices
		
		static NSMutableSet *vertexSet ;
		
		if( ! vertexSet )
			{
				vertexSet = [ [ NSMutableSet alloc ] initWithCapacity:10 ] ;
			}
			
		[ vertexSet removeAllObjects ] ;
		
		NSEnumerator *arcEnumerator = [ arcs objectEnumerator ] ;
		
		SMArc *nextArc ;
		
		double phiMid, thetaMid ;
		
		phiMid = thetaMid = 0. ;		
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				[ vertexSet addObject:[ nextArc startVertex ] ] ;
				[ vertexSet addObject:[ nextArc endVertex ] ] ;
				
				phiMid += nextArc->phiStart ;
				phiMid += nextArc->phiEnd ;
				
				thetaMid += nextArc->thetaStart ;
				thetaMid += nextArc->thetaEnd ;
				
			}
			
		// Each angle should have been added twice
		
		phiMid /= 6 ;
		thetaMid /= 6 ;
		
		MMVector3 *pReturn ;
		
		pReturn = [ [ arcs lastObject ] computePositionForTheta:thetaMid andPhi:phiMid usingMolecule:m 
			allowSelfIntersection:m->allowSelfIntersection normal:norm ] ;
			
		[ p setX:[ pReturn X ] ] ;
		[ p setY:[ pReturn Y ] ] ;
		[ p setZ:[ pReturn Z ] ] ;
		
		[ pReturn release ] ;
			
		// That should do it!
		
		return ;
	}
		
		
			

@end
