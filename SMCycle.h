//
//  SMCycle.h
//  4SMART
//
//  Created by zauhar on 6/2/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SMLimitPlane.h"

@class SMArc ;
@class SMVertex ;

// This represents a cycle of any type

@interface SMCycle : NSObject 
{

@public 

	NSMutableArray *arcs ;
	
	NSMutableArray *forward ;
	
	BOOL active ;
	
	SMCycle *parentCycle ;
	
	// Only used for saddle cycles
	
	double thetaMax, phiMax ;
	
	NSMutableSet *atoms ;
	
	SMProbe *probe ;		// For reentrant cycles
	
	// Information for limiting plane for self-intersecting surface
	
	SMLimitPlane *theLimitPlane ;
	
	BOOL selfIntersection ;
	
	BOOL reentrantRegionContact ;
	
	int subsurface ;
	
	
}

- (void) killCycle ;

- (void) addArc:(SMArc *)a forward:(BOOL)f ;

- (NSMutableArray *) arcs ;

- (NSMutableArray *) forward ;

- (NSMutableSet *) atoms ;
- (NSMutableSet *) probes ;

- (BOOL) active ;

- (void) setParentCycle:(SMCycle *)p ;
- (SMCycle *) parentCycle ;

- (void) setThetaMax:(double)tM phiMax:(double)pM ;

- (double) thetaMax ;
- (double) phiMax ;

- (void) addAtomWithIndex:(int)a ;

- (void) setAtoms:(NSSet *)s ;

- (SMProbe *) probe ;
- (void) setProbe:(SMProbe *)p ;

- (int) subsurface ;
- (void) setSubsurface:(int) i ;

- (void) contactElementMidPoint:(MMVector3 *)p andNormal:(MMVector3 *)norm forMolecule:(SMMol *)m ;
- (void) reentrantElementMidPoint:(MMVector3 *)p andNormal:(MMVector3 *)norm forMolecule:(SMMol *)m ;
- (void) saddleElementMidPoint:(MMVector3 *)p andNormal:(MMVector3 *)norm forMolecule:(SMMol *)m ;


@end
