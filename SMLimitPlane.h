//
//  SMLimitPlane.h
//  4SMART
//
//  Created by zauhar on 8/6/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SMTorus.h"


@interface SMLimitPlane : NSObject 
{
	// This small object defines a limit plane for preventing self-inersecting surface in reentrant cycles
	
@public

	NSMutableSet *cycleAtoms ;
	
	SMProbe *cycleProbe ;
	
	NSMutableArray *cycleTori ;
	
	MMVector3 *planeCenter, *planeNormal ;
	
	double bufferHeight, probeHeight ;
	
@public

	BOOL selfIntersection ;
	BOOL bordersSelfintersectingTorus ;

}

- (id) initWithReentrantArcs:(NSArray *)arcs usingMolecule:(SMMol *)mol ;

- (BOOL) adjustPosition:(MMVector3 *)p andNormal:(MMVector3 *)n ;

@end
