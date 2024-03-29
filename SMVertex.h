//
//  SMVertex.h
//  4SMART
//
//  Created by zauhar on 18/5/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "vector3.h"

#include "platform.h"
#ifdef LINUX
#define TRUE 1
#define FALSE 0
#endif


@class SMArcEnd ;
@class SMArc ;
@class SMMol ;
@class SMCycle ;

// Define a vertex object

@interface SMVertex : NSObject
{
@public
	NSMutableSet *arcEnds ;
	
	MMVector3 *position ;
	MMVector3 *normal ;
	
	int subsurface ;
	int subsurfaceIndex ;
	
	int index ;
	
	NSMutableArray *mergeVertices ;
	
	SMVertex *mergedInto ;
	
	int oldIndex ;

}

- (void) addArc:(SMArc *)a start:(BOOL)s ;

- (void) mergeVertex:(SMVertex *) v ;

- (void) computePositionForSaddleUsingMolecule:(SMMol *)m ;
- (void) computePositionForReentrantCycle:(SMCycle *)cyc usingMolecule:(SMMol *)m ;
- (void) computePositionForContact ;

- (NSSet *)arcEnds ;

- (MMVector3 *)vertexPosition ;
- (MMVector3 *)normal ;

- (int) index ;
- (void) setIndex:(int)i ;

- (int) subsurface ;
- (void) setSubsurface:(int)i ;

- (int) subsurfaceIndex ;
- (void) setSubsurfaceIndex:(int)i ;

- (NSComparisonResult) compareSubsurfaceVertices:(SMVertex *)v ;

@end
