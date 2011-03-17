//
//  SMArc.h
//  4SMART
//
//  Created by Randy Zauhar on 1/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "vector3.h"
#import "SMProbe.h"
#import "SMVertex.h"
#import "SMLimitPlane.h"

#include "platform.h"
#ifdef LINUX
#define TRUE 1
#define FALSE 0
#endif

@class SMTorus ;
@class SMCycle ;
@class SMMol ;

//#import "SMMol.h"


@interface SMArc : NSObject 
{
	// This class represents a circular arc (which could be a generalized arc or a geodesic). 
	
	// The arc has an associated atom/probe center, atom/probe radius, arc center, axis (which points away from probe/atom),
	// unit vector toward start, unit vector toward end, and arc radius. 
	
	// To handle the possibility of self-intersecting surface, the arc is also optionally associated with a "limiting plane". If the
	// arc intersects this plane, the portion of the arc that falls opposite the plane from the probe center is replaced by a line segment
	// connecting the points of intersection. We provide a special flag that indicates if both beginning and end of an arc are intersection points;
	// in that case the "arc" is simply a straight-line segment. 
	
	// The arc is oriented counterclockwise as you look DOWN the arc axis. 
	
	// It is assumed that an arc subtends less than 2 pi radians
	
@public

	MMVector3 *hostCenter ;
	double hostRadius ;
	SMProbe *hostProbe ;
	
	SMTorus *torusSection ;
	int arcType ;		// 0 = contactI, 1 = reentrantR, 2 = contactJ, 3 = reentrantL, 5 = interior saddle, 6 = interior geodesic, 7 = interior circular arc
	
	
	MMVector3 *arcCenter, *arcAxis ;
	
	double arcRadius ;
	
	MMVector3 *startU, *endU, *uPerp ;
	BOOL lineSegment ;
	
	double length ;
	double angle ;
	
	// Flag to allow skipping this arc
	
	BOOL skip ;
	
	
	BOOL startTouchLimit, endTouchLimit ;
	
	// Intersects with limit plane - anglePStart is closest to startU, angleLPEnd closest to endU
	
	double angleLPStart, angleLPEnd ;
	
	// parentArc non-nil marks an arc that has been deleted. 
	
	SMArc *parentArc ;
	
	NSMutableArray *parentCycles ;
	
	
	MMVector3 *startPosition, *endPosition ;
	
	MMVector3 *startTangent, *endTangent ;
	
	// For saddle arcs
	
	// The following angles are defined for internal saddle edges, and for contact/reentrant arcs that outline a saddle region
	double thetaStart, thetaEnd, phiStart, phiEnd ;
	
	
	SMArc *twin ;
	
	SMVertex *startVertex, *endVertex ;
	
	// I am including this new machinery to make it easier to compute vertices for the arc ends
	
	NSMutableArray *startConnections ;
	NSMutableArray *endConnections ;
	
	NSMutableArray *startConnectStart ;
	NSMutableArray *endConnectStart ;
	
	
	
}

- (id) initWithHostCenter:(MMVector3 *)hc radius:(double)hr torusSection:(SMTorus *)ts arcType:(int)at ;

- (void) initializeWithArcCenter:(MMVector3 *)ac arcRadius:(double)ar axis:(MMVector3 *)ax start:(MMVector3 *)s end:(MMVector3 *)e hostProbe:(SMProbe *)hp ;

- (id) initWithTorusSection:(SMTorus *)ts molecule:(SMMol *)mol phiStart:(double)phiS thetaStart:(double)thetaS phiEnd:(double)phiE thetaEnd:(double)thetaE ;

- (id) copyArc ;

- (void) useLimitPlaneWithPoint:(MMVector3 *)p andNormal:(MMVector3 *)n ;

- (double) length ;
- (double) angle ;

- (void) reverse ;

- (BOOL) skip ;

- (void) setSkip:(BOOL)s ;

- (MMVector3 *) hostCenter ;
- (double) hostRadius ;

- (SMProbe *) hostProbe ;

- (double) arcRadius ;

- (SMTorus *) torusSection ;

- (int) arcType ;

- (BOOL) intersectWith:(SMArc *)t ;
- (BOOL) intersectWith2:(SMArc *)arc2 ;

- (void) setStartVertex: (SMVertex *)s ;
- (void) setEndVertex: (SMVertex *)e ;

- (void) setStartU:(MMVector3 *)sU ;
- (void) setEndU:(MMVector3 *)eU ;

- (SMVertex *) startVertex ;
- (SMVertex *) endVertex ;

- (NSMutableArray *) parentCycles ;
- (void) addParentCycle:(SMCycle *)c ;
- (void) removeParentCycle:(SMCycle *)c ;

- (void) mergeParentCyclesFromArc:(SMArc *)a ;

- (MMVector3 *) startPosition ;
- (MMVector3 *) endPosition ;

- (MMVector3 *) startU ;
- (MMVector3 *) endU ;
- (MMVector3 *) uPerp ;
- (MMVector3 *) arcAxis ;

-(MMVector3 *) arcCenter ;

- (void) setStartVertexPosition:(MMVector3 *)s endVertexPosition:(MMVector3 *)e ;
- (void) computeVertexPositions ;

- (void) restoreVectorsAndAngle ;

- (int) computeThetaAndPhiAnglesUsingMolecule:(SMMol *)mol;

- (MMVector3 *) startTangent ;
- (MMVector3 *) endTangent ;

- (SMArc *) twin ;
- (void) setTwin:(SMArc *)t ;

- (double) phiStart ;
- (double) phiEnd ;

- (double) thetaStart ;
- (double) thetaEnd ;

- (void) setPhiStart:(double)ps ;
- (void) setPhiEnd:(double)pe ;

- (void) setThetaStart:(double)ts ;
- (void) setThetaEnd:(double)te ;

- (MMVector3 *) probePositionForPhi:(double) p ;

- (NSMutableArray *) startConnections ;
- (NSMutableArray *) endConnections ;
	
- (NSMutableArray *) startConnectStart ;
- (NSMutableArray *) endConnectStart ;

- (BOOL) withinArc:(MMVector3 *) pos ;

//- (void) addConnectionToArc:(SMArc *)a ;

- (void) arcPoint:(MMVector3 *)p atFraction:(double)f usingMolecule:(SMMol *)m ;

- (MMVector3 *) computePositionForTheta:(double)t andPhi:(double)p usingMolecule:(SMMol *)mol allowSelfIntersection:(BOOL)SIFlag normal:(MMVector3 *)norm ;
 
- (void) theArcPoint:(MMVector3 *)p atFraction:(double) f usingMolecule:(SMMol *)m ;

@end
