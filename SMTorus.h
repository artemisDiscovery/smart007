//
//  SMTorus.h
//  4SMART
//
//  Created by Randy Zauhar on 1/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SMMol.h"
#import "SMArc.h"
//@class SMMol ;

#include "platform.h"
#ifdef LINUX
#define TRUE 1
#define FALSE 0
#endif


@interface SMTorus : NSObject 
{

@public 

	// This class represents a torus section. It is defined by a probe position (which may represent a number of 3-atom contacts)
	// and four circular arcs. In the case of a free torus, the two of the arcs are geometrically co-incident. 
	
	// Circulation around the torus is counterclockwise (as seen from "above"). 
	
	// Contact arc I is on atom I, J on atom J
	// Reentrant arc R is oriented from I to J, while L is oriented J to I.
	
	BOOL freeTorus ;
	
	int atomI, atomJ ;
	
	MMVector3 *axis, *refR, *refL, *refPerp ;
	
	double saddleRadius, phiAngleSubtended ;
	
	MMVector3 *base ;
	
	SMProbe *probeR, *probeL ;
	
	SMArc *contactArcI, *reentrantArcR, *contactArcJ, *reentrantArcL ;
	
	BOOL skipTorus ;  // For tori that are too narrow
	
	// To make it easy to assemble saddle cycles
	
	NSMutableArray *contactIArcs, *contactJArcs, *reentrantRArcs, *reentrantLArcs ;
		
	// The next variables are associated with torus self-intersection
	// bufferI and bufferJ are 1/2 the distance from probe contact points at atoms I and J to the 
	// interatomic axis
	
	double heightI, heightJ, heightIJ ;
	double bufferIJ ;
	
	// Also compute limits on theta for self-intersection
	
	double thetaSelfLo, thetaSelfHi, thetaBufferLo, thetaBufferHi ;
	double thetaMax ;
	
	BOOL selfIntersection ;
	
	

}

- (id) initWithProbeR:(SMProbe *)pi probeL:(SMProbe *)pj atomI:(int)i atomJ:(int)j free:(BOOL)f usingMolecule:(SMMol *)m usingBase:(MMVector3 *)b 
			usingTorusAxis:(MMVector3 *)a usingSaddleRadius:(double)sr usingProbeRadius:(double)probeRad;


- (void) computeContactAndReentrantArcsUsingMolecule:(SMMol *)m andSkipWidth:(double)skipW ;

- (SMArc *) contactArcI ;
- (SMArc *) contactArcJ ;

- (SMArc *) reentrantArcR ;
- (SMArc *) reentrantArcL ;

- (MMVector3 *)base ;
- (MMVector3 *)axis ;

- (double) phiAngleSubtended ;

- (double) saddleRadius ;

- (MMVector3 *)refR ;
- (MMVector3 *)refL ;
- (MMVector3 *)refPerp ;

- (int) atomI ;
- (int) atomJ ;

- (SMProbe *) probeR ;
- (SMProbe *) probeL ;

- (BOOL) skipTorus ;

- (BOOL) freeTorus ;

- (BOOL) withinTorusPhiLimits:(MMVector3 *) disp ;

//- (void) setSkipTorus:(BOOL)s ;

- (void) registerArcs:(NSArray *)arcs parentArc:(SMArc *)p ;

- (void) registerContactIArcs:(NSArray *)cI parentArc:(SMArc *)p ;
- (void) registerContactJArcs:(NSArray *)cJ parentArc:(SMArc *)p ;
- (void) registerReentrantRArcs:(NSArray *)rR parentArc:(SMArc *)p ;
- (void) registerReentrantLArcs:(NSArray *)rL parentArc:(SMArc *)p ;

- (NSArray *)contactIArcs ;
- (NSArray *)contactJArcs ;
- (NSArray *)reentrantRArcs ;
- (NSArray *)reentrantLArcs ;

@end
