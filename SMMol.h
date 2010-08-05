//
//  SMMol.h
//  4SMART
//
//  Created by zauhar on 1/1/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "vector3.h"
//#import "SMProbe.h"
//#import "SMTorus.h"
#import "SMArc.h"

#include "platform.h"
#ifdef LINUX
#define TRUE 1
#define FALSE 0
#endif

@class SMProbe ;
@class SMTorus ;
//@class SMArc ;
@class SMCycle ;

@interface SMMol : NSObject 
{
	// This class represents a molecule in smart007
	
// This is a little clunky - I need the coordinates to be publically accessible
@public 

	double *xAtom, *yAtom, *zAtom, *radii ;
	BOOL *buriedAtom ;
	int nAtoms ;
	
	double probeRadius, maxAtomRadius ;
	
	double gridSpacing, xMin, xMax, yMin, yMax, zMin, zMax ;
	int nGridX, nGridY, nGridZ ;
	
	int ****gridToAtom, ***nAtomsForGrid, ***gridToAtomAlloc ;
	
	// Atom pairs are organized as i < j, mutual neighbors (probe touch distance) are 
	// accumulated for each pair; only unique triple (i<j<k) are used for probe generation
	
	int **atomToPairs ;
	int *nPairsForAtom ;
	
	int nAtomPairs ;
	int **atomPairs ;
	
	MMVector3 **pairAxes ;
	MMVector3 **pairBases ;
	double *torusRadii ;
	
	int **mutualNeighbors ;
	int *nMutualNeighbors ;
	
	BOOL *freeTorus ;
	int nFreeTori ;
	
	NSMutableArray *probes ;
	NSMutableArray **atomToProbes ;
	
	NSMutableDictionary *atomTripleToProbes ;
	
	// Tori
	
	SMTorus **tori ;
	int nTori, torusAlloc ;
	
	// Use hybrid structure here - regular C array points at NSMutableArray instances
	
	NSMutableArray **atomsToTori, **atomsToCycles ;
	
	NSMutableArray *saddleCycles, *contactCycles, *reentrantCycles ;
	
	int nVertices, nEdges, nElements, nVertexAlloc, nContactElements, nReentrantElements, nSaddleElements ;
	
	// Note that vertex norms are derived from position of associated atom or probe, and need not agree with element normal
	
	NSMutableArray *vertices ;
	
	// Allow self-intersecting surface?
	
	BOOL allowSelfIntersection ;
	
	int warningLevel ;
	
	int nSubsurfaces ;
	

}

- (id) initWithMOL2File:(NSString *)m andRadiiFile:(NSString *)rad allowSelfIntersection:(BOOL)asi randomizeUsing:(double)rand ;
- (void) assignAtomsToGridUsingProbeRadius:(double)r ;
- (void) findAtomNeighbors ;

- (void) assignProbePositionsUsingContactTolerance:(double)tol  ;
- (void) assignProbePositionsVer2UsingContactTolerance:(double)tol generateBlockingProbes:(BOOL)generateBlockingProbes usingSolventStart:(int)solventStartIndex ;

- (BOOL) checkForBumpAtProbePosition:(MMVector3 *)pos andTolerance:(double)tol withAtomI:(int)iAtom J:(int)jAtom K:(int)kAtom ;

- (BOOL) addProbe:(SMProbe *)p forTripleI:(int)atomI J:(int)atomJ K:(int)atomK ;

- (NSArray *) getProbesForTripleI:(int)atomI J:(int)atomJ K:(int)atomK ;

- (void) assignOrientationForProbe:(SMProbe *)p usingAtomI:(int)i J:(int)j K:(int)k ;

- (void) generateToriUsingSkipWidth:(double)skipW ;

- (void) generateContactCyclesUsingDivision:(double)div ;

- (void) combineContactCycles ;

- (void) generateReentrantCyclesWithReentrantHandling:(BOOL)rh ;

- (void) generateSaddleCycles ;

- (void) generateVerticesUsingSubsurfaces:(BOOL)sub ;

- (void) exportAsFlatsUsingPath:(NSString *)p oldStyle:(BOOL)s useSubsurfaces:(BOOL)sub ;

- (void) exportAsCubicUsingPath:(NSString *)p useSubsurfaces:(BOOL)sub ;

- (BOOL) reduceContactCyclesUsingGeodesics:(BOOL)useGeo andDivisionParameter:(double)div ;

- (BOOL) reduceReentrantCyclesUsingDivisionParameter:(double)div ;

- (BOOL) reduceSaddleCyclesUsingDivisionParameter:(double)div ;

- (double) computeTangentAngleUsingArc:(SMArc *) a start:(BOOL) s previousArc:(SMArc *) p previousForward:(BOOL)pf
			nextArc:(SMArc *)n nextForward:(BOOL)nf reentrantRegion:(BOOL)reent ;

- (SMArc *) buildCircularArcUsingSourcePre:(SMArc *)sourcePreArc sourcePreForward:(BOOL)spf 
				sourceNext:(SMArc *)sourceNextArc sourceNextForward:(BOOL)snf
				targetPre:(SMArc *)targetPreArc targetPreForward:(BOOL)tpf
				targetNext:(SMArc *)targetNextArc targetNextForward:(BOOL)tnf ;

- (SMArc *) buildCircularArcUsingHostCenter:(MMVector3 *)hc hostRadius:(double)hr start:(MMVector3 *)v1 
					end:(MMVector3 *)v2 startTangent:(MMVector3 *)t ;
					
- (SMArc *) buildGeoUsingStart:(SMArc *)sourceArc startForward:(BOOL)sf end:(SMArc *)targetArc endForward:(BOOL)ef reentrantRegion:(BOOL)reent ;

- (SMArc *) buildSaddleUsingStart:(SMArc *)sourceArc 
					startForward:(BOOL)sf 
					end:(SMArc *)targetArc
					endForward:(BOOL)ef ;

- (double) computeSaddleTangentAngleUsingArc:(SMArc *) a start:(BOOL) s previousArc:(SMArc *) p
														previousForward:(BOOL)pf
														nextArc:(SMArc *)n  nextForward:(BOOL)nf ;
														
- (NSArray *)  subdivideCycle:(SMCycle *)cyc firstIndex:(int)first
					secondIndex:(int)second usingArc:(SMArc *)arc ;
					
- (SMCycle *) combineCycle:(SMCycle *)cycleJ firstIndex:(int)startJ withCycle:(SMCycle *)cycleK secondIndex:(int)startK
				usingArc:(SMArc *)a ;	
				
- (void) subdivideArc:(SMArc *)a ;
- (void) subdivideArc:(SMArc *)a usingDivision:(double)div ;
- (NSArray *) subdivideTheArc:(SMArc *)a ;
- (NSArray *) subdivideTheArc:(SMArc *)a usingDivision:(double) div ;

- (void) probesToMOL2UsingFile:(NSString *)f ;

// For debugging

- (void) writeMOEGraphicsForArc:(SMArc *)a usingColorIndex:(unsigned int)c andGraphicsObject:(NSString *)g ;
- (void) writeMOEGraphicsForCycle:(SMCycle *)cyc usingObjectName:(NSString *)name ;
- (void) writeMOEGraphicsForHistoryOfCycle:(SMCycle *)cyc ;

@end
