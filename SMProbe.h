//
//  SMProbe.h
//  4SMART
//
//  Created by zauhar on 8/1/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "vector3.h"


@interface SMProbe : NSObject 
{
	double xP, yP, zP ;
	
	int nContacts ;
	int *contacts ; // Contact atoms 
	
	int mutualAtom ; // For probe sorting
	
	BOOL rightProbe ;
	BOOL leftProbe ;
	
	double angle ;
	
	
}

- (id) initWithPosition:(MMVector3 *)p andAtomI:(int)i J:(int)j K:(int)k ;

- (id) initWithPosition:(MMVector3 *)p ;

- (void) addContactAtom:(int)a ;

- (void) setMutualAtom:(int)a ;

- (int) mutualAtom ;

- (int) nContacts ;
- (int *) contacts ;

- (double) X ;
- (double) Y ;
- (double) Z ;

- (BOOL) leftProbe ;
- (BOOL) rightProbe ;

- (double) angle ;

- (void) setProbeTypeRight:(BOOL)r left:(BOOL)l ;

- (void) setAngle:(double)a ;

- (NSComparisonResult) compareAngles:(SMProbe *)p ;


@end
