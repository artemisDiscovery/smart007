//
//  vector.h
//  MolMon
//
//  Created by zauhar on Mon Feb 25 2002.
//  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

#include "platform.h"
#ifdef LINUX
#define TRUE 1
#define FALSE 0
#endif


// This class implements simple defs of a 3D vector and operations on vector. 
// Will be expanded as needed


@interface MMVector3 : NSObject 
{

    double X ;
    double Y ;
    double Z ;

}

+ (MMVector3 *) xAxis ;
+ (MMVector3 *) yAxis ;
+ (MMVector3 *) zAxis ;

+ (BOOL) isVector:(MMVector3 *)t betweenVector:(MMVector3 *)v1 andVector:(MMVector3 *)v2 usingNormal:(MMVector3 *)n ; 


- (double) X ;
- (double) Y ;
- (double) Z ;

- (double *) XLoc ;
- (double *) YLoc ;
- (double *) ZLoc ;

- (double) length ;

- (double) dotWith: (MMVector3 *)d ;

- (double) distWith: (MMVector3 *)b ;

- (void) setX : (double) x ;
- (void) setY : (double) y ;
- (void) setZ : (double) z ;

- (id) initX : (double)x Y: (double)y Z: (double)z ;
- (id) initAlong: (MMVector3 *) a perpTo: (MMVector3 *) p ;
- (id) initPerpTo: (MMVector3 *) p byCrossWith: (MMVector3 *) c ;

- (id) initBySubtracting:(MMVector3 *)A minus:(MMVector3 *)B ;
- (id) initByAdding:(MMVector3 *)A plus::(MMVector3 *)B ;

- (id) initUsingVector:(MMVector3 *) v ;

- (id) initByCrossing:(MMVector3 *)u and:(MMVector3 *)v ;

- (void) scaleBy: (double)s ;

- (void) addVector:(MMVector3 *)v ;

- (void) reverse ;

- (void) normalize ;
- (BOOL) normalizeWithZero: (double) z ;

- (void) coordPointersX: (double *)xp Y:(double *)yp Z:(double *) zp ;


@end
