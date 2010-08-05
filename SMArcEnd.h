//
//  SMArcEnd.h
//  4SMART
//
//  Created by zauhar on 18/5/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

#include "platform.h"
#ifdef LINUX
#define TRUE 1
#define FALSE 0
#endif

@class SMArc ;


@interface SMArcEnd : NSObject 
{
	SMArc *theArc ;
	
	BOOL start ;

}

- (id) initWithArc:(SMArc *)a atStart:(BOOL)s ;

- (SMArc *) arc ;
- (BOOL) atStart ;


@end
