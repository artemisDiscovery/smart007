//
//  SMArcEnd.h
//  4SMART
//
//  Created by zauhar on 18/5/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

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
