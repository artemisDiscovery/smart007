//
//  SMArcEnd.m
//  4SMART
//
//  Created by zauhar on 18/5/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SMArcEnd.h"


@implementation SMArcEnd

- (id) initWithArc:(SMArc *)a atStart:(BOOL)s 
	{
		self = [ super init ] ;
		
		theArc = a ;
		
		start = s ;
		
		return self ;
	}

- (BOOL) isEqual:(id)aEnd
	{
		if( theArc == [ aEnd arc ]  && start  == [ aEnd start ]  )
			{
				return YES ;
			}
		else
			{
				return NO ;
			}
			
	}


- (unsigned) hash
	{
		return [ theArc hash ] ^ [ [ NSNumber numberWithBool:start ] hash ] ;
	}

- (SMArc *) arc 
	{
		return theArc ;
	}
	
- (BOOL) start 
	{
		return start ;
	}
	
- (NSString *) description 
	{
		NSString *returnString = [ NSString stringWithFormat:@"\tArc:%p Start:%d\n", theArc, start ] ;
		
		return returnString ;
	}
	
@end
