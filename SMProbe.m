//
//  SMProbe.m
//  4SMART
//
//  Created by zauhar on 8/1/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SMProbe.h"
#include <math.h>


@implementation SMProbe

- (id) initWithPosition:(MMVector3 *)p andAtomI:(int)i J:(int)j K:(int)k 
	{
		self = [ super init ] ;
		
		xP = [ p X ] ;
		yP = [ p Y ] ;
		zP = [ p Z ] ;
		
		nContacts = 3 ;
		
		contacts = (int *) malloc( nContacts * sizeof( int ) ) ;
		
		contacts[0] = i ;
		contacts[1] = j ;
		contacts[2] = k ;
		
		angle = 0. ;
		
		// We assume i < j < k and will not check that here!
		
		return self ;
	}
	
- (id) initWithPosition:(MMVector3 *)p
	{
		self = [ super init ] ;
		
		xP = [ p X ] ;
		yP = [ p Y ] ;
		zP = [ p Z ] ;
		
		nContacts = 0 ;
		
		contacts = nil ;
		
		leftProbe = rightProbe = NO ;
		
		return self ;
	}
	

- (void) addContactAtom:(int)a 
	{
		++nContacts ;
		
		if( contacts )
			{
				contacts = (int *) realloc( contacts, nContacts * sizeof( int ) ) ;
			}
		else
			{
				contacts = (int *) malloc( nContacts * sizeof( int ) ) ;
			}
			
		contacts[nContacts - 1] = a ;
	}
	
- (void) setMutualAtom:(int)a 
	{
		mutualAtom = a ;
		
		return ;
	}

- (int) mutualAtom 
	{
		return mutualAtom ;
	}

- (int) nContacts
	{
		return nContacts ;
	}
	
- (int *) contacts
	{
		return contacts ;
	}
	
- (BOOL) leftProbe 
	{
		return leftProbe ;
	}
	
- (BOOL) rightProbe 
	{
		return rightProbe ;
	}

- (void) setProbeTypeRight:(BOOL)r left:(BOOL)l 
	{
		leftProbe = l ;
		rightProbe = r ;
		
		return ;
	}
		
- (void) setAngle:(double) a 
	{
		angle = a ; 
		
		return ;
	}
	
- (double) angle
	{
		return angle ;
	}
	
- (double) X 
	{
		return xP ;
	}
	
- (double) Y 
	{
		return yP ;
	}
	
- (double) Z 
	{
		return zP ;
	}
	
- (NSComparisonResult) compareAngles:(SMProbe *)p
	{
		// See if probes are very close
		
		if( fabs( [ self angle ] - [ p angle ] ) < 1.e-6 )
			{
				// Keep close right probe first
				
				if( rightProbe == YES && [ p leftProbe ] == YES )
					{
						return NSOrderedAscending ;
					}
				else if( leftProbe == YES && [ p rightProbe ] == YES )
					{
						return NSOrderedDescending ;
					}
					
				// else, press on...
			}
			
		if( [ self angle ] < [ p angle ] )
			{
				return NSOrderedAscending ;
			}
		else if( [ self angle ] > [ p angle ] )
			{
				return NSOrderedDescending ;
			}
		else
			{
				return NSOrderedSame ;
			}
			
		return NSOrderedSame ;
	}
	
- (NSString *) description
	{
		NSString *ret = [ NSString stringWithFormat:@"ID: %p - Angle:%f Left:%d Right:%d", self, angle, leftProbe, rightProbe ] ;
		
		return ret ;
	}

@end
