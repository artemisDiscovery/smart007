//
//  vector3.m
//  MolMon
//
//  Created by zauhar on Mon Feb 25 2002.
//  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
//

#import "vector3.h"
#include <math.h>


@implementation MMVector3

+ (MMVector3 *) xAxis 
    {
        MMVector3 *returnVector ;
        
        returnVector = [ [ MMVector3 alloc ] initX:1.0 Y:0.0 Z:0.0 ] ;
        
        [ returnVector autorelease ] ;
        
        return returnVector ;
    }

+ (MMVector3 *) yAxis 
    {
        MMVector3 *returnVector ;
        
        returnVector = [ [ MMVector3 alloc ] initX:0.0 Y:1.0 Z:0.0 ] ;
        
        [ returnVector autorelease ] ;
        
        return returnVector ;
    }

+ (MMVector3 *) zAxis 
    {
        MMVector3 *returnVector ;
        
        returnVector = [ [ MMVector3 alloc ] initX:0.0 Y:0.0 Z:1.0 ] ;
        
        [ returnVector autorelease ] ;
        
        return returnVector ;
    }

+ (BOOL) isVector:(MMVector3 *)t betweenVector:(MMVector3 *)v1 andVector:(MMVector3 *)v2 usingNormal:(MMVector3 *)n
	{
		// Test if vectors span more or less than 180 deg
		
		BOOL moreThan180 ;
		MMVector3 *testCross, *smallN ;
		
		double dot1, dot2 ;
		
		// Generate normal assuming smallest angle rotating v1 into v2
		
		smallN = [ [ MMVector3 alloc ] initByCrossing:v1 and:v2 ] ;
		
		if( [ smallN dotWith:n ] >= 0. )
			{
				moreThan180 = NO ;
			}
		else
			{
				moreThan180 = YES ;
			}
			
		testCross = [ [ MMVector3 alloc ] initByCrossing:v1 and:t ] ;
		
		dot1 = [ testCross dotWith:n ] ;
		
		[ testCross release ] ;
		
		testCross = [ [ MMVector3 alloc ] initByCrossing:t and:v2 ] ;
		
		dot2 = [ testCross dotWith:n ] ;
		
		[ testCross release ] ;
		
		[ smallN release ] ;
		
		if( dot1 >= 0. && dot2 >= 0. )
			{
				// We are inside the smallest angle subtended by v1/v2
				
				if( moreThan180 == NO )
					{
						// Inside
						
						return YES ;
					}
				else
					{
						return NO ;
					}
			}
		else if( dot1 <= 0. && dot2 <= 0. )
			{
				
				if( moreThan180 == NO )
					{
						// Outside
						
						return NO ;
					}
				else
					{
						return YES ;
					}
			}
		else
			{
				return NO ;
			}
			
		return NO ;
	}
				

- (id)initX: (double)x  Y: (double)y  Z: (double)z 
{

    if( ( self = [ super init ] ) )
        {
            X = x ;
            Y = y ;
            Z = z ;
        }
        
    return self ;
}

- (id)initAlong: (MMVector3 *)a   perpTo: (MMVector3 *)p  
{

    double dot ;

    if( ( self = [ super init ] ) )
        {
            X = [ a X ] ;
            Y = [ a Y ] ;
            Z = [ a Z ] ;
        }
    else
        {
            return self ;
        }
        
    
    dot = [ a dotWith:p ] ;
    
    X = X - dot*[p X] ;
    Y = Y - dot*[p Y] ;
    Z = Z - dot*[p Z] ;
    
    if( [ self normalizeWithZero:0.0001 ] == NO )
		{
			return nil ;
		}
    
    return self ;
}

- (id) initByCrossing:(MMVector3 *)u and:(MMVector3 *)v 
	{
		double ux, uy, uz, vx, vy, vz ;
		double rx, ry, rz ;
		
		if( ! ( self = [ super init ] ) )
			{
				return nil ;
			}
			
		ux = [ u X ] ;
		uy = [ u Y ] ;
		uz = [ u Z ] ;
		
		vx = [ v X ] ;
		vy = [ v Y ] ;
		vz = [ v Z ] ;
		
		rx = uy*vz - vy*uz ;
		ry = uz*vx - vz*ux ;
		rz = ux*vy - vx*uy ;
		
		[ self setX:rx ] ;
		[ self setY:ry ] ;
		[ self setZ:rz ] ;
		
		return self ;
	}

- (id)initPerpTo: (MMVector3 *)p   byCrossWith: (MMVector3 *)c  
{
    double xp, yp, zp, xc, yc, zc ;
    double xr, yr, zr ;

    if( ( self = [ super init ] ) )
        {
            xp = [ p X ] ;
            yp = [ p Y ] ;
            zp = [ p Z ] ;
            
            xc = [ c X ] ;
            yc = [ c Y ] ;
            zc = [ c Z ] ;
            
            xr = yp*zc - yc*zp ;
            yr = zp*xc - zc*xp ;
            zr = xp*yc - xc*yp ;
            
            [ self setX:xr ] ;
            [ self setY:yr ] ;
            [ self setZ:zr ] ;
            
            return self ;
            
        }
    else
        {
            return self ;
        }
        
}

- (id) initUsingVector:(MMVector3 *) v
	{
		self = [ super init ] ;
		
		[ self setX:[ v X ] ] ;
		[ self setY:[ v Y ] ] ;
		[ self setZ:[ v Z ] ] ;
		
		return self ;
	}

- (void) dealloc
    {
        [ super dealloc ] ;
    
        return ;
    }
    

- (double) length
{
    double l ;
    
    l = sqrt( X*X + Y*Y + Z*Z ) ;
    
    return l ;
}

- (double) X
{
    return X ;
}

- (double) Y
{
    return Y ;
}

- (double) Z
{
    return Z ;
}

- (double *) XLoc
	{
		return & X ;
	}
	
- (double *) YLoc 
	{
		return & Y ;
	}
	
- (double *) ZLoc
	{
		return & Z ;
	}

- (double) dotWith: (MMVector3 *) d
    {
        double dot ;
        
        dot = X*[d X] + Y*[d Y] + Z*[d Z] ;
        
        return dot ;
    }

- (double) distWith: (MMVector3 *)b 
	{
		double dx, dy, dz ;
		
		dx = X - [ b X ] ;
		dy = Y - [ b Y ] ;
		dz = Z - [ b Z ] ;
		
		return sqrt( dx*dx + dy*dy + dz*dz ) ;
	}

- (void) setX : (double) x
    {
        X = x ;
    }
    
- (void) setY : (double) y
    {
        Y = y ;
    }

- (void) setZ : (double) z
    {
        Z = z ;
    }
    

- (void) normalize
    {
        double SIZE ;
        
        SIZE = [ self length ] ;
		
		if( SIZE == 0. ) return ;
        
        X /= SIZE ;
        Y /= SIZE ;
        Z /= SIZE ;
        
    }
	
- (void) scaleBy: (double)s
	{
		X *= s ;
		Y *= s ;
		Z *= s ;
		
		return ;
	}
    
- (void) reverse
	{
		X = -X ;
		Y = -Y ;
		Z = -Z ;
		
		return ;
	}
	
- (BOOL) normalizeWithZero: (double) z
    {
        double SIZE ;
        
        SIZE = [ self length ] ;
        
        if( SIZE < z ) return NO ;
        
        X /= SIZE ;
        Y /= SIZE ;
        Z /=  SIZE ;
        
        return YES ;
    }
        

- (void) coordPointersX: (double *)xp Y:(double *)yp Z:(double *) zp 
	{
		xp = & X ;
		yp = & Y ;
		zp = & Z ;
		
		return ;
	}
	
			
		
    
- (void) encodeWithCoder: (NSCoder *)coder
	{
		[ coder encodeValueOfObjCType:@encode(double) at:&X ] ;
		[ coder encodeValueOfObjCType:@encode(double) at:&Y ] ;
		[ coder encodeValueOfObjCType:@encode(double) at:&Z ] ;
		
		return ;
	}
	
- (id) initWithCoder: (NSCoder *)coder
	{
		self = [ super init ] ;
		
		[ coder decodeValueOfObjCType:@encode(double) at:&X ] ;
		[ coder decodeValueOfObjCType:@encode(double) at:&Y ] ;
		[ coder decodeValueOfObjCType:@encode(double) at:&Z ] ;
		
		return self ;
	}
	
    
@end
