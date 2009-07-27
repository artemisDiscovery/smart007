//
//  SMMol.m
//  4SMART
//
//  Created by zauhar on 1/1/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SMMol.h"
#import "SMTorus.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#import "SMCycle.h"
#include <sys/types.h>
#include <unistd.h>
#import "SMArcEnd.h"

#define MAXSUBDIVIDECOUNT 10

//#define SURFDEBUG

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

@implementation SMMol

- (id) initWithMOL2File:(NSString *)m andRadiiFile:(NSString *)rad allowSelfIntersection:(BOOL)asi randomizeUsing:(double)rand
	{
		self = [ super init ] ;

		xAtom = yAtom = zAtom = radii = buriedAtom = nil ;
		nAtoms = 0 ;
		
		probeRadius = 1.4 ;
		maxAtomRadius = 0. ;
		
		gridSpacing = 0. ;
		xMin = yMin = zMin = 1.e+308 ;
		xMax = yMax = zMax = -1.e+308 ;
		
		 gridToAtom = NULL ;
		 nAtomsForGrid = NULL ;
		
		nAtomPairs = 0 ;
		atomPairs = NULL ;
		freeTorus = NULL ;
		nFreeTori = 0 ;
		
		probes = NULL ;
		atomToProbes = NULL ;
		
		allowSelfIntersection = asi ;
		
		saddleCycles = contactCycles = reentrantCycles = nil ;

		// Get radii
		
		NSMutableDictionary *mol2TypeToRadius = [ [ NSMutableDictionary alloc ] initWithCapacity:10 ] ;
		
		FILE *param ;
		
		if( !( param = fopen( [ rad cString ], "r" ) ) )
			{
				printf( "COULD NOT OPEN RADII FILE!\n" ) ;
				return nil ;
			}
			
			
		char buffer[1000], *charPtr, *mol2Type, *typeRadius ;
		
		while( fgets( buffer, 1000, param ) )
			{
				charPtr = buffer ;
				
				while( *charPtr == ' ' || *charPtr == '\t' )
					{
						++charPtr ;
					}
					
				if( *charPtr == '#' || *charPtr == '\n' )
					{
						continue ;
					}
					
				mol2Type = strtok( charPtr, " \t\n" ) ;
				typeRadius = strtok( NULL, " \t\n" ) ;
				
				[ mol2TypeToRadius setObject:[ NSNumber numberWithDouble:atof(typeRadius) ] 
					forKey:[ [ NSString stringWithCString:mol2Type ] uppercaseString ] ] ;
			}
			
		fclose( param ) ;
		
		FILE *mol2 ;
		
		if( !( mol2 = fopen( [ m cString ], "r" ) ) )
			{
				printf( "COULD NOT OPEN MOL2 FILE!\n" ) ;
				return nil ;
			}
			
		
		int atomAlloc = 10 ;
		
		xAtom = (double *) malloc( atomAlloc * sizeof( double ) ) ;
		yAtom = (double *) malloc( atomAlloc * sizeof( double ) ) ;
		zAtom = (double *) malloc( atomAlloc * sizeof( double ) ) ;
		
		radii = (double *) malloc( atomAlloc * sizeof( double ) ) ;
		
		buriedAtom = (BOOL *) malloc( atomAlloc * sizeof( BOOL ) ) ;
		
		enum parseStates { SEEK, ATOM } ;
		
		int parseState = SEEK ;
		
		char *xPtr, *yPtr, *zPtr ;
		NSNumber *radNumber ;
		double r ;
		
		if( rand > 0. )
			{
				srandom( (unsigned long ) getpid() ) ;
			}
		
		while( fgets( buffer, 1000, mol2 ) )
			{
			
				// NOTE: The code below should allow multimol2 files to be successfully read
				
				if( parseState == SEEK )
					{
						if( strstr( buffer, "@<TRIPOS>ATOM" ) )
							{
								parseState = ATOM ; 
								continue ;
							}
					}
				else
					{
						if( strstr( buffer, "@<TRIPOS>BOND" ) )
							{
								parseState = SEEK ;
								continue ;
							}
							
						// Process next atom
						
						charPtr = buffer ;
						
						while( *charPtr == ' ' || *charPtr == '\t' )
							{
								++charPtr ;
							}
						
						// Skip atom number
						
						strtok( charPtr, " \t\n" ) ;
						
						// Skip atom name 
						
						strtok( NULL, " \t\n" ) ;
						
						xPtr = strtok( NULL, " \t\n" ) ;
						yPtr = strtok( NULL, " \t\n" ) ;
						zPtr = strtok( NULL, " \t\n" ) ;
						
						mol2Type = strtok( NULL, " \t\n" ) ;
						
						radNumber = [ mol2TypeToRadius objectForKey:[ [ NSString stringWithCString:mol2Type ] uppercaseString ] ] ;
						
						if( radNumber )
							{
								r = [ radNumber doubleValue ] ;
							}
						else
							{
								r = 1.0 ;
								printf( "WARNING: DEFAULT RADIUS OF 1.0 ASSIGNED FOR ATOM TYPE %s\n", mol2Type ) ;
							}
							
						if( nAtoms == atomAlloc )
							{
								atomAlloc += 10 ;
								
								xAtom = (double *) realloc( xAtom, atomAlloc * sizeof( double ) ) ;
								yAtom = (double *) realloc( yAtom, atomAlloc * sizeof( double ) ) ;
								zAtom = (double *) realloc( zAtom, atomAlloc * sizeof( double ) ) ;
								
								radii = (double *) realloc( radii, atomAlloc * sizeof( double ) ) ;
								
								buriedAtom = (BOOL *) realloc( buriedAtom, atomAlloc * sizeof( BOOL ) ) ;
								
							}
							
						xAtom[nAtoms] = atof( xPtr ) ;
						yAtom[nAtoms] = atof( yPtr ) ;
						zAtom[nAtoms] = atof( zPtr ) ;
						
						// Randomize if needed
						
						if( rand > 0. )
							{
								double disp ;
								
								disp = 2.*( (((double)random())/pow(2,31)) - 1. ) * rand ;
								xAtom[nAtoms] += disp ;
								disp = 2.*( (((double)random())/pow(2,31)) - 1. ) * rand ;
								yAtom[nAtoms] += disp ;
								disp = 2.*( (((double)random())/pow(2,31)) - 1. ) * rand ;
								zAtom[nAtoms] += disp ;
							}
						
						if( xAtom[nAtoms] < xMin ) xMin =  xAtom[nAtoms] ;
						if( yAtom[nAtoms] < yMin ) yMin =  yAtom[nAtoms] ;
						if( zAtom[nAtoms] < zMin ) zMin =  zAtom[nAtoms] ;
						
						if( xAtom[nAtoms] > xMax ) xMax =  xAtom[nAtoms] ;
						if( yAtom[nAtoms] > yMax ) yMax =  yAtom[nAtoms] ;
						if( zAtom[nAtoms] > zMax ) zMax =  zAtom[nAtoms] ;
						
						
						radii[nAtoms] = r ;
						
						buriedAtom[nAtoms] = NO ;
						
						if( r > maxAtomRadius ) maxAtomRadius = r ;
						
						++nAtoms ;
							
					}
			}
			
		fclose( mol2 ) ;
		
		atomTripleToProbes = [ [ NSMutableDictionary alloc ] initWithCapacity:nAtoms ] ;
		
		torusAlloc = 50 ;
		
		tori = (SMTorus **) malloc( torusAlloc * sizeof( SMTorus * ) )  ;
		nTori = 0 ;
		
		atomsToTori = (NSMutableArray **) malloc( nAtoms * sizeof( NSMutableArray *) ) ;
		atomsToCycles = (NSMutableArray **) malloc( nAtoms * sizeof( NSMutableArray *) ) ;
		
		int i ;
		
		for( i = 0 ; i < nAtoms ; ++i )
			{
				atomsToTori[i] = nil ;
				atomsToCycles[i] = nil ;
			}
			
			
		nVertices = 0 ;
		nElements = 0 ;
		nContactElements = 0 ;
		nReentrantElements = 0 ;
		nSaddleElements = 0 ;
		
		vertices = [ [ NSMutableArray alloc ] initWithCapacity:10000 ] ;
		//vertexNorms = [ [ NSMutableArray alloc ] initWithCapacity:10000 ] ;
		
		// Print all warnings
		
		warningLevel = 0 ;
						
				
		return self ;
	}

						
						

						
- (void) assignAtomsToGridUsingProbeRadius:(double)r 
	{
		// Find grid spacing
		
		probeRadius = r ;
		
		// Increase slightly over "minimal" value to avoid roundoff issues
		
		gridSpacing = (2.*(maxAtomRadius + probeRadius)) + 0.1 ;
		
		nGridX = ceil( (xMax - xMin)/gridSpacing ) ;
		nGridY = ceil( (yMax - yMin)/gridSpacing ) ;
		nGridZ = ceil( (zMax - zMin)/gridSpacing ) ;
		
		if( nGridX == 0 ) nGridX = 1 ;
		if( nGridY == 0 ) nGridY = 1 ;
		if( nGridZ == 0 ) nGridZ = 1 ;
		
		gridToAtom = (int ****) malloc( nGridX * sizeof( int *** ) ) ;
		gridToAtomAlloc = (int ***) malloc( nGridX * sizeof( int ** ) ) ;
		nAtomsForGrid = (int ***) malloc( nGridX * sizeof( int ** ) ) ;
		
		int i, j, k, iG, jG, kG ;
		
		for( i = 0 ; i < nGridX ; ++i )
			{
				gridToAtom[i] = (int ***) malloc( nGridY * sizeof( int ** ) ) ;
				gridToAtomAlloc[i] = (int **) malloc( nGridY * sizeof( int * ) ) ;
				nAtomsForGrid[i] = (int **) malloc( nGridY * sizeof( int * ) ) ;
				
				for( j = 0 ; j < nGridY ; ++j )
					{
						gridToAtom[i][j] = (int **) malloc( nGridZ * sizeof( int * ) ) ;
						gridToAtomAlloc[i][j] = (int *) malloc( nGridZ * sizeof( int  ) ) ;
						nAtomsForGrid[i][j] = (int *) malloc( nGridZ * sizeof( int  ) ) ;
						
						for( k = 0 ; k < nGridZ ; ++k )
							{
								gridToAtomAlloc[i][j][k] = 10 ;
								nAtomsForGrid[i][j][k] = 0 ;
								gridToAtom[i][j][k] = (int *) malloc( 10 * sizeof( int ) ) ;
							}
					}
			}
			
		for( i = 0 ; i < nAtoms ; ++i )
			{
				iG = floor( (xAtom[i] - xMin) / gridSpacing ) ;
				jG = floor( (yAtom[i] - yMin) / gridSpacing ) ;
				kG = floor( (zAtom[i] - zMin) / gridSpacing ) ;
				
				if( nAtomsForGrid[iG][jG][kG] == gridToAtomAlloc[iG][jG][kG] )
					{
						gridToAtomAlloc[iG][jG][kG] += 10 ;
						
						gridToAtom[iG][jG][kG] = (int *) realloc( gridToAtom[iG][jG][kG], gridToAtomAlloc[iG][jG][kG] * sizeof( int ) ) ;
					}
					
				gridToAtom[iG][jG][kG][ nAtomsForGrid[iG][jG][kG] ] = i ;
				
				++nAtomsForGrid[iG][jG][kG] ;
				
			}
			
		return ;
	}
		
		


- (void) findAtomNeighbors 
	{
		int i, j, k, a, iG, jG, kG, iG2, jG2, kG2, deltaI, deltaJ, deltaK, iPair ;
		
		atomToPairs = (int **) malloc( nAtoms * sizeof( int * ) ) ;
		nPairsForAtom = (int *) malloc( nAtoms * sizeof( int ) ) ;
		
		
		int *atomToPairsAlloc = (int *) malloc( nAtoms * sizeof( int ) ) ;
		
		for( i = 0 ; i < nAtoms ; ++i )
			{
				atomToPairsAlloc[i] = 10 ;
				
				atomToPairs[i] = (int *) malloc( 10 * sizeof( int ) ) ;
				
				nPairsForAtom[i] = 0 ;
			}
			
			
		int atomPairsAlloc = 100 ;
		
		atomPairs = (int **) malloc( 100 * sizeof( int * ) ) ;
		
		freeTorus = (BOOL *) malloc( 100 * sizeof( BOOL ) ) ;
		
		for( i = 0 ; i < 100 ; ++i )
			{
				atomPairs[i] = (int *) malloc( 2 * sizeof( int ) ) ;
				freeTorus[i] = NO ;
			}
			
		nAtomPairs = 0 ;
		
		for( i = 0 ; i < nAtoms ; ++i )
			{
				iG = floor( ( xAtom[i] - xMin ) / gridSpacing ) ;
				jG = floor( ( yAtom[i] - yMin ) / gridSpacing ) ;
				kG = floor( ( zAtom[i] - zMin ) / gridSpacing ) ;
				
				for( deltaI = -1 ; deltaI <= 1 ; ++deltaI )
					{
						iG2 = iG + deltaI ;
						if( iG2 < 0 || iG2 >= nGridX ) continue ;
						
						for( deltaJ = -1 ; deltaJ <= 1 ; ++deltaJ )
							{
								jG2 = jG + deltaJ ;
								if( jG2 < 0 || jG2 >= nGridY ) continue ;
							
								for( deltaK = -1 ; deltaK <= 1 ; ++deltaK )
									{
										kG2 = kG + deltaK ;
										if( kG2 < 0 || kG2 >= nGridZ ) continue ;
										
										int *atomPtr ;
										
										atomPtr = gridToAtom[iG2][jG2][kG2] ;
										
										for( j = 0 ; j < nAtomsForGrid[iG2][jG2][kG2] ; ++j )
											{
												a = gridToAtom[iG2][jG2][kG2][j] ;
												
												if( a <= i ) continue ;
												
												double d, dx, dy, dz ;
												
												dx = xAtom[i] - xAtom[a] ;
												dy = yAtom[i] - yAtom[a] ;
												dz = zAtom[i] - zAtom[a] ;
												
												d = sqrt( dx*dx + dy*dy + dz*dz ) ;
												
												if( d >= radii[i] + radii[a] + 2.*probeRadius ) continue ;
												
												if( (d + radii[i]) <= radii[a] )
													{
														buriedAtom[i] = YES ;
														printf( "Atom %d is buried in %d \n", i, a ) ;
													}
													
												if( (d + radii[a]) <= radii[i] )
													{
														buriedAtom[a] = YES ;
														printf( "Atom %d is buried in %d \n", a, i ) ;
													}
												
												
												if( nAtomPairs == atomPairsAlloc )
													{
														atomPairsAlloc += 100 ;
									
														atomPairs = (int **) realloc(atomPairs,  atomPairsAlloc * sizeof( int * ) ) ;
														freeTorus = (BOOL *) realloc(freeTorus,  atomPairsAlloc * sizeof( BOOL ) ) ;

														
														for( k = nAtomPairs ; k < atomPairsAlloc ; ++k )
															{
																atomPairs[k] = (int *) malloc( 2 * sizeof( int ) ) ;
																freeTorus[k] = NO ;
															}
															
													}
													
												atomPairs[nAtomPairs][0] = i ;
												atomPairs[nAtomPairs][1] = a ;
												
												if( atomToPairsAlloc[i] == nPairsForAtom[i] )
													{
														atomToPairsAlloc[i] += 10 ;
														
														atomToPairs[i] = (int *) realloc( atomToPairs[i], atomToPairsAlloc[i] * sizeof( int ) ) ;
													}
												
												if( atomToPairsAlloc[a] == nPairsForAtom[a] )
													{
														atomToPairsAlloc[a] += 10 ;
														
														atomToPairs[a] = (int *) realloc(atomToPairs[a],  atomToPairsAlloc[a] * sizeof( int ) ) ;
													}
													
												atomToPairs[i][ nPairsForAtom[i] ] = nAtomPairs ;
												atomToPairs[a][ nPairsForAtom[a] ] = nAtomPairs ;
												
												++nPairsForAtom[i] ;
												++nPairsForAtom[a] ;
												
												++nAtomPairs ;
											}
									}
							}
					}
			}
			
		// Now collect all mutual neighbors of each atom pair
		
		mutualNeighbors = (int **) malloc( nAtomPairs * sizeof( int * ) ) ;
		nMutualNeighbors = (int *) malloc( nAtomPairs * sizeof( int ) ) ;
		
		NSMutableSet *collectSet1 = [ NSMutableSet setWithCapacity:10 ] ;
		NSMutableSet *collectSet2 = [ NSMutableSet setWithCapacity:10 ] ;
		
		NSAutoreleasePool *localPool = [ [ NSAutoreleasePool alloc ] init ] ;
		
		for( iPair = 0 ; iPair < nAtomPairs ; ++iPair )
			{
				if( iPair % 100 == 0 )
					{
						[ localPool release ] ;
						localPool = [ [ NSAutoreleasePool alloc ] init ] ;
					}
					
				if( buriedAtom[ atomPairs[iPair][0] ] == YES || buriedAtom[ atomPairs[iPair][1] ] == YES )
					{
						nMutualNeighbors[iPair] = 0 ;
						continue ;
					}
					
				[ collectSet1 removeAllObjects ] ;
				[ collectSet2 removeAllObjects ] ;
				
				int a1, a2 ;
				
				a1 = atomPairs[iPair][0] ;
				a2 = atomPairs[iPair][1] ;
				
				for( j = 0 ; j < nPairsForAtom[a1] ; ++j )
					{
						[ collectSet1 addObject:[ NSNumber numberWithInt:atomPairs[ atomToPairs[a1][j] ][0] ] ] ;
						[ collectSet1 addObject:[ NSNumber numberWithInt:atomPairs[ atomToPairs[a1][j] ][1] ] ] ;
					}
					
				for( j = 0 ; j < nPairsForAtom[a2] ; ++j )
					{
						[ collectSet2 addObject:[ NSNumber numberWithInt:atomPairs[ atomToPairs[a2][j] ][0] ] ] ;
						[ collectSet2 addObject:[ NSNumber numberWithInt:atomPairs[ atomToPairs[a2][j] ][1] ] ] ;
					}
					
				[ collectSet1 removeObject:[ NSNumber numberWithInt:a1 ] ] ;
				[ collectSet1 removeObject:[ NSNumber numberWithInt:a2 ] ] ;
				
				[ collectSet2 removeObject:[ NSNumber numberWithInt:a1 ] ] ;
				[ collectSet2 removeObject:[ NSNumber numberWithInt:a2 ] ] ;
				
				// Mutual neighbors only
				
				[ collectSet1 intersectSet:collectSet2 ] ;
				
				mutualNeighbors[iPair] = (int *) malloc( [ collectSet1 count ] * sizeof( int ) ) ;
				
				NSEnumerator *mutualEnumerator = [ collectSet1 objectEnumerator ] ;
				NSNumber *nextNeighbor ;
				
				j = 0 ;
				
				
				while( ( nextNeighbor = [ mutualEnumerator nextObject ] ) )
					{
						int m ;
						m = [ nextNeighbor intValue ] ;
						
						if( buriedAtom[m] == YES ) continue ;
						
						mutualNeighbors[iPair][j] = m ;
						++j ;
					}
					
				nMutualNeighbors[iPair] = j ;
				
			}
			
		if( nAtomPairs % 100 != 0 )
			{
				[ localPool release ] ;
			}
				
				
		return ;
	}
												
		


- (void) assignProbePositionsUsingContactTolerance:(double)tol 
	{
		// Probes will be generated for unique triples of atoms i < j < k
		//
		// Note that this also defines an orientation for each probe w.r.t. the 
		// normal of the plane, which is defined using the order of atoms (i -> j -> k ) and 
		// the right-hand rule. '+' probes lie above the plane, '-' below. 
		
		// A special problem is posed by atoms that are buried inside other atoms. No probes will be generated 
		// that include these as contacts, yet they are not part of free tori. We will have a separate check, after all
		// probes are generated, to determine if an atom with no probes is in fact buried. 
		
		int iPair, i, probeCount ;
		
		probeCount = 0 ;
		
		// Allocate atomToProbes array
		
		atomToProbes = (NSMutableArray **) malloc( nAtoms * sizeof( NSMutableArray * ) ) ;
		
		for( i = 0 ; i < nAtoms ; ++i )
			{
				atomToProbes[i] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
			}
		
		// Allocate room for iteratomic axis info
		
		pairAxes = (MMVector3 **) malloc( nAtomPairs * sizeof( MMVector3 *) ) ;
		pairBases = (MMVector3 **) malloc( nAtomPairs * sizeof( MMVector3 *) ) ;
		
		torusRadii = ( double * ) malloc( nAtomPairs * sizeof( double ) ) ;

		
		for( iPair = 0 ; iPair < nAtomPairs ; ++iPair )
			{
				NSAutoreleasePool * localPool = [[NSAutoreleasePool alloc] init];
				
				// Compute axis from 1 to 2
				
				double dx, dy, dz, d, cosTheta, Ri, Rj ;
				
				dx = xAtom[ atomPairs[iPair][1] ] - xAtom[ atomPairs[iPair][0] ] ;
				dy = yAtom[ atomPairs[iPair][1] ] - yAtom[ atomPairs[iPair][0] ] ;
				dz = zAtom[ atomPairs[iPair][1] ] - zAtom[ atomPairs[iPair][0] ] ;
				
				d = sqrt( dx*dx + dy*dy + dz*dz ) ;
				
				if( d == 0. || buriedAtom[atomPairs[iPair][0]] == YES ||  buriedAtom[atomPairs[iPair][1]] == YES) // Atoms co-incide, or one buried
					{
						pairAxes[iPair] = nil ;
						pairBases[iPair] = nil ;
						continue ;
					}
				
				pairAxes[iPair] = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
				[ pairAxes[iPair] scaleBy:(1./d) ] ;
				
				Ri = radii[ atomPairs[iPair][0] ] ;
				Rj = radii[ atomPairs[iPair][1] ] ;
				
				cosTheta = ( d*d + pow( probeRadius + Ri, 2) - pow( probeRadius + Rj, 2) ) / (2. * (probeRadius + Ri) * d) ;
				
				// Compute base point of probe rotation
				
				dx = xAtom[ atomPairs[iPair][0] ] + (Ri + probeRadius)*cosTheta*[ pairAxes[iPair] X ] ;
				dy = yAtom[ atomPairs[iPair][0] ] + (Ri + probeRadius)*cosTheta*[ pairAxes[iPair] Y ] ;
				dz = zAtom[ atomPairs[iPair][0] ] + (Ri + probeRadius)*cosTheta*[ pairAxes[iPair] Z ] ;
				
				pairBases[iPair] = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ]  ;
				
				torusRadii[iPair] = (Ri + probeRadius) * sqrt( 1. - pow( cosTheta, 2 ) ) ;

				
				// For all mutual probes (with acceptable index) find probe positions that do not intersect other atoms
				
				int j, mutAtom ;
				
				for( j = 0 ; j < nMutualNeighbors[iPair] ; ++j )
					{
						mutAtom = mutualNeighbors[iPair][j] ;
						
						
						if( mutAtom < atomPairs[iPair][1] ) continue ;
						
						//printf( "AT I, J, K = %d %d %d \n", atomPairs[iPair][0], atomPairs[iPair][1], mutAtom ) ;

						
						// Generate two probe positions, check for bump
						
						MMVector3 *g, *h ;  // g and h form a right-handed frame with the axis
						
						MMVector3 *axis ;
						
						axis = pairAxes[iPair] ;
						
						if( fabs( [ axis X ] ) < fabs( [ axis Y ] ) )
							{
								g = [ [ MMVector3 alloc ] initByCrossing:axis and:[ MMVector3 xAxis ] ] ;
							}
						else
							{
								g = [ [ MMVector3 alloc ] initByCrossing:axis and:[ MMVector3 yAxis ] ] ;
							}
							
						[ g normalize ] ;
						
						h = [ [ MMVector3 alloc ] initByCrossing:axis and:g ] ;
						
						[ g autorelease ] ; [ h autorelease ] ;
						
						// Probe position p = b + e*g + f*h where b is the torus base point
						// If rk is the position of the mutual probe, s the torus radius,
						// rkb the vector from b to k, and Rk = |rk|, then we demand
						//		(p - rk).(p - rk) = (Rp + Rk)^2
						// Substituting
						//		e(2.*g.rkb) + f(2.*h.rkb) = (Rp + Rk)^2 - |rkb|^2 - s^2
						//		and e^2 + f^2 = 1
						// We substitute f = +-sqrt(1 - e^2), and find
						//
						// e^2(C^2 + D^2) - 2QC+ (Q^2 -s^2D^2) = 0
						// where
						// Q = (Rp + Rk)^2 - |rkb|^2 - s^2
						// Solve for e, then find f = +-sqrt(1 - e^2)
						
						double  e1, f1, e2, f2, C, D, Q, dx, dy, dz, drkb, drk, a, b, c, disc, s2 ;
						double dx1, dy1, dz1, dx2, dy2, dz2, dk1, dk2 ;	
						
						// Solve for vector from base to mutual probe
						
						dx = -(xAtom[mutAtom] - [ pairBases[iPair] X ]) ;
						dy = -(yAtom[mutAtom] - [ pairBases[iPair] Y ]) ;
						dz = -(zAtom[mutAtom] - [ pairBases[iPair] Z ]) ;
						
						MMVector3 *rkb = [ [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] autorelease ] ;
						
						drkb = [ rkb length ] ;
						
						drk = sqrt( pow( xAtom[mutAtom], 2) +  pow( yAtom[mutAtom], 2) + pow( zAtom[mutAtom], 2) ) ;
						
						C = 2. * [ g dotWith:rkb ] ;
						D = 2. * [ h dotWith:rkb ] ;
						
						s2 = pow( torusRadii[iPair], 2 ) ;
												
						Q = pow( probeRadius + radii[mutAtom], 2 ) - pow( drkb, 2 ) - s2 ;
						
						// Corrected solution strategy
						
						// Solve in terms of e
						
						a = C*C + D*D ;
						
						b = -2.*Q*C ;
						
						c = Q*Q - s2*D*D ;
						
						disc = b*b - 4.*a*c ;
						
						if( disc <= 0.) continue ;  // Next mut neighbor
						
						e1 = (-b + sqrt(disc))/(2.*a) ;
						e2 = (-b - sqrt(disc))/(2.*a) ;
						
						f1 = sqrt(s2 - e1*e1 ) ;
						f2 = sqrt(s2 - e2*e2 ) ;
					
						
						MMVector3 *probe1Pos, *probe2Pos ;
						
						dx1 = [ pairBases[iPair] X ] + e1 * [ g X ] + f1 * [ h X ] ;
						dy1 = [ pairBases[iPair] Y ] + e1 * [ g Y ] + f1 * [ h Y ] ;
						dz1 = [ pairBases[iPair] Z ] + e1 * [ g Z ] + f1 * [ h Z ] ;
						
						dx2 = [ pairBases[iPair] X ] + e1 * [ g X ] - f1 * [ h X ] ;
						dy2 = [ pairBases[iPair] Y ] + e1 * [ g Y ] - f1 * [ h Y ] ;
						dz2 = [ pairBases[iPair] Z ] + e1 * [ g Z ] - f1 * [ h Z ] ;
						
						dk1 = sqrt( pow( dx1 - xAtom[ mutAtom ], 2 ) + pow( dy1 - yAtom[ mutAtom ], 2 ) +
							      pow( dz1 - zAtom[ mutAtom ], 2 ) ) ;
								  
						dk2 = sqrt( pow( dx2 - xAtom[ mutAtom ], 2 ) + pow( dy2 - yAtom[ mutAtom ], 2 ) +
							      pow( dz2 - zAtom[ mutAtom ], 2 ) ) ;
								  
						if( fabs( dk1 - (probeRadius + radii[mutAtom] )) < fabs( dk2 - (probeRadius + radii[mutAtom] )) )
							{
								probe1Pos = [ [ [ MMVector3 alloc ] initX:dx1 Y:dy1 Z:dz1 ] autorelease ] ;
								//dx = dx1 ;
								//dy = dy1 ; 
								//dz = dz1 ;
							}
						else
							{
								probe1Pos = [ [ [ MMVector3 alloc ] initX:dx2 Y:dy2 Z:dz2 ] autorelease ] ;
								//dx = dx2 ;
								//dy = dy2 ; 
								//dz = dz2 ;
								
							}
					
						
						// AS A CHECK - distance of probe to each atom
						/*
						double di, dj, dk ;
						
						di = sqrt( pow( dx - xAtom[ atomPairs[iPair][0] ], 2 ) + pow( dy - yAtom[ atomPairs[iPair][0] ], 2 ) +
							      pow( dz - zAtom[ atomPairs[iPair][0] ], 2 ) ) ;
						
						dj = sqrt( pow( dx - xAtom[ atomPairs[iPair][1] ], 2 ) + pow( dy - yAtom[ atomPairs[iPair][1] ], 2 ) +
							      pow( dz - zAtom[ atomPairs[iPair][1] ], 2 ) ) ;
						
						dk = sqrt( pow( dx - xAtom[ mutAtom ], 2 ) + pow( dy - yAtom[ mutAtom ], 2 ) +
							      pow( dz - zAtom[ mutAtom ], 2 ) ) ;
								  
						printf( "ERROR DIST TO I = %f\n", fabs( di - (probeRadius + Ri )) ) ;
						printf( "ERROR DIST TO J = %f\n", fabs( dj - (probeRadius + Rj )) ) ;
						printf( "ERROR DIST TO K = %f\n", fabs( dk - (probeRadius + radii[mutAtom] )) ) ;
						
		
						if( fabs( di - (probeRadius + Ri )) > 0.0001 ) printf( "WARNING: PROBE IS %f FROM ATOM I!\n", fabs( di - (probeRadius + Ri )) ) ;								  
						if( fabs( dj - (probeRadius + Rj )) > 0.0001 ) printf( "WARNING: PROBE IS %f FROM ATOM J!\n", fabs( dj - (probeRadius + Rj )) ) ;
								  						  
						if( fabs( dk - (probeRadius + radii[mutAtom] )) > 0.0001 ) 
							{
								printf( "WARNING: PROBE IS %f FROM ATOM K!\n", fabs( dk - (probeRadius + radii[mutAtom] )) ) ;
							}
						
						*/
						
						dx1 = [ pairBases[iPair] X ] + e2 * [ g X ] + f2 * [ h X ] ;
						dy1 = [ pairBases[iPair] Y ] + e2 * [ g Y ] + f2 * [ h Y ] ;
						dz1 = [ pairBases[iPair] Z ] + e2 * [ g Z ] + f2 * [ h Z ] ;
						
						dx2 = [ pairBases[iPair] X ] + e2 * [ g X ] - f2 * [ h X ] ;
						dy2 = [ pairBases[iPair] Y ] + e2 * [ g Y ] - f2 * [ h Y ] ;
						dz2 = [ pairBases[iPair] Z ] + e2 * [ g Z ] - f2 * [ h Z ] ;

						dk1 = sqrt( pow( dx1 - xAtom[ mutAtom ], 2 ) + pow( dy1 - yAtom[ mutAtom ], 2 ) +
							      pow( dz1 - zAtom[ mutAtom ], 2 ) ) ;
								  
						dk2 = sqrt( pow( dx2 - xAtom[ mutAtom ], 2 ) + pow( dy2 - yAtom[ mutAtom ], 2 ) +
							      pow( dz2 - zAtom[ mutAtom ], 2 ) ) ;
								  
						if( fabs( dk1 - (probeRadius + radii[mutAtom] )) < fabs( dk2 - (probeRadius + radii[mutAtom] )) )
							{
								probe2Pos = [ [ [ MMVector3 alloc ] initX:dx1 Y:dy1 Z:dz1 ] autorelease ] ;
								//dx = dx1 ;
								//dy = dy1 ; 
								//dz = dz1 ;
								
							}
						else
							{
								probe2Pos = [ [ [ MMVector3 alloc ] initX:dx2 Y:dy2 Z:dz2 ] autorelease ] ;
								//dx = dx2 ;
								//dy = dy2 ; 
								//dz = dz2 ;
								
							}
					
						
						/*
						di = sqrt( pow( dx - xAtom[ atomPairs[iPair][0] ], 2 ) + pow( dy - yAtom[ atomPairs[iPair][0] ], 2 ) +
							      pow( dz - zAtom[ atomPairs[iPair][0] ], 2 ) ) ;
						
						dj = sqrt( pow( dx - xAtom[ atomPairs[iPair][1] ], 2 ) + pow( dy - yAtom[ atomPairs[iPair][1] ], 2 ) +
							      pow( dz - zAtom[ atomPairs[iPair][1] ], 2 ) ) ;
						
						dk = sqrt( pow( dx - xAtom[ mutAtom ], 2 ) + pow( dy - yAtom[ mutAtom ], 2 ) +
							      pow( dz - zAtom[ mutAtom ], 2 ) ) ;
								  
						printf( "ERROR DIST TO I = %f\n", fabs( di - (probeRadius + Ri )) ) ;
						printf( "ERROR DIST TO J = %f\n", fabs( dj - (probeRadius + Rj )) ) ;
						printf( "ERROR DIST TO K = %f\n", fabs( dk - (probeRadius + radii[mutAtom] )) ) ;
						
		
						if( fabs( di - (probeRadius + Ri )) > 0.0001 ) printf( "WARNING: PROBE IS %f FROM ATOM I!\n", fabs( di - (probeRadius + Ri )) ) ;								  
						if( fabs( dj - (probeRadius + Rj )) > 0.0001 ) printf( "WARNING: PROBE IS %f FROM ATOM J!\n", fabs( dj - (probeRadius + Rj )) ) ;
								  						  
						if( fabs( dk - (probeRadius + radii[mutAtom] )) > 0.0001 ) 
							{
								printf( "WARNING: PROBE IS %f FROM ATOM K!\n", fabs( dk - (probeRadius + radii[mutAtom] )) ) ;
							}
						*/
						
						// Check for bumps, add to probe list if clear
						
						if( [ self checkForBumpAtProbePosition:probe1Pos andTolerance:tol withAtomI:atomPairs[iPair][0]
								J:atomPairs[iPair][1] K:mutAtom ] == NO )
							{
								// Add probe
								
								SMProbe *probe1 = [ [ SMProbe alloc ] initWithPosition:probe1Pos andAtomI:atomPairs[iPair][0]
									J:atomPairs[iPair][1] K:mutAtom ] ;
									
								[ atomToProbes[ atomPairs[iPair][0] ] addObject:probe1 ] ;
								[ atomToProbes[ atomPairs[iPair][1] ] addObject:probe1 ] ;
								[ atomToProbes[ mutAtom ] addObject:probe1 ] ;
								
								[ self addProbe:probe1 forTripleI:atomPairs[iPair][0] J:atomPairs[iPair][1] K:mutAtom ] ;
								
								++probeCount ;
							}
							
							
						if( [ self checkForBumpAtProbePosition:probe2Pos andTolerance:tol withAtomI:atomPairs[iPair][0]
								J:atomPairs[iPair][1] K:mutAtom ] == NO )
							{
								// Add probe
								
								SMProbe *probe2 = [ [ SMProbe alloc ] initWithPosition:probe2Pos andAtomI:atomPairs[iPair][0]
									J:atomPairs[iPair][1] K:mutAtom ] ;
									
								[ atomToProbes[ atomPairs[iPair][0] ] addObject:probe2 ] ;
								[ atomToProbes[ atomPairs[iPair][1] ] addObject:probe2 ] ;
								[ atomToProbes[ mutAtom ] addObject:probe2 ] ;
								
								[ self addProbe:probe2 forTripleI:atomPairs[iPair][0] J:atomPairs[iPair][1] K:mutAtom ] ;
								
								++probeCount ;
							}
							
					}
					
				[ localPool release ] ;
			}
			
		// Guess that's all
		
		printf( "Assigned %d probe positions", probeCount ) ;
		
		return ;
	}
			
		
- (void) assignProbePositionsVer2UsingContactTolerance:(double)tol generateBlockingProbes:(BOOL)generateBlockingProbes usingSolventStart:(int)solventStartIndex
	{
		// Probes will be generated for unique triples of atoms i < j < k
		//
		
		// This uses the algorithm from SMART - it looks like I may be missing some probes with the method above 
		
		int iPair, i, probeCount, blockingProbeCount ; ;
		
		FILE *blockingProbeFile ;
		NSMutableArray *blockingProbes ;
		
		probeCount = 0 ;
		
		blockingProbeCount = 0 ;
				
		if( generateBlockingProbes == YES )
			{
				blockingProbes = [ [ NSMutableArray alloc ] initWithCapacity:100 ] ;
			}
		
		// Allocate atomToProbes array
		
		atomToProbes = (NSMutableArray **) malloc( nAtoms * sizeof( NSMutableArray * ) ) ;
		
		for( i = 0 ; i < nAtoms ; ++i )
			{
				atomToProbes[i] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
			}
		
		// Allocate room for iteratomic axis info
		
		pairAxes = (MMVector3 **) malloc( nAtomPairs * sizeof( MMVector3 *) ) ;
		pairBases = (MMVector3 **) malloc( nAtomPairs * sizeof( MMVector3 *) ) ;
		
		torusRadii = ( double * ) malloc( nAtomPairs * sizeof( double ) ) ;

		
		for( iPair = 0 ; iPair < nAtomPairs ; ++iPair )
			{
				NSAutoreleasePool * localPool = [[NSAutoreleasePool alloc] init];
				
				// Compute axis from 1 to 2
				
				double dx, dy, dz, d, cosTheta, Ri, Rj ;
				int iAtom, jAtom ;
				
				iAtom = atomPairs[iPair][0] ;
				jAtom = atomPairs[iPair][1] ;
				
				dx = xAtom[ jAtom ] - xAtom[ iAtom ] ;
				dy = yAtom[ jAtom ] - yAtom[ iAtom ] ;
				dz = zAtom[ jAtom ] - zAtom[ iAtom ] ;
				
				d = sqrt( dx*dx + dy*dy + dz*dz ) ;
				
				if( d == 0. || buriedAtom[iAtom] == YES ||  buriedAtom[jAtom] == YES) // Atoms co-incide, or one buried
					{
						pairAxes[iPair] = nil ;
						pairBases[iPair] = nil ;
						continue ;
					}
				
				pairAxes[iPair] = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
				[ pairAxes[iPair] scaleBy:(1./d) ] ;
				
				Ri = radii[ atomPairs[iPair][0] ] ;
				Rj = radii[ atomPairs[iPair][1] ] ;
				
				cosTheta = ( d*d + pow( probeRadius + Ri, 2) - pow( probeRadius + Rj, 2) ) / (2. * (probeRadius + Ri) * d) ;
				
				// Compute base point of probe rotation
				
				dx = xAtom[ iAtom ] + (Ri + probeRadius)*cosTheta*[ pairAxes[iPair] X ] ;
				dy = yAtom[ iAtom ] + (Ri + probeRadius)*cosTheta*[ pairAxes[iPair] Y ] ;
				dz = zAtom[ iAtom ] + (Ri + probeRadius)*cosTheta*[ pairAxes[iPair] Z ] ;
				
				pairBases[iPair] = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ]  ;
				
				torusRadii[iPair] = (Ri + probeRadius) * sqrt( 1. - pow( cosTheta, 2 ) ) ;

				
				// For all mutual probes (with acceptable index) find probe positions that do not intersect other atoms
				
				int j, mutAtom ;
				
				// To check for free torus
				
				
				int nTest ;
				
				nTest = 0 ;
				
				for( j = 0 ; j < nMutualNeighbors[iPair] ; ++j )
					{
						mutAtom = mutualNeighbors[iPair][j] ;
						
						
						// NOT HERE! SCREWS UP FREE TORUS TEST!
						//if( mutAtom < atomPairs[iPair][1] ) continue ;
						
						
						//printf( "AT I, J, K = %d %d %d \n", atomPairs[iPair][0], atomPairs[iPair][1], mutAtom ) ;

						
						// Generate two probe positions, check for bump
						
						double Rk ;
						
						 Rk =  radii[ mutAtom ] ;
						 
						 MMVector3 *toK, *v, *w, *probe1Pos, *probe2Pos ;
						 
						 toK = [ [ [ MMVector3 alloc ] initX:(xAtom[mutAtom] - xAtom[iAtom])
															Y:(yAtom[mutAtom] - yAtom[iAtom])
															Z:(zAtom[mutAtom] - zAtom[iAtom]) ] autorelease ] ;
															
						if( ! ( v = [ [ MMVector3 alloc ] initAlong:toK perpTo:pairAxes[iPair] ] ) )
							{
								continue ;
							}
							
						[ v autorelease ] ;
						
						
						w = [ [ [ MMVector3 alloc ] initByCrossing:pairAxes[iPair] and:v ] autorelease ] ;
						
						// Local coordinates for mutual atom
						
						 toK = [ [ [ MMVector3 alloc ]  initX:(xAtom[mutAtom] - [ pairBases[iPair] X ])
															Y:(yAtom[mutAtom] - [ pairBases[iPair] Y ])
															Z:(zAtom[mutAtom] - [ pairBases[iPair] Z ]) ] autorelease ] ;
						
						double ku, kv ;
						
						ku = [ toK dotWith:pairAxes[iPair] ] ;
						kv = [ toK dotWith:v ] ;
						
						double dbk ;
						
						dbk = [ toK length ] ;
						
						double alpha[2], beta[2] ;
						
						alpha [0] = ( dbk*dbk + torusRadii[iPair]*torusRadii[iPair] - 
								(Rk + probeRadius)*(Rk + probeRadius) ) / (2.*torusRadii[iPair]*kv) ;
						
						 if(fabs(alpha [0]) > 1.) 
							{
								// Skip, but make test for free torus first! Compute probe pos with alpha = 1
								
								MMVector3 *t ;
								
								 t = [ [ MMVector3 alloc ] initX:( [ pairBases[iPair] X ] + torusRadii[iPair]*[ v X ] )
											Y:( [ pairBases[iPair] Y ] + torusRadii[iPair]*[ v Y ] )
											Z:( [ pairBases[iPair] Z ] + torusRadii[iPair]*[ v Z ] ) ] ;
											
								double dTest ;
								
								dTest = pow( ( [ t X ] - xAtom[mutAtom] ), 2 ) +
										pow( ( [ t Y ] - yAtom[mutAtom] ), 2 ) +
										pow( ( [ t Z ] - zAtom[mutAtom] ), 2 ) ;
										
								[ t release ] ;
										
								if( dTest < ( Rk + probeRadius )*( Rk + probeRadius ) )
									{
										++nTest ;
									}
								
								continue ;
							}
							
						++nTest ;
						
						if( mutAtom < atomPairs[iPair][1] ) continue ;
							
						 beta [0] = sqrt( 1. - alpha[0]*alpha[0] ) ;

						 alpha [1] = alpha [0] ;
						 beta [1] =  - beta [0] ;
						 
						 dx = [ pairBases[iPair] X ] + torusRadii[iPair]*( alpha[0]*[v X] + beta[0]*[w X]) ;
						 dy = [ pairBases[iPair] Y ] + torusRadii[iPair]*( alpha[0]*[v Y] + beta[0]*[w Y]) ;
						 dz = [ pairBases[iPair] Z ] + torusRadii[iPair]*( alpha[0]*[v Z] + beta[0]*[w Z]) ;
						
						probe1Pos = [ [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] autorelease ] ;
						
						 dx = [ pairBases[iPair] X ] + torusRadii[iPair]*( alpha[1]*[v X] + beta[1]*[w X]) ;
						 dy = [ pairBases[iPair] Y ] + torusRadii[iPair]*( alpha[1]*[v Y] + beta[1]*[w Y]) ;
						 dz = [ pairBases[iPair] Z ] + torusRadii[iPair]*( alpha[1]*[v Z] + beta[1]*[w Z]) ;
						
						probe2Pos = [ [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] autorelease ] ;
						
						// Check for bumps, add to probe list if clear
						
						if( [ self checkForBumpAtProbePosition:probe1Pos andTolerance:tol withAtomI:iAtom
								J:jAtom K:mutAtom ] == NO )
							{
								// Add probe
								
								SMProbe *probe1 = [ [ SMProbe alloc ] initWithPosition:probe1Pos andAtomI:iAtom
									J:jAtom K:mutAtom ] ;
									
								[ atomToProbes[ iAtom ] addObject:probe1 ] ;
								[ atomToProbes[ jAtom ] addObject:probe1 ] ;
								[ atomToProbes[ mutAtom ] addObject:probe1 ] ;
								
								[ self addProbe:probe1 forTripleI:iAtom J:jAtom K:mutAtom ] ;
								
								++probeCount ;
								
								if( generateBlockingProbes == YES )
									{
										if( iAtom >= solventStartIndex || jAtom >= solventStartIndex || mutAtom >= solventStartIndex )
											{
												[ blockingProbes addObject:probe1 ] ;
												++blockingProbeCount ;
											}
											
										
									}
								
								//printf( "%f %f %f\t%d %d %d\n", [ probe1Pos X ], [ probe1Pos Y ], [ probe1Pos Z ], iAtom, jAtom, mutAtom ) ;
							}
							
							
						if( [ self checkForBumpAtProbePosition:probe2Pos andTolerance:tol withAtomI:iAtom
								J:jAtom K:mutAtom ] == NO )
							{
								// Add probe
								
								SMProbe *probe2 = [ [ SMProbe alloc ] initWithPosition:probe2Pos andAtomI:iAtom
									J:jAtom K:mutAtom ] ;
									
								[ atomToProbes[ iAtom ] addObject:probe2 ] ;
								[ atomToProbes[ jAtom ] addObject:probe2 ] ;
								[ atomToProbes[ mutAtom ] addObject:probe2 ] ;
								
								[ self addProbe:probe2 forTripleI:iAtom J:jAtom K:mutAtom ] ;
								
								++probeCount ;
								
								if( generateBlockingProbes == YES )
									{
										if( iAtom >= solventStartIndex || jAtom >= solventStartIndex || mutAtom >= solventStartIndex )
											{
												[ blockingProbes addObject:probe2 ] ;
												++blockingProbeCount ;
											}
											
										
									}
								
								
								//printf( "%f %f %f\t%d %d %d\n", [ probe2Pos X ], [ probe2Pos Y ], [ probe2Pos Z ], iAtom, jAtom, mutAtom ) ;
							}
							
					}
					
				if( nTest == 0 )
					{
						// This is a free torus
						
						++nFreeTori ;
						
						freeTorus[iPair] = YES ;
						
					}

						
					
				[ localPool release ] ;
			}
			
		// Guess that's all
		
		printf( "\nGenerated %d probe positions and identified %d free tori\n", probeCount, nFreeTori ) ;
		
		if( generateBlockingProbes == YES )
			{
				printf( "\nGenerated %d blocking probe positions, which will be written to blockingProbes.mol2\n", blockingProbeCount ) ;
				
				blockingProbeFile = fopen( "blockingProbes.mol2", "w" ) ;
				
				fprintf( blockingProbeFile, "@<TRIPOS>MOLECULE\nBlockingProbes\n%d 0 1 0 0\nSMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n", blockingProbeCount ) ;
				
				SMProbe *nextProbe ;
				
				for( i = 0 ; i < blockingProbeCount ; ++i )
					{
						nextProbe = [ blockingProbes objectAtIndex:i ] ;
						
						fprintf( blockingProbeFile, "%7d DU      %10.5f%10.5f%10.5f Du        1 BLK1        0.000\n", i+1, [ nextProbe X ], [ nextProbe Y ], [ nextProbe Z ] ) ;
					}
				
				fprintf( blockingProbeFile, "@<TRIPOS>BOND\n@<TRIPOS>SUBSTRUCTURE\n     1 BLK1        1 GROUP\n" ) ;
				
				[ blockingProbes release ] ;
				
				fclose( blockingProbeFile ) ;
			}
		
		// We can deallocate the atom grid at this point
		
		int iG, jG, kG ;
		
		for( iG = 0 ; iG < nGridX ; ++iG )
			{
				
				for( jG = 0 ; jG < nGridY ; ++jG )
					{
						for( kG = 0 ; kG < nGridZ ; ++kG )
							{
								free( gridToAtom[iG][jG][kG]  ) ;
							}
					
						free( gridToAtom[iG][jG] ) ;
						free( gridToAtomAlloc[iG][jG] ) ;
						free( nAtomsForGrid[iG][jG] ) ;
						
					}
					
				free( gridToAtom[iG]  ) ;
				free( gridToAtomAlloc[iG]  ) ;
				free( nAtomsForGrid[iG]  ) ;
					
			}

		
		
		return ;
	}
			
						
						
- (BOOL) checkForBumpAtProbePosition:(MMVector3 *)pos andTolerance:(double)tol withAtomI:(int)iAtom J:(int)jAtom K:(int)kAtom
	{
		double xP, yP, zP ;
		int iG, jG, kG, iG2, jG2, kG2, deltaI, deltaJ, deltaK, j, a ;
		
		xP = [ pos X ] ;
		yP = [ pos Y ] ;
		zP = [ pos Z ] ;
		
		
		iG = floor( ( xP - xMin ) / gridSpacing ) ;
		jG = floor( ( yP - yMin ) / gridSpacing ) ;
		kG = floor( ( zP - zMin ) / gridSpacing ) ;
		
		for( deltaI = -1 ; deltaI <= 1 ; ++deltaI )
			{
				iG2 = iG + deltaI ;
				if( iG2 < 0 || iG2 >= nGridX ) continue ;
				
				for( deltaJ = -1 ; deltaJ <= 1 ; ++deltaJ )
					{
						jG2 = jG + deltaJ ;
						if( jG2 < 0 || jG2 >= nGridY ) continue ;
					
						for( deltaK = -1 ; deltaK <= 1 ; ++deltaK )
							{
								kG2 = kG + deltaK ;
								if( kG2 < 0 || kG2 >= nGridZ ) continue ;
								
								int *atomPtr ;
								double d, dx, dy, dz ;
								
								atomPtr = gridToAtom[iG2][jG2][kG2] ;
								
								for( j = 0 ; j < nAtomsForGrid[iG2][jG2][kG2] ; ++j )
									{
										a = gridToAtom[iG2][jG2][kG2][j] ;
										
										if( a == iAtom || a == jAtom || a == kAtom ) continue ;
										
										dx = xAtom[a] - xP ;
										dy = yAtom[a] - yP ;
										dz = zAtom[a] - zP ;
										
										d = sqrt( dx*dx + dy*dy + dz*dz ) ;
										
										if( d < (probeRadius + radii[a] - tol) )
											{
												return YES ;
											}
									}
									
							}
					}
			}
			
		// Have completed all bump tests without a bump 
											
		return NO ;
	}

						
- (BOOL) addProbe:(SMProbe *)p forTripleI:(int)atomI J:(int)atomJ K:(int)atomK
	{
		// Make hash
		
		NSString *key = [ NSString stringWithFormat:@"%d_%d_%d",atomI,atomJ,atomK ] ;
		
		NSMutableArray *probeArray ;
		
		if( ( probeArray = [ atomTripleToProbes objectForKey:key ] ) )
			{
				[ probeArray addObject:p ] ;
				
				if( [ probeArray count ] > 2 )
					{
						if( warningLevel <= 1 )
							{
								printf( "WARNING: ATOM TRIPLE %d %d %d SUPPORTS %d PROBES!\n", 
									atomI, atomJ, atomK, (int)[ probeArray count ] ) ;
							}
					}
			}
		else
			{
				probeArray = [ NSMutableArray arrayWithCapacity:2 ] ;
				
				[ probeArray addObject:p ] ;
				
				[ atomTripleToProbes setObject:probeArray forKey:key ] ;
			}
			
		return YES ;
	}
		
- (NSArray *) getProbesForTripleI:(int)atomI J:(int)atomJ K:(int)atomK
	{
		// Make hash
		
		NSString *key = [ NSString stringWithFormat:@"%d_%d_%d",atomI,atomJ,atomK ] ;
		
		return [ atomTripleToProbes objectForKey:key ] ;
	}
		

- (void) generateToriUsingSkipWidth:(double)skipW
	{
		// This method finds all torus sections. Strategy:
		//		1) For each atom pair, collect all mutual neighbors, and check for probes. This will involve sorting atom 
		//			indices to make a unique key to present to the probe dictionary
		//		2) Orient each probe as involved in a "right" or "left" collision with the mutual neighbor. Right probes rotate
		//			counterclockwise with respect to the I-J axis into left probes
		//		3) Sort the probes by angle about the interatomic axis. This is started from the first right probe in the list. 
		//			If the probe at maximum angle is also a right probe, the angles are reset with that probe as reference. 
		//		4) Starting from the first right probe, collect any adjacent right-probes at larger angle. These should all be 
		//			"close" to the starting probe. Then collect the first left probe, and any adjacent left probes. The resulting 
		//			right- and left-probe clusters delimit a torus section. 
				
		int iPair, iMut, mutAtom, iAtom, jAtom, sort1, sort2, sort3 ;
		
		NSMutableArray *probesForTorus = [ NSMutableArray arrayWithCapacity:10 ] ;
		SMProbe *nextProbe ;
		NSEnumerator *probeEnumerator ;
		
		
		for( iPair = 0 ; iPair < nAtomPairs ; ++iPair )
			{
				NSAutoreleasePool *localPool = [ [ NSAutoreleasePool alloc ] init ] ;
				
				// May already be marked as a free torus pair
				
				if( freeTorus[iPair] == YES ) continue ;
				
				iAtom = atomPairs[iPair][0] ;
				jAtom = atomPairs[iPair][1] ;
				
				if( nMutualNeighbors[iPair] == 0 )
					{
						if( buriedAtom[ iAtom ] == NO && buriedAtom[ jAtom ] == NO )
							{
								++nFreeTori ;
								
								freeTorus[iPair] = YES ;
							}
							
						continue ;
					}
					
				[ probesForTorus removeAllObjects ] ;
				
				
				for( iMut=0 ; iMut < nMutualNeighbors[iPair] ; ++iMut )
					{
						mutAtom = mutualNeighbors[iPair][iMut] ;
						
						if( mutAtom < iAtom )
							{
								sort1 = mutAtom ;
								sort2 = iAtom ;
								sort3 = jAtom ;
							}
						else if( mutAtom < jAtom )
							{
								sort1 = iAtom ;
								sort2 = mutAtom ;
								sort3 = jAtom ;
							}
						else
							{
								sort1 = iAtom ;
								sort2 = jAtom ;
								sort3 = mutAtom ;
							}
						
						NSArray *nextProbes ;
						
						if( ( nextProbes = [ self getProbesForTripleI:sort1 J:sort2 K:sort3 ] ) )
							{
								probeEnumerator = [ nextProbes objectEnumerator ] ;
								
								while( ( nextProbe = [ probeEnumerator nextObject ] ) )
									{
										[ nextProbe setMutualAtom:mutAtom ] ;
										[ probesForTorus addObject:nextProbe ] ;
									}
							}
							
					}
					
				// Check for buried torus, or free torus
				
				if( [ probesForTorus count ] == 0 ) continue ;
				
				// If only one probe, probable trapped probe
				
				if( [ probesForTorus count ] == 1 )
					{
						printf( "ONLY ONE PROBE FOR ATOMS %d - %d ; PROBABLE TRAPPED PROBE ; Continuing \n", iAtom, jAtom ) ;
						continue ;
					}
				
				// Assign probes as left or right collision
				
				
				probeEnumerator = [ probesForTorus objectEnumerator ] ;
				
				SMProbe *firstRightProbe ;
				
				firstRightProbe = nil ;
				
				while( ( nextProbe = [ probeEnumerator nextObject ] ) )
					{
						[ self assignOrientationForProbe:nextProbe usingAtomI:iAtom J:jAtom K:[ nextProbe mutualAtom ] ] ;
						
						if( [ nextProbe rightProbe ] == YES && [ nextProbe leftProbe ] == NO )
							{
								firstRightProbe = nextProbe ;
							}
					}
					
				// Find first right-probe in list
				
				if( ! firstRightProbe )
					{
						printf( "IN FUNCTION generateTori, ATOM PAIR %d - %d SUPPORTS %d PROBES BUT NO PROBE WITH RIGHT-HAND COLLISION - Exit!\n",
							iAtom, jAtom, (int)[ probesForTorus count ] ) ;
							
						exit(1) ;
					}
					
				// Set angle for each probe
				
				MMVector3 *reference ;
				
				reference = [ [ MMVector3 alloc ] initX:([ firstRightProbe X ] - [ pairBases[iPair] X ])
					Y:([ firstRightProbe Y ] - [ pairBases[iPair] Y ])
					Z:([ firstRightProbe Z ] - [ pairBases[iPair] Z ]) ] ;
					
				[ reference autorelease ] ;
					
				[ reference normalize ] ;
								
				probeEnumerator = [ probesForTorus objectEnumerator ] ;
				
				MMVector3 *probeDir ;
				double ang ;
				
				probeDir = [ [ [ MMVector3 alloc ] initX:0. Y:0. Z:0. ] autorelease ] ;
				
				while( ( nextProbe = [ probeEnumerator nextObject ] ) )
					{
						[ probeDir setX:( [ nextProbe X ] - [ pairBases[iPair] X ] ) ] ;
						[ probeDir setY:( [ nextProbe Y ] - [ pairBases[iPair] Y ] ) ] ;
						[ probeDir setZ:( [ nextProbe Z ] - [ pairBases[iPair] Z ] ) ] ;
						
						[ probeDir normalize ] ;
						
						double angDot ;
						
						angDot = [ reference dotWith:probeDir ] ;
						
						if( fabs(angDot) > 1. )
							{
								angDot = angDot/fabs(angDot) ;
							}
						
						ang = acos( angDot ) ;
						
						MMVector3 *crossPdct = [ [ [ MMVector3 alloc ] initByCrossing:probeDir and:reference ] autorelease ] ;
						
						if( [ pairAxes[iPair] dotWith:crossPdct ] < -1.e-10 )
							{
								ang = 2.*acos(-1.) - ang ;
								
								// We don't expect a "non-free" torus to have probe with angle too close to 2 Pi.
								
								if( fabs( ang - 2.*acos(-1.) ) < 1.e-6 )
									{
										ang = 0. ;
									}
							}
							
						[ nextProbe setAngle:ang ] ;
						
					}
					
				// Sort by angle
				
				[ probesForTorus sortUsingSelector:@selector(compareAngles:) ] ;
				
				// Check if last probe is right probe - if so, reset reference to that probe
				
				if( [ [ probesForTorus lastObject ] rightProbe ] == YES )
					{
						[ reference  setX:([ firstRightProbe X ] - [ pairBases[iPair] X ]) ] ;
						[ reference  setY:([ firstRightProbe Y ] - [ pairBases[iPair] Y ]) ] ;
						[ reference  setZ:([ firstRightProbe Z ] - [ pairBases[iPair] Z ]) ] ;
							
						[ reference normalize ] ;
										
						probeEnumerator = [ probesForTorus objectEnumerator ] ;
												
						while( ( nextProbe = [ probeEnumerator nextObject ] ) )
							{
								[ probeDir setX:( [ nextProbe X ] - [ pairBases[iPair] X ] ) ] ;
								[ probeDir setY:( [ nextProbe Y ] - [ pairBases[iPair] Y ] ) ] ;
								[ probeDir setZ:( [ nextProbe Z ] - [ pairBases[iPair] Z ] ) ] ;
								
								[ probeDir normalize ] ;
								
								double angDot = [ reference dotWith:probeDir ] ;
								
								if( fabs(angDot) > 1. )
									{
										angDot = angDot/fabs(angDot) ;
									}
								
								ang = acos( angDot ) ;
								
								MMVector3 *crossPdct = [ [ [ MMVector3 alloc ] initByCrossing:probeDir and:reference ] autorelease ] ;
								
								if( [ pairAxes[iPair] dotWith:crossPdct ] <  -1.e-10 )
									{
										ang = 2.*acos(-1.) - ang ;
										
										// We don't expect a "non-free" torus to have probe with angle too close to 2 Pi.
								
										if( fabs( ang - 2.*acos(-1.) ) < 1.e-6 )
											{
												ang = 0. ;
											}

									}
									
								[ nextProbe setAngle:ang ] ;
								
							}
							
						// Sort by angle
						
						[ probesForTorus sortUsingSelector:@selector(compareAngles:) ] ;
								
						if( [ [ probesForTorus lastObject ] rightProbe ] == YES )
							{
								if( warningLevel <= 1 )
									{
										printf( "WARNING: COULD NOT ORIENT PROBES FOR ATOM PAIR %d - %d - Probable trapped probe - skipping torus!\n", 
											iAtom, jAtom ) ;
									}
								break ;
							}
					}
					
				//	Have probes, can now generate torus sections for this atom pair
				
				SMProbe *rightTorusProbe, *leftTorusProbe ;
				SMTorus *nextTorusSection ;
				int iProbe ;
				
				rightTorusProbe = leftTorusProbe = nil ;
				
				for( iProbe = 0 ; iProbe < [ probesForTorus count ] ; ++iProbe )
					{
						nextProbe = [ probesForTorus objectAtIndex:iProbe ] ;
						
						if( [ nextProbe rightProbe ] == YES  )
							{
								if( ! rightTorusProbe ) rightTorusProbe = nextProbe ;
								continue ;
							}
							
						
						if( [ nextProbe leftProbe ]  == YES )
							{
								if( ! rightTorusProbe ) 
									{
										if( warningLevel <= 1 )
											{
												printf( "WARNING: PROBES OUT OF LEFT-RIGHT SEQUENCE AT ATOM PAIR %d - %d - Probable trapped probe - Skipping Torus\n",
													iAtom, jAtom ) ;
											}
											
										break ;
									}
								
								leftTorusProbe = nextProbe ;
								
								// Generate torus section 
								
								
								nextTorusSection = [ [ SMTorus alloc ] initWithProbeR:rightTorusProbe probeL:leftTorusProbe atomI:iAtom atomJ:jAtom 
									free:NO usingMolecule:self usingBase:pairBases[iPair] usingTorusAxis:pairAxes[iPair] 
									usingSaddleRadius:torusRadii[iPair]  usingProbeRadius:probeRadius ] ;
									 
								// Generate torus arcs
								
								[ nextTorusSection computeContactAndReentrantArcsUsingMolecule:self andSkipWidth:skipW ] ;
								
								// Add to global list 
								
								
								if( nTori == torusAlloc )
									{
										torusAlloc += 50 ;
										
										tori = (SMTorus **) realloc( tori, torusAlloc * sizeof( SMTorus * ) ) ;
									}
									
								tori[nTori] = nextTorusSection ;
								
								++nTori ;
								
								// Add to atom-to-tori data 
								
								if( ! atomsToTori[iAtom] )
									{
										atomsToTori[iAtom] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
										atomsToCycles[iAtom] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
									}
									
								[ atomsToTori[iAtom] addObject:nextTorusSection ] ;
								
								if( ! atomsToTori[jAtom] )
									{
										atomsToTori[jAtom] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
										atomsToCycles[jAtom] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
									}
									
								[ atomsToTori[jAtom] addObject:nextTorusSection ] ;
								
								
								// Get past subsequent leftProbes
								
								while( (iProbe < [ probesForTorus count ] - 1) && 
									[ [ probesForTorus objectAtIndex:(iProbe + 1) ] leftProbe ] == YES )
									{
												++iProbe ;
									}
									
								rightTorusProbe = nil ;
								
								
							}
							
						// Next probe 
					}
					
				// Next atom pair
				
				[ localPool release ] ;
			}
			
			
		// Process free tori
		
		// Approach here is to generate a R and L probe at the same position, assign phi = 2*Pi to te L arc, and to 
		// twin the arcs
		// 
		for( iPair = 0 ; iPair < nAtomPairs ; ++iPair )
			{
				NSAutoreleasePool *localPool = [ [ NSAutoreleasePool alloc ] init ] ;
				
				iAtom = atomPairs[iPair][0] ;
				jAtom = atomPairs[iPair][1] ;
				
				if( freeTorus[iPair] == NO ) continue ;
				
				SMProbe *rightTorusProbe, *leftTorusProbe ;
				SMTorus *nextTorusSection ;
				
				// First need to generate a direction for the probes
				
				// Take cross-product of torus axis with X and Y axes, largest will be the probe direction
				
				MMVector3 *probePosX, *probePosY, *probePos ;
				
				probePosX = [ [ MMVector3 alloc ] initByCrossing:pairAxes[iPair] and:[ MMVector3 xAxis ] ] ;
				probePosY = [ [ MMVector3 alloc ] initByCrossing:pairAxes[iPair] and:[ MMVector3 yAxis ] ] ;
				
				if( [ probePosX length ] > [ probePosY length ] )
					{
						probePos = probePosX ;
					}
				else
					{
						probePos = probePosY ;
					}
					
				// Note that pobePos initially holds the displacement vector towards the probe
				
				[ probePos normalize ] ;
				
				[ probePos setX:( [ pairBases[iPair] X ] + torusRadii[iPair] * [ probePos X ] ) ] ;
				[ probePos setY:( [ pairBases[iPair] Y ] + torusRadii[iPair] * [ probePos Y ] ) ] ;
				[ probePos setZ:( [ pairBases[iPair] Z ] + torusRadii[iPair] * [ probePos Z ] ) ] ;
					
				rightTorusProbe = [ [ SMProbe alloc ] initWithPosition:probePos andAtomI:iAtom J:jAtom K:(-1) ] ;
				[ rightTorusProbe setProbeTypeRight:YES left:NO ] ;
				[ rightTorusProbe setAngle:0. ] ;
				
				leftTorusProbe = [ [ SMProbe alloc ] initWithPosition:probePos andAtomI:iAtom J:jAtom K:(-1) ] ;
				[ leftTorusProbe setProbeTypeRight:NO left:YES ] ;
				[ leftTorusProbe setAngle:(2.*acos(-1.)) ] ;
				
				nextTorusSection = [ [ SMTorus alloc ] initWithProbeR:rightTorusProbe probeL:leftTorusProbe atomI:iAtom atomJ:jAtom 
									free:YES usingMolecule:self usingBase:pairBases[iPair] usingTorusAxis:pairAxes[iPair] 
									usingSaddleRadius:torusRadii[iPair] usingProbeRadius:probeRadius ] ;
				
				// Generate torus arcs
				
				[ nextTorusSection computeContactAndReentrantArcsUsingMolecule:self andSkipWidth:skipW ] ;
				
				// IF we have a self-intersecting (and tight) free torus, ignore it!
				
				if( nextTorusSection->selfIntersection == YES && nextTorusSection->bufferIJ < skipW ) 
					{
						continue ;
					}
				
				// Add to global list 
				
				if( nTori == torusAlloc )
					{
						torusAlloc += 50 ;
						
						tori = (SMTorus **) realloc( tori, torusAlloc * sizeof( SMTorus * ) ) ;
					}
					
				tori[nTori] = nextTorusSection ;
				
				++nTori ;
				
				// Add to atom-to-tori data 
				
				if( ! atomsToTori[iAtom] )
					{
						atomsToTori[iAtom] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
						atomsToCycles[iAtom] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
					}
					
				[ atomsToTori[iAtom] addObject:nextTorusSection ] ;
				
				if( ! atomsToTori[jAtom] )
					{
						atomsToTori[jAtom] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
						atomsToCycles[jAtom] = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
					}
					
				[ atomsToTori[jAtom] addObject:nextTorusSection ] ;
				
				// Next atom pair
				
				[ localPool release ] ;
			}
								
		
		// All done!
		
		printf( "\nGenerated %d tori\n", nTori ) ;
		
		return ;
	}
	
- (void) generateContactCycles
	{
		// Generate cycles of contact arcs 
		
		BOOL working = YES ;
		
		contactCycles = [ [ NSMutableArray alloc ] initWithCapacity:nAtoms ] ;
		
		enum { initCycleI, initCycleJ, accumulateCycle } state ;
		
		state = initCycleI ;
		
		int currentTorusIndex = 0 ;
		int currentAtom ;
		
		SMTorus *currentTorus = tori[currentTorusIndex] ;
		
		// We will handle free tori separately
		
		while( currentTorusIndex < nTori && [ currentTorus freeTorus ] == YES )
			{
				++currentTorusIndex ;
				currentTorus = tori[currentTorusIndex] ;
			}
			
		// All free tori? It can happen...
		
		if( currentTorusIndex == nTori ) goto PROCESS_FREE_TORI ;
		
		SMArc *startArc, *currentArc, *closestArc ;
		SMCycle *currentCycle ;
		SMTorus *neighborTorus ;
		double closestDist, d ;

		
		
		while( working == YES )
			{
				switch (state)
					{
						case initCycleI:
							if( [ [ currentTorus contactArcI ] parentCycles ] != nil )
								{
									state = initCycleJ ;
									break ;
								}
								
							startArc = [ currentTorus contactArcI ] ;
							currentArc = startArc ;
							
							currentAtom = [ currentTorus atomI ] ;
							
							currentCycle = [ [ SMCycle alloc ] init ] ;
							
							[ currentCycle addArc:startArc forward:YES ] ;
							
							[ currentCycle addAtomWithIndex:[ currentTorus atomI ] ] ;
							
							//[ startArc addParentCycle:currentCycle ] ;
							
							state = accumulateCycle ;
							
							break ;
							
						case initCycleJ:
							if( [ [ currentTorus contactArcJ ] parentCycles ] != nil )
								{
									++currentTorusIndex ;
									
									while( currentTorusIndex < nTori && [ tori[currentTorusIndex] freeTorus ] == YES  )
										{
											++currentTorusIndex ;
										}
									
									if( currentTorusIndex == nTori ) 
										{
											working = NO ;
											break ;
										}
										
									currentTorus = tori[currentTorusIndex] ;
									
									state = initCycleI ;
									
									break ;
								}
								
							startArc = [ currentTorus contactArcJ ] ;
							currentAtom = [ currentTorus atomJ ] ;
							
							currentArc = startArc ;
							
							currentCycle = [ [ SMCycle alloc ] init ] ;
							
							[ currentCycle addArc:startArc forward:YES ] ;
							
							[ currentCycle addAtomWithIndex:[ currentTorus atomJ ] ] ;
							
							//[ startArc addParentCycle:currentCycle ] ;
							
							state = accumulateCycle ;
							
							break ;

						case accumulateCycle: 
					
							// Cycle through neighboring tori
							
							
							
							closestArc = nil ;
							closestDist = 100000000. ;
							
							NSEnumerator *torusEnumerator ;
							
							torusEnumerator = [ atomsToTori[currentAtom] objectEnumerator ] ;
							
							while( ( neighborTorus = [ torusEnumerator nextObject ] ) )
								{
									
									if( [ [ neighborTorus contactArcI ] parentCycles ] == nil ||
											[ neighborTorus contactArcI ] == startArc )
										{
											d = [ [ currentArc endPosition ] distWith:[ [ neighborTorus contactArcI ] startPosition ] ] ;
											
											if( d < closestDist )
												{
													closestDist = d ;
													closestArc = [ neighborTorus contactArcI ] ;
												}
										}
							
									if( [ [ neighborTorus contactArcJ ] parentCycles ] == nil || 
											[ neighborTorus contactArcJ ] == startArc )
										{
											d = [ [ currentArc endPosition ] distWith:[ [ neighborTorus contactArcJ ] startPosition ] ] ;
											
											if( d < closestDist )
												{
													closestDist = d ;
													closestArc = [ neighborTorus contactArcJ ] ;
												}
										}
								}
								
							// We should have an arc
							
							if( ! closestArc || closestDist > 0.01 )
								{
									printf( "ERROR IN CLOSING CONTACT CYCLE AT ATOM %d - Exit!\n", currentAtom ) ;
									exit(1) ;
								}
								
							
							if( closestArc == startArc )
								{
									// Done with this cycle 
									
									[ contactCycles addObject:currentCycle ] ;
									
									[ atomsToCycles[ currentAtom ] addObject:currentCycle ] ;
									
									[ currentCycle release ] ;
									
														
											
									
									state = initCycleI ;
								}
							else
								{
									// Add to cycle
									
									[ currentCycle addArc: closestArc forward:YES ] ;
									
									// Assign arc to cycle
									
									//[ closestArc addParentCycle:currentCycle ] ;
									
									// New current arc
									
									currentArc = closestArc ;
									
								}
									
								
								
							break ;
					
					}
			}
			
		// Now handle free tori
		
PROCESS_FREE_TORI:
		
		for( currentTorusIndex = 0 ; currentTorusIndex < nTori ; ++currentTorusIndex )
			{
				currentTorus = tori[ currentTorusIndex ] ;
				
				if( [ currentTorus freeTorus ] == NO ) continue ;
				
				// Add two contact cycles
				
				currentCycle = [ [ SMCycle alloc ] init ] ;
				
				[ currentCycle addArc:[ currentTorus contactArcI ] forward:YES ] ;
				
				[ currentCycle addAtomWithIndex:[ currentTorus atomI ] ] ;
				
				[ contactCycles addObject:currentCycle ] ;
				
				[ atomsToCycles[ [ currentTorus atomI ] ] addObject:currentCycle ] ;
				
				[ currentCycle release ] ;
				
				// Immediately subdivide into four arcs
				
				[ self subdivideArc:[ currentTorus contactArcI ] ] ;
				
				SMArc *divide1, *divide2 ;
				
				divide1 = [ [ currentCycle arcs ] objectAtIndex:0 ] ;
				divide2 = [ [ currentCycle arcs ] objectAtIndex:1 ] ;
				
				[ self subdivideArc:divide1 ] ;
				[ self subdivideArc:divide2 ] ;
				
				currentCycle = [ [ SMCycle alloc ] init ] ;
				
				[ currentCycle addArc:[ currentTorus contactArcJ ] forward:YES ] ;
				
				[ currentCycle addAtomWithIndex:[ currentTorus atomJ ] ] ;
				
				[ contactCycles addObject:currentCycle ] ;
				
				[ atomsToCycles[ [ currentTorus atomJ ] ] addObject:currentCycle ] ;
				
				[ currentCycle release ] ;
				
				// Immediately subdivide into four arcs
				
				[ self subdivideArc:[ currentTorus contactArcJ ] ] ;
				
				divide1 = [ [ currentCycle arcs ] objectAtIndex:0 ] ;
				divide2 = [ [ currentCycle arcs ] objectAtIndex:1 ] ;
				
				[ self subdivideArc:divide1 ] ;
				[ self subdivideArc:divide2 ] ;
				
			}

		
			
		printf( "\nGenerated %d contact cycles ...\n", (int)[ contactCycles count ] ) ;
		
		return ;
	}
									
		
- (void) cullContactCycles
	{
		// This method will examine each tentative contact cycle, and check for the number of "non-skipped" arcs. If the the total is less than 
		// three, then actions will be taken:
		//
		//	non-skip == 2
		//		Subdivide each of the arcs
		//	non-skip == 1
		//		Change the non-skip torus to a skipped torus
		//		Twin the reentrant arcs of the newly-skipped torus
		//		Remove the cycle from the list
		//
		//	The process will be carried out iteratively until there are no more modified tori
		//
		
		NSEnumerator *contactCycleEnumerator, *arcEnumerator ;
		NSMutableArray *unskippedArcs ;
		
		unskippedArcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
		
		SMCycle *nextCycle ;
		SMArc *nextArc ;
		
		int MNoSkip ;
		BOOL hadAKill ;
		
		
		
		hadAKill = YES ;
		
		while( hadAKill )
			{
				hadAKill = NO ;
				
				contactCycleEnumerator = [ contactCycles objectEnumerator ] ;
				
				while( ( nextCycle = [ contactCycleEnumerator nextObject ] ) )
					{
						if( [ nextCycle active ] == NO ) continue ;
						
						arcEnumerator = [ [ nextCycle arcs ] objectEnumerator ] ;
						
						[ unskippedArcs removeAllObjects ] ;
						
						MNoSkip = 0 ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								if( [ [ nextArc torusSection ] skipTorus ] == NO )
									{
										++MNoSkip ;
										[ unskippedArcs addObject:nextArc ] ;
									}
							}
							
						if( MNoSkip == 2 )
							{
								// Subdivide both non-skip arcs 
								
								arcEnumerator = [ unskippedArcs objectEnumerator ] ;
								
								while( ( nextArc = [ arcEnumerator nextObject ] ) )
									{
										[ self subdivideArc:nextArc ] ;
									}
							}
						else if( MNoSkip == 1 )
							{
								nextArc = [ unskippedArcs lastObject ] ;
								
								SMTorus *theTorusSection = [ nextArc torusSection] ;
								
								theTorusSection->skipTorus = YES ;
								
								// Twin the reentrant arcs
								
								int M ;
								
								M = [ theTorusSection->reentrantRArcs count ] ;
								
								if( M != [ theTorusSection->reentrantLArcs count ] )
									{
										printf( "FATAL ERROR - ATTEMPT TO TWIN TORUS REENTRANT ARCS FAILS - Exit!\n" ) ;
										exit(1) ;
									}
									
								SMArc *arcR, *arcL ;
								int i ;
								
								for( i = 0 ; i < M ; ++i )
									{
										arcR = [ theTorusSection->reentrantRArcs objectAtIndex:i ] ;
										arcL = [ theTorusSection->reentrantLArcs objectAtIndex:(M - 1 - i ) ] ;
										
										[ arcR setTwin:arcL ] ;
										[ arcL setTwin:arcR ] ;
									}
									
								// Set the contact arcs as "skip" arcs
								
								NSEnumerator *contactArcEnumerator ;
								SMArc *arcI, *arcJ ;
								
								contactArcEnumerator = [ theTorusSection->contactIArcs objectEnumerator ] ;

								
								while( ( arcI = [ contactArcEnumerator nextObject ] ) )
									{
										[ arcI setSkip:YES ] ;
									}
									
								contactArcEnumerator = [ theTorusSection->contactJArcs objectEnumerator ] ;

								
								while( ( arcJ = [ contactArcEnumerator nextObject ] ) )
									{
										[ arcJ setSkip:YES ] ;
									}
									
									
								[ nextCycle killCycle ] ;
								
								hadAKill = YES ;
							}
						else if( MNoSkip == 0 )
							{
								[ nextCycle killCycle ] ;
								
								hadAKill = YES ;
							}
								
								
					}
							
			}
	}
								
								
		
		
								
- (void) generateReentrantCyclesWithReentrantHandling:(BOOL)rh
	{
		// Generate cycles of reentrant arcs 
		
		// Search tori connected up to 10 neighboring atoms
		
		//NSMutableSet *cycleProbes ;
		
		int nNeighborAtoms, neighborAtoms[10] ;
		
		BOOL working = YES ;
		
		reentrantCycles = [ [ NSMutableArray alloc ] initWithCapacity:nAtoms ] ;
		
		int nSelfIntInitial = 0 ;
		int nSelfIntFinal = 0 ;
		
		//cycleProbes = [ [ NSMutableSet alloc ] initWithCapacity:10 ] ;
		
		enum { initCycleR, initCycleL, accumulateCycle } state ;
		
		state = initCycleR ;
		
		int currentTorusIndex = 0 ;
		
		
		SMTorus *currentTorus = tori[currentTorusIndex] ;
		
		// We will watch out for free tori
		
		while( currentTorusIndex < nTori && [ currentTorus freeTorus ] == YES  )
			{
				++currentTorusIndex ;
				currentTorus = tori[currentTorusIndex] ;
			}
			
		// All free tori? It can happen...
		
		if( currentTorusIndex == nTori ) return ;
		
		
		SMArc *startArc, *currentArc, *closestArc ;
		SMCycle *currentCycle ;
		SMTorus *neighborTorus, *closestTorus ;
		double closestDist, d ;
		int i, *contacts, nContacts ; 
		
		nNeighborAtoms = 2 ;
		
		neighborAtoms[0] = [ currentTorus atomI ] ;
		neighborAtoms[1] = [ currentTorus atomJ ] ;

		
		
		while( working == YES )
			{
				switch (state)
					{
						case initCycleR:
							if( [ [ currentTorus reentrantArcR ] parentCycles ] != nil )
								{
									state = initCycleL ;
									break ;
								}
								
							startArc = [ currentTorus reentrantArcR ] ;
							currentArc = startArc ;
							
							//[ cycleProbes addObject:[ currentTorus probeR ] ] ;
							
							currentCycle = [ [ SMCycle alloc ] init ] ;
							
							[ currentCycle setProbe:[ currentTorus probeR ] ] ;
							
							[ currentCycle addArc:startArc forward:YES ] ;
							
							contacts = [ [ currentTorus probeR ] contacts ] ;
							
							for( i = 0 ; i < [ [ currentTorus probeR ] nContacts ] ; ++i )
								{
									[ currentCycle addAtomWithIndex:contacts[i] ] ;
								}
							
							//[ startArc addParentCycle:currentCycle ] ;
							
							state = accumulateCycle ;
							
							break ;
							
						case initCycleL:
							if( [ [ currentTorus reentrantArcL ] parentCycles ] != nil )
								{
									++currentTorusIndex ;
									
									while( currentTorusIndex < nTori && [ tori[currentTorusIndex] freeTorus ] == YES )
										{
											++currentTorusIndex ;
										}
									
									if( currentTorusIndex == nTori ) 
										{
											working = NO ;
											break ;
										}
										
									currentTorus = tori[currentTorusIndex] ;
									
									state = initCycleR ;
									
									nNeighborAtoms = 2 ;
		
									neighborAtoms[0] = [ currentTorus atomI ] ;
									neighborAtoms[1] = [ currentTorus atomJ ] ;
									
									break ;
								}
								
							startArc = [ currentTorus reentrantArcL ] ;
							
							//[ cycleProbes addObject:[ currentTorus probeL ] ] ;
							
							currentArc = startArc ;
							
							currentCycle = [ [ SMCycle alloc ] init ] ;
							
							[ currentCycle setProbe:[ currentTorus probeL ] ] ;
							
							[ currentCycle addArc:startArc  forward:YES ] ;
							
							contacts = [ [ currentTorus probeL ] contacts ] ;
							
							for( i = 0 ; i < [ [ currentTorus probeL ] nContacts ] ; ++i )
								{
									[ currentCycle addAtomWithIndex:contacts[i] ] ;
								}
							
							
							//[ startArc addParentCycle:currentCycle ] ;
							
							state = accumulateCycle ;
							
							break ;

						case accumulateCycle: 
					
							// Cycle through neighboring tori, both atoms
							
							closestArc = nil ;
							closestTorus = nil ;
							closestDist = 100000000. ;
							
							int iNeighborAtom ;
							
							for( iNeighborAtom = 0 ; iNeighborAtom < nNeighborAtoms ; ++iNeighborAtom )
								{
							
									NSEnumerator *torusEnumerator ;
									
									torusEnumerator = [ atomsToTori[  neighborAtoms[iNeighborAtom] ] objectEnumerator ] ;
									
									while( ( neighborTorus = [ torusEnumerator nextObject ] ) )
										{
											
											if( [ [ neighborTorus reentrantArcR ] parentCycles ] == nil ||
													[ neighborTorus reentrantArcR ] == startArc )
												{
													d = [ [ currentArc endPosition ] distWith:[ [ neighborTorus reentrantArcR ] startPosition ] ] ;
													
													if( d < closestDist )
														{
															closestDist = d ;
															closestArc = [ neighborTorus reentrantArcR ] ;
															
															contacts = [ [ neighborTorus probeR ] contacts ] ;
															nContacts = [ [ neighborTorus probeR ] nContacts  ] ;
															
															closestTorus = neighborTorus ;
														}
												}
									
											if( [ [ neighborTorus reentrantArcL ] parentCycles ] == nil || 
													[ neighborTorus reentrantArcL ] == startArc )
												{
													d = [ [ currentArc endPosition ] distWith:[ [ neighborTorus reentrantArcL ] startPosition ] ] ;
													
													if( d < closestDist )
														{
															closestDist = d ;
															closestArc = [ neighborTorus reentrantArcL ] ;
															
															contacts = [ [ neighborTorus probeL ] contacts ] ;
															nContacts = [ [ neighborTorus probeL ] nContacts  ] ;
															
															closestTorus = neighborTorus ;
														}
												}
										}
								}
								
							// Expand neighbor atom list
							
							BOOL foundI, foundJ  ;
							
							foundI = NO ;
							foundJ = NO ;
							
							for( iNeighborAtom = 0 ; iNeighborAtom < nNeighborAtoms ; ++iNeighborAtom )
								{
									if( neighborAtoms[iNeighborAtom] == [ closestTorus atomI ] ) foundI = YES ;
									if( neighborAtoms[iNeighborAtom] == [ closestTorus atomJ ] ) foundJ = YES ;
								}
								
							if( foundI == NO )
								{
									if( nNeighborAtoms == 10 )
										{
											printf( "ERROR CLOSING CONTACT CYCLE WITHIN 10 TORI - Exit!\n" ) ;
											exit(1) ;
										}
										
									neighborAtoms[nNeighborAtoms] = [ closestTorus atomI ] ;
									
									++nNeighborAtoms ;
								}
								
							if( foundJ == NO )
								{
									if( nNeighborAtoms == 10 )
										{
											printf( "ERROR CLOSING CONTACT CYCLE WITHIN 10 TORI - Exit!\n" ) ;
											exit(1) ;
										}
										
									neighborAtoms[nNeighborAtoms] = [ closestTorus atomJ ] ;
									
									++nNeighborAtoms ;
								}
							
								
							// We should have an arc
							
							if( ! closestArc || closestDist > 0.01 )
								{
									printf( "ERROR IN CLOSING REENTRANT CYCLE NEAR ATOMS %d - Exit!\n", [ currentTorus atomI ],
										[ currentTorus atomJ ] ) ;
									exit(1) ;
								}
								
							
							// Assign arc to cycle
							
							//[ closestArc addParentCycle:currentCycle ] ;
							
							if( closestArc == startArc )
								{
									// Done with this cycle 
									
									//[ currentCycle setProbes:cycleProbes ] ;
									
									//[ cycleProbes removeAllObjects ] ;
									
									[ reentrantCycles addObject:currentCycle ] ;
									
									state = initCycleR ;
									
									[ currentCycle release ] ;
									
									// Add limit plane to the cycle
									// For the moment I am going to do this for all cycles, as I may need to check for 
									// overlapping probes even for cycles that do not include self-interesecting tori
									
									currentCycle->theLimitPlane = [ [ SMLimitPlane alloc ] initWithReentrantArcs:[ currentCycle arcs ] usingMolecule:self ] ;
									
									if( currentCycle->theLimitPlane->selfIntersection == YES )
										{
											currentCycle->selfIntersection = YES ;
											++nSelfIntInitial ;
											++nSelfIntFinal ;
										}
								
								}
							else
								{
									// Add to cycle
									
									[ currentCycle addArc:closestArc forward:YES ] ;
									
									for( i = 0 ; i < nContacts ; ++i )
										{
											[ currentCycle addAtomWithIndex:contacts[i] ] ;
										}
								
								
									// New current arc
									
									currentArc = closestArc ;
								}
								
							break ;
					
					}
			}
			
		printf( "\nGenerated %d reentrant cycles ...\n", (int)[ reentrantCycles count ] ) ;
		
		if( rh == NO ) return ;
		
		
		// We need some more careful testing of reentrant sections to ascertain whether or not self-intersection assignments are correct. 
		
		// Compare each probe against all other probes. We only need to check reent cyckes that are initially non-self-intersectig, AND 
		// cycles that are self-intersecting, but NOT associated with self-int torus sections.
		
		// 1) Non-self-intersecting sections: If there is a contact with another reentrant region, then switch on self-inersection flag
		// 2) Initially self-intersecting regions: If there is NO contact with another reent region, then switch off self-interection flag
		
		int nSelfIntOverturn = 0, nNonSelfIntOverturn = 0 ;
		double dist ;
		
		NSEnumerator *cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
		SMCycle *nextCycle ;
		SMProbe *cycleProbe, *compareProbe ;
		MMVector3 *disp = [ [ MMVector3 alloc ] initX:0 Y:0 Z:0 ] ;
		MMVector3 *contactPoint = [ [ MMVector3 alloc ] initX:0 Y:0 Z:0 ] ;
		MMVector3 *compareContactPoint = [ [ MMVector3 alloc ] initX:0 Y:0 Z:0 ] ;
		SMLimitPlane *nextLimitPlane, *nextCompareLimitPlane ;
		MMVector3 *nextPlaneCenter, *nextPlaneNormal, *nextComparePlaneCenter, *nextComparePlaneNormal;
		BOOL haveContact ;
		
		while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
			{
				//if( nextCycle->selfIntersection == NO ) continue ;
				
				nextLimitPlane = nextCycle->theLimitPlane ;
				
				// If this cycle is bordered by a self-intersecting torus, we can't overturn!
				
				if( nextLimitPlane->bordersSelfintersectingTorus == YES ) continue ;
				
				nextPlaneCenter = nextLimitPlane->planeCenter ;
				nextPlaneNormal = nextLimitPlane->planeNormal ;
				
				NSEnumerator *compareCycleEnumerator = [ reentrantCycles objectEnumerator ] ;
				
				SMCycle *compareCycle ;
				
				cycleProbe = [ nextCycle probe ] ;
				
				haveContact = NO ;
				
				while( ( compareCycle = [ compareCycleEnumerator nextObject ] ) )
					{
						if( compareCycle == nextCycle ) continue ;
						
						compareProbe = [ compareCycle probe ] ;
						
						nextCompareLimitPlane = compareCycle->theLimitPlane ;
						nextComparePlaneCenter = nextCompareLimitPlane->planeCenter ;
						nextComparePlaneNormal = nextCompareLimitPlane->planeNormal ;
						
						// Get contact points on both probes						
						
						[ disp setX:( [ compareProbe X ] - [ cycleProbe X ] ) ] ;
						[ disp setY:( [ compareProbe Y ] - [ cycleProbe Y ] ) ] ;
						[ disp setZ:( [ compareProbe Z ] - [ cycleProbe Z ] ) ] ;
						
						dist = [ disp length ] ;
						
						if( dist >= 2.*probeRadius ) continue ;
						
						[ contactPoint setX:( [ cycleProbe X ] + (probeRadius/dist)*[ disp X ] ) ] ;
						[ contactPoint setY:( [ cycleProbe Y ] + (probeRadius/dist)*[ disp Y ] ) ] ;
						[ contactPoint setZ:( [ cycleProbe Z ] + (probeRadius/dist)*[ disp Z ] ) ] ;
						
						[ compareContactPoint setX:( [ compareProbe X ] - (probeRadius/dist)*[ disp X ] ) ] ;
						[ compareContactPoint setY:( [ compareProbe Y ] - (probeRadius/dist)*[ disp Y ] ) ] ;
						[ compareContactPoint setZ:( [ compareProbe Z ] - (probeRadius/dist)*[ disp Z ] ) ] ;
						
						// Check if contact point is below probe limit plane
						
						[ disp setX:( [ contactPoint X ] - [ nextPlaneCenter X ] ) ] ;
						[ disp setY:( [ contactPoint Y ] - [ nextPlaneCenter Y ] ) ] ;
						[ disp setZ:( [ contactPoint Z ] - [ nextPlaneCenter Z ] ) ] ;
						
						if( [ disp dotWith:nextPlaneNormal ] >= 0. ) continue ;
						
						// Check if compare contact point is below compare probe limit plane
						
						[ disp setX:( [ compareContactPoint X ] - [ nextComparePlaneCenter X ] ) ] ;
						[ disp setY:( [ compareContactPoint Y ] - [ nextComparePlaneCenter Y ] ) ] ;
						[ disp setZ:( [ compareContactPoint Z ] - [ nextComparePlaneCenter Z ] ) ] ;
						
						if( [ disp dotWith:nextComparePlaneNormal ] >= 0. ) continue ;
						
						// Have another reentrant regions impinging on the subject region
						
						haveContact = YES ;
						break ;
					}
					
				
				if( haveContact == YES )
					{
						if( nextCycle->selfIntersection == NO )
							{
								nextCycle->selfIntersection = YES ;
								++nNonSelfIntOverturn ;
								++nSelfIntFinal ;
								/*
								// Need to adjust limit plane position, using last COMPARE contact point
								
								[ disp setX:( [ compareContactPoint X ] - [ cycleProbe X ] ) ] ;
								[ disp setY:( [ compareContactPoint Y ] - [ cycleProbe Y ] ) ] ;
								[ disp setZ:( [ compareContactPoint Z ] - [ cycleProbe Z ] ) ] ;
								
								dist = fabs( [ disp dotWith:nextPlaneNormal ] ) ;  // This will be new height of probe center over limit plane
								
								nextLimitPlane->probeHeight = dist ;
								
								// Adjust position of plane center
								
								[ nextPlaneCenter setX:( [ nextPlaneCenter X ]  - dist * [ nextPlaneNormal X ] ) ] ;
								[ nextPlaneCenter setY:( [ nextPlaneCenter Y ]  - dist * [ nextPlaneNormal Y ] ) ] ;
								[ nextPlaneCenter setZ:( [ nextPlaneCenter Z ]  - dist * [ nextPlaneNormal Z ] ) ] ;
								*/
								
							}
					}
				else
					{
						if( nextCycle->selfIntersection == YES )
							{
								nextCycle->selfIntersection = NO ;
								++nSelfIntOverturn ;
								--nSelfIntFinal ;
							}
					}
			}
			
		printf( "Self-intersection status of %d cycles was reversed:\n", (+nNonSelfIntOverturn + nSelfIntOverturn) ) ;
		printf( "\t%d self-intersecting regions converted to non-self-intersecting, %d non-self-intersecting converted to self-intersecting\n", nSelfIntOverturn, nNonSelfIntOverturn ) ;
		printf( "Initial no. self-intersecting = %d, final count = %d\n", nSelfIntInitial, nSelfIntFinal ) ;
		
		  
		
		return ;
	}
	
				

	
- (void) generateSaddleCycles 
	{
		// One saddle cycle per unskipped torus
		
		// Method is straightforward; arcs in contactI, reentrantR, contactJ and reentrantL should immediately form a cycle. 
		// Phi and Theta are added for each cycle, and a simple representation of the cycle is constructed in the phi-theta
		// plane. 
		
		int i ;
		double thetaMax, phiMax ;
		
		saddleCycles = [ [ NSMutableArray alloc ] initWithCapacity:nTori ] ;
		
		for( i = 0 ; i < nTori ; ++i )
			{
				// Note - we need to generate a cycle for a skipped torus, even if it will not be triangulated. 
				// This is to maintain integrity of vertices involved in skipped tori
				//if( [ tori[i] skipTorus ] == YES ) continue ;
				
				
				SMCycle *nextCycle ;
				
					
				nextCycle = [ [ SMCycle alloc ] init ] ;
				
				[ nextCycle addAtomWithIndex:[ tori[i] atomI ] ] ;
				[ nextCycle addAtomWithIndex:[ tori[i] atomJ ] ] ;
				
				thetaMax = 0. ;
				phiMax = 0. ;
				
				NSEnumerator *arcEnumerator ;
				
				// Contact I
				
				arcEnumerator = [ [ tori[i] contactIArcs ] objectEnumerator ] ;
				
				SMArc *nextArc ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						[ nextCycle addArc:nextArc forward:YES ] ;
						
						if( nextArc == [ [ tori[i] contactIArcs ] objectAtIndex:0 ] )
							{
								[ nextArc setPhiStart:[ tori[i] phiAngleSubtended ] ] ;
							}
						
						if( nextArc == [ [ tori[i] contactIArcs ] lastObject ] )
							{
								[ nextArc setPhiEnd:0. ] ;
							}
					}
					
				// Reentrant R
				
				arcEnumerator = [ [ tori[i] reentrantRArcs ] objectEnumerator ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						[ nextCycle addArc:nextArc forward:YES ] ;
						
						[ nextArc setPhiStart:0. ] ;
						[ nextArc setPhiEnd:0. ] ;
					}
				
				// Contact J
				
				arcEnumerator = [ [ tori[i] contactJArcs ] objectEnumerator ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						[ nextCycle addArc:nextArc forward:YES ] ;
						
						if( nextArc == [ [ tori[i] contactJArcs ] objectAtIndex:0 ] )
							{
								[ nextArc setPhiStart:0. ] ;
							}
							
						if( nextArc == [ [ tori[i] contactJArcs ] lastObject ] )
							{
								[ nextArc setPhiEnd:[ tori[i] phiAngleSubtended ] ] ;
							}
							
					}

				// Reentrant L
				
				arcEnumerator = [ [ tori[i] reentrantLArcs ] objectEnumerator ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						[ nextCycle addArc:nextArc forward:YES ] ;
						
						[ nextArc setPhiStart:[ tori[i] phiAngleSubtended ] ] ;
						[ nextArc setPhiEnd:[ tori[i] phiAngleSubtended ] ] ;
					}
				
				
				// Compute theta and phi angles
				
				arcEnumerator = [ [ nextCycle arcs ] objectEnumerator ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						[ nextArc computeThetaAndPhiAnglesUsingMolecule:self ] ;
						
						if( [ nextArc phiStart ] > phiMax ) phiMax = [ nextArc phiStart ] ;
						if( [ nextArc phiEnd ] > phiMax ) phiMax = [ nextArc phiEnd ] ;
						
						if( [ nextArc thetaStart ] > thetaMax ) thetaMax = [ nextArc thetaStart ] ;
						if( [ nextArc thetaEnd ] > thetaMax ) thetaMax = [ nextArc thetaEnd ] ;
						
					}
					
				// Sanity check - end of each arc should match beginning of next
				
				// ALSO, all arc ends should have vertex assignment!
				
				int j ;
				
				for( j = 0 ; j < [ [ nextCycle arcs ] count ] ; ++j )
					{
						int preArc ;
						
						if( j == 0 )
							{
								preArc = [ [ nextCycle arcs ] count ] - 1 ;
							}
						else
							{
								preArc = j - 1 ;
							}
							
						SMArc *currentArc, *previousArc ;
						MMVector3 *diff ;
						
						currentArc = [ [ nextCycle arcs ] objectAtIndex:j ] ;
						previousArc = [ [ nextCycle arcs ] objectAtIndex:preArc ] ;
						
						diff = [ [ MMVector3 alloc ] initX:([ [ currentArc startPosition ] X ] - [ [ previousArc endPosition ] X ])
											Y:([ [ currentArc startPosition ] Y ] - [ [ previousArc endPosition ] Y ])
											Z:([ [ currentArc startPosition ] Z ] - [ [ previousArc endPosition ] Z ])] ;
											
						if( [ diff length ] > 0.001 )
							{
								if( warningLevel <= 1 )
									{
										printf( "WARNING: BROKEN SADDLE CYCLE!\n" ) ;
									}
							}
							
						[ diff release ] ;
						
					}
				
				[ nextCycle setThetaMax:thetaMax phiMax:phiMax ] ;
				
				if( tori[i]->selfIntersection == YES )
					{
						nextCycle->selfIntersection = YES ;
					}
				
				[ saddleCycles addObject:nextCycle ] ;
				
			}
			
		return ;
		
			
	}
	
	
- (BOOL) reduceSaddleCyclesUsingDivisionParameter:(double)div
	{
		NSEnumerator *cycleEnumerator ;
		
		SMCycle *nextCycle ;
		
		BOOL hadReduction = NO ;
		
		cycleEnumerator = [ saddleCycles objectEnumerator ] ;
		
		NSEnumerator *arcEnumerator ;
		SMArc *nextArc ;
		
		static NSMutableArray *arcsToSubdivide = nil ;
		
		if( ! arcsToSubdivide )
			{
				arcsToSubdivide = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
			}
			
		
		while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
			{
				
				
				// Note that if we have executed combine cycle operations, some cycles may already be dead
								
					
				if( [ nextCycle active ] == NO  ) continue ;
				
				// Need some extra care in this seemingly simple operation. If contact cycles have been combined into one, we will have an arc
				// AND IT'S TWIN in the same cycle! When one of the pair is encountered, an arc will be removed from the cycle. 
				
				// Initially, I will simply verify that each arc is still present in the cycle 
				
				arcEnumerator = [ [ nextCycle arcs ] objectEnumerator ] ;
				
				[ arcsToSubdivide removeAllObjects ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						// Check for no parent cycle > nextCycle 
						
						BOOL skip ;
						NSEnumerator *parentCycleEnumerator ;
						SMCycle *nextParentCycle ;
						
						skip = NO ;
						parentCycleEnumerator = [ [ nextArc parentCycles ] objectEnumerator ] ;
						
						while( ( nextParentCycle = [ parentCycleEnumerator nextObject ] ) )
							{
								if( nextParentCycle > nextCycle )
									{
										skip = YES ;
										break ;
									}
							}
							
						if( skip == YES ) continue ;
						
						// Make sure the arc is still there! (Not removed owing to subdivision of a twin in the same cycle!)
						
						if( [ [ nextCycle arcs ] indexOfObject:nextArc ] == NSNotFound ) continue ;
						
						[ arcsToSubdivide addObject:nextArc ] ;
						
						
					}
			}
		
		arcEnumerator = [ arcsToSubdivide objectEnumerator ] ;
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				// Check for a short arc that subtends too big an angle
				
				if( floor([ nextArc length ] / div ) < 2 )
					{
						if( fabs( [ nextArc phiStart ] - [ nextArc phiEnd ] ) > acos(-1.)/2. )
							{
								[ self subdivideArc:nextArc ] ;
							}
					}
				else
					{
						[ self subdivideArc:nextArc usingDivision:div ] ;
					}
			}
				
		NSMutableArray *tryArcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
		
#ifdef SURFDEBUG
		NSMutableArray *newCycles = [ [ NSMutableArray alloc ] initWithCapacity:[ saddleCycles count ] ] ;
#else
		NSMutableArray *newCycles ;
#endif
		
		while( TRUE )
			{
				NSEnumerator *cycleEnumerator ;
				int M ;
				
				SMCycle *nextCycle ;
				
#ifdef SURFDEBUG
				[ newCycles removeAllObjects ] ;
#else
				newCycles = [ [ NSMutableArray alloc ] initWithCapacity:[ saddleCycles count ] ] ;
#endif
				
				
				BOOL newCycleCreated = NO ;
				
				// For now I am going to hardwire the minimal angle in the tangent test. I will make this 1 degrees
				
				double minTangentAngle = 1. * acos(-1.)/180. ;
				
				//double minTangentAngle = 0. ;
				
				cycleEnumerator = [ saddleCycles objectEnumerator ] ;
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
					
					
						if( [ nextCycle active ] == NO  ) continue ;
						
												
						 
						NSArray *cycleArcs ;
						NSArray *cycleForward ;
						
						cycleArcs = [ nextCycle arcs ] ;
						cycleForward = [ nextCycle forward ] ;
					
						NSEnumerator *arcEnumerator ;
						SMArc *nextArc ;
							
										
						M = [ cycleArcs count ] ;
						
						//if( M == 3 ) continue ;
						if( M == 3 ) 
							{
#ifndef SURFDEBUG
								[ newCycles addObject:nextCycle ] ;
#endif
								continue ;
							}
						
						// DEBUG - Print out cycle
						/*
						printf( "******CYCLE: %p\n", nextCycle ) ;
						
						printf( "%s\n", [ [ cycleArcs description ] cString ] ) ; 
						*/
						
						// Check for skipped torus (these are still active)
						
						if( [ [ [ cycleArcs lastObject ] torusSection ] skipTorus ] == YES ) 
							{
								// We need to keep the skipped torus for vertex generation (so there will not be "holes" in the surface)
#ifndef SURFDEBUG
								[ newCycles addObject:nextCycle ] ;
#endif
								  
								
								continue ;
							}
					
						BOOL success ;
						
						success = NO ;
						
						int subdivideCount = 0 ;
						
						while( success == NO )
							{
								[ tryArcs removeAllObjects ] ;
								
								int j, k, preArc ;
								SMArc *buildArc ;
								
								
								// Need to recompute size of cycle, as we may have subdivided
								
								M = [ cycleArcs count ] ;
								
								for( j = 0 ; j <= M - 3 ; ++j )
									{
											
										for( k = j + 2 ; k <= M - 1 ; ++k )
											{													
												// Build arc from start of j to start of k
												
												// This should really be a class method of SMArc, but I only use it in this class
												
												if( j == 0 )
													{
														preArc = M - 1 ;
													}
												else
													{
														preArc = j - 1 ;
													}
													
												// Check if arc is disallowed, based on arc types
												
												SMArc *jArc, *kArc ;
												
												jArc = [ cycleArcs objectAtIndex:j ] ;
												kArc = [ cycleArcs objectAtIndex:k ] ;
												
												// This test needs to be more sophisticated - I will include additional test in buildSaddleUsingStart: (etc)
												
												if( [ jArc arcType ] != 5 && [ jArc arcType ] == [ kArc arcType ] ) continue ;
												
												
												buildArc = [ self buildSaddleUsingStart:[ cycleArcs objectAtIndex:j ] 
													startForward:[ [ cycleForward objectAtIndex:j ] boolValue ] 
													end:[ cycleArcs objectAtIndex:k ]
													endForward:[ [ cycleForward objectAtIndex:k ] boolValue ]  ] ;
												
												if( ! buildArc ) continue ;
												
												double tA, smallestTA ;
												
												smallestTA = acos(-1.) ;
												
													
												tA = [ self computeSaddleTangentAngleUsingArc:buildArc start:YES previousArc:[ cycleArcs objectAtIndex:preArc ]
														previousForward:[ [ cycleForward objectAtIndex:preArc ] boolValue ]
														nextArc:[ cycleArcs objectAtIndex:j ]  nextForward:[ [ cycleForward objectAtIndex:j ] boolValue ]
														 ];
												
												if( tA < smallestTA ) smallestTA = tA ;
												
												tA = [ self computeSaddleTangentAngleUsingArc:buildArc start:NO previousArc:[ cycleArcs objectAtIndex:(k - 1) ]
														previousForward:[ [ cycleForward objectAtIndex:(k - 1) ] boolValue ]
														nextArc:[ cycleArcs objectAtIndex:k ] nextForward:[ [ cycleForward objectAtIndex:k ] boolValue ]	
														 ] ;
												
												if( tA < smallestTA ) smallestTA = tA ;
										
														
												if( smallestTA < minTangentAngle )
													{
														[ buildArc release ] ;
														
														continue ;
													}
													
													
												NSArray *arcInfo ;
												NSNumber *arcIndexJ, *arcIndexK, *worstTangent ;
												
												arcIndexJ = [ [ NSNumber alloc ] initWithInt:j ] ;
												arcIndexK = [ [ NSNumber alloc ] initWithInt:k ] ;
												worstTangent = [ [ NSNumber alloc ] initWithDouble:smallestTA ] ;
												
												arcInfo = [ [ NSArray alloc ] initWithObjects:buildArc, arcIndexJ, arcIndexK, worstTangent, nil] ;
												
												[ tryArcs addObject:arcInfo ] ;
												
												[ buildArc release ] ;
												[ arcIndexJ release ] ; [ arcIndexK release ] ; [ worstTangent release ] ;
												
											}
									}
									
								// Now choose the best arc
								
								// We expect that subdivision should NOT be necessary - but it is!
								
								if( [ tryArcs count ] == 0 )
									{
										if( subdivideCount == MAXSUBDIVIDECOUNT )
											{
												printf( "SADDLE CYCLE COULD NOT BE REDUCED AFTER %d SUBDIVISIONS\n", MAXSUBDIVIDECOUNT ) ;
												break ;
											}
										// Find longest arc
										
										SMArc *theLongest ;
										double maxLength ;
										
										maxLength = 0. ;
										
										arcEnumerator = [ cycleArcs objectEnumerator ] ;
										
										while( ( nextArc = [ arcEnumerator nextObject ] ) )
											{
												double l ;
												
												l = [ nextArc length ] ;
												
												if( l > maxLength )
													{
														maxLength = l ;
														theLongest = nextArc ;
													}
													
											}
											
										[ self subdivideArc:theLongest ] ;
										
										++subdivideCount ;
										
										continue ;
									
									}
									
								
								/*
								if( [ tryArcs count ] == 0 )
									{
										if( subdivideCount == MAXSUBDIVIDECOUNT )
											{
												printf( "CONTACT CYCLE COULD NOT BE REDUCED AFTER %d SUBDIVISIONS\n", MAXSUBDIVIDECOUNT ) ;
												break ;
											}
											
										// Find longest arc
										
										SMArc *theLongest ;
										double maxLength ;
										
										maxLength = 0. ;
										
										arcEnumerator = [ cycleArcs objectEnumerator ] ;
										
										while( ( nextArc = [ arcEnumerator nextObject ] ) )
											{
												double l ;
												
												l = [ nextArc length ] ;
												
												if( l > maxLength )
													{
														maxLength = l ;
														theLongest = nextArc ;
													}
													
											}
											
										[ self subdivideArc:theLongest ] ;
										
										++subdivideCount ;
										
										continue ;
									}
									
									
								*/
								
								// OK, we have at least one valid subdivision
									
								// Select best-quality subdivision
								
								// Due to continuing problems with this portion of the algorithm when handling large proteins, I am going to beef up the 
								// evalution of subdivision quality by including this factor: min(edge length) / max(edge length) where the edge length will be appox., and 
								// where both subcycles are considered when arriving at the minimum value. 
								
								double q, ratio, bestQuality ;
								
								NSEnumerator *trialArcEnumerator ;
								NSArray *nextTrial, *bestTrial ;
								
								bestQuality = 0. ;
								bestTrial = nil ;
								
								double perimeterJ, perimeterK, minJ, maxJ, minK, maxK, ratioJ, ratioK ;
								
								
								trialArcEnumerator = [ tryArcs objectEnumerator ] ;
								
								while( ( nextTrial = [ trialArcEnumerator nextObject ] ) )
									{
										j = [ [ nextTrial objectAtIndex:1 ] intValue ] ;
										k = [ [ nextTrial objectAtIndex:2 ] intValue ] ;
										
										perimeterJ = perimeterK = 0. ;
										
										minJ = minK = 1.e12 ;
										maxJ = maxK = -1.e12 ;
										
										int arcIndex ;
										
										arcIndex = 0 ;
										
										arcEnumerator = [ cycleArcs objectEnumerator ] ;
										
										while( ( nextArc = [ arcEnumerator nextObject ] ) )
											{
												//deltaTheta = [ nextArc thetaEnd ] - [ nextArc thetaStart ] ;
												//deltaPhi = [ nextArc phiEnd ] - [ nextArc phiStart ] ;
												
												if( arcIndex < j || arcIndex >= k )
													{
														perimeterK += [ nextArc length ] ;
														//perimeterK += sqrt( deltaTheta*deltaTheta + deltaPhi*deltaPhi ) ;
														
														if( [ nextArc length ] < minK ) minK = [ nextArc length ] ;
														if( [ nextArc length ] > maxK ) maxK = [ nextArc length ] ;
													
													}
												else
													{
														perimeterJ += [ nextArc length ] ;
														
														if( [ nextArc length ] < minJ ) minJ = [ nextArc length ] ;
														if( [ nextArc length ] > maxJ ) maxJ = [ nextArc length ] ;
														
														//perimeterJ += sqrt( deltaTheta*deltaTheta + deltaPhi*deltaPhi ) ;
													}
													
												++arcIndex ;
											}
											
										ratio = (MIN(perimeterJ,perimeterK)/MAX(perimeterJ,perimeterK)) ;
										
										ratioJ = minJ/maxJ ;
										ratioK = minK/maxK ;
										
										// NOTE - due to problems that seem to arise with the perimeter test (can end up with very narrow
										// elements), I am only including tangent test in assessing division quality
										
										q = ratio * sin( [ [ nextTrial objectAtIndex:3 ] doubleValue ] ) * ratioJ * ratioK ;
										//q = sin( [ [ nextTrial objectAtIndex:3 ] doubleValue ] ) ;
										
										if( q > bestQuality )
											{
												bestQuality = q ;
												bestTrial = nextTrial ;
											}
									}
									
								// Subdivide cycle using best quality division
								/*
								printf( "**BEST SUBDIVIDING ARC: %s\n", [ [ [ bestTrial objectAtIndex:0 ] description ] cString ] ) ;
								printf( "**WORST TANGENT: %f\n", [ [ bestTrial objectAtIndex:3 ] doubleValue ] ) ;
								printf( "**QUALITY: %f\n", bestQuality ) ;
								*/
								
								NSArray *subCycles ;
								
								subCycles = [ self subdivideCycle:nextCycle firstIndex:[ [ bestTrial objectAtIndex:1 ] intValue ]
									secondIndex:[ [ bestTrial objectAtIndex:2 ] intValue ] usingArc:[ bestTrial objectAtIndex:0 ] ] ;
									
								/*
								printf( "***BEFORE SUBDIVISION ------\n" ) ;
								printf( "\n**NEW CYCLE 1: %p\n%s\n", [ subCycles objectAtIndex:0 ], [ [ [ [ subCycles objectAtIndex:0 ] arcs ] description ] cString ] ) ;
								printf( "\n**NEW CYCLE 2: %p\n%s\n", [ subCycles objectAtIndex:1 ], [ [ [ [ subCycles objectAtIndex:1 ] arcs ] description ] cString ] ) ;
								*/
									
								//DEBUGING - did we generate a cycle with all theta = 0 ?
								
								BOOL thetaNotZero ;
								
								
								/*
								subArcEnumerator = [ [ [ subCycles objectAtIndex:0 ] arcs ] objectEnumerator ] ;
								
								thetaNotZero = NO ;
								
								while( ( subArc = [ subArcEnumerator nextObject ] ) )
									{
										if( [ subArc thetaStart ] > 0.001 || [ subArc thetaEnd ]  > 0.001 )
											{
												thetaNotZero = YES ;
											}
									}
									
								if( thetaNotZero == NO )
									{
										if( warningLevel <= 1 )
											{
												printf( "WARNING - MADE SADDLE SUBCYCLE WITH ALL THETA ZERO!\n" ) ;
											}
									}
									
								subArcEnumerator = [ [ [ subCycles objectAtIndex:1 ] arcs ] objectEnumerator ] ;
								
								thetaNotZero = NO ;
								
								while( ( subArc = [ subArcEnumerator nextObject ] ) )
									{
										if( [ subArc thetaStart ]  > 0.001 || [ subArc thetaEnd ]  > 0.001 )
											{
												thetaNotZero = YES ;
											}
									}
									
								if( thetaNotZero == NO )
									{
										if( warningLevel <= 1 )
											{
												printf( "WARNING - MADE SADDLE SUBCYCLE WITH ALL THETA ZERO!\n" ) ;
											}
									}
								*/
									
								// Now that the arc is embedded in new cycles, subdivide it!
								
								double x, rDiv ;
								SMArc *theBest ;
								
								theBest = [ bestTrial objectAtIndex:0 ] ;
								
								x = [ theBest length ] ;
								
								rDiv = floor( x / div ) ;
								
								if( rDiv < 2. )
									{
										double dPhi ;
										dPhi = fabs( [ [ bestTrial objectAtIndex:0 ] phiStart ] - [ [ bestTrial objectAtIndex:0 ] phiEnd ] ) ;
										
										if(  dPhi > acos(-1.)/2. )
											{
												[ self subdivideArc:[ bestTrial objectAtIndex:0 ] ] ;
											}
									}
								else
									{
										[ self subdivideArc:[ bestTrial objectAtIndex:0 ] usingDivision:div ] ;
									}
								
								//[ self subdivideArc:[ bestTrial objectAtIndex:0 ] usingDivision:div ] ;
								
								/*
								printf( "***AFTER SUBDIVISION ------\n" ) ;
								printf( "\n**NEW CYCLE 1: %p\n%s\n", [ subCycles objectAtIndex:0 ], [ [ [ [ subCycles objectAtIndex:0 ] arcs ] description ] cString ] ) ;
								printf( "\n**NEW CYCLE 2: %p\n%s\n", [ subCycles objectAtIndex:1 ], [ [ [ [ subCycles objectAtIndex:1 ] arcs ] description ] cString ] ) ;
								*/
								
																
								[ nextCycle killCycle ] ;
								
								[ newCycles addObjectsFromArray:subCycles ] ;
								newCycleCreated = YES ;
								
								[ subCycles release ] ;
								
								success = YES ;
								
								hadReduction = YES ;
								
							}
							
						if( success == NO )
							{
								printf( "COULD NOT REDUCE SADDLE CYCLE - Exit!\n" ) ;
								exit(1) ;
							}
							
						// Next cycle 
					}
					
				// Did we reduce any cycles?
				
#ifdef SURFDEBUG
				if( [ newCycles count ] == 0 ) break ;
				
				[ saddleCycles addObjectsFromArray:newCycles ] ;
#else
				if( newCycleCreated == NO )
				//if( [ newCycles count ] == [ saddleCycles count ] ) 
					{
						[ newCycles release ] ;
						break ;
					}
				
				[ saddleCycles release ] ;
				saddleCycles = newCycles ;
#endif
			}
			
			
		return hadReduction ;
	}
		
		
- (void)  generateVerticesUsingSubsurfaces:(BOOL)sub
	{
		// This is my 5th go at this!
		
		// My brain has broken on this problem, which is frankly a byproduct of my initial prejudice that the problem was "simple". It is not simple, owing to the
		// large degree of independance that I have assigned to arcs in the formal development
		
		// Define a vertex object. A vertex object points back at arc ends using a small SMArcEnd object
		
		// Algorithm:
		//
		// Enumerate all live cycles. As each corner is discovered, add a vertex object if one does not already exist
		//
		// If it is discovered that a a corner is already associated with a vertex, add the arc ends. 
		//
		// If it is discovered that two different vertices are already associated with the same corner (this can happen, depending on how the arcs are enumerated),
		// then these are merged into a single vertex
		// 
		// If subsurfaces are selected, we will identify them, with 0 being the largest (in number of vertices), 1 the next largest, etc
		//
		// An output file will be generated for each subsurface, using the name <basename.flats>.# , where # = 0, 1, 2, ...
		//
		// Algorithm to identify subsurfaces:
		//
		// Prelim: Each vertex holds a subsurface index, initially set to -1 (no assignment), and likewise each element
		// nSubsurface is # of subsurfaces identified. 
		//
		// Find first unassigned vertex, vU
		//
		// Initialize boundaryArcs set from vertex, initialize boundaryElements set from boundaryArcs
		//
		// WHILE boundaryElements not Empty:
		//		Assign boundaryElements to Subsurface, add to accumElements set
		//		Generate newBoundaryArce from boundaryElements
		//		Set boundaryArcs = newBoundaryArcs - boundaryArcs
		//		Empty boundaryElements
		//		Fill boundaryElements using boundaryArcs, using only UNassigned elements
		//	ENDWHILE
		//
		//	Find next unassigned vertex, etc. 
		
		//NSArray *cycleTypes = [ [ NSArray alloc ] initWithObjects:contactCycles,reentrantCycles,saddleCycles,nil ] ;
		
		
		NSArray *cycleTypes = [ [ NSArray alloc ] initWithObjects:saddleCycles,reentrantCycles,contactCycles,nil ] ;
		NSArray *nextCycleType ;
		NSEnumerator *cycleTypeEnumerator ;
		NSMutableArray *twinnedArcs ;
		
		
		
		twinnedArcs = [ [ NSMutableArray alloc ] initWithCapacity:1000 ] ;
		
		nElements = 0 ;
		nContactElements = 0 ;
		nReentrantElements = 0 ;
		nSaddleElements = 0 ;
		
		cycleTypeEnumerator = [ cycleTypes objectEnumerator ] ;
		
		while( ( nextCycleType = [ cycleTypeEnumerator nextObject ] ) )
			{
				NSEnumerator *cycleEnumerator ;
				SMCycle *nextCycle ;
				
				cycleEnumerator = [ nextCycleType objectEnumerator ] ;
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
				
						int iArc, iPreArc ;
						NSArray *theArcs ;
						SMArc *nextArc, *preArc ;
						BOOL nextForward, preForward ;
						
							
						if( [ nextCycle active ] == NO ) continue ;
						
						//printf( "AT CYCLE: %p \n", nextCycle ) ;
						
						theArcs = [ nextCycle arcs ] ;
						
						
						for( iArc = 0 ; iArc < [ theArcs count ] ; ++iArc )
							{
								BOOL newVertex ;
								
								newVertex = NO ;
								
								iPreArc = iArc - 1 ;
								
								if( iPreArc < 0 )
									{
										iPreArc = [ theArcs count ] - 1 ;
									}
									
								nextArc = [ theArcs objectAtIndex:iArc ] ;
								preArc = [ theArcs objectAtIndex:iPreArc ] ;
								
								if( [ nextArc twin ] ) [ twinnedArcs addObject:nextArc ] ;
								
								nextForward = [ [ [ nextCycle forward ] objectAtIndex:iArc ] boolValue ] ;
								preForward = [ [ [ nextCycle forward ] objectAtIndex:iPreArc ] boolValue ] ;
								
								SMVertex *nextVertex, *preVertex, *keepVertex ;
								
								if( nextForward == YES )
									{
										nextVertex = [ nextArc startVertex ] ;
									}
								else
									{
										nextVertex = [ nextArc endVertex ] ;
									}
									
								if( preForward == YES )
									{
										preVertex = [ preArc endVertex ] ;
									}
								else
									{
										preVertex = [ preArc startVertex ] ;
									}
								
								if( nextVertex )
									{
										if( preVertex )
											{
												if( nextVertex != preVertex )
													{
														// Keep the nextVertex
														
														[ nextVertex mergeVertex:preVertex ] ;
														[ vertices removeObject:preVertex ] ;
														
														keepVertex = nextVertex ;
													}
												else
													{
														keepVertex = nextVertex ;
													}
											}
										else
											{
												keepVertex = nextVertex ;
											}
									}
								else
									{
										if( preVertex )
											{
												keepVertex = preVertex ;
											}
										else
											{
												keepVertex = nil ;
											}
									}
									
								if( ! keepVertex )
									{
										// Make a new vertex
										
										newVertex = YES ;
										
										keepVertex = [ [ SMVertex alloc ] init ] ;
										[ vertices addObject:keepVertex ] ;
										[ keepVertex release ] ;
										
										// Compute vertex position below
										
									}
										
								if( nextForward == YES )
									{
										[ keepVertex addArc:nextArc start:YES ] ;
										[ nextArc setStartVertex:keepVertex ] ;
									}
								else
									{
										[ keepVertex addArc:nextArc start:NO ] ;
										[ nextArc setEndVertex:keepVertex ] ;
									}
								
								if( preForward == YES )
									{
										[ keepVertex addArc:preArc start:NO ] ;
										[ preArc setEndVertex:keepVertex ] ;
									}
								else
									{
										[ keepVertex addArc:preArc start:YES ] ;
										[ preArc setStartVertex:keepVertex ] ;
									}
									
								if( newVertex == YES )
									{
										if( nextCycleType == saddleCycles )
											{
												[ keepVertex computePositionForSaddleUsingMolecule:self ] ;
												
											}
										else if( nextCycleType == reentrantCycles )
											{
													
												[ keepVertex computePositionForReentrantCycle:nextCycle usingMolecule:self ] ;
											}
										else
											{
												[ keepVertex computePositionForContact ] ;
											}
									}
									
									
								// Next arc in cycle 
							}
							
						// Next cycle
						
						++nElements ;
						
						if( nextCycleType == saddleCycles)
							{
								++nSaddleElements ;
							}
						else if( nextCycleType == reentrantCycles )
							{
								++nReentrantElements ;
							}
						else
							{
								++nContactElements ;
							}
					}
				// Next cycle type
			}
			
		// Want to check here that every cycle has valid vertices
		/*
		cycleTypeEnumerator = [ cycleTypes objectEnumerator ] ;
		
		while( ( nextCycleType = [ cycleTypeEnumerator nextObject ] ) )
			{
				NSEnumerator *cycleEnumerator ;
				SMCycle *nextCycle ;
				
				cycleEnumerator = [ nextCycleType objectEnumerator ] ;
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						if( [ nextCycle active ] == NO ) continue ;
						
						if( nextCycleType == contactCycles )
							{
								printf( "CYCLE: %p (CONTACT)\n", nextCycle ) ;
							}
						else if( nextCycleType == reentrantCycles )
							{
								printf( "CYCLE: %p (REENTRANT)\n", nextCycle ) ;
							}
						else
							{
								printf( "CYCLE: %p (SADDLE)\n", nextCycle ) ;
							}
						
						
						NSArray *theArcs = [ nextCycle arcs ] ;
						
						NSEnumerator *arcEnumerator ;
						SMArc *nextArc ;
						
						arcEnumerator = [ theArcs objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								printf( "\tArc:%p\n", nextArc ) ;
								printf( "\tStart Vertex:\n" ) ;
								printf( "%s", [ [ [ nextArc startVertex ] description ] cString ] ) ;
								printf( "\tEnd Vertex:\n" ) ;
								printf( "%s", [ [ [ nextArc endVertex ] description ] cString ] ) ;
							}
					}
			}
						
		*/
		
		// Check integrity of vertices
		
		//BOOL integrity ;
		
		//integrity = [ self checkVertexIntegrity ] ;
		
		//if( integrity == NO )
		//	{
		//		printf( "WARNING: VERTEX INTEGRITY COMPROMISED BEFORE MERGE!\n" ) ;
		//	}
		
		
		// Now need to merge vertices that are twinned
		
		// DEBUGGING - give each vertex a current index
		
		/*
		NSEnumerator *vertexEnumerator ;
		SMVertex *nextVertex ;
		
		vertexEnumerator = [ vertices objectEnumerator ] ;
		
		int oldIndex = 0 ;
		
		while( ( nextVertex = [ vertexEnumerator nextObject ] ) )
			{
				[ nextVertex setOldIndex:oldIndex ] ;
				[ nextVertex setIndex:oldIndex ] ;
				
				++oldIndex ;
			}
			
		
		nVertices = oldIndex ;
		*/
		
		// For the purposes of debugging, write out the surface prior to vertex merge so that we have valid vertex indices
		
		//[ self exportAsFlatsUsingPath:@"/Users/zauhar/SVN/smart007Testing/indinavirPriorMerge.flats" oldStyle:YES useSubsurfaces:NO ] ;
		
		// Subsurface generation requires that twinned arcs be appropriately associated with cycles. That means that cycles that 
		// lie on opposite sides of a deleted torus need to be associated with both twinned arcs
		
		// NOTE that our design with twinned arcs ASSUMES that they are always associated with deleted (skipped) tori!
		
		NSEnumerator *arcEnumerator ;
		SMArc *nextArc ;
		BOOL didAMerge ;
		
		didAMerge = YES ;
		
		
		while( didAMerge == YES )
			{
				didAMerge = NO ;
				
				
				arcEnumerator = [ twinnedArcs objectEnumerator ] ;
				
				// Note that twinned arcs have opposite orientation!
				
				SMVertex *toRemove ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						// Make sure that twinned arcs point two elements
						
						if( [ [ nextArc parentCycles ] count ] != 3 )
							{
								// Why is the test for 3, not 2? 
								// Because the four-cycle that represents the torus-to-be-deleted also needs to be
								// included - it will be delted in a step below. 
								
								// Add parent cycle of twin 
								
								[ nextArc mergeParentCyclesFromArc:[ nextArc twin ] ] ;
								
								[ [ nextArc twin ] mergeParentCyclesFromArc:nextArc ] ;
							}
							
						//printf( "Process arc %p with twin %p\n", nextArc, [ nextArc twin ] ) ;
						if( [ [ nextArc twin ] endVertex ] != [ nextArc startVertex ] )
							{
								didAMerge = YES ;
								
								toRemove = [ [ nextArc twin ] endVertex ] ;
								
								[ [ nextArc startVertex ] mergeVertex:[ [ nextArc twin ] endVertex ] ] ;
								//printf( "Merge twin end vertex %p with start vertex %p\n", toRemove, [ nextArc startVertex ] ) ;
								
								[ vertices removeObject:toRemove ] ;
								//printf( "Remove twin end vertex %p\n", toRemove ) ;
							}
							
						if( [ [ nextArc twin ] startVertex ] != [ nextArc endVertex ] )
							{
								didAMerge = YES ;
								
								toRemove = [ [ nextArc twin ] startVertex ] ;
								
								[ [ nextArc endVertex ] mergeVertex:[ [ nextArc twin ] startVertex ] ] ;
								//printf( "Merge twin start vertex %p with end vertex %p\n", toRemove, [ nextArc endVertex ] ) ;
								
								[ vertices removeObject:toRemove ] ;
								//printf( "Remove twin start vertex %p\n", toRemove ) ;
							}
					}
			}
			
		//integrity = [ self checkVertexIntegrity ] ;
		
		//if( integrity == NO )
		//	{
		//		printf( "WARNING: VERTEX INTEGRITY COMPROMISED AFTER MERGE, BEFORE CYCLE DELETION!\n" ) ;
		//	}
		

		// All vertices should be assigned and merges complete!
		
		// Now generate vertices (easy - hah hah)
		
		
		NSEnumerator *vertexEnumerator = [ vertices objectEnumerator ] ;
		SMVertex *nextVertex ;
		
		nVertices = 0 ; 
		
		while( ( nextVertex = [ vertexEnumerator nextObject ] ) )
			{
			
				// Check if any vertex has illegal position - bail out if so
				
				MMVector3 *nextVertexPos ;
				
				nextVertexPos = [ nextVertex vertexPosition ] ;
				
				if( isnan( [ nextVertexPos X ] ) || isnan( [ nextVertexPos Y ] ) || isnan( [ nextVertexPos Z ] ) )
					{
						printf( "ILLEGAL VERTEX POSITION (NaN) - Exit!\n" ) ;
						exit(1) ;
					}
					
				[ nextVertex setIndex:nVertices ] ;
				
					
				++nVertices ;
			}
			
		// Check for collapsed cycles, which can arise from skipped arcs (and some that end up effectively being skipped owing to vertex merging)
		
		// We will also remove saddle cycles assoc with skipped tori
		
		cycleTypeEnumerator = [ cycleTypes objectEnumerator ] ;
		
		while( ( nextCycleType = [ cycleTypeEnumerator nextObject ] ) )
			{
				NSEnumerator *cycleEnumerator ;
				SMCycle *nextCycle ;
				
				cycleEnumerator = [ nextCycleType objectEnumerator ] ;
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						int M ;
						NSEnumerator *arcEnumerator ;
						SMArc *nextArc ;
						
						if( [ nextCycle active ] == NO ) continue ;
						
						if( nextCycleType == saddleCycles )
							{
								if( [ [ [ [ nextCycle arcs ] lastObject ] torusSection ] skipTorus ] == YES )
									{
										[ nextCycle killCycle ] ;
										--nSaddleElements ;
										
										--nElements ;
										
										continue ;
									}
							}
						
						arcEnumerator = [ [ nextCycle arcs ] objectEnumerator ] ;
						
						M = 0 ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								// The following commented-out code is meant to address an issue arising from free tori, where in principle a
								// arc can have the same vertex at start and finish. I have introduced a constraint that saddle arcs must always
								// subtend an angle less than 90 deg, so this should no longer be required.
								
								/*
								if( [ nextArc phiStart ] >= 0 && [ nextArc skip ] == NO )
									{
										if( [ nextArc phiStart ] != [ nextArc phiEnd ] || [ nextArc thetaStart ] != [ nextArc thetaEnd ] )
											{
												++M ;
											}
											
										if( [ nextArc startVertex ] == [ nextArc endVertex ] )
											{
												printf( "WARNING: VERTICES THE SAME ON ONE EDGE!\n" ) ;
											}
									}
								else if( [ nextArc startVertex ] != [ nextArc endVertex ] )
								*/
								if( [ nextArc skip ] == NO && [ nextArc startVertex ] != [ nextArc endVertex ] )
									{
										++M ;
									}
							}
							
						if( M < 3 )
							{
								[ nextCycle killCycle ] ;
								
								if( nextCycleType == contactCycles )
									{
										--nContactElements ;
										--nElements ;
									}
								else if( nextCycleType == reentrantCycles )
									{
										--nReentrantElements ;
										--nElements ;
									}
								else
									{
										--nSaddleElements ;
										--nElements ;
									}
							}
							
						
					}
					
				// Next cycle type 
			}
			
		//integrity = [ self checkVertexIntegrity ] ;
		
		//if( integrity == NO )
		//	{
		//		printf( "WARNING: VERTEX INTEGRITY COMPROMISED AFTER CYCLE DELETION!\n" ) ;
		//	}

		printf( "\nSurface includes %d vertices, %d contact elements, %d reentrant elements, %d saddle elements (%d total)\n",
			nVertices, nContactElements, nReentrantElements, nSaddleElements, nElements ) ;
		
		// Now generate subsurfaces, if needed
		
		if( sub == NO ) return ;
		
		nSubsurfaces = 0 ;
		
		NSMutableSet *boundaryArcs = [ [ NSMutableSet alloc ] initWithCapacity:( (int)floor( 1.5 * nElements ) ) ] ;
		NSMutableSet *newBoundaryArcs = [ [ NSMutableSet alloc ] initWithCapacity:( (int)floor( 1.5 * nElements ) ) ] ;
		
		NSMutableSet *boundaryElements = [ [ NSMutableSet alloc ] initWithCapacity:nElements  ] ;
		
		
		// First vertex should be unassigned
		
		SMVertex *nextUnassignedVertex = [ vertices objectAtIndex:0 ] ;
		
		
		
		
		while( nextUnassignedVertex )
			{
			
				int nVertexLocal = 0 ;
				
				// Use arcEnds of vertex to seed boundary arcs array
				
				[ boundaryArcs removeAllObjects ] ;
				[ boundaryElements removeAllObjects ] ;
				
				NSEnumerator *arcEndEnumerator = [ [ nextUnassignedVertex arcEnds ] objectEnumerator ] ;
				SMArcEnd *nextArcEnd ;
				
				while( ( nextArcEnd = [ arcEndEnumerator nextObject ] ) )
					{
						[ boundaryArcs addObject:[ nextArcEnd arc ] ] ;
					}
				
				NSEnumerator *arcEnumerator ;
				SMArc *nextArc ;
				
				while( [ boundaryArcs count ] > 0 )
					{
						arcEnumerator = [ boundaryArcs objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								[ boundaryElements addObjectsFromArray:[ nextArc parentCycles ] ] ;
								
								if( [ [ nextArc startVertex ] subsurface ] <  0 ) 
									{
										[ [ nextArc startVertex ] setSubsurface:nSubsurfaces ] ;
										[ [ nextArc startVertex ] setSubsurfaceIndex:nVertexLocal ] ;
										++nVertexLocal ;
									}
									
								if( [ [ nextArc endVertex ] subsurface ] <  0 ) 
									{
										[ [ nextArc endVertex ] setSubsurface:nSubsurfaces ] ;
										[ [ nextArc endVertex ] setSubsurfaceIndex:nVertexLocal ] ;
										++nVertexLocal ;
									}
									
							}
							
						// Set current collection of boundary elements as belonging to this subsurface - also, collect all arcs
						
						[ newBoundaryArcs removeAllObjects ] ;
						
						NSEnumerator *cycleEnumerator = [ boundaryElements objectEnumerator ] ;
						SMCycle *nextCycle ;
						
						while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
							{
								if( [ nextCycle active ] == NO ) continue ;
								
								if( [ nextCycle subsurface ] >= 0 ) continue ;
								
								[ nextCycle setSubsurface:nSubsurfaces ] ;
								
								[ newBoundaryArcs addObjectsFromArray:[ nextCycle arcs ] ] ;
							}
							
						[ newBoundaryArcs minusSet:boundaryArcs ] ;
						
						[ boundaryArcs setSet:newBoundaryArcs ] ;
					}
					
				// Increment number of subsurfaces
				
				++nSubsurfaces ;
					
				// Now see if there are unassigned vertices
				
				NSEnumerator *vertexEnumerator = [ vertices objectEnumerator ] ;
				SMVertex *nextVertex ;
				
				nextUnassignedVertex = nil ;
				
				while( ( nextVertex = [ vertexEnumerator nextObject ]  ) )
					{
						if( [ nextVertex subsurface ] < 0 )
							{
								nextUnassignedVertex = nextVertex ;
								break ;
							}
					}
					
			}
				
				
				
	   printf( "\nGenerated %d subsurfaces\n", nSubsurfaces ) ;
				
				
	
				
		
		
		return ;
	}
							
		
- (BOOL) checkVertexIntegrity
	{
		// This utility function checks that 
		//	1) Every cycle corner is associated with a vertex, and both arc ends in the corner point at the same vertex
		//	2) Every vertex collects the same arc ends 
				
		BOOL hadError ;
		NSArray *cycleTypes = [ [ NSArray alloc ] initWithObjects:saddleCycles,reentrantCycles,contactCycles,nil ] ;
		NSArray *nextCycleType ;
		NSEnumerator *cycleTypeEnumerator ;
		
		NSMutableArray *checkVertices = [ [ NSMutableArray  alloc ] initWithArray:vertices ] ;
		
		cycleTypeEnumerator = [ cycleTypes  objectEnumerator ] ;
		
		hadError = NO ;
		
		while( ( nextCycleType = [ cycleTypeEnumerator nextObject ] ) )
			{
				
				NSEnumerator *cycleEnumerator ;
				SMCycle *nextCycle ;
				
				cycleEnumerator = [ nextCycleType objectEnumerator ] ;
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
				
						int iArc, iPreArc ;
						NSArray *theArcs ;
						SMArc *nextArc, *preArc ;
						BOOL nextForward, preForward ;
						SMVertex *nextVertex, *preVertex ;
						
						if( [ nextCycle active ] == NO ) continue ;
						
						//printf( "AT CYCLE: %p \n", nextCycle ) ;
						
						theArcs = [ nextCycle arcs ] ;
						
						
						for( iArc = 0 ; iArc < [ theArcs count ] ; ++iArc )
							{
								
								iPreArc = iArc - 1 ;
								
								if( iPreArc < 0 )
									{
										iPreArc = [ theArcs count ] - 1 ;
									}
									
								nextArc = [ theArcs objectAtIndex:iArc ] ;
								preArc = [ theArcs objectAtIndex:iPreArc ] ;
								
								if( [ [ nextArc parentCycles ] count ] != 2 )
									{
										if( ![ nextArc torusSection ] || [ [ nextArc torusSection ] freeTorus ] != YES || [ nextArc twin ] == nil )
											{
												printf( "ARC %p HAS %d PARENT CYCLES!\n", nextArc, (int)[ [ nextArc parentCycles ] count ] ) ;
												hadError = YES ;
											}
									}
												
		
								nextForward = [ [ [ nextCycle forward ] objectAtIndex:iArc ] boolValue ] ;
								preForward = [ [ [ nextCycle forward ] objectAtIndex:iPreArc ] boolValue ] ;
								
								if( nextForward )
									{
										nextVertex = [ nextArc startVertex ] ;
									}
								else
									{
										nextVertex = [ nextArc endVertex ]  ;
									}
								
								if( preForward )
									{
										preVertex = [ preArc endVertex ] ;
									}
								else
									{
										preVertex = [ preArc startVertex ]  ;
									}
									
								if( preVertex != nextVertex )
									{
										printf( "UNMATCHED VERTICES: %p %p AT CORNER IN CYCLE %p \n", preArc, nextArc, nextCycle ) ;
										hadError = YES ;
										continue ;
									}
									
								[ checkVertices removeObject:nextVertex ] ;
								
								// Next arc
							}
						// Next cycle
					}
					
				// Next cycle type
			}
			
		// Any vertices unused?
		
		if( [ checkVertices count ] > 0 )
			{
				printf( "WARNING: %d UNUSED VERTICES!\n", (int)[ checkVertices count ] ) ;
				hadError = YES ;
			}
			
		[ checkVertices release ] ;
									
		
		// Check all vertices for agreement with arcends (all arcs pointed to by vertex should be referenced by the arcs)
		
		NSEnumerator *vertexEnumerator ;
		SMVertex *nextVertex ;
		NSEnumerator *arcEndEnumerator ;
		SMArcEnd *nextArcEnd ;
		SMArc *theArc ;
		
		vertexEnumerator = [ vertices objectEnumerator ] ;
		
		while( ( nextVertex = [ vertexEnumerator nextObject ] ) )
			{
				arcEndEnumerator = [ [ nextVertex arcEnds ] objectEnumerator ] ;
				
				while( ( nextArcEnd = [ arcEndEnumerator nextObject ] ) )
					{
						theArc = [ nextArcEnd arc ] ;
						
						if( [ nextArcEnd atStart ] )
							{
								if( [ theArc startVertex ] != nextVertex )
									{
										printf( "BROKEN REFERENCE: VERTEX %p POINTS TO ARC %p WHICH REFERS TO VERTEX %p\n",
											nextVertex, theArc, [ theArc startVertex ] ) ;
										hadError = YES ;
									}
							}
						else
							{
								if( [ theArc endVertex ] != nextVertex )
									{
										printf( "BROKEN REFERENCE: VERTEX %p POINTS TO ARC %p WHICH REFERS TO VERTEX %p\n",
											nextVertex, theArc, [ theArc endVertex ] ) ;
										hadError = YES ;
									}
							}
					}
			}
			
		// Finally check that each edge is part of two and only two elements
		
		
			
		if( hadError )
			return NO ;
		else
			return YES ;
		
	}
								
								
									
								
	
- (void) combineContactCycles 
	{
		// This method combines contact cycles at each atom, if possible. It does this by considering all unique pairs of cycles, and attempting to 
		// build a geodesic arc from the start of arc 0 in the first cycle to any arc in the second of the pair. If success, the cycles are joined into one. 
		// This is only possible if one cydle is interior to the other. 
		
		int iAtom, cycleCount, j, k ;
		
		double minTangentAngle = 0. * acos(-1.)/180. ;
		
		NSMutableArray *tryArcs = [ NSMutableArray arrayWithCapacity:10 ] ;
		
		for( iAtom = 0 ; iAtom < nAtoms ; ++iAtom )
			{
				if( ! atomsToCycles[iAtom] ) continue ;
				
				// While loop is executed while we still need to check for cycle combination
				
				while( TRUE )
					{
				
						int N ;
						
						NSEnumerator *cycleEnumerator ;
						SMCycle *nextCycle ;
						
						cycleEnumerator = [ atomsToCycles[iAtom] objectEnumerator ] ;
						
						N = 0 ;
						
						while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
							{
								if( [ nextCycle active ] == YES )
									{
										++N ;
									}
							}
						
						if( N == 1 ) break ;
						
						cycleCount = [ atomsToCycles[iAtom] count ] ;
						
						SMCycle *cycleJ ;
						
						BOOL combinedCycles ;
						
						combinedCycles = NO ;
						
						for( j = 0 ; j < cycleCount - 1 ; ++j )
							{
								cycleJ = [ atomsToCycles[iAtom] objectAtIndex:j ] ;
								
								if( [ cycleJ active ] == NO ) continue ;
								
								SMArc *sourceArc, *preSourceArc ;
								
								
										for( k = j + 1 ; k < cycleCount ; ++k )
											{
												// Try to build arc from arcs of cycle j to any of the arcs
												// of cycle k. Tangent and intersection tests are applied.
												
												int M ;
										
												SMCycle *cycleK ;
											
												cycleK = [ atomsToCycles[iAtom] objectAtIndex:k ] ;
												M = [ [ cycleK arcs ] count ] ;
												
												if( [ cycleK active ] == NO ) continue ;
												
												SMArc *targetArc, *buildArc ;
												int kArc, preArc ;
												
												[ tryArcs removeAllObjects ] ;
												
												int jArc, jPreArc ;
												
												for( jArc = 0 ; jArc < [ [ cycleJ arcs ] count ] ; ++jArc )
													{
												
														sourceArc = [ [ cycleJ arcs ] objectAtIndex:jArc ] ;
														
														jPreArc = jArc - 1 ;
														
														if( jPreArc < 0 )
															{
																jPreArc = [ [ cycleJ arcs ] count ] - 1 ;
															}
															
														preSourceArc = [ [ cycleJ arcs ] objectAtIndex:jPreArc ] ;
												
												
														for( kArc = 0 ; kArc < M ; ++kArc )
															{
																if( kArc == 0 )
																	{
																		preArc = M - 1 ;
																	}
																else
																	{
																		preArc = kArc - 1 ;
																	}
																	
																targetArc = [ [ cycleK arcs ] objectAtIndex:kArc ] ;
																
																buildArc = [ self buildGeoUsingStart:sourceArc startForward:YES end:targetArc endForward:YES reentrantRegion:NO ] ;
																
																if( ! buildArc ) continue ;
																
																double tA, smallestTA ;
																
																smallestTA = acos(-1.) ;
																													
																tA = [ self computeTangentAngleUsingArc:buildArc start:YES previousArc:preSourceArc
																		previousForward:YES
																		nextArc:sourceArc nextForward:YES reentrantRegion:NO ] ;
																
																if( tA < smallestTA ) smallestTA = tA ;
																
																tA = [ self computeTangentAngleUsingArc:buildArc start:NO previousArc:[ [ cycleK arcs ] objectAtIndex:preArc ]
																		previousForward:YES
																		nextArc:targetArc nextForward:YES reentrantRegion:NO	] ;
																
																if( tA < smallestTA ) smallestTA = tA ;
														
																		
																if( smallestTA < minTangentAngle )
																	{
																		[ buildArc release ] ;
																		
																		continue ;
																	}
																	
																// Compute intersections
																
																// It's time-consuming, but this test must be performed on all arcs of all cycles
																
																SMArc *testArc ;
																NSEnumerator *arcEnumerator ;
																int arcIndex ;
																BOOL haveIntersection ;
																
																haveIntersection = NO ;
																
																arcEnumerator = [ [ cycleK arcs ] objectEnumerator ] ;
																
																arcIndex = -1 ;
																
																while( ( testArc = [ arcEnumerator nextObject ] ) )
																	{
																		++arcIndex ;
																		
																		if(  arcIndex == kArc - 1 || arcIndex == kArc )
																			{
																				continue ;
																			}
																			
																		if( [ buildArc intersectWith2:testArc ] == YES )
																			{
																				haveIntersection = YES ;
																				break ;
																			}
																	}
																	
																if( haveIntersection == YES )
																	{
																		[ buildArc release ] ;
																		
																		continue ;
																	}
																	
																arcEnumerator = [ [ cycleJ arcs ] objectEnumerator ] ;
																
																arcIndex = -1 ;
																
																while( ( testArc = [ arcEnumerator nextObject ] ) )
																	{
																		++arcIndex ;
																		
																		if(  arcIndex == jArc || arcIndex == jPreArc )
																			{
																				continue ;
																			}
																			
																		if( [ buildArc intersectWith2:testArc ] == YES )
																			{
																				haveIntersection = YES ;
																				break ;
																			}
																	}
																	
																if( haveIntersection == YES )
																	{
																		[ buildArc release ] ;
																		
																		continue ;
																	}
																	
																// Now test against all arcs in other cycles
																	
																NSEnumerator *testCycles = [ atomsToCycles[iAtom] objectEnumerator ] ;
																SMCycle *nextTestCycle ;
																
																while( ( nextTestCycle = [ testCycles nextObject ] ) )
																	{
																		if( [ nextTestCycle active ] == NO ) continue ;
																		
																		if( nextTestCycle == cycleK || nextTestCycle == cycleJ ) continue ;
																	
																		arcEnumerator = [ [ nextTestCycle arcs ] objectEnumerator ] ;
																		
																		while( ( testArc = [ arcEnumerator nextObject ] ) )
																			{
																				if( [ buildArc intersectWith2:testArc ] == YES )
																					{
																						haveIntersection = YES ;
																						break ;
																					}
																			}
																	
																		if( haveIntersection == YES )
																			{
																				break ;
																			}
																	}
																	
																	
																if( haveIntersection == YES )
																	{
																		[ buildArc release ] ;
																		
																		continue ;
																	}
																	
																NSArray *arcInfo ;
																NSNumber *arcIndexJ, *arcIndexK, *worstTangent ;
																
																arcIndexJ = [ [ NSNumber alloc ] initWithInt:jArc ] ;
																arcIndexK = [ [ NSNumber alloc ] initWithInt:kArc ] ;
																worstTangent = [ [ NSNumber alloc ] initWithDouble:smallestTA ] ;
																
																arcInfo = [ [ NSArray alloc ] initWithObjects:buildArc, arcIndexJ, arcIndexK, worstTangent, nil] ;
																
																[ tryArcs addObject:arcInfo ] ;
																
																[ buildArc release ] ;
																[ arcIndexJ release ] ; [ arcIndexK release ] ; [ worstTangent release ] ;
																	
																// Next kArc
															}
														// Next jArc
													}
											
												if( [ tryArcs count ] > 0 )
													{
														combinedCycles = YES ;
														
														SMCycle *newCycle ;
														
														double q, bestQuality ;
														
														NSEnumerator *trialArcEnumerator ;
														NSArray *nextTrial, *bestTrial ;
														
														bestQuality = 0. ;
														bestTrial = nil ;
																										
														trialArcEnumerator = [ tryArcs objectEnumerator ] ;
														
														while( ( nextTrial = [ trialArcEnumerator nextObject ] ) )
															{
																
																	
																q = sin( [ [ nextTrial objectAtIndex:3 ] doubleValue ] ) ;
																
																if( q > bestQuality )
																	{
																		bestQuality = q ;
																		bestTrial = nextTrial ;
																	}
															}
															
														
														newCycle = [ self combineCycle:cycleJ firstIndex:[ [ bestTrial objectAtIndex:1 ] intValue ]
																		withCycle:cycleK secondIndex:[ [ bestTrial objectAtIndex:2 ] intValue ]
																		usingArc:[ bestTrial objectAtIndex:0 ] ] ;
																		
														[ contactCycles addObject:newCycle ] ;
														
														[ cycleJ killCycle ] ;
														[ cycleK killCycle ] ;
														
														[ atomsToCycles[iAtom] addObject:newCycle ] ;
																		
														// Break out of target cycle loop
														
														break ;
													}
											
												// Next target cycle
											}
									
										// If we combined a cycle, need to restart
										
										if( combinedCycles == YES ) break ;
										
									// Otherwise, next source cycle
									
								}
							
						// If NO cycles combined, exit
						
						if( combinedCycles == NO ) break ;
						
					}
					
				// Next atom
				
			}
					
		return ;
										
										
	}
										
/*			

// WARNING - The following algorithm would not really work, and on top of that should be obviated by the "twinned" arc approach implemented.
- (void) combineReentrantCycles
	{
		// This method combines reentrant cycles that share a skipped torus.
		

		SMTorus *nextTorus ;
		int i ;
		
		NSMutableDictionary *cycleToSkipTorus ;
		NSMutableArray *toriToProcess ;
		
		cycleToSkipTorus = [ [ NSMutableDictionary alloc ] initWithCapacity:100 ] ;
		toriToProcess = [ [ NSMutableArray alloc ] initWithCapacity:100 ] ;
		
		
		NSMutableSet *skipToriSet ;
		
		skipToriSet = [ [ NSMutableSet alloc ] initWithCapacity:nTori ] ;
		
		for( i = 0 ; i < nTori ; ++i )
			{
				nextTorus = tori[i] ;
				
				if( [ nextTorus skipTorus ] == NO ) continue ;
				
				[ skipToriSet addObject:nextTorus ] ;
				
				// Get cycles, see if collision with other skipped tori
				
				SMCycle *cycleR, *cycleL ;
				
				if( [ [ [ nextTorus reentrantArcR ] parentCycles ] count ] != 1 )
					{
						printf( "IMPOSSIBLE CYCLE COUNT %d FOR TORUS RIGHT ARC - Exit!\n", 
							[ [ [ nextTorus reentrantArcR ] parentCycles ] count ] ) ;
						exit(1) ;
					}
					
				cycleR = [ [ [ nextTorus reentrantArcR ] parentCycles ] objectAtIndex:0 ] ;
				
				if( [ [ [ nextTorus reentrantArcL ] parentCycles ] count ] != 1 )
					{
						printf( "IMPOSSIBLE CYCLE COUNT %d FOR TORUS LEFT ARC - Exit!\n", 
							[ [ [ nextTorus reentrantArcL ] parentCycles ] count ] ) ;
						exit(1) ;
					}
					
				cycleL = [ [ [ nextTorus reentrantArcL ] parentCycles ] objectAtIndex:0 ] ;
				
				SMTorus *collideTorus ;
				
				collideTorus = [ cycleToSkipTorus objectForKey:[ NSNumber numberWithUnsignedLong:(unsigned long)cycleR ] ] ;
				
				if( collideTorus )
					{
						// Remove collideTorus from tori to be processed - skip this one as well!
						
						[ toriToProcess removeObject:collideTorus ] ;
						
						// Skip current torus
						
						continue ;
					}
					
				
				collideTorus = [ cycleToSkipTorus objectForKey:[ NSNumber numberWithUnsignedLong:(unsigned long)cycleL ] ] ;
				
				if( collideTorus )
					{
						// Remove collideTorus from tori to be processed - skip this one as well!
						
						[ toriToProcess removeObject:collideTorus ] ;
						
						// Skip current torus
						
						continue ;
					}
				
				[ toriToProcess addObject:nextTorus ] ;
				
				[ cycleToSkipTorus setObject:nextTorus forKey:[ NSNumber numberWithUnsignedLong:(unsigned long)cycleR ] ] ;
				[ cycleToSkipTorus setObject:nextTorus forKey:[ NSNumber numberWithUnsignedLong:(unsigned long)cycleL ] ] ;
				
			}
			
		// Now, combine reetrant cycles for all surviving narrow tori
		
		// BE CAREFUL - if a narrow torus has been included in reentrant-cycle generation, it must ALSO be
		// included in the contact cycles!
		
		NSMutableSet *skipToriToKeep ;
		
		skipToriToKeep = [ [ NSMutableSet alloc ] initWithArray:toriToProcess ] ;
		
		[ skipToriSet minusSet:skipToriToKeep ] ;
		
		
		NSEnumerator *torusEnumerator ;
		
		// Anything left in skipToriSet should be set skip = NO 
		
		torusEnumerator = [ skipToriSet objectEnumerator ] ;
		
		while( ( nextTorus = [ torusEnumerator nextObject ] ) )
			{
				[ nextTorus setSkipTorus:NO ] ;
			}
			
		[ skipToriSet release ] ;
		[ skipToriToKeep release ] ;


		torusEnumerator = [ toriToProcess objectEnumerator ] ;
		
		while( ( nextTorus = [ torusEnumerator nextObject ] ) )
			{
				SMCycle *newCycle, *cycleR, *cycleL ;
				
				cycleR = [ [ [ nextTorus reentrantArcR ] parentCycles ] objectAtIndex:0 ] ;
				cycleL = [ [ [ nextTorus reentrantArcL ] parentCycles ] objectAtIndex:0 ] ;
				
				// Combine R and L using the contact arcs of the skipped torus
				
				newCycle = [ [ SMCycle alloc ] init ] ;
				
				// Begin new cycle using arc in cycleL that succeeds reentrantArcL
				
				SMArc *startArc ; int iNext, iL, iR ;
				
				iL = [ [ cycleL arcs ] indexOfObject:[ nextTorus reentrantArcL ] ] ;
				
				if( iL == NSNotFound )
					{
						printf( "REENTRANT CYCLE BROKEN - Exit!\n" ) ;
						exit(1) ;
					}
					
				iNext = iL + 1 ;
				
				if( iNext == [ [ cycleL arcs ] count ] )
					{
						iNext = 0 ;
					}
					
				// Keep adding arcs until reentrantArcL encountered
				
				while( TRUE )
					{
						SMArc *nextArc ;
						
						nextArc = [ [ cycleL arcs ] objectAtIndex:iNext ] ;
						
						if( nextArc == [ nextTorus reentrantArcL ] )
							{
								break ;
							}
							
						[ newCycle addArc:nextArc forward:[ [ cycleL forward ] objectAtIndex:iNext ] ] ;
						
						++iNext ;
						
						if( iNext == [ [ cycleL arcs ] count ] )
							{
								iNext = 0 ;
							}
					}
					
				// Add contactArcJ with reverse orientation
				
				[ newCycle addArc:[ nextTorus contactArcJ ] forward:[ NSNumber numberWithBool:NO ] ] ;
				
				// Add arcs from cycleR
				
				iR = [ [ cycleR arcs ] indexOfObject:[ nextTorus reentrantArcR ] ] ;
				
				if( iR == NSNotFound )
					{
						printf( "REENTRANT CYCLE BROKEN - Exit!\n" ) ;
						exit(1) ;
					}
					
				iNext = iR + 1 ;
				
				if( iNext == [ [ cycleR arcs ] count ] )
					{
						iNext = 0 ;
					}
					
				// Keep adding arcs until reentrantArcR encountered
				
				while( TRUE )
					{
						SMArc *nextArc ;
						
						nextArc = [ [ cycleR arcs ] objectAtIndex:iNext ] ;
						
						if( nextArc == [ nextTorus reentrantArcR ] )
							{
								break ;
							}
							
						[ newCycle addArc:nextArc forward:[ [ cycleR forward ] objectAtIndex:iNext ] ] ;
						
						++iNext ;
						
						if( iNext == [ [ cycleR arcs ] count ] )
							{
								iNext = 0 ;
							}
						
					}
					
				// Add contactArcI with reverse orientation
				
				[ newCycle addArc:[ nextTorus contactArcI ] forward:[ NSNumber numberWithBool:NO ] ] ;
				
				// Add new cycle
				
				[ reentrantCycles addObject:newCycle ] ;
				
				// Kill cycles cycleL and cycleR
				
				[ cycleR killCycle ] ;
				[ cycleL killCycle ] ;
				
				// Next torus 
				
			}
			
		// That's it
		
		[ cycleToSkipTorus release ] ;
		[ toriToProcess release ] ;
		
		return ;
		
	}
	
*/
				
		
	
	
- (BOOL) reduceContactCyclesUsingGeodesics:(BOOL)useGeo andDivisionParameter:(double)div
	{
		// This method reduces all cycles to three cycles, generating new arcs in the process
		
		// Outline of algorithm (in sort-of pseudocode) :
		//
		// NEXT:
		// find next cycle (currentCycle) > 3
		// START:
		//		M = # arcs in cycle
		//	
		//		tryArcs = {} 
		//		
		//		for j in 0 to M-3
		//			for k = j+2, M-1
		//				construct g = geo(j,k) from start of j to start of k
		//				
		//				if( j == 0 ) preArc = M - 1 else preArc = j - 1
		//				tangetOK = testTangent(g, j) && testTangent( j, preArc ) && testTangent( g, k ) && testTangent( g, k-1 )
		//				if( ! tangentOK ) continue 
		//				
		//				{intersectTest} = {M} \ {j, preArc, k, k-1 } 
		//				intersectOK = TRUE 
		//
		//				foreach a{intersectTest} intersectOK = intersectOK && testIntersect( g, a )
		//
		//				if( tangentOK && intersectOK )
		//					{tryArcs} = {tryArcs} U {g,j,k}
		//
		//			endfor(k)
		//		endfor(k)
		//
		//		if( {tryArcs} EMPTY )
		//			bigArc = longestArc({M})
		//			{arc1,arc2} = split( bigArc )
		//			replace( bigArc, {arc1,arc2} )
		//			goto START
		//		else
		//			qs = {quality({tryArcs})} (where quality depends on ratio of lengths of two sides of divided cycle, min tangent)
		//			index = best({qs})
		//			{newCycle1, newCycle2} = divideCycle( currentCycle, {tryArcs}[index]
		//			addCycles({newCycle1, newCycle2})
		//			killCycle(currentCycle)
		//			goto NEXT
		//		endif
		//			
		
		NSEnumerator *cycleEnumerator ;
		
		SMCycle *nextCycle ;
		
		BOOL hadReduction ;
		
 		hadReduction = NO ;
		
		
		// Start by subdividing all existing arcs
		
		NSEnumerator *arcEnumerator ;
		SMArc *nextArc ;
		
		static NSMutableArray *arcsToSubdivide = nil ;
		
		if( ! arcsToSubdivide )
			{
				arcsToSubdivide = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
			}
		
		cycleEnumerator = [ contactCycles objectEnumerator ] ;
		
		while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
			{
				
				
				// Note that if we have executed combine cycle operations, some cycles may already be dead
				
				if( [ nextCycle active ] == NO ) continue ;
				
				// Need some extra care in this seemingly simple operation. If contact cycles have been combined into one, we will have an arc
				// AND IT'S TWIN in the same cycle! When one of the pair is encountered, an arc will be removed from the cycle. 
				
				// Initially, I will simply verify that each arc is still present in the cycle 
				
				arcEnumerator = [ [ nextCycle arcs ] objectEnumerator ] ;
				
				[ arcsToSubdivide removeAllObjects ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						// Check for no parent cycle > nextCycle 
						
						BOOL skip ;
						NSEnumerator *parentCycleEnumerator ;
						SMCycle *nextParentCycle ;
						
						skip = NO ;
						parentCycleEnumerator = [ [ nextArc parentCycles ] objectEnumerator ] ;
						
						while( ( nextParentCycle = [ parentCycleEnumerator nextObject ] ) )
							{
								if( nextParentCycle > nextCycle )
									{
										skip = YES ;
										break ;
									}
							}
							
						if( skip == YES ) continue ;
						
						// Make sure the arc is still there! (Not removed owing to subdivision of a twin in the same cycle!)
						
						if( [ [ nextCycle arcs ] indexOfObject:nextArc ] == NSNotFound ) continue ;
						
						[ arcsToSubdivide addObject:nextArc ] ;
						
						//[ self subdivideArc:nextArc usingDivision:div ] ;
					}
			}
		
		arcEnumerator = [ arcsToSubdivide objectEnumerator ] ;
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				[ self subdivideArc:nextArc usingDivision:div ] ;
			}
		
		NSMutableArray *tryArcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
		
		
#ifdef SURFDEBUG
		NSMutableArray *newCycles = [ [ NSMutableArray alloc ] initWithCapacity:[ contactCycles count ] ] ;
#else
		NSMutableArray *newCycles ;
#endif
		
		int subdivideCount = 0 ;
		
		while( TRUE )
			{
				NSEnumerator *cycleEnumerator ;
				int M ;
				
				SMCycle *nextCycle ;
				
#ifdef SURFDEBUG
				[ newCycles removeAllObjects ] ;
#else
				newCycles = [ [ NSMutableArray alloc ] initWithCapacity:[ contactCycles count ] ] ;
#endif


				BOOL newCycleCreated = NO ;
				
				// For now I am going to hardwire the minimal angle in the tangent test. I will make this 1 degrees
				
				double minTangentAngle = 1. * acos(-1.)/180. ;
				
				cycleEnumerator = [ contactCycles objectEnumerator ] ;
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						if( [ nextCycle active ] == NO ) continue ;
						
						
						 
						NSArray *cycleArcs ;
						NSArray *cycleForward ;
						
						cycleArcs = [ nextCycle arcs ] ;
						cycleForward = [ nextCycle forward ] ;
					
						NSEnumerator *arcEnumerator ;
						SMArc *nextArc ;
							
						int MNoSkip ;
						
						MNoSkip = 0 ;
						
						arcEnumerator = [ cycleArcs objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								if( [ nextArc skip ] == NO  )
									{
										++MNoSkip ;
									}
							}
		
					
						if( MNoSkip == 3 ) 
							{
#ifndef SURFDEBUG
								[ newCycles addObject:nextCycle ] ;
#endif
								continue ;
							}
						
						int iArc ;
						
						if( MNoSkip == 2 )
							{
								// Subdivide both available arcs
								
								//arcEnumerator = [ cycleArcs objectEnumerator ] ;
								
								int initialCount = [ cycleArcs count ] ;
								
								for( iArc = 0 ; iArc < initialCount ; ++iArc )
								//while( ( nextArc = [ arcEnumerator nextObject ] ) )
									{
										nextArc = [ cycleArcs objectAtIndex:iArc ] ;
										
										if( [ nextArc skip ] == NO )
											{
												[ self subdivideArc:nextArc ] ;
											}
									}
								
							}
						else if( MNoSkip < 2 )
							{
								printf( "FATAL ERROR - CONTACT CYCLE HAS %d UNSKIPPED ARCS - Exit!\n", MNoSkip ) ;
								exit(1) ;
							}
								
							
										
						M = [ cycleArcs count ] ;
						
					
						BOOL success ;
						
						success = NO ;
						subdivideCount = 0 ;
						
						while( success == NO )
							{
					
																
								[ tryArcs removeAllObjects ] ;
								
								int j, k, preArc ;
								SMArc *buildArc ;
								
								
								// Need to recompute size of cycle, as we may have subdivided
								
								M = [ cycleArcs count ] ;
						
						
								for( j = 0 ; j <= M - 3 ; ++j )
									{
										if( [ [ cycleArcs objectAtIndex:j ] skip ]  == YES )
											{
												continue ;
											}
											
										for( k = j + 2 ; k <= M - 1 ; ++k )
											{
												if( [ [ cycleArcs objectAtIndex:k ] skip ] == YES )
													{
														continue ;
													}
													
												// Build arc from start of j to start of k
												
												// This should really be a class method of SMArc, but I only use it in this class
												
												if( j == 0 )
													{
														preArc = M - 1 ;
													}
												else
													{
														preArc = j - 1 ;
													}
												
												if( useGeo == YES )
													{
														buildArc = [ self buildGeoUsingStart:[ cycleArcs objectAtIndex:j ] 
															startForward:[ [ cycleForward objectAtIndex:j ] boolValue ] 
															end:[ cycleArcs objectAtIndex:k ]
															endForward:[ [ cycleForward objectAtIndex:k ] boolValue ] reentrantRegion:NO ] ;
													}
												else
													{
														buildArc = [ self buildCircularArcUsingSourcePre:[ cycleArcs objectAtIndex:preArc ] 
															sourcePreForward:[ [ cycleForward objectAtIndex:preArc ] boolValue ] 
															sourceNext:[ cycleArcs objectAtIndex:j ] 
															sourceNextForward:[ [ cycleForward objectAtIndex:j ] boolValue ]
															targetPre:[ cycleArcs objectAtIndex:(k - 1) ]
															targetPreForward:[ [ cycleForward objectAtIndex:(k - 1) ] boolValue ]
															targetNext:[ cycleArcs objectAtIndex:k ] 
															targetNextForward:[ [ cycleForward objectAtIndex:k ] boolValue ] ] ;
													}
													
												// All arcs in the cycle should have same host atom 
												
												if( ! buildArc ) continue ;
												
												
												// Reject an arc with length identically zero (can be generated with combined contact cycles - adjacent arcs with 
												// opposite orientations)
												
												if( [ buildArc length ] == 0. )
													{
														[ buildArc release ] ;
														continue ;
													}
												
												double tA, smallestTA ;
												
												smallestTA = acos(-1.) ;
												
													
												tA = [ self computeTangentAngleUsingArc:buildArc start:YES previousArc:[ cycleArcs objectAtIndex:preArc ]
														previousForward:[ [ cycleForward objectAtIndex:preArc ] boolValue ]
														nextArc:[ cycleArcs objectAtIndex:j ]  nextForward:[ [ cycleForward objectAtIndex:j ] boolValue ]
														reentrantRegion:NO ];
												
												if( tA < smallestTA ) smallestTA = tA ;
												
												tA = [ self computeTangentAngleUsingArc:buildArc start:NO previousArc:[ cycleArcs objectAtIndex:(k - 1) ]
														previousForward:[ [ cycleForward objectAtIndex:(k - 1) ] boolValue ]
														nextArc:[ cycleArcs objectAtIndex:k ] nextForward:[ [ cycleForward objectAtIndex:k ] boolValue ]	
														reentrantRegion:NO ] ;
												
												if( tA < smallestTA ) smallestTA = tA ;
										
														
												if( smallestTA < minTangentAngle )
													{
														[ buildArc release ] ;
														
														continue ;
													}
													
												// Compute intersections
												
												
												SMArc *testArc ;
												int arcIndex ;
												BOOL haveIntersection ;
												
												haveIntersection = NO ;
												
												arcEnumerator = [ cycleArcs objectEnumerator ] ;
												
												arcIndex = -1 ;
												
												while( ( testArc = [ arcEnumerator nextObject ] ) )
													{
														++arcIndex ;
														
														if( arcIndex == preArc || arcIndex == j || arcIndex == k - 1 || arcIndex == k )
															{
																continue ;
															}
															
														if( [ buildArc intersectWith2:testArc ] == YES )
															{
																haveIntersection = YES ;
																break ;
															}
													}
													
												if( haveIntersection == YES )
													{
														[ buildArc release ] ;
														
														continue ;
													}
													
												NSArray *arcInfo ;
												NSNumber *arcIndexJ, *arcIndexK, *worstTangent ;
												
												arcIndexJ = [ [ NSNumber alloc ] initWithInt:j ] ;
												arcIndexK = [ [ NSNumber alloc ] initWithInt:k ] ;
												worstTangent = [ [ NSNumber alloc ] initWithDouble:smallestTA ] ;
												
												arcInfo = [ [ NSArray alloc ] initWithObjects:buildArc, arcIndexJ, arcIndexK, worstTangent, nil] ;
												
												[ tryArcs addObject:arcInfo ] ;
												
												[ buildArc release ] ;
												[ arcIndexJ release ] ; [ arcIndexK release ] ; [ worstTangent release ] ;
												
											}
									}
									
								// Now choose the best arc
								
								if( [ tryArcs count ] == 0 )
									{
										if( subdivideCount == MAXSUBDIVIDECOUNT )
											{
												printf( "CONTACT CYCLE COULD NOT BE REDUCED AFTER %d SUBDIVISIONS\n", MAXSUBDIVIDECOUNT ) ;
												break ;
											}
											
										// Find longest arc
										
										SMArc *theLongest ;
										double maxLength ;
										
										maxLength = 0. ;
										
										arcEnumerator = [ cycleArcs objectEnumerator ] ;
										
										while( ( nextArc = [ arcEnumerator nextObject ] ) )
											{
												double l ;
												
												l = [ nextArc length ] ;
												
												if( l > maxLength )
													{
														maxLength = l ;
														theLongest = nextArc ;
													}
													
											}
											
										[ self subdivideArc:theLongest ] ;
										
										++subdivideCount ;
										
										continue ;
									}
									
									
								// OK, we have at least one valid subdivision
									
								// Select best-quality subdivision
								
								double q, ratio, bestQuality ;
								
								NSEnumerator *trialArcEnumerator ;
								NSArray *nextTrial, *bestTrial ;
								
								bestQuality = 0. ;
								bestTrial = nil ;
								
								double perimeterJ, perimeterK ;
								
								trialArcEnumerator = [ tryArcs objectEnumerator ] ;
								
								while( ( nextTrial = [ trialArcEnumerator nextObject ] ) )
									{
										j = [ [ nextTrial objectAtIndex:1 ] intValue ] ;
										k = [ [ nextTrial objectAtIndex:2 ] intValue ] ;
										
										perimeterJ = perimeterK = 0. ;
										int arcIndex ;
										
										arcIndex = 0 ;
										
										arcEnumerator = [ cycleArcs objectEnumerator ] ;
										
										while( ( nextArc = [ arcEnumerator nextObject ] ) )
											{
												if( arcIndex < j || arcIndex >= k )
													{
														perimeterK += [ nextArc length ] ;
													}
												else
													{
														perimeterJ += [ nextArc length ] ;
													}
													
												++arcIndex ;
											}
											
										ratio = MIN(perimeterJ,perimeterK)/MAX(perimeterJ,perimeterK) ;
										
										q = ratio + sin( [ [ nextTrial objectAtIndex:3 ] doubleValue ] ) ;
										
										if( q > bestQuality )
											{
												bestQuality = q ;
												bestTrial = nextTrial ;
											}
									}
									
								// Subdivide cycle using best quality division
								
								NSArray *subCycles ;
								
								subCycles = [ self subdivideCycle:nextCycle firstIndex:[ [ bestTrial objectAtIndex:1 ] intValue ]
									secondIndex:[ [ bestTrial objectAtIndex:2 ] intValue ] usingArc:[ bestTrial objectAtIndex:0 ] ] ;
									
								// Now that the arc is embedded in new cycles, subdivide it!
								
								[ self subdivideArc:[ bestTrial objectAtIndex:0 ] usingDivision:div ] ;
																
								[ nextCycle killCycle ] ;
								
								[ newCycles addObjectsFromArray:subCycles ] ;
								newCycleCreated = YES ;
								
								[ subCycles release ] ;
								
								success = YES ;
								
								hadReduction = YES ;
								
							}
							
						if( success == NO )
							{
								printf( "COULD NOT REDUCE CONTACT CYCLE - Exit!\n" ) ;
								exit(1) ;
							}
							
						// Next cycle 
					}
					
				// Did we reduce any cycles?
				
#ifdef SURFDEBUG
				if( [ newCycles count ] == 0 ) break ;
				 
				[ contactCycles addObjectsFromArray:newCycles ] ;
#else
				if( newCycleCreated == NO )
				//if( [ newCycles count ] == [ contactCycles count ] ) 
					{
						[ newCycles release ] ;
						break ;
					}
				
				[ contactCycles release ] ;
				contactCycles = newCycles ;
#endif
			}
			
			
		return hadReduction ;
	}
								
								
- (BOOL) reduceReentrantCyclesUsingDivisionParameter:(double)div
	{
		// This method reduces all cycles to three cycles, generating new arcs in the process
		
		// Outline of algorithm (in sort-of pseudocode) :
		//
		// NEXT:
		// find next cycle (currentCycle) > 3
		// START:
		//		M = # arcs in cycle
		//	
		//		tryArcs = {} 
		//		
		//		for j in 0 to M-3
		//			for k = j+2, M-1
		//				construct g = geo(j,k) from start of j to start of k
		//				
		//				if( j == 0 ) preArc = M - 1 else preArc = j - 1
		//				tangetOK = testTangent(g, j) && testTangent( j, preArc ) && testTangent( g, k ) && testTangent( g, k-1 )
		//				if( ! tangentOK ) continue 
		//				
		//				{intersectTest} = {M} \ {j, preArc, k, k-1 } 
		//				intersectOK = TRUE 
		//
		//				foreach a{intersectTest} intersectOK = intersectOK && testIntersect( g, a )
		//
		//				if( tangentOK && intersectOK )
		//					{tryArcs} = {tryArcs} U {g,j,k}
		//
		//			endfor(k)
		//		endfor(k)
		//
		//		if( {tryArcs} EMPTY )
		//			bigArc = longestArc({M})
		//			{arc1,arc2} = split( bigArc )
		//			replace( bigArc, {arc1,arc2} )
		//			goto START
		//		else
		//			qs = {quality({tryArcs})} (where quality depends on ratio of lengths of two sides of divided cycle, min tangent)
		//			index = best({qs})
		//			{newCycle1, newCycle2} = divideCycle( currentCycle, {tryArcs}[index]
		//			addCycles({newCycle1, newCycle2})
		//			killCycle(currentCycle)
		//			goto NEXT
		//		endif
		//			
		
		NSEnumerator *cycleEnumerator ;
		
		SMCycle *nextCycle ;
		
		BOOL hadReduction = NO ;
		
		
		// Start by subdividing all existing arcs
		
		cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
		
		static NSMutableArray *arcsToSubdivide = nil ;
		
		if( ! arcsToSubdivide )
			{
				arcsToSubdivide = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
			}
		
		NSEnumerator *arcEnumerator ;
		SMArc *nextArc ;
		
		while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
			{
				
				[ arcsToSubdivide removeAllObjects ] ;
				
				arcEnumerator = [ [ nextCycle arcs ] objectEnumerator ] ;
				
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						// Check for no parent cycle > nextCycle 
						
						BOOL skip ;
						NSEnumerator *parentCycleEnumerator ;
						SMCycle *nextParentCycle ;
						
						skip = NO ;
						parentCycleEnumerator = [ [ nextArc parentCycles ] objectEnumerator ] ;
						
						while( ( nextParentCycle = [ parentCycleEnumerator nextObject ] ) )
							{
								if( nextParentCycle > nextCycle )
									{
										skip = YES ;
										break ;
									}
							}
							
						if( skip == YES ) continue ;
						
						[ arcsToSubdivide addObject:nextArc ] ;
						
						//[ self subdivideArc:nextArc usingDivision:div ] ;
					}
			}
		
		arcEnumerator = [ arcsToSubdivide objectEnumerator ] ;
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				[ self subdivideArc:nextArc usingDivision:div ] ;
			}
		
		NSMutableArray *tryArcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
		
#ifdef SURFDEBUG
		NSMutableArray *newCycles = [ [ NSMutableArray alloc ] initWithCapacity:[ reentrantCycles count ] ] ;
#else
		NSMutableArray *newCycles ;
#endif
		
		 
		 
		int subdivideCount = 0 ;
		
		while( TRUE )
			{
				NSEnumerator *cycleEnumerator ;
				int M ;
				
				SMCycle *nextCycle ;
				
#ifdef SURFDEBUG
				[ newCycles removeAllObjects ] ;
#else
				newCycles = [ [ NSMutableArray alloc ] initWithCapacity:[ reentrantCycles count ] ] ;
#endif
				
				
				BOOL newCycleCreated = NO ;
				
				// For now I am going to hardwire the minimal angle in the tangent test. I will make this 1 degrees
				
				double minTangentAngle = 1. * acos(-1.)/180. ;
				
				cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						if( [ nextCycle active ] == NO ) continue ;
						
						 
						NSArray *cycleArcs ;
						NSArray *cycleForward ;
						
						cycleArcs = [ nextCycle arcs ] ;
						cycleForward = [ nextCycle forward ] ;
					
						NSEnumerator *arcEnumerator ;
						SMArc *nextArc ;
							
						int MNoSkip ;
						
						MNoSkip = 0 ;
						
						arcEnumerator = [ cycleArcs objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								if( [ nextArc skip ] == NO )
									{
										++MNoSkip ;
									}
									
							}
		
					
						if( MNoSkip <= 3 ) 
							{
								
#ifndef SURFDEBUG
								[ newCycles addObject:nextCycle ] ;
#endif
								
								continue ;
							}
										
						M = [ cycleArcs count ] ;
						
					
						BOOL success ;
						
						success = NO ;
						subdivideCount = 0 ;
						
						while( success == NO )
							{
					
																
								[ tryArcs removeAllObjects ] ;
								
								int j, k, preArc ;
								SMArc *buildArc ;
								
								
								// Need to recompute size of cycle, as we may have subdivided
								
								M = [ cycleArcs count ] ;
								
								SMTorus *theTorus ;
								SMArc *cI, *cJ ;
						
						
								for( j = 0 ; j <= M - 3 ; ++j )
									{
										if( [ [ cycleArcs objectAtIndex:j ] torusSection ] != nil )
											{
												theTorus = [ [ cycleArcs objectAtIndex:j ] torusSection ] ;
												cI = [ theTorus contactArcI ] ;
												cJ = [ theTorus contactArcJ ] ;
												
												if( [ cycleArcs objectAtIndex:j ] == cI ||
													[ cycleArcs objectAtIndex:j ] == cJ ) 
													{
														continue ;
													}
											}
											
										for( k = j + 2 ; k <= M - 1 ; ++k )
											{
												if( [ [ cycleArcs objectAtIndex:k ] torusSection ] != nil )
													{
														theTorus = [ [ cycleArcs objectAtIndex:k ] torusSection ] ;
														cI = [ theTorus contactArcI ] ;
														cJ = [ theTorus contactArcJ ] ;
														
														if( [ cycleArcs objectAtIndex:k ] == cI ||
															[ cycleArcs objectAtIndex:k ] == cJ ) 
															{
																continue ;
															}
													}
													
												// Build arc from start of j to start of k
												
												// This should really be a class method of SMArc, but I only use it in this class
												
												if( j == 0 )
													{
														preArc = M - 1 ;
													}
												else
													{
														preArc = j - 1 ;
													}
												
												buildArc = [ self buildGeoUsingStart:[ cycleArcs objectAtIndex:j ] 
													startForward:[ [ cycleForward objectAtIndex:j ] boolValue ] 
													end:[ cycleArcs objectAtIndex:k ]
													endForward:[ [ cycleForward objectAtIndex:k ] boolValue ] reentrantRegion:YES ] ;
												
												if( ! buildArc ) continue ;
												
												double tA, smallestTA ;
												
												smallestTA = acos(-1.) ;
												
													
												tA = [ self computeTangentAngleUsingArc:buildArc start:YES previousArc:[ cycleArcs objectAtIndex:preArc ]
														previousForward:[ [ cycleForward objectAtIndex:preArc ] boolValue ]
														nextArc:[ cycleArcs objectAtIndex:j ]  nextForward:[ [ cycleForward objectAtIndex:j ] boolValue ] 
														reentrantRegion:YES ];
												
												if( tA < smallestTA ) smallestTA = tA ;
												
												tA = [ self computeTangentAngleUsingArc:buildArc start:NO previousArc:[ cycleArcs objectAtIndex:(k - 1) ]
														previousForward:[ [ cycleForward objectAtIndex:(k - 1) ] boolValue ]
														nextArc:[ cycleArcs objectAtIndex:k ] nextForward:[ [ cycleForward objectAtIndex:k ] boolValue ]	
														reentrantRegion:YES ] ;
												
												if( tA < smallestTA ) smallestTA = tA ;
										
														
												if( smallestTA < minTangentAngle )
													{
														[ buildArc release ] ;
														
														continue ;
													}
													
												// Compute intersections
												
												
												SMArc *testArc ;
												int arcIndex ;
												BOOL haveIntersection ;
												
												haveIntersection = NO ;
												
												arcEnumerator = [ cycleArcs objectEnumerator ] ;
												
												arcIndex = -1 ;
												
												while( ( testArc = [ arcEnumerator nextObject ] ) )
													{
														++arcIndex ;
														
														if( arcIndex == preArc || arcIndex == j || arcIndex == k - 1 || arcIndex == k )
															{
																continue ;
															}
															
														if( [ buildArc intersectWith2:testArc ] == YES )
															{
																haveIntersection = YES ;
																break ;
															}
													}
													
												if( haveIntersection == YES )
													{
														[ buildArc release ] ;
														
														continue ;
													}
													
												NSArray *arcInfo ;
												NSNumber *arcIndexJ, *arcIndexK, *worstTangent ;
												
												arcIndexJ = [ [ NSNumber alloc ] initWithInt:j ] ;
												arcIndexK = [ [ NSNumber alloc ] initWithInt:k ] ;
												worstTangent = [ [ NSNumber alloc ] initWithDouble:smallestTA ] ;
												
												arcInfo = [ [ NSArray alloc ] initWithObjects:buildArc, arcIndexJ, arcIndexK, worstTangent, nil] ;
												
												[ tryArcs addObject:arcInfo ] ;
												
												[ buildArc release ] ;
												[ arcIndexJ release ] ; [ arcIndexK release ] ; [ worstTangent release ] ;
												
											}
									}
									
								// Now choose the best arc
								
								if( [ tryArcs count ] == 0 )
									{
										if( subdivideCount == MAXSUBDIVIDECOUNT )
											{
												printf( "REENTRANT CYCLE COULD NOT BE REDUCED AFTER %d SUBDIVISIONS\n", MAXSUBDIVIDECOUNT ) ;
												break ;
											}
											
										// Find longest arc
										
										SMArc *theLongest ;
										double maxLength ;
										
										maxLength = 0. ;
										
										arcEnumerator = [ cycleArcs objectEnumerator ] ;
										
										while( ( nextArc = [ arcEnumerator nextObject ] ) )
											{
												double l ;
												
												l = [ nextArc length ] ;
												
												if( l > maxLength )
													{
														maxLength = l ;
														theLongest = nextArc ;
													}
													
											}
											
										[ self subdivideArc:theLongest ] ;
										
										++subdivideCount ;
										
										continue ;
									}
									
									
								// OK, we have at least one valid subdivision
									
								// Select best-quality subdivision
								
								double q, ratio, bestQuality ;
								
								NSEnumerator *trialArcEnumerator ;
								NSArray *nextTrial, *bestTrial ;
								
								bestQuality = 0. ;
								bestTrial = nil ;
								
								double perimeterJ, perimeterK ;
								
								trialArcEnumerator = [ tryArcs objectEnumerator ] ;
								
								while( ( nextTrial = [ trialArcEnumerator nextObject ] ) )
									{
										j = [ [ nextTrial objectAtIndex:1 ] intValue ] ;
										k = [ [ nextTrial objectAtIndex:2 ] intValue ] ;
										
										perimeterJ = perimeterK = 0. ;
										int arcIndex ;
										
										arcIndex = 0 ;
										
										arcEnumerator = [ cycleArcs objectEnumerator ] ;
										
										while( ( nextArc = [ arcEnumerator nextObject ] ) )
											{
												if( arcIndex < j || arcIndex >= k )
													{
														perimeterK += [ nextArc length ] ;
													}
												else
													{
														perimeterJ += [ nextArc length ] ;
													}
													
												++arcIndex ;
											}
											
										ratio = MIN(perimeterJ,perimeterK)/MAX(perimeterJ,perimeterK) ;
										
										q = ratio + sin( [ [ nextTrial objectAtIndex:3 ] doubleValue ] ) ;
										
										if( q > bestQuality )
											{
												bestQuality = q ;
												bestTrial = nextTrial ;
											}
									}
									
								// Subdivide cycle using best quality division
								
								NSArray *subCycles ;
								
								subCycles = [ self subdivideCycle:nextCycle firstIndex:[ [ bestTrial objectAtIndex:1 ] intValue ]
									secondIndex:[ [ bestTrial objectAtIndex:2 ] intValue ] usingArc:[ bestTrial objectAtIndex:0 ] ] ;
									
								// Now that the arc is embedded in new cycles, subdivide it!
								
								[ self subdivideArc:[ bestTrial objectAtIndex:0 ] usingDivision:div ] ;
																
								[ nextCycle killCycle ] ;
								
								[ newCycles addObjectsFromArray:subCycles ] ;
								newCycleCreated = YES ;
								
								[ subCycles release ] ;
								
								success = YES ;
								
								hadReduction = YES ;
								
							}
							
						if( success == NO )
							{
								printf( "COULD NOT REDUCE REENTRANT CYCLE - Exit!\n" ) ;
								exit(1) ;
							}
							
						// Next cycle 
					}
					
				// Did we reduce any cycles?
#ifdef SURFDEBUG
				if( [ newCycles count ] == 0 ) break ;
				
				[ reentrantCycles addObjectsFromArray:newCycles ] ;
#else
				
				//if( [ newCycles count ] == [ reentrantCycles count ] ) 
				if( newCycleCreated == NO )

					{
						[ newCycles release ] ;
						break ;
					}
				
				[ reentrantCycles release ] ;
				reentrantCycles = newCycles ;
#endif
				
			}
			
		
		return hadReduction ;
	}
								
								

- (double) computeTangentAngleUsingArc:(SMArc *) a start:(BOOL) s previousArc:(SMArc *) p previousForward:(BOOL)pf
				nextArc:(SMArc *) n nextForward:(BOOL)nf reentrantRegion:(BOOL)reent
	{
		double angle, minAngle ;
		
		MMVector3 *normal, *center, *vertexPosition ;
		
		center = [ n hostCenter ] ;
		
		if( nf == YES )
			{
				vertexPosition = [ n startPosition ] ;
			}
		else
			{
				vertexPosition = [ n endPosition ] ;
			}
		
		normal = [ [ MMVector3 alloc ] initX:([ vertexPosition X ] - [ center X ]) 
			Y:([ vertexPosition Y ] - [ center Y ])
			Z:([ vertexPosition Z ] - [ center Z ]) ] ;
			
		[ normal normalize ] ;
		
		minAngle = acos(-1.) ;
		
		
		MMVector3 *testCross ;
		double dot, dot2 ;
		
		
		// To simplify the code...
		
		MMVector3 *useTangent ;
		
		if( s == YES )
			{
				useTangent = [ [ MMVector3 alloc ] initX:[ [ a startTangent ] X ] Y:[ [ a startTangent ] Y ] Z:[ [ a startTangent ] Z ] ] ;
			}
		else
			{
				useTangent = [ [ MMVector3 alloc ] initX:(-[ [ a endTangent ] X ]) Y:(-[ [ a endTangent ] Y ]) Z:(-[ [ a endTangent ] Z ]) ] ;
			}
		
		
		if( nf == YES )
			{
				testCross = [ [ MMVector3 alloc ] initByCrossing:useTangent and:[ n startTangent ] ] ;
			}
		else
			{
				testCross = [ [ MMVector3 alloc ] initByCrossing:useTangent and:[ n endTangent ] ] ;
				[ testCross reverse ] ;
			}
		
		dot = [ normal dotWith:testCross ] ;
		
		if( reent == YES )
			{
				dot = -dot ;
			}
		
		if( nf == YES )
			{
				dot2 = [ useTangent dotWith:[ n startTangent ] ] ;
			}
		else
			{
				dot2 = -[ useTangent dotWith:[ n endTangent ] ] ;
			}
		
		if( fabs(dot2) > 1.0 ) dot2 = dot2/fabs(dot2) ;
		
		angle = acos( dot2 ) ;
		
		
		if( dot < 0. )
			{
				angle = -angle ;
			}
			
		if( angle < minAngle ) minAngle = angle ;
		
		// Other side ---
		
		// A little lame to repeat this code, but it is short...
		
		[ testCross release ] ;
		
		if( pf == YES )
			{
				testCross = [ [ MMVector3 alloc ] initByCrossing:useTangent and:[ p endTangent ] ] ;
			}
		else
			{
				testCross = [ [ MMVector3 alloc ] initByCrossing:useTangent and:[ p startTangent ] ] ;
				[ testCross reverse ] ;
			}
		
		dot = [ normal dotWith:testCross ] ;
		
		if( reent == YES )
			{
				dot = -dot ;
			}

		
		if( pf == YES )
			{
				dot2 = -[ useTangent dotWith:[ p endTangent ] ] ;
			}
		else
			{
				dot2 = [ useTangent dotWith:[ p startTangent ] ] ;
			}
		
		
		if( fabs(dot2) > 1.0 ) dot2 = dot2/fabs(dot2) ;
		
		angle = acos( dot2 ) ;
		
		
		if( dot < 0. )
			{
				angle = -angle ;
			}
		
			
		if( angle < minAngle ) minAngle = angle ;
		
		[ testCross release ] ;
		[ useTangent release ] ;
		
		return minAngle ;
	}
		
- (double) computeSaddleTangentAngleUsingArc:(SMArc *) a start:(BOOL) s previousArc:(SMArc *) p
														previousForward:(BOOL)pf
														nextArc:(SMArc *)n  nextForward:(BOOL)nf 
	{
		// This method computes the tangent in the phi-theta plane. This is not the angle that would be found in R3
		
		typedef struct { double tComp ; double pComp ; } saddleVec ;
		
		saddleVec aTangent ;
		double size ;
		
		if( s == YES )
			{
				aTangent.tComp = [ a thetaEnd ] - [ a thetaStart ] ;
				aTangent.pComp = [ a phiEnd ] - [ a phiStart ] ;
			}
		else
			{
				aTangent.tComp = [ a thetaStart ] - [ a thetaEnd ] ;
				aTangent.pComp = [ a phiStart ] - [ a phiEnd] ;
			}
			
		size = sqrt( aTangent.tComp*aTangent.tComp + aTangent.pComp*aTangent.pComp ) ;
		
		aTangent.tComp /= size ;
		aTangent.pComp /= size ;
		
		saddleVec preTangent ;
		
		if( pf == NO )
			{
				preTangent.tComp = [ p thetaEnd ] - [ p thetaStart ] ;
				preTangent.pComp = [ p phiEnd ] - [ p phiStart ] ;
			}
		else
			{
				preTangent.tComp = [ p thetaStart ] - [ p thetaEnd ] ;
				preTangent.pComp = [ p phiStart ] - [ p phiEnd] ;
			}
		
		size = sqrt( preTangent.tComp*preTangent.tComp + preTangent.pComp*preTangent.pComp ) ;
		
		preTangent.tComp /= size ;
		preTangent.pComp /= size ;
		
		saddleVec nextTangent ;

		if( nf == YES )
			{
				nextTangent.tComp = [ n thetaEnd ] - [ n thetaStart ] ;
				nextTangent.pComp = [ n phiEnd ] - [ n phiStart ] ;
			}
		else
			{
				nextTangent.tComp = [ n thetaStart ] - [ n thetaEnd ] ;
				nextTangent.pComp = [ n phiStart ] - [ n phiEnd] ;
			}
			
		size = sqrt( nextTangent.tComp*nextTangent.tComp + nextTangent.pComp*nextTangent.pComp ) ;
		
		nextTangent.tComp /= size ;
		nextTangent.pComp /= size ;
			
		double dot, minAngle, angle ;
		
		minAngle = acos(-1.) ;
		
		// With next arc 
		
		dot = aTangent.tComp*nextTangent.tComp + aTangent.pComp*nextTangent.pComp ;
		
		if( fabs(dot) > 1. ) dot = dot/fabs(dot) ;
		
		angle = acos(dot) ;
		
		if( angle < minAngle ) minAngle = angle ;
		
		// With previous arc 
		
		dot = aTangent.tComp*preTangent.tComp + aTangent.pComp*preTangent.pComp ;
		
		if( fabs(dot) > 1. ) dot = dot/fabs(dot) ;
		
		angle = acos(dot) ;
		
		if( angle < minAngle ) minAngle = angle ;
		
		return minAngle ;
	}
		
		
		
		
- (SMArc *) buildGeoUsingStart:(SMArc *)sourceArc startForward:(BOOL)sf end:(SMArc *)targetArc endForward:(BOOL)ef reentrantRegion:(BOOL)reent
	{
		// Build a geodesic arc from beginning of start arc to beginning of end arc
		
		MMVector3 *center, *start, *end, *startVertexPosition, *endVertexPosition ;
		
		SMArc *returnArc ;
		
		center = [ sourceArc hostCenter ] ;
		
		MMVector3 *startTangent, *endTangent ;
		
		if( sf == YES )
			{
				start = [ [ MMVector3 alloc ] initX:([ [sourceArc startPosition ] X ] - [ center X ])
								Y:([ [sourceArc startPosition ] Y ] - [ center Y ])
								Z:([ [sourceArc startPosition ] Z ] - [ center Z ]) ] ;
								
				startVertexPosition = [ [ MMVector3 alloc ] initX:[ [sourceArc startPosition ] X ]
											Y:[ [sourceArc startPosition ] Y ]
											Z:[ [sourceArc startPosition ] Z ] ] ; 
											
				startTangent = [ [ MMVector3 alloc ] initUsingVector:[ sourceArc startTangent ] ] ;
								
			}
		else
			{
				start = [ [ MMVector3 alloc ] initX:([ [sourceArc endPosition ] X ] - [ center X ])
								Y:([ [sourceArc endPosition ] Y ] - [ center Y ])
								Z:([ [sourceArc endPosition ] Z ] - [ center Z ]) ] ;
								
				startVertexPosition = [ [ MMVector3 alloc ] initX:[ [sourceArc endPosition ] X ]
											Y:[ [sourceArc endPosition ] Y ]
											Z:[ [sourceArc endPosition ] Z ] ] ; 
											
				startTangent = [ [ MMVector3 alloc ] initUsingVector:[ sourceArc endTangent ] ] ;
				
				[ startTangent reverse ] ;
								
			}
						
		[ start normalize ] ;
		
		if( ef == YES )
			{
				end = [ [ MMVector3 alloc ] initX:([ [targetArc startPosition ] X ] - [ center X ])
								Y:([ [targetArc startPosition ] Y ] - [ center Y ])
								Z:([ [targetArc startPosition ] Z ] - [ center Z ]) ] ;
								
				endVertexPosition = [ [ MMVector3 alloc ] initX:[ [targetArc startPosition ] X ]
											Y:[ [targetArc startPosition ] Y ]
											Z:[ [targetArc startPosition ] Z ] ] ; 
											
				endTangent = [ [ MMVector3 alloc ] initUsingVector:[ targetArc startTangent ] ] ;
			}
		else
			{
				end = [ [ MMVector3 alloc ] initX:([ [targetArc endPosition ] X ] - [ center X ])
								Y:([ [targetArc endPosition ] Y ] - [ center Y ])
								Z:([ [targetArc endPosition ] Z ] - [ center Z ]) ] ;
								
				endVertexPosition = [ [ MMVector3 alloc ] initX:[ [targetArc endPosition ] X ]
											Y:[ [targetArc endPosition ] Y ]
											Z:[ [targetArc endPosition ] Z ] ] ;
											
				endTangent = [ [ MMVector3 alloc ] initUsingVector:[ targetArc endTangent ] ] ;
				
				[ endTangent reverse ] ;
			}
		
		[ end normalize ] ;
		
		
		// Need arc axis - have to be CAREFUL here, as this needs to be consistent with 
		// the orientations of the tangents of start and end arcs
		
		// Create tangent initially by taking start ^ end
		
		MMVector3 *tryAxis ;
		
		tryAxis = [ [ MMVector3 alloc ] initByCrossing:start and:end ] ;
		
		[ tryAxis normalize ] ;
				
		MMVector3 *tryStartTangent, *tryEndTangent ;
		
		tryStartTangent = [ [ MMVector3 alloc ] initByCrossing:tryAxis and:start ] ;
		
		tryEndTangent = [ [ MMVector3 alloc ] initByCrossing:tryAxis and:end ] ;
		
		// Check orientation at ends
		
		MMVector3 *testCross ;
		double dotS, dotE ;
		
		testCross = [ [ MMVector3 alloc ] initByCrossing:tryStartTangent and:startTangent ] ;
		
		dotS = [ testCross dotWith:start ] ;
		
		[ testCross release ] ;
		
		testCross = [ [ MMVector3 alloc ] initByCrossing:tryEndTangent and:endTangent ] ;
		
		dotE = [ testCross dotWith:end ] ;
		
		// If reentrant region, circulation is reversed - need to reverse these dot products
		
		if( reent == YES )
			{
				dotS = -dotS ;
				dotE = -dotE ;
			}
		
		if( dotS > 0. && dotE < 0. )
			{
				// All is right with the world!
						
				returnArc = [ [ SMArc alloc ] initWithHostCenter:center radius:[ sourceArc hostRadius ] torusSection:nil arcType:6 ] ;
				
				[ returnArc initializeWithArcCenter:center arcRadius:[ sourceArc hostRadius ] axis:tryAxis start:start end:end hostProbe:nil  ] ;
				
				[ returnArc setStartVertexPosition:startVertexPosition endVertexPosition:endVertexPosition ] ;
			}
		else if( dotS < 0. && dotE > 0. )
			{
				// Need to reverse axis
				
				[ tryAxis reverse ] ;
				
				returnArc = [ [ SMArc alloc ] initWithHostCenter:center radius:[ sourceArc hostRadius ] torusSection:nil arcType:6 ] ;
				
				[ returnArc initializeWithArcCenter:center arcRadius:[ sourceArc hostRadius ] axis:tryAxis start:start end:end hostProbe:nil  ] ;
				
				[ returnArc setStartVertexPosition:startVertexPosition endVertexPosition:endVertexPosition ] ;
			}
		else
			{
				// The world is a messed up place indeed
				
				//printf( "BUILD OF GEODESIC ARC FAILED!\n" ) ;
				
				return nil ;
			}
			
		[ tryAxis release ] ;
		[ testCross release ] ;
		[ start release ] ;
		[ end release ] ;
		[ startVertexPosition release ] ;
		[ endVertexPosition release ] ;
		[ startTangent release ] ;
		[ endTangent release ] ;
		
		return returnArc ;
	}
	
	
- (SMArc *) buildSaddleUsingStart:(SMArc *)sourceArc 
					startForward:(BOOL)sf 
					end:(SMArc *)targetArc
					endForward:(BOOL)ef
	{
	
		double theta1, phi1, theta2, phi2 ;
		
		if( sf == YES )
			{
				theta1 = [ sourceArc thetaStart ] ;
				phi1 = [ sourceArc phiStart ] ;
			}
		else
			{	
				theta1 = [ sourceArc thetaEnd ] ;
				phi1 = [ sourceArc phiEnd ] ;
			}
			
		if( ef == YES )
			{
				theta2 = [ targetArc thetaStart ] ;
				phi2 = [ targetArc phiStart ] ;
			}
		else
			{	
				theta2 = [ targetArc thetaEnd ] ;
				phi2 = [ targetArc phiEnd ] ;
			}
				
		
		SMArc *returnArc = [ [ SMArc alloc ] initWithTorusSection:[ sourceArc torusSection ]
			molecule:self phiStart:phi1 thetaStart:theta1 phiEnd:phi2 thetaEnd:theta2 ] ;
			
		return returnArc ;
	}
				
				
				
				
- (SMArc *) buildCircularArcUsingSourcePre:(SMArc *)sourcePreArc sourcePreForward:(BOOL)spf 
				sourceNext:(SMArc *)sourceNextArc sourceNextForward:(BOOL)snf
				targetPre:(SMArc *)targetPreArc targetPreForward:(BOOL)tpf
				targetNext:(SMArc *)targetNextArc targetNextForward:(BOOL)tnf
	{
		// This builds a general circular arc, using the tangent window strategy of SMART
		
		// Collect vectors
		
		MMVector3 *tangent11, *tangent12, *tangent21, *tangent22 ;
		MMVector3 *vertex1, *vertex2, *hostCenter, *normal1 ;
		
		double hostRadius ;
		
		hostCenter = [ sourcePreArc hostCenter ] ;
		hostRadius = [ sourcePreArc hostRadius ] ;
		
		if( spf == YES )
			{
				tangent11 = [ [ MMVector3 alloc ] initUsingVector:[ sourcePreArc endTangent ] ] ;
				[ tangent11 reverse ] ;
				
				vertex1 = [ sourcePreArc endPosition ] ;
			}
		else
			{
				tangent11 = [ [ MMVector3 alloc ] initUsingVector:[ sourcePreArc startTangent ] ] ;
				
				vertex1 = [ sourcePreArc startPosition ] ;
			}
			
		if( snf == YES )
			{
				tangent12 = [ [ MMVector3 alloc ] initUsingVector:[ sourceNextArc startTangent ] ] ;
			}
		else
			{
				tangent12 = [ [ MMVector3 alloc ] initUsingVector:[ sourceNextArc endTangent ] ] ;
				[ tangent12 reverse ] ;
			}
		
		if( tpf == YES )
			{
				tangent21 = [ [ MMVector3 alloc ] initUsingVector:[ targetPreArc endTangent ] ] ;
				[ tangent21 reverse ] ;
				
				vertex2 = [ targetPreArc endPosition ] ;
			}
		else
			{
				tangent21 = [ [ MMVector3 alloc ] initUsingVector:[ targetPreArc startTangent ] ] ;
				
				vertex2 = [ targetPreArc startPosition ] ;
			}
			
		if( tnf == YES )
			{
				tangent22 = [ [ MMVector3 alloc ] initUsingVector:[ targetNextArc startTangent ] ] ;
			}
		else
			{
				tangent22 = [ [ MMVector3 alloc ] initUsingVector:[ targetNextArc endTangent ] ] ;
				[ tangent22 reverse ] ;
			}
			
		
		normal1 = [ [ MMVector3 alloc ] initX:([ vertex1 X ] - [ hostCenter X ])
						Y:([ vertex1 Y ] - [ hostCenter Y ])
						Z:([ vertex1 Z ] - [ hostCenter Z ]) ] ;
						
		[ normal1 normalize ] ;
				
				
		SMArc *arc21To1, *arc22To1 ;
		MMVector3 *tangent21At1, *tangent22At1 ;
		
		arc21To1 = [ self buildCircularArcUsingHostCenter:hostCenter hostRadius:hostRadius start:vertex2 end:vertex1 startTangent:tangent21 ] ;
		arc22To1 = [ self buildCircularArcUsingHostCenter:hostCenter hostRadius:hostRadius start:vertex2 end:vertex1 startTangent:tangent22 ] ;
		
		tangent21At1 = [ arc21To1 endTangent ] ;
		tangent22At1 = [ arc22To1 endTangent ] ;
		
		// Note that we are messing up the original vectors, not copies!
		
		[ tangent21At1 reverse ] ;
		[ tangent22At1 reverse ] ;

		// Sanity check - the crossproduct tangent22At1 ^ tangent21At1 should point along normal
		
		MMVector3 *testCross ;
		
		testCross = [ [ MMVector3 alloc ] initByCrossing:tangent22At1 and:tangent21At1 ] ;
		
		if( [ testCross dotWith:normal1 ] < 0. )
			{
				return nil ;
			}
			
		// Determine limits of tangent window at vertex 1
		
		BOOL tangent22At1Inside, tangent21At1Inside ;
		
		tangent21At1Inside = [ MMVector3 isVector:tangent21At1 betweenVector:tangent11 andVector:tangent12 usingNormal:normal1 ] ;
		tangent22At1Inside = [ MMVector3 isVector:tangent22At1 betweenVector:tangent11 andVector:tangent12 usingNormal:normal1 ] ;
		
		MMVector3 *windowTangent1, *windowTangent2, *aveTangent ;
		
		// Four possibilities
		
		if( tangent21At1Inside == NO && tangent22At1Inside == NO )
			{
				// Sanity check - insure that tangent11 is inside the pair
				
				if( [ MMVector3 isVector:tangent11 betweenVector:tangent22At1 andVector:tangent21At1 usingNormal:normal1 ] == NO )
					{
						return nil ;
					}
					
				windowTangent1 = tangent11 ;
				windowTangent2 = tangent12 ;
			}
		else if( tangent21At1Inside == NO && tangent22At1Inside == YES )
			{
				// Sanity check - tangent12 should be between tangent21At1 and tangent22At1
				
				if( [ MMVector3 isVector:tangent12 betweenVector:tangent22At1 andVector:tangent21At1 usingNormal:normal1 ] == NO )
					{
						return nil ;
					}
		
				windowTangent1 = tangent22At1 ;
				windowTangent2 = tangent12 ;
			}
		else if( tangent21At1Inside == YES && tangent22At1Inside == NO )
			{
				// Sanity check - tangent11 should be between tangent21At1 and tangent22At1
				
				if( [ MMVector3 isVector:tangent11 betweenVector:tangent22At1 andVector:tangent21At1 usingNormal:normal1 ] == NO )
					{
						return nil ;
					}
		
				windowTangent1 = tangent11 ;
				windowTangent2 = tangent21At1 ;
			}
		else
			{
				// tangent21At1Inside == YES && tangent22At1Inside == YES
				
				windowTangent1 = tangent22At1 ;
				windowTangent2 = tangent21At1 ;
			}
			
		// Find average
		
		aveTangent = [ [ MMVector3 alloc ] initX:([ windowTangent1 X ] + [ windowTangent2 X ])
							Y:([ windowTangent1 Y ] + [ windowTangent2 Y ])
							Z:([ windowTangent1 Z ] + [ windowTangent2 Z ]) ] ;
							
		if( [ aveTangent length ] < 1e-04 )
			{
				// Tangents oppose - find average by cross product
				
				[ aveTangent release ] ;
				
				aveTangent = [ [ MMVector3 alloc ] initByCrossing:normal1 and:windowTangent1 ] ;
				
			}
				
				
							
		[ aveTangent normalize ] ;
		
		
		SMArc *returnArc ;
		
		returnArc = [ self buildCircularArcUsingHostCenter:hostCenter hostRadius:hostRadius start:vertex1 end:vertex2 startTangent:aveTangent ] ;
		
		[ aveTangent release ] ;
		[ tangent11 release ] ;
		[ tangent12 release ] ;
		[ arc21To1 release ] ;
		[ arc22To1 release ] ;
		[ normal1 release ] ;
		
		return returnArc ;
	}
			
		
- (SMArc *) buildCircularArcUsingHostCenter:(MMVector3 *)hc hostRadius:(double)hr start:(MMVector3 *)v1 
					end:(MMVector3 *)v2 startTangent:(MMVector3 *)t
	{
		MMVector3 *span, *axis, *arcCenter, *diff, *start, *end ;
		
		span = [ [ MMVector3 alloc ] initX:([ v2 X ] - [ v1 X ]) 
			Y:([ v2 Y ] - [ v1 Y ])
			Z:([ v2 Z ] - [ v1 Z ]) ] ;
			
		axis = [ [ MMVector3 alloc ] initByCrossing:t and:span ] ;
		
		if( [ axis length ] < 1e-06 )
			{
				printf( "ERROR, AXIS IN CIRCULAR ARC TOO SMALL!\n" ) ;
			}
		
		[ axis normalize ] ;
		
		// Get arc center 
		
		diff = [ [ MMVector3 alloc ] initX:([ v1 X ] - [ hc X ]) 
			Y:([ v1 Y ] - [ hc Y ])
			Z:([ v1 Z ] - [ hc Z ]) ] ;
		
		double lambda ;
		
		lambda = [ diff dotWith:axis ] ;
		
		arcCenter = [ [ MMVector3 alloc ] initX:( [ hc X ] + lambda*[ axis X ] )
						Y:( [ hc Y ] + lambda*[ axis Y ] )
						Z:( [ hc Z ] + lambda*[ axis Z ] ) ] ;
						
		start = [ [ MMVector3 alloc ] initX:([ v1 X ] - [ arcCenter X ])
					Y:([ v1 Y ] - [ arcCenter Y ])
					Z:([ v1 Z ] - [ arcCenter Z ]) ] ;
					
		double arcRadius ;
		
		arcRadius = [ start length ] ;
		
		[ start normalize ] ;
		
		
		end = [ [ MMVector3 alloc ] initX:([ v2 X ] - [ arcCenter X ])
					Y:([ v2 Y ] - [ arcCenter Y ])
					Z:([ v2 Z ] - [ arcCenter Z ]) ] ;
					
		[ end normalize ] ;
			
		SMArc *returnArc ;
		
		returnArc = [ [ SMArc alloc ] initWithHostCenter:hc radius:hr torusSection:nil arcType:7 ] ;
		
		[ returnArc initializeWithArcCenter:arcCenter arcRadius:arcRadius axis:axis start:start end:end hostProbe:nil ] ;
		
		return returnArc ;
	}
			
		
	
- (void) subdivideArc:(SMArc *)a
	{
		
		// This method divides an arc into two equal parts, while maintaining reference to its parent cycles
		// This method takes into account that the argument arc might be twinned, and calls a helper method
		
		NSArray *newArcs, *newTwinArcs ;
		SMArc *theTwin ;
		
		if( ( theTwin = [ a twin ] ) )
			{
				newArcs = [ self subdivideTheArc:a ] ;
				newTwinArcs = [ self subdivideTheArc:theTwin ] ;
				
				// These are expected to have opposite (geometrical) orientations
				
				int i, M ; 
				
				// Sanity check
				
				if( ( M = [ newArcs count ] ) != [ newTwinArcs count ] )
					{
						printf( "TWINNED ARCS SUBDIVIDE UNEQUALLY - Exit!\n" ) ;
						exit(1) ;
					}
				
				for( i = 0 ; i < M ; ++i )
					{
						[ [ newArcs objectAtIndex:i ] setTwin:[ newTwinArcs objectAtIndex:(M - 1 - i) ] ] ;
						[ [ newTwinArcs objectAtIndex:(M - 1 - i) ] setTwin:[ newArcs objectAtIndex:i ] ] ;
					}
					
				[ newArcs release ] ;
				[ newTwinArcs release ] ;
			}
		else
			{
				newArcs = [ self subdivideTheArc:a ] ;
				[ newArcs release ] ;
			}
			
		return ;
	}
	
- (NSArray *) subdivideTheArc:(SMArc *)a
	{
		
		SMArc *arc1, *arc2 ;
		
		if( [ a arcType ] != 5 )
			{
		

				// Find vector to divide arc
				
				double divideAngle ;
				
				divideAngle = ((double)[ a angle ]) / 2. ;
				
				double dx, dy, dz ;
				
				dx = cos(divideAngle)*[ [ a startU ] X ] + sin(divideAngle)*[ [ a uPerp ] X ] ;
				dy = cos(divideAngle)*[ [ a startU ] Y ] + sin(divideAngle)*[ [ a uPerp ] Y ] ;
				dz = cos(divideAngle)*[ [ a startU ] Z ] + sin(divideAngle)*[ [ a uPerp ] Z ] ;
				
				MMVector3 *divideVector ;
				
				divideVector = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
				
				[ divideVector normalize ] ;
				
				arc1 = [ a copyArc ] ;
				arc2 = [ a copyArc ] ;
				
				[ arc1 setEndVertex:nil ] ;
				[ arc1 setEndU:divideVector ] ;
				
				[ arc2 setStartVertex:nil ] ;
				[ arc2 setStartU:divideVector ] ;
				
				// Adjust angles 
				
				[ arc1 restoreVectorsAndAngle ] ;
				[ arc2 restoreVectorsAndAngle ] ;
				
				[ divideVector release ] ;
				
				// Careful! Even if not an interior saddle arc, we may need to adjust theta, phi angles
				
				if( [ a torusSection ] )
					{
						double theta1, phi1, theta2, phi2 ;
						
						theta1 = [ a thetaStart ] ;
						phi1 = [ a phiStart ] ;
						theta2 = [ a thetaEnd ] ;
						phi2 = [ a phiEnd ] ;
						
						[ arc1 setPhiEnd:(phi2 + phi1)/2. ] ;
						[ arc1 setThetaEnd:(theta2 + theta1)/2. ] ;
						
						[ arc2 setPhiStart:(phi2 + phi1)/2. ] ;
						[ arc2 setThetaStart:(theta2 + theta1)/2. ] ;
					}
			}
		else
			{	
				double deltaTheta, deltaPhi ;
				
				deltaTheta = ([ a thetaEnd ] - [ a thetaStart ])/2  ;
				deltaPhi =   ([ a phiEnd ] - [ a phiStart ])/2 ;
				
				int iDiv ;
				
				for( iDiv = 0 ; iDiv < 2 ; ++iDiv )
					{
						double phi1, theta1, phi2, theta2 ;
						SMArc *newSaddleArc ;
						
						phi1 = [ a phiStart ] + iDiv*deltaPhi ;
						theta1 = [ a thetaStart ] + iDiv*deltaTheta ;
						
						phi2 = phi1 + deltaPhi ;
						theta2 = theta1 + deltaTheta ;
						
						newSaddleArc = [ [ SMArc alloc ] initWithTorusSection:[ a torusSection ]
							molecule:self phiStart:phi1 thetaStart:theta1 phiEnd:phi2 thetaEnd:theta2 ] ;
							
						if( iDiv == 0 )
							{
								arc1 = newSaddleArc ;
							}
						else
							{
								arc2 = newSaddleArc ;
							}
					}
			}
		
		// Adjust cycles. I will do this the obvious way, by removing and inserting arcs
		
		int iCycle ;
		
		for( iCycle = 0 ; iCycle < [ [ a parentCycles ] count ] ; ++iCycle )
			{
				SMCycle *theParent ;
				
				theParent = [ [ a parentCycles ] objectAtIndex:iCycle ] ;
				
				int aIndex = [ [ theParent arcs ] indexOfObject:a ] ;
				
				BOOL forward ;
				
				forward = [ [ [ theParent forward ] objectAtIndex:aIndex ] boolValue ] ;
				
				if( forward == YES )
					{
						[ [ theParent arcs ] removeObjectAtIndex:aIndex ] ;
						[ [ theParent forward ] removeObjectAtIndex:aIndex ] ;
						
						[ [ theParent arcs ] insertObject:arc1 atIndex:aIndex ] ;
						[ [ theParent arcs ] insertObject:arc2 atIndex:(aIndex + 1) ] ;
						
						[ [ theParent forward ] insertObject:[ NSNumber numberWithBool:YES ] atIndex:aIndex ] ;
						[ [ theParent forward ] insertObject:[ NSNumber numberWithBool:YES ] atIndex:(aIndex + 1) ] ;
						
					}
				else
					{
						
						[ [ theParent arcs ] removeObjectAtIndex:aIndex ] ;
						[ [ theParent forward ] removeObjectAtIndex:aIndex ] ;
						
						[ [ theParent arcs ] insertObject:arc2 atIndex:aIndex ] ;
						[ [ theParent arcs ] insertObject:arc1 atIndex:(aIndex + 1) ] ;
						
						[ [ theParent forward ] insertObject:[ NSNumber numberWithBool:NO ] atIndex:aIndex ] ;
						[ [ theParent forward ] insertObject:[ NSNumber numberWithBool:NO ] atIndex:(aIndex + 1) ] ;
					}
					
			}
			
		// Check if we have torus to update
		
		SMTorus *theTorus ;
		
		if( ( theTorus = [ a torusSection ] ) != nil  && [ a arcType ] != 5 )
			{
				NSArray *arcs ;
				
				arcs = [ NSArray arrayWithObjects:arc1,arc2,nil ] ;
				
				[ theTorus registerArcs:arcs parentArc:a ] ;
			}
			
		NSArray *returnArcs = [ [ NSArray alloc ] initWithObjects:arc1,arc2,nil ] ;
		
		return returnArcs ;
	}
		
- (void) subdivideArc:(SMArc *)a usingDivision:(double)div
	{
		
		// This method divides an arc into multiple parts, while maintaining reference to its parent cycles
		// This method takes into account that the argument arc might be twinned, and calls a helper method
		
		// It also provides special handling for saddle arcs
		
		NSArray *newArcs, *newTwinArcs ;
		SMArc *theTwin ;
		
		if( ( theTwin = [ a twin ] ) )
			{
				newArcs = [ self subdivideTheArc:a usingDivision:div ] ;
				newTwinArcs = [ self subdivideTheArc:theTwin usingDivision:div ] ;
				
				if( ! newArcs ) return ;
				
				// These are expected to have opposite (geometrical) orientations
				
				int i, M ; 
				
				// Sanity check
				
				if( ( M = [ newArcs count ] ) != [ newTwinArcs count ] )
					{
						printf( "TWINNED ARCS SUBDIVIDE UNEQUALLY - Exit!\n" ) ;
						exit(1) ;
					}
				
				for( i = 0 ; i < M ; ++i )
					{
						[ [ newArcs objectAtIndex:i ] setTwin:[ newTwinArcs objectAtIndex:(M - 1 - i) ] ] ;
						[ [ newTwinArcs objectAtIndex:(M - 1 - i) ] setTwin:[ newArcs objectAtIndex:i ] ] ;
					}
					
				[ newArcs release ] ;
				[ newTwinArcs release ] ;
			}
		else
			{
				newArcs = [ self subdivideTheArc:a usingDivision:div ] ;
				[ newArcs release ] ;
			}
			
		return ;
	}
		
		
- (NSArray *) subdivideTheArc:(SMArc *)a usingDivision:(double) div
	{
		NSMutableArray *newArcs ;
		
		double deltaAngle, endAngle ;
		
		int nDiv, iDiv ;
		
		MMVector3 *endVector, *startU, *endU ;
		
		SMVertex *startVertex, *endVertex ;
		
		double phi1, theta1, phi2, theta2, deltaTheta, deltaPhi ;
		
		nDiv = (int) floor( [ a length ] / div ) ;
		
		if( nDiv < 2 ) return nil ;
		
		newArcs = [ [ NSMutableArray alloc ] initWithCapacity:10 ] ;
		
		// Note - use saddle-type subdivision on reentrant arcs (reentrantR and reentrantL) to handle self-intersecting surface
		
		if( [ a arcType ] != 5  )
			{
		
				deltaAngle = [ a angle ] / nDiv ;
				
				if( [ a torusSection ] )
					{
						deltaPhi = ([ a phiEnd ] - [ a phiStart ])/nDiv ;
						deltaTheta = ([ a thetaEnd ] - [ a thetaStart ])/nDiv ;
					}
				
				
				for( iDiv = 0 ; iDiv < nDiv ; ++iDiv )
					{
						endAngle = ( iDiv + 1 ) * deltaAngle ;
						
						double dx, dy, dz ;
						
						dx = cos(endAngle)*[ [ a startU ] X ] + sin(endAngle)*[ [ a uPerp ] X ] ;
						dy = cos(endAngle)*[ [ a startU ] Y ] + sin(endAngle)*[ [ a uPerp ] Y ] ;
						dz = cos(endAngle)*[ [ a startU ] Z ] + sin(endAngle)*[ [ a uPerp ] Z ] ;
						
						if( iDiv < (nDiv - 1) )
							{
								endVector = [ [ MMVector3 alloc ] initX:dx Y:dy Z:dz ] ;
						
								[ endVector normalize ] ;
							}
							
						if( iDiv == 0 )
							{
								startU = [ a startU ] ;
								endU = endVector ;
								
								startVertex = [ a startVertex ] ;
								endVertex = nil ;
							}
						else if( iDiv < nDiv - 1 )
							{
								startU = [ [ newArcs lastObject ] endU ] ;
								endU = endVector ;
								
								startVertex = nil ;
								endVertex = nil ;
							}
						else	
							{
								startU = [ [ newArcs lastObject ] endU ] ;
								endU = [ a endU ] ;
								
								startVertex = nil ;
								endVertex = [ a endVertex ] ;
							}
							
						SMArc *copyArc ;
						
						copyArc = [ a copyArc ] ;
						
						[ copyArc setStartU:startU ] ;
						[ copyArc setEndU:endU ] ;
						
						[ copyArc setStartVertex:startVertex ] ;
						[ copyArc setEndVertex:endVertex ] ;
						
						[ copyArc restoreVectorsAndAngle ] ;
						
						if( [ a torusSection ] )
							{
								phi1 = [ a phiStart ] + iDiv*deltaPhi ;
								theta1 = [ a thetaStart ] + iDiv*deltaTheta ;
								
								phi2 = phi1 + deltaPhi ;
								theta2 = theta1 + deltaTheta ;
								
								[ copyArc setPhiStart:phi1 ] ;
								[ copyArc setThetaStart:theta1 ] ;
								
								[ copyArc setPhiEnd:phi2 ] ;
								[ copyArc setThetaEnd:theta2 ] ;
								
							}
						
					
						[ newArcs addObject:copyArc ] ;
						
						if( iDiv < (nDiv - 1) )
							{
								[ endVector release ] ;
							}
							
						
					}
			}
		else
			{	
				deltaTheta = ([ a thetaEnd ] - [ a thetaStart ])/nDiv  ;
				deltaPhi =   ([ a phiEnd ] - [ a phiStart ])/nDiv ;
				
						for( iDiv = 0 ; iDiv < nDiv ; ++iDiv )
							{
								SMArc *newSaddleArc ;
								
								phi1 = [ a phiStart ] + iDiv*deltaPhi ;
								theta1 = [ a thetaStart ] + iDiv*deltaTheta ;
								
								phi2 = phi1 + deltaPhi ;
								theta2 = theta1 + deltaTheta ;
								
								newSaddleArc = [ [ SMArc alloc ] initWithTorusSection:[ a torusSection ]
									molecule:self phiStart:phi1 thetaStart:theta1 phiEnd:phi2 thetaEnd:theta2 ] ;
									
								[ newArcs addObject:newSaddleArc ] ;
							}
						
			}
						
						
				
		// Adjust cycles. I will do this the obvious way, by removing and inserting arcs
		
		int iCycle ;
		
		for( iCycle = 0 ; iCycle < [ [ a parentCycles ] count ] ; ++iCycle )
			{
				SMCycle *theParent ;
				
				theParent = [ [ a parentCycles ] objectAtIndex:iCycle ] ;
				
				int aIndex = [ [ theParent arcs ] indexOfObject:a ] ;
				
				BOOL forward ; 
				
				forward = [ [ [ theParent forward ] objectAtIndex:aIndex ] boolValue ] ;
				
				if( forward == YES )
					{
						[ [ theParent arcs ] removeObjectAtIndex:aIndex ] ;
						[ [ theParent forward ] removeObjectAtIndex:aIndex ] ;
						
						for( iDiv = 0 ; iDiv < nDiv ; ++iDiv )
							{
								[ [ theParent arcs ] insertObject:[ newArcs objectAtIndex:iDiv ] atIndex:(aIndex + iDiv) ] ;
								[ [ theParent forward ] insertObject:[ NSNumber numberWithBool:YES ] atIndex:(aIndex + iDiv) ] ;
							}
					}
				else
					{
						
						[ [ theParent arcs ] removeObjectAtIndex:aIndex ] ;
						[ [ theParent forward ] removeObjectAtIndex:aIndex ] ;
						
						for( iDiv = 0 ; iDiv < nDiv ; ++iDiv )
							{
								[ [ theParent arcs ] insertObject:[ newArcs objectAtIndex:(nDiv - iDiv - 1) ] atIndex:(aIndex + iDiv) ] ;
								[ [ theParent forward ] insertObject:[ NSNumber numberWithBool:NO ] atIndex:(aIndex + iDiv) ] ;
							}
						
					}
					
			}
			
		// Check if we have torus to update
		
		SMTorus *theTorus ;
		
		if( ( theTorus = [ a torusSection ] ) != nil && [ a arcType ] != 5  )
			{				
				[ theTorus registerArcs:newArcs parentArc:a ] ;
			}
			
		//[ newArcs release ] ;
			
		return newArcs ;
	}
	
- (NSArray *)  subdivideCycle:(SMCycle *)cyc firstIndex:(int)start
					secondIndex:(int)second usingArc:(SMArc *)arc
	{
		// Subdivide a cycle into two - return an array of the new cycles
		
		SMCycle *cycle1, *cycle2 ;
		
		cycle1 = [ [ SMCycle alloc ] init ] ;
		cycle2 = [ [ SMCycle alloc ] init ] ;
		
		[ cycle1 addArc:arc forward:YES ] ;
		
		//[ [ cycle1 arcs ] addObject:arc ] ;
		//[ [ cycle1 forward ] addObject:[ NSNumber numberWithBool:YES ] ] ;
		
		int arcIndex ;
	
		int M = [ [ cyc arcs ] count ] ;
		
		arcIndex = second ;
		
		while( TRUE )
			{
				//[ [ cycle1 arcs ] addObject:[ [ cyc arcs ] objectAtIndex:arcIndex ] ] ;
				//[ [ cycle1 forward ] addObject:[ [ cyc forward ] objectAtIndex:arcIndex ] ] ;
				
				SMArc *cycArc ; BOOL isForward ;
				
				cycArc = [ [ cyc arcs ] objectAtIndex:arcIndex ] ;
				isForward = [ [ [ cyc forward ] objectAtIndex:arcIndex ] boolValue ] ;
				
				[ cycArc removeParentCycle:cyc ] ;
				
				[ cycle1 addArc:cycArc forward:isForward ] ;
				
				++arcIndex ;
				
				if( arcIndex == M )
					{
						arcIndex = 0 ;
					}
					
				if( arcIndex == start ) break ;
			}
			
		// Other way around
		
		[ cycle2 addArc:arc forward:NO ] ;
		
		//[ [ cycle2 arcs ] addObject:arc ] ;
		//[ [ cycle2 forward ] addObject:[ NSNumber numberWithBool:NO ] ] ;
		
		arcIndex = start ;
		
		while( TRUE )
			{
				//[ [ cycle2 arcs ] addObject:[ [ cyc arcs ] objectAtIndex:arcIndex ] ] ;
				//[ [ cycle2 forward ] addObject:[ [ cyc forward ] objectAtIndex:arcIndex ] ] ;
				
				SMArc *cycArc ; BOOL isForward ;
				
				cycArc = [ [ cyc arcs ] objectAtIndex:arcIndex ] ;
				isForward = [ [ [ cyc forward ] objectAtIndex:arcIndex ] boolValue ] ;
				
				[ cycArc removeParentCycle:cyc ] ;
				
				[ cycle2 addArc:cycArc forward:isForward ] ;
				
				
				++arcIndex ;
									
				if( arcIndex == second ) break ;
			}
		
		// DEBUGING - check if any of the arcs match
		/*
		int i, j ;
		
		for( i = 0 ; i < [ [ cycle1 arcs ] count ] - 1 ; ++i )
			{
				for( j = i + 1 ; j < [ [ cycle1 arcs ] count ] ; ++j )
					{
						if( [ [ cycle1 arcs ] objectAtIndex:i ] == [ [ cycle1 arcs ] objectAtIndex:j ] )
							{
								printf( "WARNING - Cycle has duplicate arcs\n" ) ;
							}
					}
					
			}
			
		for( i = 0 ; i < [ [ cycle2 arcs ] count ] - 1 ; ++i )
			{
				for( j = i + 1 ; j < [ [ cycle2 arcs ] count ] ; ++j )
					{
						if( [ [ cycle2 arcs ] objectAtIndex:i ] == [ [ cycle2 arcs ] objectAtIndex:j ] )
							{
								printf( "WARNING - Cycle has duplicate arcs\n" ) ;
							}
					}
					
			}
		*/
		
		[ cycle1 setParentCycle:cyc ] ;
		[ cycle2 setParentCycle:cyc ] ;
		
		[ cycle1 setAtoms:[ cyc atoms ] ] ;
		[ cycle2 setAtoms:[ cyc atoms ] ] ;
		
		if( [ cyc probe ] )
			{
				[ cycle1 setProbe:[ cyc probe ] ] ;
				[ cycle2 setProbe:[ cyc probe ] ] ;
			}
		
		// Limit planes (only needed for reentrant cycles)
		
		cycle1->theLimitPlane = cyc->theLimitPlane ;
		cycle2->theLimitPlane = cyc->theLimitPlane ;
		
		cycle1->selfIntersection = cyc->selfIntersection ;
		cycle2->selfIntersection = cyc->selfIntersection ;
		
		NSArray *returnArray = [ [ NSArray alloc ] initWithObjects:cycle1,cycle2,nil ] ;
		
		return returnArray ;
	}
				
		
		
- (SMCycle *) combineCycle:(SMCycle *)cycleJ firstIndex:(int)startJ withCycle:(SMCycle *)cycleK secondIndex:(int)startK
				usingArc:(SMArc *)a
	{
		// This function makes the union of two cycles
		
		// Begin at startJ of cycle J, go around to preArcJ, traverse linking arc, traverse cycleK from startK to preArcK, and traverse 
		// linking arc in reverse direction
		
		int i, arcIndex ;
		
		SMCycle *returnCycle = [ [ SMCycle alloc ] init ] ;
		
		arcIndex = startJ ;
		
		for( i = 0 ; i < [ [ cycleJ arcs ] count ] ; ++i )
			{
				[ returnCycle addArc:[ [ cycleJ arcs ] objectAtIndex:arcIndex ] forward:[ [ [ cycleJ forward ] objectAtIndex:arcIndex ] boolValue ] ] ;
				
				//[ [ returnCycle arcs ] addObject:[ [ cycleJ arcs ] objectAtIndex:arcIndex ] ] ;
				//[ [ returnCycle forward ] addObject:[ [ cycleJ forward ] objectAtIndex:arcIndex ] ] ;
				
				[ [ [ cycleJ arcs ] objectAtIndex:arcIndex ] removeParentCycle:cycleJ ] ;
				
				++arcIndex ;
				
				if( arcIndex == [ [ cycleJ arcs ] count ] )
					{
						arcIndex = 0 ;
					}
			}
			
		// Traverse new arc - first, make sure of orientation
		// (That may not be needed, I think it is 'forward' by construction)
		
		MMVector3 *testVector ;
		SMArc *source ;
		
		source = [ [ cycleJ arcs ] objectAtIndex:startJ ] ;
		BOOL aForward ;
		
		testVector = [ [ MMVector3 alloc ] 
							initX:( [ [ a startPosition ] X ] - [ [ source startPosition ] X ] )
							Y:( [ [ a startPosition ] Y ] - [ [ source startPosition ] Y ] )
							Z:( [ [ a startPosition ] Z ] - [ [ source startPosition ] Z ] ) ] ;
							
		if( [ testVector length ] == 0. )
			{
				aForward = YES ;
			}
		else
			{
				aForward = NO ;
			}
			
		[ returnCycle addArc:a forward:aForward ] ;
		
		//[ [ returnCycle arcs ] addObject:a ] ;
		//[ [ returnCycle forward ] addObject:[ NSNumber numberWithBool:aForward ] ] ;
		
		// Traverse other cycle 
		
		arcIndex = startK ;
		
		for( i = 0 ; i < [ [ cycleK arcs ] count ] ; ++i )
			{
				[ returnCycle addArc:[ [ cycleK arcs ] objectAtIndex:arcIndex ] forward:[ [ [ cycleK forward ] objectAtIndex:arcIndex ] boolValue ] ] ;
				
				//[ [ returnCycle arcs ] addObject:[ [ cycleK arcs ] objectAtIndex:arcIndex ] ] ;
				//[ [ returnCycle forward ] addObject:[ [ cycleK forward ] objectAtIndex:arcIndex ] ] ;
				
				++arcIndex ;
				
				if( arcIndex == [ [ cycleK arcs ] count ] )
					{
						arcIndex = 0 ;
					}
			}
		
		// Now, COPY of new arc again, traversed opposite way
		
		SMArc *aCopy ;
		
		aCopy = [ a copyArc ] ;
		
		// Reverse the copy so we can twin it!
		
		[ aCopy reverse ] ;
		
		// Add with same orientation as original arc 
		/*
		if( aForward == YES )
			{
				aForward = NO ;
			}
		else
			{
				aForward = YES ;
			}
		*/
			
		[ returnCycle addArc:aCopy forward:aForward ] ;
		
		// Have two copies of parent cycle - NO WE DON'T. I updated the addParentCycle method to preclude this
		
		// [ [ aCopy parentCycles ] removeLastObject ] ;
		
		
		
		[ aCopy setTwin:a ] ;
		[ a setTwin:aCopy ] ;
		
		//[ [ returnCycle arcs ] addObject:a ] ;
		//[ [ returnCycle forward ] addObject:[ NSNumber numberWithBool:aForward ] ] ;
		
		// Add atoms from cycleJ - should be the same as cycleK
		
		[ returnCycle setAtoms:[ cycleJ atoms ] ] ;
		
		// I think that's it!
		
		return returnCycle ;
	}
			
		
		
								
							
						
- (void) assignOrientationForProbe:(SMProbe *)p usingAtomI:(int)i J:(int)j K:(int)k
	{
		// Is atom K above or below the I-P-J plane? (Right-hand rule) 
		// A right probe is below the plane, left is above. If in the plane, no assignments are made (that should not happen anyway)
		
		double dxIP, dyIP, dzIP, dxIJ, dyIJ, dzIJ, dxIK, dyIK, dzIK ;
		
		dxIP = [ p X ] - xAtom[i] ;
		dyIP = [ p Y ] - yAtom[i] ;
		dzIP = [ p Z ] - zAtom[i] ;
		
		dxIJ = xAtom[j] - xAtom[i] ;
		dyIJ = yAtom[j] - yAtom[i] ;
		dzIJ = zAtom[j] - zAtom[i] ;
		
		dxIK = xAtom[k] - xAtom[i] ;
		dyIK = yAtom[k] - yAtom[i] ;
		dzIK = zAtom[k] - zAtom[i] ;
		
		
		// Cross product, IP X IJ
		
		double cx, cy, cz ;
		
		cx = dyIP * dzIJ - dyIJ * dzIP ;
		cy = dzIP * dxIJ - dzIJ * dxIP ;
		cz = dxIP * dyIJ - dxIJ * dyIP ;
		
		double dot ;
		BOOL left, right ;
		
		dot = cx * dxIK + cy * dyIK + cz * dzIK ;
		
		if( dot > 0. )
			{
				right = NO ;
				left = YES ;
			}
		else if( dot < 0. )
			{
				right = YES ;
				left = NO ;
			}
		else 
			{
				right = NO ;
				left = NO ;
			}
			
		
		[ p setProbeTypeRight:right left:left ] ;
		
		return ;
	}
		
		
								
- (void) exportAsFlatsUsingPath:(NSString *)p oldStyle:(BOOL)s useSubsurfaces:(BOOL)sub
	{
		// Write flats-style output file
		
		// Two formats are available - "new style" and "old"
		
		// New style
		
		//	#vertices (%d) N(ew)
		// x y z xn yn zn (%f10.6 %f10.6 %f10.6 %f10.6 %f10.6 %f10.6)
		// (etc)
		// #elements #contact #reentrant #saddle
		// v1 v2 v3 enx eny enz (%d %d %d %f10.6 %f10.6 %f10.6 ) [0|1|2] #atoms (%d) atom1 atom2 ... (%d %d ...) [0|1 / 1=self-intersecting surface ]
		// (etc - NOTE element type code: 0=contact, 1=reentrant, 2=saddle )
		
		// Old style:
		
		//	#vertices (%d) O(ld)
		// x y z xn yn zn (%f10.6 %f10.6 %f10.6 %f10.6 %f10.6 %f10.6)
		// (etc)
		// #elements #contact #reentrant #saddle
		// v1 v2 v3 enx eny enz (%d %d %d %f10.6 %f10.6 %f10.6 ) [0|1|2] atom1 atom2 atom3 (%d %d %d)
		// (etc - NOTE element type code: 0=contact, 1=reentrant, 2=saddle; note that atom2 and atom3 may = -1 )
		
		FILE *output ;
			if( !( output = fopen( [ p cString ], "w" ) ) )
			{
				printf( "ERROR: COULD NOT OPEN OUTPUT FILE %s - Exit!\n", [ p cString ] ) ;
				exit(1) ;
			}
			
		if( s == YES )
			{
				// Old Style
				
				fprintf( output, "%d OLD\n", nVertices ) ;
			}
		else
			{
				// New style
					
				fprintf( output, "%d NEW\n", nVertices ) ;
			}
			
		int i, j ;
		
		for( i = 0 ; i < nVertices ; ++i )
			{
				fprintf( output, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", 
					[ [ [ vertices objectAtIndex:i ] vertexPosition ] X ], [ [ [ vertices objectAtIndex:i ] vertexPosition ] Y ], [ [ [ vertices objectAtIndex:i ] vertexPosition ] Z ],
					[ [ [ vertices objectAtIndex:i ] normal ] X ], [ [ [ vertices objectAtIndex:i ] normal ] Y ], [ [ [ vertices objectAtIndex:i ] normal ] Z ] ) ;
			}
			
		NSEnumerator *cycleEnumerator, *arcEnumerator, *forwardEnumerator ;
		NSArray *arcs, *forward ;
		SMCycle *nextCycle ;
		SMArc *nextArc ;
		NSNumber *nextForward ;
		
		MMVector3 *normal  ;
		
		MMVector3 *vs[3], *v12, *v13 ;
		int ivs[3], vCount ;
		
		double dx12, dy12, dz12, dx13, dy13, dz13 ;
		
			
		if( s == YES )
			{
				// Old style
				
				fprintf( output, "%d %d = N_ELEMENT, N_CONTACT\n", nElements, nContactElements ) ;
			}
		else
			{
				// New style
				
				fprintf( output, "%d %d %d %d = nEle, nCon, nReen, nSadd\n", nElements, nContactElements, nReentrantElements,
					nSaddleElements ) ;
			}
			
		int iPhase ;
		
		for( iPhase = 0 ; iPhase <= 2 ; ++iPhase )
			{
			
				if( iPhase == 0 )
					{
						cycleEnumerator = [ contactCycles objectEnumerator ] ;
					}
				else if( iPhase == 1 )
					{
						cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
					}
				else
					{	
						cycleEnumerator = [ saddleCycles objectEnumerator ] ;
					}
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						if( [ nextCycle active ] == NO ) continue ;
						
						arcs = [ nextCycle arcs ] ;
						forward = [ nextCycle forward ] ;
						
						vCount = 0 ;
						
						
						// Note - some arcs may be skipped!!
						
						// BEWARE! Saddle arcs are ordered counterclockwise (looking down), reentrant and contact are clockwise!
						
						arcEnumerator = [ arcs objectEnumerator ] ;
						forwardEnumerator = [ forward objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								nextForward = [ forwardEnumerator nextObject ] ;
								
								if( [ nextArc skip ] == YES ) continue ;
										
								/* Don't need this check here 
								if( vCount == 3 )
									{
										printf( "TOO MANY ARCS IN CYCLE OF TYPE %d - Exit!\n", iPhase ) ;
										exit(1) ;
									}
								*/
										
								if( [ nextForward boolValue ] == YES )
									{
										vs[vCount] = [ [ nextArc startVertex ] vertexPosition ] ;
										ivs[vCount] = [ [ nextArc startVertex] index ] ;
									}
								else
									{
										vs[vCount] = [ [ nextArc endVertex ] vertexPosition ] ;
										ivs[vCount] = [ [ nextArc endVertex ] index ] ;
									}
									
								++vCount ;
							}
							

						if( iPhase != 2  || ( iPhase == 2 && [ [ [ [ nextCycle arcs ] lastObject ] torusSection ] skipTorus ] == NO ) )
							{
								if( vCount != 3 )
									{
										printf( "VERTEX COUNT %d IN CYCLE OF TYPE %d - Exit!\n", vCount, iPhase ) ;
										exit(1) ;
									}
							}
							
				
						dx12 = [ vs[1] X ] - [ vs[0] X ] ;
						dy12 = [ vs[1] Y ] - [ vs[0] Y ] ;
						dz12 = [ vs[1] Z ] - [ vs[0] Z ] ;
						
						dx13 = [ vs[2] X ] - [ vs[0] X ] ;
						dy13 = [ vs[2] Y ] - [ vs[0] Y ] ;
						dz13 = [ vs[2] Z ] - [ vs[0] Z ] ;
						
						v12 = [ [ MMVector3 alloc ] initX:dx12 Y:dy12 Z:dz12 ] ;
						v13 = [ [ MMVector3 alloc ] initX:dx13 Y:dy13 Z:dz13 ] ;
						
						
						normal = [ [ MMVector3 alloc ] initByCrossing:v13 and:v12 ] ;
						
						[ normal normalize ] ;
						
						if( iPhase == 2 ) [ normal reverse ] ;
						
						if( iPhase == 2 )
							{
								fprintf( output, "%d %d %d %10.6f %10.6f %10.6f %d ", ivs[0], ivs[1], ivs[2], 
									[ normal X ], [ normal Y ], [ normal Z ], iPhase ) ;
							}
						else
							{
								fprintf( output, "%d %d %d %10.6f %10.6f %10.6f %d ", ivs[2], ivs[1], ivs[0], 
									[ normal X ], [ normal Y ], [ normal Z ], iPhase ) ;
							}
							
							
						[ normal release ] ; 
						[ v12 release ] ;
						[ v13 release ] ;
						
						NSArray *elemAtoms ;
						int nElemAtoms ;
						
						elemAtoms = [ [ nextCycle atoms ] allObjects ] ;
						nElemAtoms = [ elemAtoms count ] ;
						
						if( s == YES )
							{
								// Old style
								
								for( i = 0 ; i < 3 ; ++i )
									{
										if( i >= nElemAtoms ) 
											{
												fprintf( output, "-1 " ) ;
											}
										else
											{
												fprintf( output, "%d ", [ [ elemAtoms objectAtIndex:i ] intValue ] ) ;
											}
									}
							}
						else
							{
								// New Style
								
								fprintf( output, "%d ", nElemAtoms ) ;
								
								for( i = 0 ; i < nElemAtoms ; ++i )
									{
										fprintf( output, "%d ", [ [ elemAtoms objectAtIndex:i ] intValue ] ) ;
									}
									
								if( nextCycle->selfIntersection == YES )
									{
										fprintf( output, "S" ) ;
									}
								else
									{
										fprintf( output, "N" ) ;
									}
							}
								
						fprintf( output, "\n" ) ;
							
						// Next cycle
					}
					
				// Next phase
			}
			
		// All done with complete surface
		
		fclose( output ) ;
		
		if( sub == NO ) return ;
		
		// We will also write out subsurface files
		
		// First, sort subsurfaces by size, which we will do crudely as # vertices
		
		typedef struct { int subsurfaceName ; int subsurfaceSize ; int contactElems ; int reentrantElems ; int saddleElems ; FILE *outFile ; } subsurfaceData ;

		subsurfaceData *subsurfaceInfo = (subsurfaceData *) malloc( nSubsurfaces * sizeof( subsurfaceData ) ) ;
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				subsurfaceInfo[i].subsurfaceName = i ;
				subsurfaceInfo[i].contactElems = 0 ;
				subsurfaceInfo[i].reentrantElems = 0 ;
				subsurfaceInfo[i].saddleElems = 0 ;
				subsurfaceInfo[i].subsurfaceSize = 0 ;
			}
		
		for( i = 0 ; i < nVertices ; ++i )
			{
				int nextSubsurface = [ [ vertices objectAtIndex:i ] subsurface ] ;
				
				++subsurfaceInfo[nextSubsurface].subsurfaceSize ;
				
			}
	
		subsurfaceData tempSubsurfaceInfo ;
		
		// Bubble sort
		
		int *subsurfaceNameToRank = (int *) malloc( nSubsurfaces * sizeof( int ) ) ;
		
		for( i = 0 ; i < nSubsurfaces - 1 ; ++i )
			{
				for( j = i + 1 ; j < nSubsurfaces ; ++j )
					{
						if( subsurfaceInfo[j].subsurfaceSize > subsurfaceInfo[i].subsurfaceSize )
							{
								tempSubsurfaceInfo = subsurfaceInfo[i] ;
								subsurfaceInfo[i] = subsurfaceInfo[j] ;
								subsurfaceInfo[j] = tempSubsurfaceInfo ;
							}
					}
			}
			
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				subsurfaceNameToRank[ subsurfaceInfo[i].subsurfaceName ] = i ;
			}
			
			
		// Elements counts
				
		for( iPhase = 0 ; iPhase <= 2 ; ++iPhase )
			{
			
				if( iPhase == 0 )
					{
						cycleEnumerator = [ contactCycles objectEnumerator ] ;
					}
				else if( iPhase == 1 )
					{
						cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
					}
				else
					{	
						cycleEnumerator = [ saddleCycles objectEnumerator ] ;
					}
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						if( [ nextCycle active ] == NO ) continue ;
						
						int rank = subsurfaceNameToRank[  [ nextCycle subsurface ] ] ;
						
						if( iPhase == 0 )
							{
								++subsurfaceInfo[rank].contactElems ;
							}
						else if( iPhase == 1 )
							{
								++subsurfaceInfo[rank].reentrantElems ;
							}
						else
							{
								++subsurfaceInfo[rank].saddleElems ;
							}
					}
			}
			
		// We need to sort the vertices in the subsurfaces in accord with local (in-subsurface) index
		
		NSMutableArray *subsurfaceRankToVertices = [ [ NSMutableArray alloc ] initWithCapacity:nSubsurfaces ] ;
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				[ subsurfaceRankToVertices addObject:[ [ NSMutableArray alloc ] initWithCapacity:subsurfaceInfo[i].subsurfaceSize ] ] ;
			}
			
		for( i = 0 ; i < nVertices ; ++i )
			{
				int rank ;
				SMVertex *nextVertex = [ vertices objectAtIndex:i ] ;
				
				rank = subsurfaceNameToRank[ [ nextVertex subsurface ] ] ;
				
				[ [ subsurfaceRankToVertices objectAtIndex:rank ] addObject:nextVertex ] ;
			}
			
		// Sort vertices 
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				[ [ subsurfaceRankToVertices objectAtIndex:i ] sortUsingSelector:@selector(compareSubsurfaceVertices:) ] ;
			}
				
	
		// Open output files
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				FILE *nextFile ;
				
				char outName[500] ;
				sprintf( outName, "%s.%d", [ p cString ], i ) ;
				
				if( !( nextFile = fopen( outName, "w" ) ) )
					{
						printf( "COULD NOT OPEN OUTPUT FILE %s - Exit!\n", outName ) ;
						exit(1) ;
					}
					
				subsurfaceInfo[i].outFile = nextFile ;
			}
			
		
		if( s == YES )
			{
				// Old Style
				
				for( i = 0 ; i < nSubsurfaces ; ++i )
					{
						fprintf( subsurfaceInfo[i].outFile, "%d OLD\n", subsurfaceInfo[i].subsurfaceSize ) ;
					}
			}
		else
			{
				// New style
				
				for( i = 0 ; i < nSubsurfaces ; ++i )
					{
						fprintf( subsurfaceInfo[i].outFile, "%d NEW\n", subsurfaceInfo[i].subsurfaceSize ) ;
					}
			}
			
		for( j = 0 ; j < nSubsurfaces ; ++j )
			{
				NSArray *nextVertices = [ subsurfaceRankToVertices objectAtIndex:j ] ;
				
				for( i = 0 ; i < [ nextVertices count ] ; ++i )
					{
						fprintf( subsurfaceInfo[j].outFile, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", 
							[ [ [ nextVertices objectAtIndex:i ] vertexPosition ] X ], [ [ [ nextVertices objectAtIndex:i ] vertexPosition ] Y ], [ [ [ nextVertices objectAtIndex:i ] vertexPosition ] Z ],
							[ [ [ nextVertices objectAtIndex:i ] normal ] X ], [ [ [ nextVertices objectAtIndex:i ] normal ] Y ], [ [ [ nextVertices objectAtIndex:i ] normal ] Z ] ) ;
					}
			}
			
		if( s == YES )
			{
				// Old style
				
				for( i = 0 ; i < nSubsurfaces ; ++i )
					{
					 	fprintf( subsurfaceInfo[i].outFile, "%d %d = N_ELEMENT, N_CONTACT\n", 
							(subsurfaceInfo[i].contactElems + subsurfaceInfo[i].reentrantElems + subsurfaceInfo[i].saddleElems),
							subsurfaceInfo[i].contactElems ) ;
					}
			}
		else
			{
				// New style
				
				for( i = 0 ; i < nSubsurfaces ; ++i )
					{
					 	fprintf( subsurfaceInfo[i].outFile, "%d %d %d %d = nEle, nCon, nReen, nSadd\n", 
							(subsurfaceInfo[i].contactElems + subsurfaceInfo[i].reentrantElems + subsurfaceInfo[i].saddleElems),
							subsurfaceInfo[i].contactElems, subsurfaceInfo[i].reentrantElems, subsurfaceInfo[i].saddleElems ) ;
					}
				
			}
			
		
		for( iPhase = 0 ; iPhase <= 2 ; ++iPhase )
			{
			
				if( iPhase == 0 )
					{
						cycleEnumerator = [ contactCycles objectEnumerator ] ;
					}
				else if( iPhase == 1 )
					{
						cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
					}
				else
					{	
						cycleEnumerator = [ saddleCycles objectEnumerator ] ;
					}
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						FILE *outFile ;
						int rank ;
						
						if( [ nextCycle active ] == NO ) continue ;
						
						rank = subsurfaceNameToRank[ [ nextCycle subsurface ] ] ;
						
						outFile = subsurfaceInfo[rank].outFile ;
						
						arcs = [ nextCycle arcs ] ;
						forward = [ nextCycle forward ] ;
						
						vCount = 0 ;
						
						/*
						if( [ arcs count ] != 3 )
							{
								printf( "CYCLE OF TYPE %d HAS %d ARCS - Exit!\n", iPhase, [ arcs count ] ) ;
								exit(1) ;
							}
						*/
						
						// Note - some arcs may be skipped!!
						
						// BEWARE! Saddle arcs are ordered counterclockwise (looking down), reentrant and contact are clockwise!
						
						arcEnumerator = [ arcs objectEnumerator ] ;
						forwardEnumerator = [ forward objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								nextForward = [ forwardEnumerator nextObject ] ;
								
								if( [ nextArc skip ] == YES ) continue ;
										
								/* Don't need this check here 
								if( vCount == 3 )
									{
										printf( "TOO MANY ARCS IN CYCLE OF TYPE %d - Exit!\n", iPhase ) ;
										exit(1) ;
									}
								*/
										
								if( [ nextForward boolValue ] == YES )
									{
										vs[vCount] = [ [ nextArc startVertex ] vertexPosition ] ;
										ivs[vCount] = [ [ nextArc startVertex] subsurfaceIndex] ;
									}
								else
									{
										vs[vCount] = [ [ nextArc endVertex ] vertexPosition ] ;
										ivs[vCount] = [ [ nextArc endVertex ] subsurfaceIndex ] ;
									}
									
								++vCount ;
							}
							
						if( iPhase != 2  || ( iPhase == 2 && [ [ [ [ nextCycle arcs ] lastObject ] torusSection ] skipTorus ] == NO ) )
							{
								if( vCount != 3 )
									{
										printf( "VERTEX COUNT %d IN CYCLE OF TYPE %d - Exit!\n", vCount, iPhase ) ;
										exit(1) ;
									}
							}
				
						dx12 = [ vs[1] X ] - [ vs[0] X ] ;
						dy12 = [ vs[1] Y ] - [ vs[0] Y ] ;
						dz12 = [ vs[1] Z ] - [ vs[0] Z ] ;
						
						dx13 = [ vs[2] X ] - [ vs[0] X ] ;
						dy13 = [ vs[2] Y ] - [ vs[0] Y ] ;
						dz13 = [ vs[2] Z ] - [ vs[0] Z ] ;
						
						v12 = [ [ MMVector3 alloc ] initX:dx12 Y:dy12 Z:dz12 ] ;
						v13 = [ [ MMVector3 alloc ] initX:dx13 Y:dy13 Z:dz13 ] ;
						
						
						normal = [ [ MMVector3 alloc ] initByCrossing:v13 and:v12 ] ;
						
						[ normal normalize ] ;
						
						if( iPhase == 2 ) [ normal reverse ] ;
						
						if( iPhase == 2 )
							{
								fprintf( outFile, "%d %d %d %10.6f %10.6f %10.6f %d ", ivs[0], ivs[1], ivs[2], 
									[ normal X ], [ normal Y ], [ normal Z ], iPhase ) ;
							}
						else
							{
								fprintf( outFile, "%d %d %d %10.6f %10.6f %10.6f %d ", ivs[2], ivs[1], ivs[0], 
									[ normal X ], [ normal Y ], [ normal Z ], iPhase ) ;
							}
							
							
						[ normal release ] ; 
						[ v12 release ] ;
						[ v13 release ] ;
						
						NSArray *elemAtoms ;
						int nElemAtoms ;
						
						elemAtoms = [ [ nextCycle atoms ] allObjects ] ;
						nElemAtoms = [ elemAtoms count ] ;
						
						if( s == YES )
							{
								// Old style
								
								for( i = 0 ; i < 3 ; ++i )
									{
										if( i >= nElemAtoms ) 
											{
												fprintf( outFile, "-1 " ) ;
											}
										else
											{
												fprintf( outFile, "%d ", [ [ elemAtoms objectAtIndex:i ] intValue ] ) ;
											}
									}
							}
						else
							{
								// New Style
								
								fprintf( outFile, "%d ", nElemAtoms ) ;
								
								for( i = 0 ; i < nElemAtoms ; ++i )
									{
										fprintf( outFile, "%d ", [ [ elemAtoms objectAtIndex:i ] intValue ] ) ;
									}
									
								if( nextCycle->selfIntersection == YES )
									{
										fprintf( outFile, "S" ) ;
									}
								else
									{
										fprintf( outFile, "N" ) ;
									}
							}
								
						fprintf( outFile, "\n" ) ;
							
						// Next cycle
					}
					
				// Next phase
			}
			
		// Close output files
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{					
				fclose( subsurfaceInfo[i].outFile ) ;
			}
	
		return ;
	}
						
		
- (void) exportAsCubicUsingPath:(NSString *)p useSubsurfaces:(BOOL)sub
	{
		// Write cubic-style output file
		
		
		//	#vertices (%d) N(ew)
		// x y z xn yn zn (%f10.6 %f10.6 %f10.6 %f10.6 %f10.6 %f10.6)
		// (etc)
		// #elements #contact #reentrant #saddle
		// v1 v2 v3 enx eny enz (%d %d %d %f10.6 %f10.6 %f10.6 ) [0|1|2] #atoms (%d) atom1 atom2 ... (%d %d ...) 
		// (etc - NOTE element type code: 0=contact, 1=reentrant, 2=saddle )
		// Midpoint, normal (%f %f %f %f %f %f MP,NM )
		// v1 - v2 points, 1/3 and 2/3 ( %f %f %f %f %f %f P11 P12 )
		// v2 - v3 points, 1/3 and 2/3 ( %f %f %f %f %f %f P21 P22 )
		// v1 - v3 points, 1/3 and 2/3 ( %f %f %f %f %f %f P31 P32 )
		// Note that we maintain counterclockwise orientation throughout (i.e. following diagram):
		//	
		//	|\
		//	| \
		//	31 \
		//	|	22
		//	|	  \
		//	32	   \
		//	|		21
		//	|		  \
		//  |_11__12___\
		
		FILE *output ;
			if( !( output = fopen( [ p cString ], "w" ) ) )
			{
				printf( "ERROR: COULD NOT OPEN OUTPUT FILE %s - Exit!\n", [ p cString ] ) ;
				exit(1) ;
			}
			
		fprintf( output, "%d CUBIC\n", nVertices ) ;
		
		int i, j ;
		
		for( i = 0 ; i < nVertices ; ++i )
			{
				fprintf( output, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", 
					[ [ [ vertices objectAtIndex:i ] vertexPosition ] X ], [ [ [ vertices objectAtIndex:i ] vertexPosition ] Y ], [ [ [ vertices objectAtIndex:i ] vertexPosition ] Z ],
					[ [ [ vertices objectAtIndex:i ] normal ] X ], [ [ [ vertices objectAtIndex:i ] normal ] Y ], [ [ [ vertices objectAtIndex:i ] normal ] Z ] ) ;
			}
			
		NSEnumerator *cycleEnumerator, *arcEnumerator, *forwardEnumerator ;
		NSArray *arcs, *forward ;
		SMCycle *nextCycle ;
		SMArc *nextArc ;
		NSNumber *nextForward ;
		
		MMVector3 *pMid, *norm ;

		
		MMVector3 *vs[3] ;
		int ivs[3], vCount ;
		
		
			
				
		fprintf( output, "%d %d %d %d = nEle, nCon, nReen, nSadd\n", nElements, nContactElements, nReentrantElements,
			nSaddleElements ) ;
					
		int iPhase ;
		
		MMVector3 *interiorPoint[3][2] ;
		
		int vEdge, jPoint ;
		
		for( vEdge = 0 ; vEdge < 3 ; ++vEdge )
			{
				for( jPoint = 0 ; jPoint < 2 ; ++jPoint )
					{
						interiorPoint[vEdge][jPoint] = [ [ MMVector3 alloc ] initX:0. Y:0. Z:0. ] ;
					}
			}
		
		
		for( iPhase = 0 ; iPhase <= 2 ; ++iPhase )
			{
			
				if( iPhase == 0 )
					{
						cycleEnumerator = [ contactCycles objectEnumerator ] ;
					}
				else if( iPhase == 1 )
					{
						cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
					}
				else
					{	
						cycleEnumerator = [ saddleCycles objectEnumerator ] ;
					}
				
				
				pMid = [ [ MMVector3 alloc ] initX:0. Y:0. Z:0 ] ;
				norm = [ [ MMVector3 alloc ] initX:0. Y:0. Z:0 ] ;
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						if( [ nextCycle active ] == NO ) continue ;
						
						// Get midpoint and normal
						
						switch( iPhase )
							{
								case 0:
									[ nextCycle contactElementMidPoint:pMid andNormal:norm forMolecule:self ] ;
									break ;
									
								case 1:
									[ nextCycle reentrantElementMidPoint:pMid andNormal:norm forMolecule:self ] ;
									break ;
									
								default:
									[ nextCycle saddleElementMidPoint:pMid andNormal:norm forMolecule:self ] ;
									
							}
						
						arcs = [ nextCycle arcs ] ;
						forward = [ nextCycle forward ] ;
						
						vCount = 0 ;
						
						/*
						if( [ arcs count ] != 3 )
							{
								printf( "CYCLE OF TYPE %d HAS %d ARCS - Exit!\n", iPhase, [ arcs count ] ) ;
								exit(1) ;
							}
						*/
						
						// Note - some arcs may be skipped!!
						
						// BEWARE! Saddle arcs are ordered counterclockwise (looking down), reentrant and contact are clockwise!
						
						arcEnumerator = [ arcs objectEnumerator ] ;
						forwardEnumerator = [ forward objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								nextForward = [ forwardEnumerator nextObject ] ;
								
								if( [ nextArc skip ] == YES ) continue ;
										
								/* Don't need this check here 
								if( vCount == 3 )
									{
										printf( "TOO MANY ARCS IN CYCLE OF TYPE %d - Exit!\n", iPhase ) ;
										exit(1) ;
									}
								*/
										
								if( [ nextForward boolValue ] == YES )
									{
										vs[vCount] = [ [ nextArc startVertex ] vertexPosition ] ;
										ivs[vCount] = [ [ nextArc startVertex] index ] ;
										
										 [ nextArc arcPoint:interiorPoint[vCount][0] atFraction:(1./3.) usingMolecule:self ] ;
										 [ nextArc arcPoint:interiorPoint[vCount][1] atFraction:(2./3.) usingMolecule:self ] ;
									}
								else
									{
										vs[vCount] = [ [ nextArc endVertex ] vertexPosition ] ;
										ivs[vCount] = [ [ nextArc endVertex ] index ] ;
										
										[ nextArc arcPoint:interiorPoint[vCount][0] atFraction:(2./3.) usingMolecule:self ] ;
										[ nextArc arcPoint:interiorPoint[vCount][1] atFraction:(1./3.) usingMolecule:self ] ;
									}
									
								++vCount ;
							}
							
						if( iPhase != 2  || ( iPhase == 2 && [ [ [ [ nextCycle arcs ] lastObject ] torusSection ] skipTorus ] == NO ) )
							{
								if( vCount != 3 )
									{
										printf( "VERTEX COUNT %d IN CYCLE OF TYPE %d - Exit!\n", vCount, iPhase ) ;
										exit(1) ;
									}
							}
				
						
						if( iPhase == 2 )
							{
								fprintf( output, "%d %d %d %d ", ivs[0], ivs[1], ivs[2], 
									 iPhase ) ;
							}
						else
							{
								fprintf( output, "%d %d %d %d ", ivs[2], ivs[1], ivs[0], 
									 iPhase ) ;
							}
							
							
						
						NSArray *elemAtoms ;
						int nElemAtoms ;
						
						elemAtoms = [ [ nextCycle atoms ] allObjects ] ;
						nElemAtoms = [ elemAtoms count ] ;
						
							
						// New Style
						
						fprintf( output, "%d ", nElemAtoms ) ;
						
						for( i = 0 ; i < nElemAtoms ; ++i )
							{
								fprintf( output, "%d ", [ [ elemAtoms objectAtIndex:i ] intValue ] ) ;
							}
							
						if( nextCycle->selfIntersection == YES )
							{
								fprintf( output, "S" ) ;
							}
						else
							{
								fprintf( output, "N" ) ;
							}
							
						fprintf( output, "\n" ) ;
								
						fprintf( output, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f MPNT NORM\n",
							[ pMid X ], [ pMid Y ], [ pMid Z ], [ norm X ], [ norm Y ], [ norm Z ] ) ;
							
						if( iPhase == 2 )
							{
								fprintf( output, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P11 P12\n", 
									[ interiorPoint[0][0] X ], [ interiorPoint[0][0] Y ], [ interiorPoint[0][0] Z ],
									[ interiorPoint[0][1] X ], [ interiorPoint[0][1] Y ], [ interiorPoint[0][1] Z ] ) ;
								fprintf( output, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P21 P22\n", 
									[ interiorPoint[1][0] X ], [ interiorPoint[1][0] Y ], [ interiorPoint[1][0] Z ],
									[ interiorPoint[1][1] X ], [ interiorPoint[1][1] Y ], [ interiorPoint[1][1] Z ] ) ;
								fprintf( output, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P31 P32\n", 
									[ interiorPoint[2][0] X ], [ interiorPoint[2][0] Y ], [ interiorPoint[2][0] Z ],
									[ interiorPoint[2][1] X ], [ interiorPoint[2][1] Y ], [ interiorPoint[2][1] Z ] ) ;
							}
						else
							{
								fprintf( output, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P11 P12\n", 
									[ interiorPoint[1][1] X ], [ interiorPoint[1][1] Y ], [ interiorPoint[1][1] Z ],
									[ interiorPoint[1][0] X ], [ interiorPoint[1][0] Y ], [ interiorPoint[1][0] Z ] ) ;
								fprintf( output, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P21 P22\n", 
									[ interiorPoint[0][1] X ], [ interiorPoint[0][1] Y ], [ interiorPoint[0][1] Z ],
									[ interiorPoint[0][0] X ], [ interiorPoint[0][0] Y ], [ interiorPoint[0][0] Z ] ) ;
								fprintf( output, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P31 P32\n", 
									[ interiorPoint[2][1] X ], [ interiorPoint[2][1] Y ], [ interiorPoint[2][1] Z ],
									[ interiorPoint[2][0] X ], [ interiorPoint[2][0] Y ], [ interiorPoint[2][0] Z ] ) ;
							}
							
									
									
							
						// Next cycle
					}
					
				// Next phase
			}
			
		// All done with complete surface
		
		fclose( output ) ;
		
		if( sub == NO ) return ;
		
		// We will also write out subsurface files
		
		// First, sort subsurfaces by size, which we will do crudely as # vertices
		
		typedef struct { int subsurfaceName ; int subsurfaceSize ; int contactElems ; int reentrantElems ; int saddleElems ; FILE *outFile } subsurfaceData ;

		subsurfaceData *subsurfaceInfo = (subsurfaceData *) malloc( nSubsurfaces * sizeof( subsurfaceData ) ) ;
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				subsurfaceInfo[i].subsurfaceName = i ;
				subsurfaceInfo[i].contactElems = 0 ;
				subsurfaceInfo[i].reentrantElems = 0 ;
				subsurfaceInfo[i].saddleElems = 0 ;
				subsurfaceInfo[i].subsurfaceSize = 0 ;
			}
		
		for( i = 0 ; i < nVertices ; ++i )
			{
				int nextSubsurface = [ [ vertices objectAtIndex:i ] subsurface ] ;
				
				++subsurfaceInfo[nextSubsurface].subsurfaceSize ;
				
			}
	
		subsurfaceData tempSubsurfaceInfo ;
		
		// Bubble sort
		
		int *subsurfaceNameToRank = (int *) malloc( nSubsurfaces * sizeof( int ) ) ;
		
		for( i = 0 ; i < nSubsurfaces - 1 ; ++i )
			{
				for( j = i + 1 ; j < nSubsurfaces ; ++j )
					{
						if( subsurfaceInfo[j].subsurfaceSize > subsurfaceInfo[i].subsurfaceSize )
							{
								tempSubsurfaceInfo = subsurfaceInfo[i] ;
								subsurfaceInfo[i] = subsurfaceInfo[j] ;
								subsurfaceInfo[j] = tempSubsurfaceInfo ;
							}
					}
			}
			
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				subsurfaceNameToRank[ subsurfaceInfo[i].subsurfaceName ] = i ;
			}
			
			
		// Elements counts
				
		for( iPhase = 0 ; iPhase <= 2 ; ++iPhase )
			{
			
				if( iPhase == 0 )
					{
						cycleEnumerator = [ contactCycles objectEnumerator ] ;
					}
				else if( iPhase == 1 )
					{
						cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
					}
				else
					{	
						cycleEnumerator = [ saddleCycles objectEnumerator ] ;
					}
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						if( [ nextCycle active ] == NO ) continue ;
						
						int rank = subsurfaceNameToRank[  [ nextCycle subsurface ] ] ;
						
						if( iPhase == 0 )
							{
								++subsurfaceInfo[rank].contactElems ;
							}
						else if( iPhase == 1 )
							{
								++subsurfaceInfo[rank].reentrantElems ;
							}
						else
							{
								++subsurfaceInfo[rank].saddleElems ;
							}
					}
			}
			
		// We need to sort the vertices in the subsurfaces in accord with local (in-subsurface) index
		
		NSMutableArray *subsurfaceRankToVertices = [ [ NSMutableArray alloc ] initWithCapacity:nSubsurfaces ] ;
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				[ subsurfaceRankToVertices addObject:[ [ NSMutableArray alloc ] initWithCapacity:subsurfaceInfo[i].subsurfaceSize ] ] ;
			}
			
		for( i = 0 ; i < nVertices ; ++i )
			{
				int rank ;
				SMVertex *nextVertex = [ vertices objectAtIndex:i ] ;
				
				rank = subsurfaceNameToRank[ [ nextVertex subsurface ] ] ;
				
				[ [ subsurfaceRankToVertices objectAtIndex:rank ] addObject:nextVertex ] ;
			}
			
		// Sort vertices 
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				[ [ subsurfaceRankToVertices objectAtIndex:i ] sortUsingSelector:@selector(compareSubsurfaceVertices:) ] ;
			}
				
	
		// Open output files
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				FILE *nextFile ;
				
				char outName[500] ;
				sprintf( outName, "%s.%d", [ p cString ], i ) ;
				
				if( !( nextFile = fopen( outName, "w" ) ) )
					{
						printf( "COULD NOT OPEN OUTPUT FILE %s - Exit!\n", outName ) ;
						exit(1) ;
					}
					
				subsurfaceInfo[i].outFile = nextFile ;
			}
			
				

			
		for( j = 0 ; j < nSubsurfaces ; ++j )
			{
				NSArray *nextVertices = [ subsurfaceRankToVertices objectAtIndex:j ] ;
				
				fprintf( subsurfaceInfo[j].outFile, "%d CUBIC\n", subsurfaceInfo[j].subsurfaceSize ) ;
				
				for( i = 0 ; i < [ nextVertices count ] ; ++i )
					{
						fprintf( subsurfaceInfo[j].outFile, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", 
							[ [ [ nextVertices objectAtIndex:i ] vertexPosition ] X ], [ [ [ nextVertices objectAtIndex:i ] vertexPosition ] Y ], [ [ [ nextVertices objectAtIndex:i ] vertexPosition ] Z ],
							[ [ [ nextVertices objectAtIndex:i ] normal ] X ], [ [ [ nextVertices objectAtIndex:i ] normal ] Y ], [ [ [ nextVertices objectAtIndex:i ] normal ] Z ] ) ;
					}
			}
						
				
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{
				fprintf( subsurfaceInfo[i].outFile, "%d %d %d %d = nEle, nCon, nReen, nSadd\n", 
					(subsurfaceInfo[i].contactElems + subsurfaceInfo[i].reentrantElems + subsurfaceInfo[i].saddleElems),
					subsurfaceInfo[i].contactElems, subsurfaceInfo[i].reentrantElems, subsurfaceInfo[i].saddleElems ) ;
			}
				
		
		for( iPhase = 0 ; iPhase <= 2 ; ++iPhase )
			{
			
				if( iPhase == 0 )
					{
						cycleEnumerator = [ contactCycles objectEnumerator ] ;
					}
				else if( iPhase == 1 )
					{
						cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
					}
				else
					{	
						cycleEnumerator = [ saddleCycles objectEnumerator ] ;
					}
				
				while( ( nextCycle = [ cycleEnumerator nextObject ] ) )
					{
						FILE *outFile ;
						int rank ;
						
						if( [ nextCycle active ] == NO ) continue ;
						
						rank = subsurfaceNameToRank[ [ nextCycle subsurface ] ] ;
						
						outFile = subsurfaceInfo[rank].outFile ;
						
						arcs = [ nextCycle arcs ] ;
						forward = [ nextCycle forward ] ;
						
						vCount = 0 ;
						
						switch( iPhase )
							{
								case 0:
									[ nextCycle contactElementMidPoint:pMid andNormal:norm forMolecule:self ] ;
									break ;
									
								case 1:
									[ nextCycle reentrantElementMidPoint:pMid andNormal:norm forMolecule:self ] ;
									break ;
									
								default:
									[ nextCycle saddleElementMidPoint:pMid andNormal:norm forMolecule:self ] ;
									
							}
						
						
						/*
						if( [ arcs count ] != 3 )
							{
								printf( "CYCLE OF TYPE %d HAS %d ARCS - Exit!\n", iPhase, [ arcs count ] ) ;
								exit(1) ;
							}
						*/
						
						// Note - some arcs may be skipped!!
						
						// BEWARE! Saddle arcs are ordered counterclockwise (looking down), reentrant and contact are clockwise!
						
						arcEnumerator = [ arcs objectEnumerator ] ;
						forwardEnumerator = [ forward objectEnumerator ] ;
						
						while( ( nextArc = [ arcEnumerator nextObject ] ) )
							{
								nextForward = [ forwardEnumerator nextObject ] ;
								
								if( [ nextArc skip ] == YES ) continue ;
										
								/* Don't need this check here 
								if( vCount == 3 )
									{
										printf( "TOO MANY ARCS IN CYCLE OF TYPE %d - Exit!\n", iPhase ) ;
										exit(1) ;
									}
								*/
										
								if( [ nextForward boolValue ] == YES )
									{
										vs[vCount] = [ [ nextArc startVertex ] vertexPosition ] ;
										ivs[vCount] = [ [ nextArc startVertex] subsurfaceIndex] ;
										
										[ nextArc arcPoint:interiorPoint[vCount][0] atFraction:(1./3.) usingMolecule:self ] ;
										[ nextArc arcPoint:interiorPoint[vCount][1] atFraction:(2./3.) usingMolecule:self ] ;

									}
								else
									{
										vs[vCount] = [ [ nextArc endVertex ] vertexPosition ] ;
										ivs[vCount] = [ [ nextArc endVertex ] subsurfaceIndex ] ;
										
										[ nextArc arcPoint:interiorPoint[vCount][0] atFraction:(2./3.) usingMolecule:self ] ;
										[ nextArc arcPoint:interiorPoint[vCount][1] atFraction:(1./3.) usingMolecule:self ] ;

									}
									
								++vCount ;
							}
							
						if( iPhase != 2  || ( iPhase == 2 && [ [ [ [ nextCycle arcs ] lastObject ] torusSection ] skipTorus ] == NO ) )
							{
								if( vCount != 3 )
									{
										printf( "VERTEX COUNT %d IN CYCLE OF TYPE %d - Exit!\n", vCount, iPhase ) ;
										exit(1) ;
									}
							}
							
							
						if( iPhase == 2 )
							{
								fprintf( outFile, "%d %d %d %d ", ivs[0], ivs[1], ivs[2], 
									 iPhase ) ;
							}
						else
							{
								fprintf( outFile, "%d %d %d %d ", ivs[2], ivs[1], ivs[0], 
									 iPhase ) ;
							}
							
							
						
						NSArray *elemAtoms ;
						int nElemAtoms ;
						
						elemAtoms = [ [ nextCycle atoms ] allObjects ] ;
						nElemAtoms = [ elemAtoms count ] ;
						
							
						// New Style
						
						fprintf( outFile, "%d ", nElemAtoms ) ;
						
						for( i = 0 ; i < nElemAtoms ; ++i )
							{
								fprintf( outFile, "%d ", [ [ elemAtoms objectAtIndex:i ] intValue ] ) ;
							}
							
						if( nextCycle->selfIntersection == YES )
							{
								fprintf( outFile, "S" ) ;
							}
						else
							{
								fprintf( outFile, "N" ) ;
							}
							
						fprintf( outFile, "\n" ) ;
								
						fprintf( outFile, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f MPNT NORM\n",
							[ pMid X ], [ pMid Y ], [ pMid Z ], [ norm X ], [ norm Y ], [ norm Z ] ) ;
							
						fprintf( outFile, " MP NM\n" ) ;
						
						if( iPhase == 2 )
							{
								fprintf( outFile, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P11 P12\n", 
									[ interiorPoint[0][0] X ], [ interiorPoint[0][0] Y ], [ interiorPoint[0][0] Z ],
									[ interiorPoint[0][1] X ], [ interiorPoint[0][1] Y ], [ interiorPoint[0][1] Z ] ) ;
								fprintf( outFile, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P21 P22\n", 
									[ interiorPoint[1][0] X ], [ interiorPoint[1][0] Y ], [ interiorPoint[1][0] Z ],
									[ interiorPoint[1][1] X ], [ interiorPoint[1][1] Y ], [ interiorPoint[1][1] Z ] ) ;
								fprintf( outFile, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P31 P32\n", 
									[ interiorPoint[2][0] X ], [ interiorPoint[2][0] Y ], [ interiorPoint[2][0] Z ],
									[ interiorPoint[2][1] X ], [ interiorPoint[2][1] Y ], [ interiorPoint[2][1] Z ] ) ;
							}
						else
							{
								fprintf( outFile, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P11 P12\n", 
									[ interiorPoint[1][1] X ], [ interiorPoint[1][1] Y ], [ interiorPoint[1][1] Z ],
									[ interiorPoint[1][0] X ], [ interiorPoint[1][0] Y ], [ interiorPoint[1][0] Z ] ) ;
								fprintf( outFile, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P21 P22\n", 
									[ interiorPoint[0][1] X ], [ interiorPoint[0][1] Y ], [ interiorPoint[0][1] Z ],
									[ interiorPoint[0][0] X ], [ interiorPoint[0][0] Y ], [ interiorPoint[0][0] Z ] ) ;
								fprintf( outFile, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f P31 P32\n", 
									[ interiorPoint[2][1] X ], [ interiorPoint[2][1] Y ], [ interiorPoint[2][1] Z ],
									[ interiorPoint[2][0] X ], [ interiorPoint[2][0] Y ], [ interiorPoint[2][0] Z ] ) ;
							}
							
				
						
						
						// Next cycle
					}
					
				// Next phase
			}
			
		// Close output files
		
		for( i = 0 ; i < nSubsurfaces ; ++i )
			{					
				fclose( subsurfaceInfo[i].outFile ) ;
			}
	
		return ;
	}
						
		
							
					
- (void) probesToMOL2UsingFile:(NSString *)f
	{
		// Debugging function - report the probes in mol2 format to file
		
		FILE *outFile ;
		
		NSArray *tripleProbes = [ atomTripleToProbes allValues ] ;
		NSMutableArray *allProbes = [ NSMutableArray arrayWithCapacity:100 ] ;
		
		NSEnumerator *probeArrayEnumerator = [ tripleProbes objectEnumerator ] ;
		NSArray *nextArray ;
		
		if( ! ( outFile = fopen( [ f cString ], "w" ) ) )
			{
				printf( "COULD NOT OPEN PROBE OUTPUT FILE - Continuing ...\n" ) ;
				return ;
			}
		
		while( ( nextArray = [ probeArrayEnumerator nextObject ] ) )
			{
				[ allProbes addObjectsFromArray:nextArray ] ;
			}
		
		
		 fprintf( outFile, "# Probe positions\n\n@<TRIPOS>MOLECULE\nPROBES\n%4d%4d%4d%4d\nSMALL\nNO_CHARGES\n\n\n@<TRIPOS>ATOM\n",
				(int)[ allProbes count ], 0, 0, 0 ) ;
				
		NSEnumerator *probeEnumerator = [ allProbes objectEnumerator ] ;
		SMProbe *nextProbe ;
		int probeCount = 0 ; 
		
		while( ( nextProbe = [ probeEnumerator nextObject ] ) )
			{
				++probeCount ;
				
				fprintf( outFile,  "%4d %10s %10.4f %10.4f %10.4f O.3%8d%9s\n", probeCount, "O", [ nextProbe X ], 
   				 [ nextProbe Y ],  [ nextProbe Z ], probeCount, "WAT" ) ;
			}
			
		fprintf( outFile, "@<TRIPOS>BOND\n@end\n" ) ;
		
		fclose( outFile ) ;
		
		return ;
	}
				
				
// For debugging

- (void) writeMOEGraphicsForAllToriUsingObjectName:(NSString *)name 
	{
		int i ; 
		
		unsigned int red, blue ;
		
		red = 0xFF0000 ;
		blue = 0x0000FF ;
		
		printf( "%s = GCreate '%s' ;\n", [ name cString ], [ name cString ] ) ;
		
		for( i = 0 ; i < nTori ; ++i )
			{
				[ self writeMOEGraphicsForArc:[ tori[i] contactArcI ] usingColorIndex:blue andGraphicsObject:name ] ;
				[ self writeMOEGraphicsForArc:[ tori[i] reentrantArcR ] usingColorIndex:red andGraphicsObject:name ] ;
				[ self writeMOEGraphicsForArc:[ tori[i] contactArcJ ] usingColorIndex:blue andGraphicsObject:name ] ;
				[ self writeMOEGraphicsForArc:[ tori[i] reentrantArcL ] usingColorIndex:red andGraphicsObject:name ] ;
			}
			
		return ;
	}
	
- (void) writeMOEGraphicsForAllSaddles
	{
		
		unsigned int yellow ;
		
		yellow = 0xFFFF00 ;
		
		printf( "saddleCycles = GCreate 'saddleCycles' ;\n" ) ;
		
		NSEnumerator *cycleEnumerator ;
		
		cycleEnumerator = [ saddleCycles objectEnumerator ] ;
		
		SMCycle *nextCycle ;
		
		while( ( nextCycle = [ cycleEnumerator nextObject ] ) ) 
			{
				// Skip dead cycles
				
				if( [ nextCycle active ] == NO ) continue ;
				
				[ self writeMOEGraphicsForSaddleCycle:nextCycle usingObjectName:@"saddleCycles" ] ;
											
				// Next Cycle 
			}
			
		return ;
	}

				

- (void) writeMOEGraphicsForCycle:(SMCycle *)cyc usingObjectName:(NSString *)name
	{
		// Alternate the colors red, green, blue 
		
		NSEnumerator *arcEnumerator ;
		
		printf( "%s = GCreate '%s' ;\n", [ name cString ], [ name cString ] ) ;
		
		arcEnumerator = [ [ cyc arcs ] objectEnumerator ] ;
		
		SMArc *nextArc ;
		unsigned int color ;
		 
		color = 0xFF0000 ;
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				[ self writeMOEGraphicsForArc:nextArc usingColorIndex:color andGraphicsObject:name ] ;
				
				if( color == 0xFF0000 )
					{
						color = 0x00FF00 ;
					}
				else if( color == 0x00FF00 )
					{
						color = 0x0000FF ;
					}
				else
					{
						color = 0xFF0000 ;
					}
			}
			
		return ;
	}
		 
		 
		
		
		
		
- (void) writeMOEGraphicsForArc:(SMArc *)a usingColorIndex:(unsigned int)c andGraphicsObject:(NSString *)g
	{
		// This prints commands to standard out, which are then introduced to MOE
		
		// Use ten angular divisions
		
		int DIV = 9 ;
		
		double angle = [ a angle ] ;
		double deltaAngle = angle/DIV ;
		double radius = [ a arcRadius ] ;
		
		MMVector3 *startU = [ a startU ] ;
		MMVector3 *uPerp = [ a uPerp ] ;
		MMVector3 *center = [ a arcCenter ] ;
		
		int i ; 
		
		for( i = 0 ; i < DIV; ++i )
			{
				double startA, endA ;
				
				startA = i*deltaAngle ;
				endA = (i + 1)*deltaAngle ;
				
				double x1, y1, z1, x2, y2, z2 ;
				
				x1 = [ center X ] + radius*(cos(startA)*[ startU X ] + sin(startA)*[ uPerp X ]) ;
				y1 = [ center Y ] + radius*(cos(startA)*[ startU Y ] + sin(startA)*[ uPerp Y ]) ;
				z1 = [ center Z ] + radius*(cos(startA)*[ startU Z ] + sin(startA)*[ uPerp Z ]) ;
				
				
				x2 = [ center X ] + radius*(cos(endA)*[ startU X ] + sin(endA)*[ uPerp X ]) ;
				y2 = [ center Y ] + radius*(cos(endA)*[ startU Y ] + sin(endA)*[ uPerp Y ]) ;
				z2 = [ center Z ] + radius*(cos(endA)*[ startU Z ] + sin(endA)*[ uPerp Z ]) ;
				
				if( i < 3 )
					{
						c = 0xFF0000 ;
					}
				else if( i < 6 )
					{
						c = 0x00FF00 ;
					}
				else
					{
						c = 0x0000FF ;
					}

				
				printf( "GLine[ %s, %#x, %f, %f, %f, %f, %f, %f ]; \n", [ g cString ], c, x1, y1, z1, x2, y2, z2 ) ;
			}
			
		return ;
	}
				
			
- (void) writeMOEGraphicsForSaddleArc:(SMArc *)a usingColorIndex:(unsigned int)c andGraphicsObject:(NSString *)g
	{
		// This prints commands to standard out, which are then introduced to MOE
		
		// Use ten divisions
		
		int DIV = 9 ;
		
		if( [ a phiStart ] < 0. ) return ;
		
		double phi1 = [ a phiStart ] ;
		double phi2 = [ a phiEnd ] ;
		
		double theta1 = [ a thetaStart ] ;
		double theta2 = [ a thetaEnd ] ;
		
		double deltaPhi = (phi2 - phi1)/DIV ;
		double deltaTheta = (theta2 - theta1)/DIV ;
		
		
		int i ; 
		
		MMVector3 *prePos = [ a computePositionForTheta:theta1 andPhi:phi1 usingMolecule:self allowSelfIntersection:allowSelfIntersection normal:nil ] ;
		
		for( i = 1 ; i <= DIV; ++i )
			{
			
				double phi, theta ;
				
				phi = phi1 + i*deltaPhi ;
				theta = theta1 + i*deltaTheta ;
				
				MMVector3 *pos ;
				
				pos = [ a computePositionForTheta:theta andPhi:phi usingMolecule:self allowSelfIntersection:allowSelfIntersection normal:nil ] ;
				
				double x1, y1, z1, x2, y2, z2 ;
				
				x1 = [ prePos X ] ;
				y1 = [ prePos Y ] ;
				z1 = [ prePos Z ] ;
				
				
				x2 = [ pos X ] ;
				y2 = [ pos Y ] ;
				z2 = [ pos Z ] ;
				
				if( i <= 3 )
					{
						c = 0xFF0000 ;
					}
				else if( i <= 6 )
					{
						c = 0x00FF00 ;
					}
				else
					{
						c = 0x0000FF ;
					}
				
				printf( "GLine[ %s, %#x, %f, %f, %f, %f, %f, %f ]; \n", [ g cString ], c, x1, y1, z1, x2, y2, z2 ) ;
				
				if( i > 1 ) [ prePos release ] ;
				
				prePos = pos ;
			}
			
		return ;
	}
	
- (void) writeMOEGraphicsForSaddleCycle:(SMCycle *)cyc usingObjectName:(NSString *)name
	{
		// Write all arcs in yellow, using phi, theta. Show ref vector, and vector to R probe position from base
	
		NSEnumerator *arcEnumerator ;
		
		arcEnumerator = [ [ cyc arcs ] objectEnumerator ] ;
		
		SMArc *nextArc ;
		
		printf( "%s = GCreate '%s' ;\n", [ name cString ], [ name cString ] ) ;
		
		while( ( nextArc = [ arcEnumerator nextObject ] ) )
			{
				[ self writeMOEGraphicsForSaddleArc:nextArc usingColorIndex:0xFFFF00 andGraphicsObject:name ] ;
			}
			
		// Write green vector from torus base to probeR position
		
		double baseX, baseY, baseZ, probeRX, probeRY, probeRZ, refX, refY, refZ, x2, y2, z2 ;
		SMProbe *RP ;
		
		nextArc = [ [ cyc arcs ] objectAtIndex:0 ] ;
		
		baseX = [ [ [ nextArc torusSection ] base ] X ] ;
		baseY = [ [ [ nextArc torusSection ] base ] Y ] ;
		baseZ = [ [ [ nextArc torusSection ] base ] Z ] ;
		
		RP = [ nextArc torusSection ]->probeR ;
		
		probeRX = [ RP X ] ;
		probeRY = [ RP Y ] ;
		probeRZ = [ RP Z ] ;
		
		//printf( "GLine[ %s, %#x, %f, %f, %f, %f, %f, %f ]; \n", [ name cString ], 0x00FF00, baseX, baseY, baseZ, probeRX, probeRY, probeRZ ) ;
		
		refX = [ [ [ nextArc torusSection ] refR ] X ] ;
		refY = [ [ [ nextArc torusSection ] refR ] Y ] ;
		refZ = [ [ [ nextArc torusSection ] refR ] Z ] ;
		
		x2 = baseX + 2.*refX ;
		y2 = baseY + 2.*refY ;
		z2 = baseZ + 2.*refZ ;
		
		//printf( "GLine[ %s, %#x, %f, %f, %f, %f, %f, %f ]; \n", [ name cString ], 0xFF0000, baseX, baseY, baseZ, x2, y2, z2 ) ;
		
		// Make axis blue
		
		
		x2 = baseX + 2.*[ [ [ nextArc torusSection ] axis ] X ] ;
		y2 = baseY + 2.*[ [ [ nextArc torusSection ] axis ] Y ] ;
		z2 = baseZ + 2.*[ [ [ nextArc torusSection ] axis ] Z ] ;

		//printf( "GLine[ %s, %#x, %f, %f, %f, %f, %f, %f ]; \n", [ name cString ], 0x0000FF, baseX, baseY, baseZ, x2, y2, z2 ) ;
		
		return ;
	}
		
					
- (void) writeMOEGraphicsForHistoryOfCycle:(SMCycle *)cyc
	{
		[ self writeMOEGraphicsForCycle:cyc usingObjectName:@"Child" ] ;
		
		SMCycle *theParent ;
		
		theParent = [ cyc parentCycle ] ;
		
		int index = 0 ;
		
		while( theParent )
			{
				NSString *name ;
				
				name = [ NSString stringWithFormat:@"Parent%d",index ] ;
				
				[ self writeMOEGraphicsForCycle:theParent usingObjectName:name ] ;
				
				++index ;
				
				theParent = [ theParent parentCycle ] ;
			}
			
		return ;
	}
				
- (void) writeMOEGraphicsForHistoryOfSaddleCycle:(SMCycle *)cyc
	{
		[ self writeMOEGraphicsForSaddleCycle:cyc usingObjectName:@"Child" ] ;
		
		SMCycle *theParent ;
		
		theParent = [ cyc parentCycle ] ;
		
		int index = 0 ;
		
		while( theParent )
			{
				NSString *name ;
				
				name = [ NSString stringWithFormat:@"Parent%d",index ] ;
				
				[ self writeMOEGraphicsForSaddleCycle:theParent usingObjectName:name ] ;
				
				++index ;
				
				theParent = [ theParent parentCycle ] ;
			}
			
		return ;
	}
				
- (void) writeMOEGraphicsForContactCycles
	{
		printf( "contactCycles = GCreate 'contactCycles' ;\n" ) ;
		
		NSEnumerator *cycleEnumerator ;
		
		cycleEnumerator = [ contactCycles objectEnumerator ] ;
		
		SMCycle *nextCycle ;
		
		while( ( nextCycle = [ cycleEnumerator nextObject ] ) ) 
			{
				// Skip dead cycles
				
				if( [ nextCycle active ] == NO ) continue ;
				
				NSEnumerator *arcEnumerator ;
				
				arcEnumerator = [ [ nextCycle arcs ] objectEnumerator ] ;
				
				SMArc *nextArc ;
				unsigned int color = 0x0000FF ;
				 
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						// Check to see if arc has a parent cycle with index < current cycle
						
						NSEnumerator *parentCycleEnumerator ;
						SMCycle *nextParentCycle ;
						BOOL skip ;
						
						skip = NO ;
						
						parentCycleEnumerator = [ [ nextArc parentCycles ] objectEnumerator ] ;
						
						while( ( nextParentCycle = [ parentCycleEnumerator nextObject ] ) )
							{
								if( nextParentCycle < nextCycle )
									{
										skip = YES ;
										break ;
									}
									
							}
							
						if( skip == YES ) continue ;
						
						[ self writeMOEGraphicsForArc:nextArc usingColorIndex:color andGraphicsObject:@"contactCycles" ] ;
						
						// Next arc
						
					}
					
				// Next Cycle 
			}
			
		return ;
	}
						
- (void) writeMOEGraphicsForSaddleCycles
	{
		printf( "saddleCycles = GCreate 'saddleCycles' ;\n" ) ;
		
		NSEnumerator *cycleEnumerator ;
		
		cycleEnumerator = [ saddleCycles objectEnumerator ] ;
		
		SMCycle *nextCycle ;
		
		while( ( nextCycle = [ cycleEnumerator nextObject ] ) ) 
			{
				// Skip dead cycles
				
				if( [ nextCycle active ] == NO ) continue ;
				
				NSEnumerator *arcEnumerator ;
				
				arcEnumerator = [ [ nextCycle arcs ] objectEnumerator ] ;
				
				SMArc *nextArc ;
				unsigned int color = 0xFFFF00 ;
				 
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						// Check to see if arc has a parent cycle with index < current cycle
						
						NSEnumerator *parentCycleEnumerator ;
						SMCycle *nextParentCycle ;
						BOOL skip ;
						
						skip = NO ;
						
						parentCycleEnumerator = [ [ nextArc parentCycles ] objectEnumerator ] ;
						
						while( ( nextParentCycle = [ parentCycleEnumerator nextObject ] ) )
							{
								if( nextParentCycle < nextCycle )
									{
										skip = YES ;
										break ;
									}
									
							}
							
						if( skip == YES ) continue ;
						
						[ self writeMOEGraphicsForSaddleArc:nextArc usingColorIndex:color andGraphicsObject:@"saddleCycles" ] ;
						
						// Next arc
						
					}
					
				// Next Cycle 
			}
			
		return ;
	}

- (void) writeMOEGraphicsForReentrantCycles
	{
		printf( "reentrantCycles = GCreate 'reentrantCycles' ;\n" ) ;
		
		NSEnumerator *cycleEnumerator ;
		
		cycleEnumerator = [ reentrantCycles objectEnumerator ] ;
		
		SMCycle *nextCycle ;
		
		while( ( nextCycle = [ cycleEnumerator nextObject ] ) ) 
			{
				// Skip dead cycles
				
				if( [ nextCycle active ] == NO ) continue ;
				
				NSEnumerator *arcEnumerator ;
				
				arcEnumerator = [ [ nextCycle arcs ] objectEnumerator ] ;
				
				SMArc *nextArc ;
				unsigned int color = 0xFF0000 ;
				 
				while( ( nextArc = [ arcEnumerator nextObject ] ) )
					{
						// Check to see if arc has a parent cycle with index < current cycle
						
						NSEnumerator *parentCycleEnumerator ;
						SMCycle *nextParentCycle ;
						BOOL skip ;
						
						skip = NO ;
						
						parentCycleEnumerator = [ [ nextArc parentCycles ] objectEnumerator ] ;
						
						while( ( nextParentCycle = [ parentCycleEnumerator nextObject ] ) )
							{
								if( nextParentCycle < nextCycle )
									{
										skip = YES ;
										break ;
									}
									
							}
							
						if( skip == YES ) continue ;
						
						[ self writeMOEGraphicsForArc:nextArc usingColorIndex:color andGraphicsObject:@"reentrantCycles" ] ;
						
						// Next arc
						
					}
					
				// Next Cycle 
			}
			
		return ;
	}

		
		


@end
