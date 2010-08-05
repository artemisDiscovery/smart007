#import <Foundation/Foundation.h>
#import "SMMol.h"
#include <stdlib.h>

#include "platform.h"
#ifdef LINUX
#define TRUE 1
#define FALSE 0
#endif

int main (int argc, const char * argv[]) {
    NSAutoreleasePool * pool = [[NSAutoreleasePool alloc] init];
	
	// This is the long-awaited update of SMART
	
	// Even though SMART needs to be superseded by a physics-based surface, there is still need
	// for a relatively fast geometrical algorithm
	
	//SYNTAX:
	//smart007 [ flags ] <mol2 file> <radii file> <output file> 
	//	where legal flags are:
		
	//	-d <density length param (used by default; default value = 0.5>
	// -da <density angle param>
	// -pr <probe radius; default = 1.4>
	
	// We will add support for free atom triangulation!
	
	BOOL useLengthParam = YES, useGeodesics = YES, useOldFormat = YES, generateSubsurfaces = NO, generateBlockingProbes = NO  ;
	BOOL handleReentrantCarefully = NO ;
	BOOL writeCubicElements = NO ;
	double probeRadius = 1.4, contactTolerance = 0. ;
	double lengthParam = 0.5,  angleParam = 0. ;
	double arcSkipWidth = 0.05, randomization = 0. ;
	int warningLevel = 0 ;
	
	// Index of first solvent molecule - solvents are expected at END of the file
	int solventStartIndex = -1 ;
	
	// Allow self-intersection?
	
	BOOL SIFlag = NO ;
	
	SMMol *theMolecule ;
	enum parseStates { GETTOKEN, GETFLAG }  ;
	enum flagTypes { LENGTHDENSITY, ANGLEDENSITY, PROBERADIUS, CONTACTTOLERANCE, USEGEODESICS, SKIPWIDTH, ALLOWSELFINT, WARNINGLEVEL, RANDOMIZE, NEWFORMAT,
						SUBSURFACES, BLOCKINGPROBES, REENTRANTWITHCARE, CUBICELEMENTS, PROBEFILE } ;

	if( argc < 4 )
		{
			printf( "USAGE: smart007 [ flags ] <mol2 file> <radii file> <output file> \n" ) ;
			printf( "Flags:\n" ) ;
			printf( "\t-de <density length param (used by default; default value = 1.0)>\n" ) ;
			//printf( "\t-da <density angle param (degrees)>\n" ) ;
			printf( "\t-pr <probe radius; default = 1.4>\n" ) ;
			printf( "\t-c <contact tolerance; default = 0.>\n" ) ;
			printf( "\t-s <narrow arc skip-width; default = 0.05>\n" ) ;
			printf( "\t-asi <allow self-intersection; default = NO>\n" ) ;
			printf( "\t-ra <randomization for coodinates; default = 0., suggested = 0.001>\n" ) ;
			printf( "\t-nf <new file formatl default = NO>\n" ) ;
			printf( "\t-ss <generate subsurface files; default = NO>\n" ) ;
			printf( "\t-w <warning level; default = 0 (print all warnings)>\n" ) ;
			printf( "\t-b <atom index; generate blocking probes for solvent starting at this index>\n" ) ;
			printf( "\t-rc <handle reentrant regions with care; default = NO>\n" ) ;
			printf( "\t-cubic <generate cubic elements; default = NO>\n" ) ;
			printf( "\t-pf <file name> <generate probe file, mol2 format, atom type = O.3>\n" ) ;
			
			exit(1) ;
		}
		
	int parseState, flagType ;
	NSString *molFile = nil, *paramFile = nil, *outFile = nil, *probeFile = nil ;
	
	parseState = GETTOKEN ;
	int i ;
	
	for( i = 1 ; i < argc ; ++i )
		{
			if( parseState == GETTOKEN )
				{
					if( argv[i][0] == '-' )
						{
							parseState = GETFLAG ;
							
							if( strcasecmp( &argv[i][1], "de" ) == 0 )
								{
									flagType = LENGTHDENSITY ;
								}
							else if( strcasecmp( &argv[i][1], "da" ) == 0 )
								{
									flagType = ANGLEDENSITY ;
								}	
							else if( strcasecmp( &argv[i][1], "asi" ) == 0 )
								{
									flagType = ALLOWSELFINT ;
								}																							
							else if( strcasecmp( &argv[i][1], "pr" ) == 0 )
								{
									flagType = PROBERADIUS ;
								}
							else if( strcasecmp( &argv[i][1], "pf" ) == 0 )
								{
									flagType = PROBEFILE ;
								}								
							else if( strcasecmp( &argv[i][1], "ra" ) == 0 )
								{
									flagType = RANDOMIZE ;
								}								
							else if( strcasecmp( &argv[i][1], "c" ) == 0 )
								{
									flagType = CONTACTTOLERANCE ;
								}			
							else if( strcasecmp( &argv[i][1], "g" ) == 0 )
								{
									flagType = USEGEODESICS ;
								}
							else if( strcasecmp( &argv[i][1], "s" ) == 0 )
								{
									flagType = SKIPWIDTH ;
								}	
							else if( strcasecmp( &argv[i][1], "nf" ) == 0 )
								{
									flagType = NEWFORMAT ;
								}		
							else if( strcasecmp( &argv[i][1], "ss" ) == 0 )
								{
									flagType = SUBSURFACES ;
								}																								
							else if( strcasecmp( &argv[i][1], "w" ) == 0 )
								{
									flagType = WARNINGLEVEL ;
								}		
							else if( strcasecmp( &argv[i][1], "b" ) == 0 )
								{
									flagType = BLOCKINGPROBES ;
								}					
							else if( strcasecmp( &argv[i][1], "rc" ) == 0 )
								{
									flagType = REENTRANTWITHCARE ;
								}		
							else if( strcasecmp( &argv[i][1], "cubic" ) == 0 )
								{
									flagType = CUBICELEMENTS ;
								}										
							else
								{
									printf( "CAN'T INTERPRET FLAG - Exit!\n" ) ;
									exit(1) ;
								}
								
							continue ;
						}
					else
						{
							if( ! molFile )
								{
									molFile = [ NSString stringWithCString:argv[i] ] ;
								}
							else if( ! paramFile )
								{
									paramFile = [ NSString stringWithCString:argv[i] ] ;
								}
							else if( ! outFile )
								{
									outFile = [ NSString stringWithCString:argv[i] ] ;
								}
							else
								{
									printf( "MORE THAN THREE FILE NAMES - Exit!\n" ) ;
									exit(1) ;
								}
								
							continue ;
						}
				}
			else
				{
					// In GETFLAG state 
					
					parseState = GETTOKEN ;
					
					switch( flagType )
						{
							case LENGTHDENSITY:
								lengthParam = atof( argv[i] ) ;
								useLengthParam = YES ;
								break ;
								
							case ANGLEDENSITY:
								angleParam = atof( argv[i] ) ;
								useLengthParam = NO ;
								break ;
								
							case PROBERADIUS:
								probeRadius = atof( argv[i] ) ;
								break ;
								
							case PROBEFILE:
								probeFile = [ [ NSString alloc ] initWithCString:argv[i] encoding:NSASCIIStringEncoding ] ;
								break ;								
								
							case RANDOMIZE:
								randomization = atof( argv[i] ) ;
								break ;								
								
							case CONTACTTOLERANCE:
								contactTolerance = atof( argv[i] ) ;
								break ;

							case SKIPWIDTH:
								arcSkipWidth = atof( argv[i] ) ;
								break ;
								
							case WARNINGLEVEL:
								warningLevel = atoi( argv[i] ) ;
								break ;
								
							case BLOCKINGPROBES:
								generateBlockingProbes = YES ;
								// NOTE: WE ASSUME USER INDEXES ATOMS STARTING AT 1
								solventStartIndex = atoi( argv[i] ) - 1 ;
								break ;
								
								
							case USEGEODESICS:
								if( argv[i][0] == 'y' || argv[i][0] == 'Y' )
									{
										useGeodesics = YES ;
									}
								else
									{
										useGeodesics = NO ;
									}
									
								break ;
								
							case NEWFORMAT:
								if( argv[i][0] == 'y' || argv[i][0] == 'Y' )
									{
										useOldFormat = NO ;
									}
								else
									{
										useOldFormat = YES ;
									}
									
								break ;
								
									
							case ALLOWSELFINT:
								if( argv[i][0] == 'y' || argv[i][0] == 'Y' )
									{
										SIFlag = YES ;
									}
								else
									{
										SIFlag = NO ;
									}
									
								break ;
									
							case SUBSURFACES:
								if( argv[i][0] == 'y' || argv[i][0] == 'Y' )
									{
										generateSubsurfaces = YES ;
									}
								else
									{
										generateSubsurfaces = NO ;
									}
									
									
								break ;
								
							case REENTRANTWITHCARE:
								if( argv[i][0] == 'y' || argv[i][0] == 'Y' )
									{
										handleReentrantCarefully = YES ;
									}
								else
									{
										handleReentrantCarefully = NO ;
									}
									
									
								break ;
								
							case CUBICELEMENTS:
								if( argv[i][0] == 'y' || argv[i][0] == 'Y' )
									{
										writeCubicElements = YES ;
									}
								else
									{
										writeCubicElements = NO ;
									}
									
									
								break ;
	
								
								
						}
				}
				
			// Get next argument
		}
						
	if( ! molFile || ! paramFile || ! outFile )
		{
			printf( "USAGE: smart007 [ flags ] <mol2 file> <radii file> <output file> \n" ) ;
			printf( "Flags:\n" ) ;
			printf( "\t-d <density length param (used by default; default value = 0.5)>\n" ) ;
			//printf( "\t-da <density angle param (degrees)>\n" ) ;
			printf( "\t-pr <probe radius; default = 1.4>\n" ) ;
			printf( "\t-c <contact tolerance; default = 0.>\n" ) ;
			printf( "\t-s <narrow torus arc-width; default = 0.05>\n" ) ;
			printf( "\t-asi <allow self-intersection; default = NO>\n" ) ;
			printf( "\t-ra <randomization for coodinates; default = 0., suggested = 0.001>\n" ) ;
			printf( "\t-nf <new file formatl default = NO>\n" ) ;
			printf( "\t-ss <generate subsurface files; default = NO>\n" ) ;
			printf( "\t-w <warning level; default = 0 (print all warnings)>\n" ) ;
			printf( "\t-b <atom index; generate blocking probes for solvent starting at this index>\n" ) ;
			printf( "\t-rc <handle reentrant regions with care; default = NO>\n" ) ;
			printf( "\t-cubic <generate cubic elements; default = NO>\n" ) ;
			printf( "\t-pf <file name> <generate probe file, mol2 format, atom type = O.3>\n" ) ;

			
			
			exit(1) ;
		}

	// Import mol2 file
	
	theMolecule = [ [ SMMol alloc ] initWithMOL2File:molFile andRadiiFile:paramFile allowSelfIntersection:SIFlag randomizeUsing:randomization ] ;
	
	if( ! theMolecule )
		{
			printf( "MOLECULE IMPORT FAILED - Exit!\n" ) ;
			exit(2) ;
		}
		
	// Set warning level
	
	theMolecule->warningLevel = warningLevel ;
		
	// Create atom grid
	
	[ theMolecule assignAtomsToGridUsingProbeRadius:probeRadius ] ;
	
	// Atom pairs...
	
	[ theMolecule findAtomNeighbors ] ;
	
	// Assign probes
	
	//[ theMolecule assignProbePositionsUsingContactTolerance:contactTolerance ] ;
	[ theMolecule assignProbePositionsVer2UsingContactTolerance:contactTolerance generateBlockingProbes:(BOOL)generateBlockingProbes usingSolventStart:(int)solventStartIndex ] ;
	
	if( probeFile )
		{
			[ theMolecule probesToMOL2UsingFile:probeFile ] ;
		}
		
	
	// Torus sections
	
	[ theMolecule generateToriUsingSkipWidth:arcSkipWidth ] ;
	
	// Initial contact cycles
	
	[ theMolecule generateContactCyclesUsingDivision:lengthParam ] ;
	
	//As a test
	//[ theMolecule cullContactCycles ] ;

	// Initial reentrant cycles
	
	[ theMolecule generateReentrantCyclesWithReentrantHandling:handleReentrantCarefully ] ;
	
	//[ theMolecule combineReentrantCycles ] ;
	
	[ theMolecule reduceReentrantCyclesUsingDivisionParameter:lengthParam ] ;	
	
	[ theMolecule combineContactCycles ] ;
	
	[ theMolecule reduceContactCyclesUsingGeodesics:useGeodesics andDivisionParameter:lengthParam ] ;
	
	[ theMolecule generateSaddleCycles ] ;
	
	[ theMolecule reduceSaddleCyclesUsingDivisionParameter:lengthParam ] ;
	
	// Just in case any cycles messed with in preceding step...
	
	// NOTE that I may not need the initial invocation of these methods above - however, there may be dependencies I have forgotten
	// I am taking the safer route...
	
	BOOL hadReduction = YES ;
	
	while( hadReduction == YES )
		{
			hadReduction = NO ;
	
			if( [ theMolecule reduceContactCyclesUsingGeodesics:useGeodesics andDivisionParameter:lengthParam ] == YES )
				{
					hadReduction = YES ;
				}
				
			if( [ theMolecule reduceReentrantCyclesUsingDivisionParameter:lengthParam ] == YES )
				{
					hadReduction = YES ;
				}
				
			if( [ theMolecule reduceSaddleCyclesUsingDivisionParameter:lengthParam ] == YES )
				{
					hadReduction = YES ;
				}
		}
	
	[ theMolecule generateVerticesUsingSubsurfaces:generateSubsurfaces ] ;
	
	if( writeCubicElements == NO )
		{
			[ theMolecule exportAsFlatsUsingPath:outFile 
					oldStyle:useOldFormat useSubsurfaces:generateSubsurfaces ] ;
		}
	else
		{
			[ theMolecule exportAsCubicUsingPath:outFile 
					useSubsurfaces:generateSubsurfaces ] ;
		}

	


    [pool release];
    return 0;
}
