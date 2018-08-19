#ifndef OSTN15_h
#define OSTN15_h

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {
	int num;
	int emin;
	int nmin;
	int emax;
	int nmax;
	char nameUTF8[77];
	char sheetUTF8[46];
} OSMap;

typedef struct {
	int deg;  // usual range 0 - 180 (S or W)
	int min;  // usual range 0 - 60
	double sec;
	bool westOrSouth;
} DegMinSec;

typedef struct {
	double e;
	double n;
	double elevation;
	unsigned char geoid;
} EastingNorthing;

typedef struct {
	char letters[2];
	int e;
	int n;
	int resolution;  // 1/10/100/1000m
} GridRefComponents;

typedef struct {
	double lat;
	double lon;
	double elevation;
	unsigned char geoid;
} LatLonDecimal;

typedef struct {
	DegMinSec lat;
	DegMinSec lon;
	double elevation;
	unsigned char geoid;
} LatLonDegMinSec;

typedef struct {
	double semiMajorAxis;
	double semiMinorAxis;
} Ellipsoid;

typedef struct {
	double centralMeridianScale;
	LatLonDecimal trueOriginLatLon;
	EastingNorthing trueOriginEastingNorthing;
} MapProjection;

typedef struct {
	unsigned int eMin : 10;
	unsigned int eCount : 10;
	unsigned int offset : 20;
}  OSTN02Index;

typedef struct {
	unsigned int eShift : 15;
	unsigned int nShift : 15;
	unsigned int gShift : 14;
	unsigned int gFlag : 4;
}  OSTN02Datum;

typedef struct LTF
{
	char Name[2];
	int Xblock, Yblock;
	double **data;
	struct LTF *next;
} LIDARTerrainFile;

#define dbl(n) ((n)*(n))
#define tri(n) ((n)*(n)*(n))


EastingNorthing   *eastingNorthingFromLatLon(LatLonDecimal *latLon, Ellipsoid *ellipsoid, MapProjection *projection);
EastingNorthing   *ETRS89EastingNorthingFromETRS89LatLon(LatLonDecimal *latLon);
LatLonDecimal     *latLonFromEastingNorthing(EastingNorthing *en, Ellipsoid *ellipsoid, MapProjection *projection);
LatLonDecimal     *ETRS89LatLonFromETRS89EastingNorthing(EastingNorthing *en);
EastingNorthing   *OSTN02ShiftsForIndices(int eIndex, int nIndex);
EastingNorthing   *shiftsForEastingNorthing(EastingNorthing *en);
EastingNorthing   *OSGB36EastingNorthingFromETRS89EastingNorthing(EastingNorthing *en);
EastingNorthing   *ETRS89EastingNorthingFromOSGB36EastingNorthing(EastingNorthing *en);
DegMinSec         degMinSecFromDecimal(double dec);
double            decimalFromDegMinSec(DegMinSec *dms);
LatLonDegMinSec   latLonDegMinSecFromLatLonDecimal(LatLonDecimal dec);
LatLonDecimal     latLonDecimalFromLatLonDegMinSec(LatLonDegMinSec dms);
bool              test(bool noisily);

double findHeight(EastingNorthing *en, char *blk, LIDARTerrainFile *ter);
int nextOSExplorerMap(EastingNorthing *en, int prevMap);
int checkGridPosition(EastingNorthing *en, int prevMap);
LIDARTerrainFile * findFile(char * BLOCK, int fileE, int fileN, LIDARTerrainFile * ter);

#endif

