//
//  OSTN15.c
//  OSTN15
//
//  Created by Tony Chen on 21/7/2018.
//  Copyright (c) 2018 Tony Chen. All rights reserved.
//

#include "OSTN15.h"
#include "crc32.h"
#include "fancyOut.h"
#include "constants.data"
#include "shifts.index.data"
#include "shifts.data"
#include "geoids.data"
#include "gridRef.data"
#include "testCoords.data"
#include "testConvergences.data"
#include "explorerMaps.data"
#include "gridPosition.data"
#include "math.h"
#include <stdio.h>
#include "D:\\Program Files\\MATLAB\\R2018a\\extern\\include\mex.h"

#define LENGTH_OF(x) (sizeof (x) / sizeof *(x))
#define originalIndicesCRC 244629328  // these won't be robust against differing endianness or compiler packing of bit-structs
#define originalDataCRC    790474494

static char blkNow[2];
static int blke=0, blkn=0;

double piOver180 = 0.017453292519943295769236907684886127134428718885417254560971914401710091146034494436822415696345094822;
double oneEightyOverPi = 57.29577951308232087679815481410517033240547246656432154916024386120284714832155263244096899585111094418;


EastingNorthing *eastingNorthingFromLatLon(LatLonDecimal *latLon, Ellipsoid *ellipsoid, MapProjection *projection) 
{
	//Constants
	double a = ellipsoid->semiMajorAxis;
	double b = ellipsoid->semiMinorAxis;
	double e2 = 1.0 - (b * b) / (a * a);

	double F0 = projection->centralMeridianScale;
	double N0 = projection->trueOriginEastingNorthing.n;
	double E0 = projection->trueOriginEastingNorthing.e;

	double phi = piOver180 * latLon->lat;
	double lambda = piOver180 * latLon->lon;
	double PHI0 = piOver180 * projection->trueOriginLatLon.lat;
	double LAMBDA0 = piOver180 * projection->trueOriginLatLon.lon;

	//Annexe B Realisation
	double n = (a - b) / (a + b); //B2
	double v = a * F0 * pow(1.0 - e2 * pow(sin(phi), 2.0), -0.5); //B3
	double rho = a * F0 * (1.0 - e2) * pow(1.0 - e2 * pow(sin(phi), 2.0), -1.5); //B4
	double eta2 = v / rho - 1.0; //B5

	double M = b * F0 * (
						(1.0 + n + 1.25 * dbl(n) + 1.25 * tri(n))*(phi - PHI0) -
						(3.0 * n + 3.0 * dbl(n) + 21.0 / 8 * tri(n))*sin(phi - PHI0)*cos(phi +PHI0) +
						(1.875 * dbl(n) + 1.875 * tri(n))*sin(2*(phi - PHI0))*cos(2*(phi + PHI0)) -
						35.0 / 24.0 * tri(n) * sin(3*(phi - PHI0))*cos(3*(phi + PHI0))
						); //B6
	
	double I = M + N0;
	double II = v / 2.0 * sin(phi) * cos(phi);
	double III = v / 24.0 * sin(phi) * tri(cos(phi)) * (5.0 - dbl(tan(phi)) + 9.0 * eta2);
	double IIIA = v / 720.0 * sin(phi) * dbl(cos(phi)) * tri(cos(phi)) * (61.0 - 58.0 * dbl(tan(phi)) + dbl(dbl(tan(phi))));
	double IV = v * cos(phi);
	double V = v / 6.0 * tri(cos(phi)) * (v / rho - dbl(tan(phi)));
	double VI = v / 120.0 * dbl(cos(phi)) * tri(cos(phi)) * (5.0 - 18.0 * dbl(tan(phi)) + dbl(dbl(tan(phi))) + 14.0 * eta2 - 58.0 * eta2 * dbl(tan(phi)));
	double N = I + II * dbl(lambda - LAMBDA0) + III * dbl(dbl(lambda - LAMBDA0)) + IIIA * dbl(tri(lambda - LAMBDA0)); //B7
	double E = E0 + IV * (lambda - LAMBDA0) + V * tri(lambda - LAMBDA0) + VI * dbl(lambda - LAMBDA0) * tri(lambda - LAMBDA0); //B8

	EastingNorthing *en = (EastingNorthing*)malloc(sizeof(EastingNorthing));
	en->n = N;
	en->e = E;
	en->elevation = latLon->elevation;
	en->geoid = latLon->geoid;

	return en;
}

EastingNorthing *ETRS89EastingNorthingFromETRS89LatLon(LatLonDecimal *latLon)
{
	return eastingNorthingFromLatLon(latLon, &GRS80Ellipsoid, &NationalGridProj);
}

LatLonDecimal *latLonFromEastingNorthing(EastingNorthing *en, Ellipsoid *ellipsoid, MapProjection *projection)
{
	//Constants
	double a = ellipsoid->semiMajorAxis;
	double b = ellipsoid->semiMinorAxis;
	double e2 = 1 - (b * b) / (a * a);

	double F0 = projection->centralMeridianScale;
	double N0 = projection->trueOriginEastingNorthing.n;
	double E0 = projection->trueOriginEastingNorthing.e;

	double PHI0 = piOver180 * projection->trueOriginLatLon.lat;
	double LAMBDA0 = piOver180 * projection->trueOriginLatLon.lon;

	double N = en->n;
	double E = en->e;

	//First compute: C1
	double n = (a - b) / (a + b); //B2
	double n2 = n * n;
	double n3 = n2 * n;

	double phiN = PHI0;  // this is phi' in the OS docs
	double M = 0.0;
	double deltaPhi, sumPhi;
	do
	{
		phiN = (en->n - N0 - M) / (a * F0) + phiN;
		deltaPhi = phiN - PHI0;
		sumPhi = phiN + PHI0;
		M = b * F0 * (((1.0) + n + ((5.0) / (4.0)) * n2 + ((5.0) / (4.0)) * n3) * deltaPhi
			- ((3.0) * n + (3.0) * n2 + ((21.0) / (8.0)) * n3) * sin(deltaPhi) * cos(sumPhi)
			+ ((15.0 / 8.0) * n2 + (15.0 / 8.0) * n3) * sin(2.0 * deltaPhi) * cos((2.0) * sumPhi)
			- ((35.0) / (24.0)) * n3 * sin(3.0 * deltaPhi) * cos(3.0 * sumPhi)
			);

	} while (fabs(en->n - N0 - M) >= 0.00001);
	
	double v = a * F0 * pow(1 - e2 * pow(sin(phiN), 2), -0.5); //B3
	double rho = a * F0 * (1 - e2) * pow(1 - e2 * pow(sin(phiN), 2), -1.5); //B4
	double eta2 = v / rho - 1; //B5

	double VII = tan(phiN) / (2.0 * rho * v);
	double VIII = tan(phiN) / (24.0 * rho * tri(v)) * (5 + 3 * dbl(tan(phiN)) + eta2 - 9 * dbl(tan(phiN)) * eta2);
	double IX = tan(phiN) / (24.0 * rho * tri(v) * dbl(v)) * (61 + 90 * dbl(tan(phiN)) + 45 * dbl(dbl(tan(phiN))));
	double X = 1.0 / (v * cos(phiN));
	double XI = 1.0 / (6 * cos(phiN) * tri(v)) * (v/rho + 2 * dbl(tan(phiN)));
	double XII = 1.0 / (120 * cos(phiN) * tri(v) * dbl(v)) * (5.0 + 28.0 * dbl(tan(phiN)) + 24 * dbl(dbl(tan(phiN))));
	double XIIA = 1.0 / (5040 * cos(phiN) * tri(v) * dbl(v) * dbl(v)) * (61.0 + 662.0 * dbl(tan(phiN)) + 1320 * dbl(dbl(tan(phiN))) + 720 * dbl(tri(tan(phiN))));

	double phi = phiN - VII * dbl(E-E0) + VIII*pow(E - E0,4) - IX * pow(E-E0, 6);
	double lambda = LAMBDA0 + X*(E-E0) - XI * tri(E-E0) + XII * pow(E - E0, 5) - XIIA * pow(E - E0, 7);

	LatLonDecimal *latLon = (LatLonDecimal*)malloc(sizeof(LatLonDecimal));
	latLon->lat = oneEightyOverPi * phi;
	latLon->lon = oneEightyOverPi * lambda;
	latLon->elevation = en->elevation;
	latLon->geoid = en->geoid;

	return latLon;
}

LatLonDecimal *ETRS89LatLonFromETRS89EastingNorthing(EastingNorthing *en)
{
	return latLonFromEastingNorthing(en, &GRS80Ellipsoid, &NationalGridProj);
}

EastingNorthing *OSTN02ShiftsForIndices(int eIndex, int nIndex)
{
	EastingNorthing *shifts = (EastingNorthing*)malloc(sizeof(EastingNorthing));
	shifts->e = shifts->n = shifts->elevation = shifts->geoid = 0;
	if (nIndex < 0 || nIndex > 1250) return shifts;

	const OSTN02Index dataIndex = OSTN02Indices[nIndex];
	if (eIndex < dataIndex.eMin || eIndex >= dataIndex.eMin + dataIndex.eCount) return shifts;

	const unsigned int dataOffset = dataIndex.offset + (eIndex - dataIndex.eMin);
	const OSTN02Datum record = OSTN02Data[dataOffset];
	if (record.gFlag == 0) return shifts;

	shifts->e = (((double)record.eShift) / (1000.0)) + (86.0);
	shifts->n = (((double)record.nShift) / (1000.0)) - (82.0);
	shifts->elevation = (((double)record.gShift) / (1000.0)) + (43.0);
	shifts->geoid = record.gFlag;
	return shifts;
}

EastingNorthing *shiftsForEastingNorthing(EastingNorthing *en)
{
	EastingNorthing *shifts = (EastingNorthing*)malloc(sizeof(EastingNorthing));
	shifts->e = shifts->n = shifts->elevation = shifts->geoid = 0;

	const int e0 = (int)(en->e / (1000.0));
	const int n0 = (int)(en->n / (1000.0));
	double      dx = en->e - (double)(e0 * 1000);
	double      dy = en->n - (double)(n0 * 1000);
	double      t = dx / (1000.0);
	double      u = dy / (1000.0);

	EastingNorthing *shifts0 = OSTN02ShiftsForIndices(e0, n0);
	if (shifts0->geoid == 0) return shifts;

	EastingNorthing *shifts1 = OSTN02ShiftsForIndices(e0 + 1, n0);
	if (shifts1->geoid == 0) return shifts;

	EastingNorthing *shifts2 = OSTN02ShiftsForIndices(e0 + 1, n0 + 1);
	if (shifts2->geoid == 0) return shifts;

	EastingNorthing *shifts3 = OSTN02ShiftsForIndices(e0, n0 + 1);
	if (shifts3->geoid == 0) return shifts;

	double weight0 = ((1.0) - t) * ((1.0) - u);
	double weight1 = t * ((1.0) - u);
	double weight2 = t * u;
	double weight3 = ((1.0) - t) *           u;

	shifts->e = weight0 * shifts0->e
		+ weight1 * shifts1->e
		+ weight2 * shifts2->e
		+ weight3 * shifts3->e;

	shifts->n = weight0 * shifts0->n
		+ weight1 * shifts1->n
		+ weight2 * shifts2->n
		+ weight3 * shifts3->n;

	shifts->elevation = weight0 * shifts0->elevation
		+ weight1 * shifts1->elevation
		+ weight2 * shifts2->elevation
		+ weight3 * shifts3->elevation;

	bool left = dx < (500.0);
	bool bottom = dy < (500.0);
	EastingNorthing *nearestShift = left ? (bottom ? shifts0 : shifts3) : (bottom ? shifts1 : shifts2);
	shifts->geoid = nearestShift->geoid;

	return shifts;
}

EastingNorthing *OSGB36EastingNorthingFromETRS89EastingNorthing(EastingNorthing *en)
{
	EastingNorthing *shifts = shiftsForEastingNorthing(en);
	if (shifts->geoid == 0) return shifts;

	EastingNorthing *shifted = (EastingNorthing*)malloc(sizeof(EastingNorthing));
	shifted->e = en->e + shifts->e;
	shifted->n = en->n + shifts->n;
	shifted->elevation = en->elevation - shifts->elevation;
	shifted->geoid = shifts->geoid;

	return shifted;
}

EastingNorthing *ETRS89EastingNorthingFromOSGB36EastingNorthing(EastingNorthing *en)
{
	EastingNorthing 
		*shifts = (EastingNorthing *)malloc(sizeof(EastingNorthing)), 
		*prevShifts = (EastingNorthing *)malloc(sizeof(EastingNorthing)), 
		*shifted = (EastingNorthing *)malloc(sizeof(EastingNorthing));
	shifts->e = shifts->n = shifted->elevation = shifted->geoid = 0;  // initialising .elevation and .geoid just avoids warnings
	do
	{
		prevShifts->e = shifts->e;
		prevShifts->n = shifts->n;
		shifted->e = en->e - shifts->e;
		shifted->n = en->n - shifts->n;
		shifts = shiftsForEastingNorthing(shifted);
		if (shifts->geoid == 0) return shifts;
	} while (abs(shifts->e - prevShifts->e) > (0.0001) || abs(shifts->n - prevShifts->n) > (0.0001));

	shifted->elevation = en->elevation + shifts->elevation;
	shifted->geoid = shifts->geoid;  // tells us which geoid datum was used in conversion
	return shifted;
}


DegMinSec degMinSecFromDecimal(double dec)
{
	DegMinSec dms;
	dms.westOrSouth = dec < (0.0) ? true : false;
	dec = abs(dec);
	dms.deg = (int)floor(dec);
	dec -= dms.deg;
	dms.min = (int)floor((60.0) * dec);
	dec -= dms.min / (60.0);
	dms.sec = dec * (3600.0);
	return dms;
}

double decimalFromDegMinSec(DegMinSec dms)
{
	return (dms.westOrSouth ? (-1.0) : (1.0)) * (((double)dms.deg) + ((double)dms.min) / (60.0) + dms.sec / (3600.0));
}

LatLonDecimal latLonDecimalFromLatLonDegMinSec(LatLonDegMinSec dms)
{
	LatLonDecimal dec;
	dec.lat = decimalFromDegMinSec(dms.lat);
	dec.lon = decimalFromDegMinSec(dms.lon);
	dec.elevation = dms.elevation;
	dec.geoid = dms.geoid;
	return dec;
}

LatLonDegMinSec latLonDegMinSecFromLatLonDecimal(LatLonDecimal dec)
{
	LatLonDegMinSec dms;
	dms.lat = degMinSecFromDecimal(dec.lat);
	dms.lon = degMinSecFromDecimal(dec.lon);
	dms.elevation = dec.elevation;
	dms.geoid = dec.geoid;
	return dms;
}

bool test(bool noisily)
{
	short numTested = 0;
	short numPassed = 0;
	bool  testPassed;

	// check data integrity

	const unsigned long dataCRC = crc32((unsigned char *)OSTN02Data, sizeof(OSTN02Data));
	testPassed = dataCRC == originalDataCRC;
	if (noisily) {
		mexPrintf("\nOriginal CRC32 (data):  %li\n", originalDataCRC);
		mexPrintf("Computed CRC32 (data):  %li, %s\n\n", dataCRC, (testPassed ? "Passed" : "Not Passed"));
	}

	const unsigned long indicesCRC = crc32((unsigned char *)OSTN02Indices, sizeof(OSTN02Indices));
	testPassed = indicesCRC == originalIndicesCRC;
	if (noisily) {
		mexPrintf("Original CRC32 (index): %li\n", originalIndicesCRC);
		mexPrintf("Computed CRC32 (index): %li, %s\n\n", indicesCRC, (testPassed ? "Passed" : "Not Passed"));
	}

	// check test conversions against known good results

	int len;

	// OSTN02 conversions

	LatLonDecimal actualLatLon, *computedLatLon;
	EastingNorthing actualEN, *computedEN;
	char *actualLatLonStr = (char *)malloc(500), *computedLatLonStr = (char *)malloc(500), *actualENStr = (char *)malloc(500), *computedENStr = (char *)malloc(500);

	len = LENGTH_OF(testETRSCoords);
	for (int i = 0; i < len; i++)
	{

		// actual coords

		actualLatLon = latLonDecimalFromLatLonDegMinSec(testETRSCoords[i]);
		actualEN = testOSGB36Coords[i];

		sprintf(actualLatLonStr, "%.6lf, %.6lf: %.6lf", actualLatLon.lat, actualLatLon.lon, actualLatLon.elevation);
		sprintf(actualENStr, "%.6lf, %.6lf: %.6lf at %s - %s", actualEN.e, actualEN.n, actualEN.elevation, OSGB36GeoidRegions[actualEN.geoid], OSGB36GeoidNames[actualEN.geoid]);

		if (noisily) {
			mexPrintf("ETRS89 actual   %s\n", actualLatLonStr);
			mexPrintf("OSGB36 actual   %s\n", actualENStr);
		}

		// computed coords

		computedLatLon = ETRS89LatLonFromETRS89EastingNorthing(ETRS89EastingNorthingFromOSGB36EastingNorthing(&actualEN));
		computedEN = OSGB36EastingNorthingFromETRS89EastingNorthing(ETRS89EastingNorthingFromETRS89LatLon(&actualLatLon));

		sprintf(computedLatLonStr, "%.6lf, %.6lf: %.6lf", computedLatLon->lat, computedLatLon->lon, computedLatLon->elevation);
		sprintf(computedENStr, "%.6lf, %.6lf: %.6lf at %s - %s", computedEN->e, computedEN->n, computedEN->elevation, OSGB36GeoidRegions[computedEN->geoid], OSGB36GeoidNames[computedEN->geoid]);

		if (actualEN.geoid != 0) {  // i.e. these coordinates are not zeroes, and can be converted sensibly
			testPassed = strcmp(actualLatLonStr, computedLatLonStr) == 0;
			numTested++;
			if (testPassed) numPassed++;
			if (noisily) mexPrintf("ETRS89 computed %s, %s\n", computedLatLonStr, (testPassed ? "Passed" : "Not Passed"));
		}

		testPassed = strcmp(actualENStr, computedENStr) == 0;
		numTested++;
		if (testPassed) numPassed++;
		if (noisily) mexPrintf("OSGB36 computed %s, %s\n\n", computedENStr, (testPassed ? "Passed" : "Not Passed"));

		//free(actualLatLonStr);
		//free(computedLatLonStr);
		//free(actualENStr);
		//free(computedENStr);
	}

	bool allPassed = numTested == numPassed;
	if (noisily) mexPrintf("%i tests; %i passed; %i failed\n\n", numTested, numPassed, numTested - numPassed);
	return allPassed;
}

LIDARTerrainFile *findFile(char *BLOCK, int fileE, int fileN, LIDARTerrainFile *ter)
{
	char *name = (char *)malloc(500);
	//mexPrintf("\n        Loading file .\\LIDAR\\%c%c%02d%02d_DSM_1M.asc\n", BLOCK[0], BLOCK[1], fileE, fileN);
	sprintf(name, ".\\LIDAR\\%c%c%02d%02d_DSM_1M.asc", BLOCK[0] - 'A' + 'a', BLOCK[1] - 'A' + 'a', fileE, fileN);
	LIDARTerrainFile *tmp = ter;
	FILE *thisFile;
	if (!(thisFile = fopen(name, "r")))
	{
		char *msg[100];
		sprintf(msg, "Error missing LIDAR Elevation file %c%c %d%d\n", BLOCK[0], BLOCK[1], fileE, fileN);
		mexWarnMsgIdAndTxt("MATLAB:OSTN15_Matlab:LackingLIDARElevationFile",
			msg);
		return NULL;
	}
	ter->next = (LIDARTerrainFile*)malloc(sizeof(LIDARTerrainFile));
	tmp = ter->next;
	tmp->next = NULL;
	strncpy(tmp->Name, BLOCK, 2);
	int sizeE, sizeN, LLX, LLY, cellsize, neg;
	tmp->Xblock = fileE;
	tmp->Yblock = fileN;
	char temp[255];
	fgets(temp, 255, thisFile);
	fgets(temp, 255, thisFile);
	fgets(temp, 255, thisFile);
	fgets(temp, 255, thisFile);
	fgets(temp, 255, thisFile);
	fgets(temp, 255, thisFile);
	double **res = (double **)malloc(1001 * sizeof(double*)); 
	for (int i = 999; i>=0; i--)
	{
		res[i] = (double *)malloc(1001 * sizeof(double));
		for (int j = 0; j < 1000; j++)
		{
			fscanf(thisFile, "%lf", &res[i][j]);
		}
	}
	tmp->data = res;
	return tmp;
}

double findHeight(EastingNorthing *en, char *blk, LIDARTerrainFile *ter)
{
	int E = lrint(en->e), N = lrint(en->n);
	int fileE = (E % 100000) / 1000, fileN = (N % 100000) / 1000;
	int offsetE = E % 1000, offsetN = N % 1000;
	if (strcmp(blk, blkNow)!=0 && blke == (fileE / 10) && blkn == (fileN / 10) && ter->next)
	{
		LIDARTerrainFile *temp = ter->next, *temp2 = temp;
		while (temp)
		{
			temp = temp2->next;
			free(temp2->data);
			free(temp2);
			temp2 = temp;
		}
	}
	strncpy(blkNow, blk, 2);
	char BLOCK[2];
	strncpy(BLOCK, blk, 2);
	LIDARTerrainFile *current = ter->next, *last = ter;
	while (current)
	{
		if (!(current->Xblock == fileE && current->Yblock == fileN && !strncmp(BLOCK, current->Name, 2)))
		{
			current = current->next;
			last = last->next;
		}
		else break;
	}
	if (!current) current = findFile(BLOCK, fileE, fileN, last);
	if (!current) return -9999;
	//mexPrintf("%s%d%d=%s%d%d", current->Name, current->Xblock, current->Yblock, BLOCK, fileE, fileN);
	return (current->data)[offsetN][offsetE];
}

int checkGridPosition(EastingNorthing *en, int prevMap)
{  // start with prevMap = -1
	for (int i = prevMap + 1, len = LENGTH_OF(gridPosition); i < len; i++)
	{
		OSMap map = gridPosition[i];
		if (en->e >= (double)map.emin
			&& en->e <= (double)map.emax
			&& en->n >= (double)map.nmin
			&& en->n <= (double)map.nmax) return i;
	}
	return -1;
}