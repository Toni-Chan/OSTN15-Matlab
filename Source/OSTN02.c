//
//  OSTN02.c
//  OSTN02
//
//  Created by George MacKerron on 24/12/2011.
//  Copyright (c) 2011 George MacKerron. All rights reserved.
//

#include "OSTN02.h"
#include "dblRelated.h"
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
#include "D:\\Program Files\\MATLAB\\R2018a\\extern\\include\matrix.h"

#define LENGTH_OF(x) (sizeof (x) / sizeof *(x))
#define originalIndicesCRC 244629328  // these won't be robust against differing endianness or compiler packing of bit-structs
#define originalDataCRC    790474494

double piOver180 = 0.017453292519943295769236907684886127134428718885417254560971914401710091146034494436822415696345094822;
double oneEightyOverPi = 57.29577951308232087679815481410517033240547246656432154916024386120284714832155263244096899585111094418;


double gridConvergenceDegreesFromLatLon(LatLonDecimal *latLon, Ellipsoid *ellipsoid, MapProjection *projection)
{

	double phi = piOver180 * latLon->lat;
	double lambda = piOver180 * latLon->lon;
	double lambda0 = piOver180 * projection->trueOriginLatLon.lon;

	double a = ellipsoid->semiMajorAxis;
	double b = ellipsoid->semiMinorAxis;
	double f0 = projection->centralMeridianScale;

	double deltaLambda = lambda - lambda0;
	double deltaLambda2 = deltaLambda  * deltaLambda;
	double deltaLambda3 = deltaLambda2 * deltaLambda;
	double deltaLambda5 = deltaLambda3 * deltaLambda2;

	double sinPhi = sin(phi);
	double sinPhi2 = sinPhi * sinPhi;
	double cosPhi = cos(phi);
	double cosPhi2 = cosPhi * cosPhi;
	double cosPhi4 = cosPhi2 * cosPhi2;
	double tanPhi = tan(phi);
	double tanPhi2 = tanPhi * tanPhi;

	double af0 = a * f0;
	double af02 = af0 * af0;
	double bf0 = b * f0;
	double bf02 = bf0 * bf0;

	double e2 = (af02 - bf02) / af02;
	double nu = af0 / sqrt(L(1.0) - (e2 * sinPhi2));
	double rho = (nu * (L(1.0) - e2)) / (L(1.0) - (e2 * sinPhi2));
	double eta2 = (nu / rho) - L(1.0);
	double xiv = ((sinPhi * cosPhi2) / L(3.0))  * (L(1.0) + L(3.0) * eta2 + L(2.0) * eta2 * eta2);
	double xv = ((sinPhi * cosPhi4) / L(15.0)) * (L(2.0) - tanPhi2);

	double cRads = (deltaLambda * sinPhi) + (deltaLambda3 * xiv) + (deltaLambda5 * xv);
	return oneEightyOverPi * cRads;
}

double gridConvergenceDegreesFromEastingNorthing(EastingNorthing *en, Ellipsoid *ellipsoid, MapProjection *projection) 
{
	double n0 = projection->trueOriginEastingNorthing.n;
	double e0 = projection->trueOriginEastingNorthing.e;
	double f0 = projection->centralMeridianScale;
	double af0 = ellipsoid->semiMajorAxis * f0;
	double bf0 = ellipsoid->semiMinorAxis * f0;

	double af02 = af0 * af0;
	double bf02 = bf0 * bf0;

	double n = (af0 - bf0) / (af0 + bf0);
	double n2 = n * n;
	double n3 = n2 * n;
	double e2 = (af02 - bf02) / af02;

	double phi0 = piOver180 * projection->trueOriginLatLon.lat;
	double phi = (en->n - n0) / af0 + phi0;

	while (true) 
	{
		double deltaPhi = phi - phi0;
		double sumPhi = phi + phi0;
		double j3 = (L(1.0) + n + L(5.0) / L(4.0) * n2 + L(5.0) / L(4.0) * n3) * deltaPhi;
		double j4 = (L(3.0) * n + L(3.0) * n2 + L(21.0) / L(8.0) * n3) * sin(deltaPhi) * cos(sumPhi);
		double j5 = (L(15.0) / L(8.0) * n2 + L(15.0) / L(8.0) * n3) * sin(L(2.0) * deltaPhi) * cos(L(2.0) * sumPhi);
		double j6 = L(35.0) / L(24.0) * n3 * sin(L(3.0) * deltaPhi) * cos(L(3.0) * sumPhi);
		double M = bf0 * (j3 - j4 + j5 - j6);
		if (ABS(en->n - n0 - M) < L(0.001)) break;
		phi += (en->n - n0 - M) / af0;
	};

	double sinPhi = sin(phi);
	double sinPhi2 = sinPhi * sinPhi;
	double nu = af0 / sqrt(L(1.0) - e2 * sinPhi2);
	double rho = nu * (L(1.0) - e2) / (L(1.0) - e2 * sinPhi2);
	double eta2 = nu / rho - L(1.0);
	double eta4 = eta2 * eta2;

	double nu2 = nu * nu;
	double nu3 = nu2 * nu;
	double nu5 = nu3 * nu2;
	double tanPhi = tan(phi);
	double tanPhi2 = tanPhi * tanPhi;
	double tanPhi4 = tanPhi2 * tanPhi2;

	double y = en->e - e0;
	double y3 = y * y * y;
	double y5 = y3 * y * y;

	double j31 = tanPhi / nu;
	double j41 = tanPhi / (L(3.0) * nu3) * (L(1.0) + tanPhi2 - eta2 - L(2.0) * eta4);
	double j51 = tanPhi / (L(15.0) * nu5) * (L(2.0) + L(5.0) * tanPhi2 + L(3.0) * tanPhi4);

	double c = y * j31 - y3 * j41 + y5 * j51;

	return oneEightyOverPi * c;
}

double gridConvergenceDegreesFromOSGB36EastingNorthing(EastingNorthing *en) 
{
	return gridConvergenceDegreesFromEastingNorthing(en, &Airy1830Ellipsoid, &NationalGridProj);
}

double gridConvergenceDegreesFromETRS89LatLon(LatLonDecimal *latLon) 
{
	return gridConvergenceDegreesFromLatLon(latLon, &GRS80Ellipsoid, &NationalGridProj);
}

int nextOSExplorerMap(EastingNorthing *en, int prevMap)
{  // start with prevMap = -1
	for (int i = prevMap + 1, len = LENGTH_OF(OSExplorerMaps); i < len; i++)
	{
		OSMap map = OSExplorerMaps[i];
		if (en->e >= (double)map.emin
			&& en->e <= (double)map.emax
			&& en->n >= (double)map.nmin
			&& en->n <= (double)map.nmax) return i;
	}
	return -1;
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

EastingNorthing *eastingNorthingFromLatLon(LatLonDecimal *latLon, Ellipsoid *ellipsoid, MapProjection *projection) {
	double a = ellipsoid->semiMajorAxis;
	double b = ellipsoid->semiMinorAxis;
	double f0 = projection->centralMeridianScale;
	double e0 = projection->trueOriginEastingNorthing.e;
	double n0 = projection->trueOriginEastingNorthing.n;

	double phi = piOver180 * latLon->lat;
	double lambda = piOver180 * latLon->lon;
	double phi0 = piOver180 * projection->trueOriginLatLon.lat;
	double lambda0 = piOver180 * projection->trueOriginLatLon.lon;

	double deltaPhi = phi - phi0;
	double sumPhi = phi + phi0;
	double sinPhi = sin(phi);
	double sinPhi2 = sinPhi * sinPhi;
	double cosPhi = cos(phi);
	double cosPhi2 = cosPhi * cosPhi;
	double cosPhi3 = cosPhi2 * cosPhi;
	double cosPhi5 = cosPhi3 * cosPhi2;
	double tanPhi = tan(phi);
	double tanPhi2 = tanPhi * tanPhi;
	double tanPhi4 = tanPhi2 * tanPhi2;

	double a2 = a * a;
	double e2 = (a2 - b * b) / a2;
	double n = (a - b) / (a + b);
	double n2 = n * n;
	double n3 = n2 * n;
	double oneMinusE2SinPhi2 = L(1.0) - e2 * sinPhi2;
	double sqrtOneMinusE2SinPhi2 = sqrt(oneMinusE2SinPhi2);
	double v = a * f0 / sqrtOneMinusE2SinPhi2;
	double rho = a * f0 * (L(1.0) - e2) / (oneMinusE2SinPhi2 * sqrtOneMinusE2SinPhi2);
	double eta2 = v / rho - L(1.0);
	double m = b * f0 * ((L(1.0) + n + (L(5.0) / L(4.0)) * n2 + (L(5.0) / L(4.0)) * n3) * deltaPhi
		- (L(3.0) * n + L(3.0) * n2 + (L(21.0) / L(8.0)) * n3) * sin(deltaPhi) * cos(sumPhi)
		+ ((L(15.0) / L(8.0)) * n2 + (L(15.0) / L(8.0)) * n3) * sin(L(2.0) * deltaPhi) * cos(L(2.0) * sumPhi)
		- (L(35.0) / L(24.0)) * n3 * sin(L(3.0) * deltaPhi) * cos(L(3.0) * sumPhi)
		);

	double one = m + n0;
	double two = (v / L(2.0)) * sinPhi * cosPhi;
	double three = (v / L(24.0)) * sinPhi * cosPhi3 * (L(5.0) - tanPhi2 + L(9.0) * eta2);
	double threeA = (v / L(720.0)) * sinPhi * cosPhi5 * (L(61.0) - L(58.0) * tanPhi2 + tanPhi4);
	double four = v * cosPhi;
	double five = (v / L(6.0)) * cosPhi3 * (v / rho - tanPhi2);
	double six = (v / L(120.0)) * cosPhi5 * (L(5.0) - L(18.0) * tanPhi2 + tanPhi4 + L(14.0) * eta2 - L(58.0) * tanPhi2 * eta2);

	double deltaLambda = lambda - lambda0;
	double deltaLambda2 = deltaLambda  * deltaLambda;
	double deltaLambda3 = deltaLambda2 * deltaLambda;
	double deltaLambda4 = deltaLambda3 * deltaLambda;
	double deltaLambda5 = deltaLambda4 * deltaLambda;
	double deltaLambda6 = deltaLambda5 * deltaLambda;

	EastingNorthing *en = (EastingNorthing*)malloc(sizeof(EastingNorthing));
	en->n = one + two * deltaLambda2 + three * deltaLambda4 + threeA * deltaLambda6;
	en->e = e0 + four * deltaLambda + five * deltaLambda3 + six * deltaLambda5;
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
	double a = ellipsoid->semiMajorAxis;
	double b = ellipsoid->semiMinorAxis;
	double f0 = projection->centralMeridianScale;
	double e0 = projection->trueOriginEastingNorthing.e;
	double n0 = projection->trueOriginEastingNorthing.n;
	double n = (a - b) / (a + b);
	double n2 = n * n;
	double n3 = n2 * n;
	double bf0 = b * f0;

	double phi0 = piOver180 * projection->trueOriginLatLon.lat;
	double lambda0 = piOver180 * projection->trueOriginLatLon.lon;

	double phi = phi0;  // this is phi' in the OS docs
	double m = L(0.0);
	double deltaPhi, sumPhi;
	do
	{
		phi = (en->n - n0 - m) / (a * f0) + phi;
		deltaPhi = phi - phi0;
		sumPhi = phi + phi0;
		m = bf0 * ((L(1.0) + n + (L(5.0) / L(4.0)) * n2 + (L(5.0) / L(4.0)) * n3) * deltaPhi
			- (L(3.0) * n + L(3.0) * n2 + (L(21.0) / L(8.0)) * n3) * sin(deltaPhi) * cos(sumPhi)
			+ ((L(15.0) / L(8.0)) * n2 + (L(15.0) / L(8.0)) * n3) * sin(L(2.0) * deltaPhi) * cos(L(2.0) * sumPhi)
			- (L(35.0) / L(24.0)) * n3 * sin(L(3.0) * deltaPhi) * cos(L(3.0) * sumPhi)
			);
	} while (ABS(en->n - n0 - m) >= L(0.00001));

	double sinPhi = sin(phi);
	double sinPhi2 = sinPhi * sinPhi;
	double cosPhi = cos(phi);
	double secPhi = L(1.0) / cosPhi;
	double tanPhi = tan(phi);
	double tanPhi2 = tanPhi * tanPhi;
	double tanPhi4 = tanPhi2 * tanPhi2;
	double tanPhi6 = tanPhi4 * tanPhi2;

	double a2 = a * a;
	double e2 = (a2 - b * b) / a2;
	double oneMinusE2SinPhi2 = L(1.0) - e2 * sinPhi2;
	double sqrtOneMinusE2SinPhi2 = sqrt(oneMinusE2SinPhi2);
	double v = a * f0 / sqrtOneMinusE2SinPhi2;
	double rho = a * f0 * (L(1.0) - e2) / (oneMinusE2SinPhi2 * sqrtOneMinusE2SinPhi2);
	double eta2 = v / rho - L(1.0);

	double v2 = v  * v;
	double v3 = v2 * v;
	double v5 = v3 * v2;
	double v7 = v5 * v2;

	double seven = tanPhi / (L(2.0) * rho * v);
	double eight = (tanPhi * (L(5.0) + L(3.0) * tanPhi2 + eta2 - L(9.0) * tanPhi2 * eta2)) / (L(24.0) * rho * v3);
	double nine = (tanPhi * (L(61.0) + L(90.0) * tanPhi2 + L(45.0) * tanPhi4)) / (L(720.0) * rho * v5);
	double ten = secPhi / v;
	double eleven = (secPhi * ((v / rho) + L(2.0) * tanPhi2)) / (L(6.0) * v3);
	double twelve = (secPhi * (L(5.0) + L(28.0) * tanPhi2 + L(24.0) * tanPhi4)) / (L(120.0) * v5);
	double twelveA = (secPhi * (L(61.0) + L(662.0) * tanPhi2 + L(1320.0) * tanPhi4 + L(720.0) * tanPhi6)) / (L(5040.0) * v7);

	double deltaE = en->e - e0;
	double deltaE2 = deltaE  * deltaE;
	double deltaE3 = deltaE2 * deltaE;
	double deltaE4 = deltaE2 * deltaE2;
	double deltaE5 = deltaE3 * deltaE2;
	double deltaE6 = deltaE3 * deltaE3;
	double deltaE7 = deltaE4 * deltaE3;

	LatLonDecimal *latLon = (LatLonDecimal*)malloc(sizeof(LatLonDecimal));
	latLon->lat = oneEightyOverPi * (phi - seven * deltaE2 + eight * deltaE4 - nine * deltaE6);
	latLon->lon = oneEightyOverPi * (lambda0 + ten * deltaE - eleven * deltaE3 + twelve * deltaE5 - twelveA * deltaE7);
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

	shifts->e = (((double)record.eShift) / L(1000.0)) + L(86.0);
	shifts->n = (((double)record.nShift) / L(1000.0)) - L(82.0);
	shifts->elevation = (((double)record.gShift) / L(1000.0)) + L(43.0);
	shifts->geoid = record.gFlag;
	return shifts;
}

EastingNorthing *shiftsForEastingNorthing(EastingNorthing *en)
{
	EastingNorthing *shifts = (EastingNorthing*)malloc(sizeof(EastingNorthing));
	shifts->e = shifts->n = shifts->elevation = shifts->geoid = 0;

	const int e0 = (int)(en->e / L(1000.0));
	const int n0 = (int)(en->n / L(1000.0));
	double      dx = en->e - (double)(e0 * 1000);
	double      dy = en->n - (double)(n0 * 1000);
	double      t = dx / L(1000.0);
	double      u = dy / L(1000.0);

	EastingNorthing *shifts0 = OSTN02ShiftsForIndices(e0, n0);
	if (shifts0->geoid == 0) return shifts;

	EastingNorthing *shifts1 = OSTN02ShiftsForIndices(e0 + 1, n0);
	if (shifts1->geoid == 0) return shifts;

	EastingNorthing *shifts2 = OSTN02ShiftsForIndices(e0 + 1, n0 + 1);
	if (shifts2->geoid == 0) return shifts;

	EastingNorthing *shifts3 = OSTN02ShiftsForIndices(e0, n0 + 1);
	if (shifts3->geoid == 0) return shifts;

	double weight0 = (L(1.0) - t) * (L(1.0) - u);
	double weight1 = t  * (L(1.0) - u);
	double weight2 = t  *           u;
	double weight3 = (L(1.0) - t) *           u;

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

	bool left = dx < L(500.0);
	bool bottom = dy < L(500.0);
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
	EastingNorthing *shifts = (EastingNorthing *)malloc(sizeof(EastingNorthing)), *prevShifts = (EastingNorthing *)malloc(sizeof(EastingNorthing)), *shifted = (EastingNorthing *)malloc(sizeof(EastingNorthing));
	shifts->e = shifts->n = shifted->elevation = shifted->geoid = 0;  // initialising .elevation and .geoid just avoids warnings
	do
	{
		prevShifts->e = shifts->e;
		prevShifts->n = shifts->n;
		shifted->e = en->e - shifts->e;
		shifted->n = en->n - shifts->n;
		shifts = shiftsForEastingNorthing(shifted);
		if (shifts->geoid == 0) return shifts;
	} while (ABS(shifts->e - prevShifts->e) > L(0.0001) || ABS(shifts->n - prevShifts->n) > L(0.0001));

	shifted->elevation = en->elevation + shifts->elevation;
	shifted->geoid = shifts->geoid;  // tells us which geoid datum was used in conversion
	return shifted;
}

char *tetradFromOSGB36EastingNorthing(EastingNorthing *en) 
{
	// see http://www.bto.org/volunteer-surveys/birdatlas/taking-part/correct-grid-references/know-your-place
	const int  eTrunc = (int)en->e;  // note: unlike for a grid ref, we never round -- we (implictly) truncate
	const int  nTrunc = (int)en->n;
	const int  firstEIndex = eTrunc / 500000;
	const int  firstNIndex = nTrunc / 500000;
	const int  secondEIndex = (eTrunc % 500000) / 100000;
	const int  secondNIndex = (nTrunc % 500000) / 100000;
	const char sq0 = firstLetters[firstNIndex][firstEIndex];
	const char sq1 = secondLetters[secondNIndex][secondEIndex];
	const int  eMinusSq = eTrunc - (500000 * firstEIndex) - (100000 * secondEIndex);
	const int  nMinusSq = nTrunc - (500000 * firstNIndex) - (100000 * secondNIndex);
	const int  eDigit = eMinusSq / 10000;
	const int  nDigit = nMinusSq / 10000;
	const int  tetradEIndex = (eMinusSq % 10000) / 2000;
	const int  tetradNIndex = (nMinusSq % 10000) / 2000;
	char tetradLetter = tetradLetters[tetradNIndex][tetradEIndex];
	char *tetrad = (char *)malloc(2560);
	sprintf(tetrad, "%c%c%d%d%c", sq0, sq1, eDigit, nDigit, tetradLetter);
	return tetrad;
}

DegMinSec degMinSecFromDecimal(double dec) 
{
	DegMinSec dms;
	dms.westOrSouth = dec < L(0.0) ? true : false;
	dec = ABS(dec);
	dms.deg = (int)FLOOR(dec);
	dec -= dms.deg;
	dms.min = (int)FLOOR(L(60.0) * dec);
	dec -= dms.min / L(60.0);
	dms.sec = dec * L(3600.0);
	return dms;
}

double decimalFromDegMinSec(DegMinSec dms) 
{
	return (dms.westOrSouth ? L(-1.0) : L(1.0)) * (((double)dms.deg) + ((double)dms.min) / L(60.0) + dms.sec / L(3600.0));
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
	FILE *out = fopen("checkRecord.txt","w");

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
	char *actualLatLonStr=(char *)malloc(500), *computedLatLonStr = (char *)malloc(500), *actualENStr = (char *)malloc(500), *computedENStr = (char *)malloc(500);

	len = LENGTH_OF(testETRSCoords);
	for (int i = 0; i < len; i++) 
	{

		// actual coords

		actualLatLon = latLonDecimalFromLatLonDegMinSec(testETRSCoords[i]);
		actualEN = testOSGB36Coords[i];

		sprintf(actualLatLonStr, LLFMT, actualLatLon.lat, actualLatLon.lon, actualLatLon.elevation);
		sprintf(actualENStr, ENFMT, actualEN.e, actualEN.n, actualEN.elevation, OSGB36GeoidRegions[actualEN.geoid], OSGB36GeoidNames[actualEN.geoid]);

		if (noisily) {
			mexPrintf("ETRS89 actual   %s\n", actualLatLonStr);
			mexPrintf("OSGB36 actual   %s\n", actualENStr);
		}

		// computed coords

		computedLatLon = ETRS89LatLonFromETRS89EastingNorthing(ETRS89EastingNorthingFromOSGB36EastingNorthing(&actualEN));
		computedEN = OSGB36EastingNorthingFromETRS89EastingNorthing(ETRS89EastingNorthingFromETRS89LatLon(&actualLatLon));

		sprintf(computedLatLonStr, LLFMT, computedLatLon->lat, computedLatLon->lon, computedLatLon->elevation);
		sprintf(computedENStr, ENFMT, computedEN->e, computedEN->n, computedEN->elevation, OSGB36GeoidRegions[computedEN->geoid], OSGB36GeoidNames[computedEN->geoid]);

		if (actualEN.geoid != 0) {  // i.e. these coordinates are not zeroes, and can be converted sensibly
			testPassed = strcmp(actualLatLonStr, computedLatLonStr) == 0;
			numTested++;
			if (testPassed) numPassed++;
			if (noisily) mexPrintf("ETRS89 computed %s, %s\n", computedLatLonStr, (testPassed ? "Passed" : "Not Passed"));
		}

		testPassed = strcmp(actualENStr, computedENStr) == 0;
		numTested++;
		if (testPassed) numPassed++;
		if (noisily) mexPrintf("OSGB36 computed %s, %s\n\n",  computedENStr, (testPassed ? "Passed" : "Not Passed"));

		//free(actualLatLonStr);
		//free(computedLatLonStr);
		//free(actualENStr);
		//free(computedENStr);
	}

	// convergences by Easting/Northing

	EastingNorthing convergenceEN;
	DegMinSec actualC, computedC;
	char *actualCStr = (char *)malloc(500), *computedCStr = (char *)malloc(500);

	len = LENGTH_OF(testConvergenceOSGB36Coords);
	for (int i = 0; i < len; i++) {

		// actual coords

		convergenceEN = testConvergenceOSGB36Coords[i];
		actualC = testConvergencesFromOSGB36Coords[i];
		sprintf(actualCStr, "%ideg  %02imin %.4fsec %c", actualC.deg, actualC.min, actualC.sec, actualC.westOrSouth ? 'W' : 'E');
		if (noisily) mexPrintf("Convergence reference from E/N  %s\n", actualCStr);

		// computed coords

		computedC = degMinSecFromDecimal(gridConvergenceDegreesFromOSGB36EastingNorthing(&convergenceEN));
		sprintf(computedCStr, "%ideg  %02imin %.4fsec %c", computedC.deg, computedC.min, computedC.sec, computedC.westOrSouth ? 'W' : 'E');

		testPassed = strcmp(actualCStr, computedCStr) == 0;
		numTested++;
		if (testPassed) numPassed++;
		if (noisily) mexPrintf("Convergence computed from E/N   %s, %s\n\n",  computedCStr, (testPassed ? "Passed" : "Not Passed"));

		//free(actualCStr);
		//free(computedCStr);
	}

	// convergences by lat/lon
	//jjj

	LatLonDecimal convergenceLatLon;

	len = LENGTH_OF(testConvergenceLatLons);
	for (int i = 0; i < len; i++) 
	{

		// actual coords

		convergenceLatLon = latLonDecimalFromLatLonDegMinSec(testConvergenceLatLons[i]);
		actualC = testConvergencesFromLatLons[i];
		sprintf(actualCStr, "%ideg  %02imin %.4fsec %c", actualC.deg, actualC.min, actualC.sec, actualC.westOrSouth ? 'W' : 'E');
		if (noisily) mexPrintf("Convergence reference from lat/lon  %s\n", actualCStr);

		// computed coords

		computedC = degMinSecFromDecimal(gridConvergenceDegreesFromLatLon(&convergenceLatLon, &Airy1830Ellipsoid, &NationalGridProj));
		sprintf(computedCStr, "%ideg  %02imin %.4fsec %c", computedC.deg, computedC.min, computedC.sec, computedC.westOrSouth ? 'W' : 'E');

		testPassed = strcmp(actualCStr, computedCStr) == 0;
		numTested++;
		if (testPassed) numPassed++;
		if (noisily) mexPrintf("Convergence computed from lat/lon   %s ,%s\n",  computedCStr, (testPassed ? "Passed" : "Not Passed"));

		//free(actualCStr);
		//free(computedCStr);
	}

	bool allPassed = numTested == numPassed;
	if (noisily) mexPrintf("%i tests; %i passed; %i failed\n\n", numTested, numPassed, numTested - numPassed);
	return allPassed;
}

LIDARTerrainFile *FindFile(char *BLOCK, int fileE, int fileN, LIDARTerrainFile *ter)
{
	char *name = (char *)malloc(500);
	sprintf(name, "U:\\Cambridge\\LIDAR\\%s%02d%02d_DSM_1M.asc", BLOCK, fileE, fileN);
	FILE *thisFile;
	if (!(thisFile = fopen(name, "r")))
	{
		thisFile = getMissingMap(BLOCK, fileE, fileN);
	}
	if (!ter) ter = (LIDARTerrainFile*)malloc(sizeof(LIDARTerrainFile));
	else
	{
		ter->next = (LIDARTerrainFile*)malloc(sizeof(LIDARTerrainFile));
		ter = ter->next;
	}
	ter->Name = BLOCK;
	char temp[25];
	int sizeE, sizeN, LLX, LLY, cellsize, neg;
	scanf("%s%d%s%d%s%d%s%d%s%d%s%d", temp, &sizeE, temp, &sizeN, temp, &LLX, temp, &LLY, temp, &cellsize, temp, &neg);
	ter->Xblock = LLX;
	ter->Yblock = LLY;
	double **res = (double **)malloc(1001*sizeof(double*));
	for (int i = sizeN - 1; i >= 0; i--)
	{
		res[i] = (double *)malloc(1001 * sizeof(double));
			for (int j = 0; j < sizeE; j++)
				scanf("%lf", &res[i][j]);
	}
	return ter;
}

double findHeight(EastingNorthing *en, char *blk, LIDARTerrainFile *ter)
{
	char BLOCK[2];
	strncpy(BLOCK, blk, 2);
	int E = lrint(en->e), N = lrint(en->n);
	int fileE = E / 1000, fileN = N / 1000;
	int offsetE = E % 1000, offsetN = N % 1000;
	
}