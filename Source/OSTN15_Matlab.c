//
//  main.c
//  OSTN15
//
//  Created by George MacKerron on 24/12/2011.
//  Copyright (c) 2011 George MacKerron. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "OSTN15.h"
#include "fancyOut.h"
#include "geoids.data"
#include "explorerMaps.data"
#include "gridPosition.data"

#define LENGTH_OF(x) (sizeof (x) / sizeof *(x))
#include <math.h>
#include "D:\\Program Files\\MATLAB\\R2018a\\extern\\include\mex.h"
#include "D:\\Program Files\\MATLAB\\R2018a\\extern\\include\matrix.h"

/* Input Arguments */

#define TYP		prhs[0]
#define LAT_IN	prhs[1]
#define LNG_IN	prhs[2]
#define NOR_DEC_IN  prhs[2]
#define EAS_DEC_IN	prhs[1]
#define OS36BLK_IN  prhs[3]
#define OS36NOR_IN  prhs[2]
#define OS36EAS_IN	prhs[1]
#define Location	prhs[1]

/* Output Arguments */

#define	LAT_OUT	plhs[0]
#define LNG_OUT plhs[1]
#define OS36NOR_OUT plhs[1]
#define OS36EAS_OUT plhs[0]
#define OS36Height plhs[2]
#define LATLNGHeight plhs[2]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/*		mexWarnMsgIdAndTxt("MATLAB:yprime:divideByZero",
			"Division by zero!\n");*/

int checkCommand(mxArray *input, char *words)
{
	char *input_buf = mxArrayToString(input);
	return strcmp(input_buf, words);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[])
{
	/* Check for proper number of arguments */
	LIDARTerrainFile *LTFile = (LIDARTerrainFile*)malloc(sizeof(LIDARTerrainFile));
	LTFile->next = NULL;

	if (nrhs < 1)
	{
		mexErrMsgIdAndTxt("MATLAB:OSTN15_Matlab:invalidNumInputs",
			"input arguments required. test, Lat-Lng, OSGB36 Decimal, OSGB36 Grided");
	}
	else if (nlhs > 3)
	{
		mexErrMsgIdAndTxt("MATLAB:OSTN15_Matlab:maxlhs",
			"Too many output arguments.");
	}

	if (nrhs == 1 && checkCommand(TYP, "test") == 0)
	{
		bool passed = test(true);
		return;
	}
	else if (nrhs == 1 && checkCommand(TYP, "list-geoids") == 0)
	{
		mexPrintf("\nFlag  Datum          Region");
		int len = LENGTH_OF(OSGB36GeoidNames);
		for (int i = 0; i < len; i++)
		mexPrintf("  %2d  %-13s  %s\n", i, OSGB36GeoidNames[i], OSGB36GeoidRegions[i]);
		mexPrintf("\n");
		return;
	}
	else if ((nrhs == 2 || nrhs == 3) && checkCommand(TYP, "gps-to-grid") == 0)
	{
		LatLonDecimal *latLon = (LatLonDecimal*)malloc(sizeof(LatLonDecimal));
		EastingNorthing *en;
				
		if (nrhs == 2)
		{
			char *filename;
			if (mxIsChar(Location) != 1)
				mexErrMsgIdAndTxt("MATLAB:OSTN15_Matlab:inputNotString",
					"Input file name must be a string.");
			filename = mxArrayToString(Location);
			mexPrintf("%s\n", filename);
			FILE *inputSequence = fopen(filename, "r");

			char *outname = (char*)malloc(60);
			strncpy(outname, filename, strlen(filename) - 4);
			strncpy(outname + strlen(filename) - 3, "\0", 1);
			strcat(outname, "-Output.csv\0");
			FILE *outputSequence = fopen(outname, "w");

			char *outname2 = (char*)malloc(60);
			strncpy(outname2, filename, strlen(filename) - 4);
			strncpy(outname2 + strlen(filename) - 3, "\0", 1);
			strcat(outname2, "-Raw.csv\0");
			FILE *outputSequence2 = fopen(outname2, "w");

			int scanResult;
			while ((scanResult = fscanf(inputSequence, "%lf,%lf", &latLon->lat, &latLon->lon) == 2))
			{
				//mexPrintf("Dealing with Lat %lf , Lon %lf ", latLon->lat, latLon->lon);
				latLon->elevation = 0;
				en = OSGB36EastingNorthingFromETRS89EastingNorthing(ETRS89EastingNorthingFromETRS89LatLon(latLon));
				char *txt = mxCalloc(2, sizeof(char)); //should be 2

				int mapIndex = checkGridPosition(en, -1);
				while (mapIndex >= 0)
				{
					OSMap map = gridPosition[mapIndex];

					sprintf(txt + strlen(txt), "%s", map.nameUTF8);
					mapIndex = checkGridPosition(en, mapIndex);
				};

				double NorRemnant = en->n - floor(en->n / 100000) * 100000;
				double EasRemnant = en->e - floor(en->e / 100000) * 100000;
				int EasFile = floor(EasRemnant) / 10000;
				int NorFile = floor(NorRemnant) / 10000;
				int EFileN = floor(EasRemnant - 10000 * EasFile) / 1000;
				int NFileN = floor(NorRemnant - 10000 * NorFile) / 1000;

				char fileT[2];
				if (EFileN >= 5) fileT[1] = 'E'; else fileT[1] = 'W';
				if (NFileN >= 5) fileT[0] = 'N'; else fileT[0] = 'S';

				//mexPrintf("as %s%01d%01d%s:%.6lf,%.6lf", txt, EasFile, NorFile, fileT, EasRemnant, NorRemnant);
				double height = findHeight(en, txt, LTFile);
				//if (height == -9999) continue;

				//mexPrintf("as %s%01d%01d%s:%.6lf,%.6lf\n", txt, EasFile, NorFile, fileT, EasRemnant, NorRemnant);
				//fprintf(outputSequence, "%s%01d%01d%s: %.6lf %.6lf\n", txt, EasFile, NorFile, fileT, EasRemnant, NorRemnant);
				//mexPrintf(", height %.2lf\n", height);
				fprintf(outputSequence, "%s%01d%01d%c%c,%.6lf,%.6lf,%.6lf\n", txt, EasFile, NorFile, fileT[0], fileT[1], EasRemnant, NorRemnant, height);
				fprintf(outputSequence2, "%.6lf,%.6lf,%.6lf\n", en->e, en->n, height);
			}
			fclose(outputSequence);
			fclose(inputSequence);
			fclose(outputSequence2);
			return;
		}
		else 
		{
			double *LatIN, *LngIN, *nDisplay, *eDisplay, *height;
			char *blkDisplay;
			
			latLon->lat = latLon->lon = latLon->elevation = 0;
			
			LatIN = mxGetPr(LAT_IN);
			LngIN = mxGetPr(LNG_IN);

			int sizeM = mxGetM(LAT_IN), sizeN = mxGetN(LAT_IN);
			int checkM = mxGetM(LNG_IN), checkN = mxGetN(LNG_IN);
			if (sizeM != checkM || sizeN != checkN)
			{
				mexErrMsgIdAndTxt("MATLAB:OSTN15_Matlab:InputMismatch",
					"Lat and Long input rows does not match.");
				return;
			}

			OS36NOR_OUT = mxCreateDoubleMatrix(sizeM, sizeN, mxREAL);
			OS36EAS_OUT = mxCreateDoubleMatrix(sizeM, sizeN, mxREAL);
			OS36Height = mxCreateDoubleMatrix(sizeM, sizeN, mxREAL);

			nDisplay = mxGetPr(OS36NOR_OUT);
			eDisplay = mxGetPr(OS36EAS_OUT);
			height = mxGetPr(OS36Height);

			for (int i = 0; i < sizeM*sizeN; i++)
			{
				latLon->lat = LatIN[i];
				latLon->lon = LngIN[i];

				en = OSGB36EastingNorthingFromETRS89EastingNorthing(ETRS89EastingNorthingFromETRS89LatLon(latLon));
				if (en->geoid == 0)
				{
					mexErrMsgIdAndTxt("MATLAB:OSTN15_Matlab:OutofOSGB36",
						"Input Position out of OSGB36 Range.");
				}

				char *txt = mxCalloc(5, sizeof(char));

				int mapIndex = checkGridPosition(en, -1);
				while (mapIndex >= 0)
				{
					OSMap map = gridPosition[mapIndex];

					//mexPrintf("%d, %s: %s", map.num, map.nameUTF8, map.sheetUTF8);
					sprintf(txt + strlen(txt), "%s", map.nameUTF8);
					if (strlen(map.sheetUTF8) > 0) sprintf(txt + strlen(txt), ":%s", map.sheetUTF8);
					mapIndex = checkGridPosition(en, mapIndex);
				};
				//mexPrintf("\n%lf,%lf\n", en->e, en->n);

				nDisplay[i] = en->n;
				eDisplay[i] = en->e;
				height[i] = findHeight(en, txt, LTFile);
			}
			return;
		}
	}
	else if ((nrhs == 2 || nrhs == 3) && checkCommand(TYP, "grid-to-gps") == 0)
	{
		EastingNorthing *en = (EastingNorthing*)malloc(sizeof(EastingNorthing)), *enETRS89;
		LatLonDecimal *latLon;

		// with no additional args, convert CSV (or similar) on STDIN to CSV on STDOUT
		if (nrhs == 2)
		{
			char *filename;
			if (mxIsChar(Location) != 1)
				mexErrMsgIdAndTxt("MATLAB:OSTN15_Matlab:inputNotString",
					"Input file name must be a string.");
			filename = mxArrayToString(Location);
			FILE *inputSequence = fopen(filename, "r");
			char *outname = (char*)malloc(60);
			strncpy(outname, filename, strlen(filename) - 4);
			strncpy(outname + strlen(filename) - 3, "\0", 1);
			strcat(outname, "-Output.csv\0");
			FILE *outputSequence = fopen(outname, "w");
			int scanResult;
			while ((scanResult = fscanf(inputSequence, "%lf,%lf", &en->e, &en->n)) == 2)
			{
				//mexPrintf("Dealing with Easting %lf , Northing %lf ", en->e, en->n);
				en->elevation = 0;
				enETRS89 = ETRS89EastingNorthingFromOSGB36EastingNorthing(en);
				if (enETRS89->geoid == 0) latLon->lat = latLon->lon = latLon->elevation = latLon->geoid = 0;
				else latLon = ETRS89LatLonFromETRS89EastingNorthing(enETRS89);
				//mexPrintf("as %.6lf,%.6lf\n", latLon->lat, latLon->lon);
				fprintf(outputSequence, "%.8lf,%.8lf\n", latLon->lat, latLon->lon);
			}
			fclose(outputSequence);
			fclose(inputSequence);
			return;
			
		}
		else 
		{
			double *EstIN, *NorIN, *latDisplay, *lngDisplay;
			
			EstIN = mxGetPr(OS36EAS_IN);
			NorIN = mxGetPr(OS36NOR_IN);

			int sizeM = mxGetM(OS36EAS_IN), sizeN = mxGetN(OS36EAS_IN);
			int checkM = mxGetM(OS36NOR_IN), checkN = mxGetN(OS36NOR_IN);
			if (sizeM != checkM || sizeN != checkN)
			{
				mexErrMsgIdAndTxt("MATLAB:OSTN15_Matlab:InputMismatch",
					"Lat and Long input rows does not match.");
				return;
			}

			LAT_OUT = mxCreateDoubleMatrix(sizeM, sizeN, mxREAL);
			LNG_OUT = mxCreateDoubleMatrix(sizeM, sizeN, mxREAL);

			latDisplay = mxGetPr(LAT_OUT);
			lngDisplay = mxGetPr(LNG_OUT);

			for (int i = 0; i < sizeM*sizeN; i++)
			{
				en->elevation = 0;
				en->e = EstIN[i];
				en->n = NorIN[i];

				enETRS89 = ETRS89EastingNorthingFromOSGB36EastingNorthing(en);
				if (enETRS89->geoid == 0)
				{
					mexErrMsgIdAndTxt("MATLAB:OSTN15_Matlab:inputExcess",
						"Input out of OSTN15 defined conversion limits.");
				}
				latLon = ETRS89LatLonFromETRS89EastingNorthing(enETRS89);

				latDisplay[i] = latLon->lat;
				lngDisplay[i] = latLon->lon;
				
			}
			return;
		}
	}
	else
	{
		mexPrintf("\n"  "OSTN15-Matlab -- \n"
			"\n"
			"This tool converts coordinates between ETRS89 (WGS84, GPS) lat/lon/elevation and OSGB36 easting/northing/elevation.\n"
			"Conversions make use of Ordnance Survey's OSTN15 and OSGM15 transformations, and should thus be accurate to within 1m.\n"
			"\n"
			"Usage: "  "\"[eas,nor,hei] = OSTN15_Matlab('gps-to-grid',lat,lon)\" converts ETRS89 to OSGB36 and checks LIDAR DSM heights data in ./LIDAR . If no data available hei will be -9999.\n"
			"       "  "\"[lat,lon] = OSTN15_Matlab('grid-to-gps',eas,nor)\" converts OSGB36 to ETRS89.\n"
			"       "  "\"OSTN15_Matlab('gps-to-grid',filename)\" converts a batch of data in filename from ETRS89 to OSGB36 and checks LIDAR DSM heights data in ./LIDAR . If no data available hei will be -9999, outputing to filename-output.csv and filename-raw.csv .\n"
			"       "  "\"OSTN15_Matlab('grid-to-gps',filename)\" converts a batch of data in filename OSGB36 to ETRS89, outputing to filename-output.csv.\n"
			"       "  "\"OSTN15_Matlab('test')\" checks embedded data integrity and runs conversion tests with known coordinates.\n"
			"       "  "\"OSTN15_Matlab\" displays this message, and runs a test.\n"
			"\n"
			"The conversion commands "  "gps-to-grid"  " and "  "grid-to-gps"  " can be used in two ways:\n"
			"* Given a set of coordinates as command-line arguments, they will convert this set with user-friendly output\n"
			"* Given no command-line arguments, they will convert batches of coordinates"
			"\n"
			"In the batch conversion case:\n"
			"* Input rows must have 2 columns -- lat or easting, lon or northing, with comma separator\n"
			"* Output rows have 3 columns -- easting or lat, northing or lon, elevation(if from gps to grid), separated with commas\n"
			"* In case of out-of-range input coordinates, all output columns will be zero\n"
			"* Malformatted input terminates processing and results in a non-zero exit code\n"
			"\n"
			"Software copyright (c) Tony Chen 2018\n"
			"Released under the MIT licence (http://www.opensource.org/licenses/mit-license.php)\n\n"
			"OSTN15 and OSGM15 are trademarks of Ordnance Survey\n"
			"Embedded OSTN15 and OSGM15 data (c) Crown copyright 2015. All rights reserved.\n");
		bool passed = test(true);
		return;
	}
}