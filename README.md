# OSTN15-Matlab
An implementation of Ordnance Survey conversion from OSGB36 to WGS84. For MATLAB use, but can be changed into general use

This tool converts coordinates between ETRS89 (WGS84, GPS) lat/lon/elevation and OSGB36 easting/northing/elevation.

Conversions make use of Ordnance Survey's OSTN15 and OSGM15 transformations, and should thus be accurate to within 1m.
Usage:   
 - [eas,nor,hei] = OSTN15_Matlab('gps-to-grid',lat,lon)\ converts matrices of ETRS89 to OSGB36 and checks LIDAR DSM heights data in ./LIDAR . If no data available hei will be -9999.
 - [lat,lon] = OSTN15_Matlab('grid-to-gps',eas,nor)\ converts matrices of OSGB36 to ETRS89.
 - OSTN15_Matlab('gps-to-grid',filename)\ converts a batch of data in filename from ETRS89 to OSGB36 and checks LIDAR DSM heights data in ./LIDAR . If no data available hei will be -9999, outputing to filename-output.csv and filename-raw.csv .
 - OSTN15_Matlab('grid-to-gps',filename)\ converts a batch of data in filename OSGB36 to ETRS89, outputing to filename-output.csv.
 - OSTN15_Matlab('test')\ checks embedded data integrity and runs conversion tests with known coordinates.
 - OSTN15_Matlab\ displays this message, and runs a test.
			
In the batch conversion case:
- Input rows must have 2 columns -- lat or easting, lon or northing, with comma separator
- Output rows have 3 columns -- easting or lat, northing or lon, elevation(if from gps to grid), separated with commas. If no elevation data available, elevation returns -9999. At eas-nor mode, a heading line tells which block the point is in. This helps to decide which LIDAR data file to download
- In case of out-of-range input coordinates, all output columns will be zero
- Malformatted input leads unexpected output

Software copyright (c) Tony Chen 2018
Released under the MIT licence (http://www.opensource.org/licenses/mit-license.php)
OSTN15 and OSGM15 are trademarks of Ordnance Survey
Embedded OSTN15 and OSGM15 data (c) Crown copyright 2015. All rights reserved.)
LIDAR Data uses 1m DSM data at best. LIDAR data use under the Open Government Licence v3.0 of Environmental Agency, GB Govt.
[LIDAR DSM-1m](http://environment.data.gov.uk/ds/survey/#/survey)
