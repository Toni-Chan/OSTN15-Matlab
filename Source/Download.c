#include"OSTN02.h"
#include<windows.h>
#include <tchar.h>
#include <urlmon.h>
#pragma comment(lib, "urlmon.lib")
#pragma comment(lib,"wininet.lib")
#include "D:\\Program Files\\MATLAB\\R2018a\\extern\\include\mex.h"

FILE* getMissingMap(char *BLOCK, int fileE, int fileN)
{
	HRESULT hr;
	char *name = (char *)malloc(500);
	sprintf(name, ".\\LIDAR\\%s%02d%02d_DSM_1M.asc", BLOCK, fileE, fileN);
	LPCTSTR Url = _T(""), File = _T(name);
	hr = URLDownloadToFile(0, Url, File, 0, 0);
	if (hr != S_OK)
	{
		switch (hr)
		{
		case E_OUTOFMEMORY:
			mexErrMsgIdAndTxt("MATLAB:OSTN02_Matlab:fileNotFound",
				"Error while downloading a missing local file: No remaining memory");
			break;
		case INET_E_DOWNLOAD_FAILURE:
			mexErrMsgIdAndTxt("MATLAB:OSTN02_Matlab:fileNotFound",
				"Error while downloading a missing local file: Unable to fetch server file");
			break;
		default:
			mexErrMsgIdAndTxt("MATLAB:OSTN02_Matlab:fileNotFound",
				"Error while downloading a missing local file: Unknown error"); 
			break;
		}
	}
	FILE *f = fopen(name, "r");
	return f;
}