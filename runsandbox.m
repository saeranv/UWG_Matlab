%epwDir = "C:\\saeran\\master\\git\\UWG_Matlab\\data";
%uwgParamDir = "C:\\saeran\\master\\git\\UWG_Matlab\\data";
epwDir = ".\\data";
epwFileName = "CAN_ON_Toronto.716240_CWEC.epw";%"SGP_Singapore.486980_IWEC.epw";
uwgParamDir = ".\\data";
%uwgParamFileName = "initialize_toronto.m";
uwgParamFileName = "initialize.m"; %_2_feb16.m";
UWG(epwDir, epwFileName, uwgParamDir, uwgParamFileName);


