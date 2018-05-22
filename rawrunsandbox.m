% raw matlab
epwDir = ".\\data";
uwgParamDir = ".\\data";
%uwgParamDir = "C:\\ladybug\\sl_uwg_fatal_1\\UWG\";
epwFileName = "CAN_ON_Toronto.716240_CWEC.epw";%"SGP_Singapore.486980_IWEC.epw";

%epwFileName = "CHN_Beijing.Beijing.545110_IWEC.epw";
%uwgParamFileName = "initialize_beijing.m";

uwgParamFileName = "initialize.m";
%uwgParamFileName = "initialize_toronto.m";
%uwgParamFileName = "sl_uwg_fatal_1.xml";

newDir = ".\\data";
%newFileName = "CAN_ON_Toronto.716240_CWEC_heatdemand_UWG_Matlab.epw";
newFileName = "sl_test";
UWG(epwDir, epwFileName, uwgParamDir, uwgParamFileName, newDir, newFileName);
