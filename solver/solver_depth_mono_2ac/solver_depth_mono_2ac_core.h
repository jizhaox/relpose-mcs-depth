#include <Eigen/Dense>
#include "mex.h"

using namespace Eigen;

MatrixXcd solver_depth_mono_2ac_core(const VectorXd& data)
{
	// Compute coefficients
    VectorXd coeffs = data;

	// Setup elimination template
	static const int coeffs0_ind[] = { 0,35,70,105,140,175,315,350,385,0,1,35,36,70,71,105,106,140,141,175,176,315,316,350,351,385,386,1,2,36,37,71,72,106,107,141,142,176,177,316,317,351,352,386,387,2,3,37,38,72,73,107,108,142,143,177,178,317,318,352,353,387,388,3,4,38,39,73,74,108,109,143,144,178,179,318,319,353,354,388,389,4,39,74,109,144,179,319,354,389,0,5,35,40,70,75,105,110,140,145,175,180,315,320,350,355,385,390,1,5,6,36,40,41,71,75,76,106,110,111,141,145,146,176,180,181,316,320,321,351,355,356,386,390,391,2,6,7,37,41,42,72,76,77,107,111,112,142,146,147,177,181,182,317,321,322,352,356,357,387,391,392,3,7,8,38,42,43,73,77,78,108,112,113,143,147,148,178,182,183,318,322,323,353,357,358,388,392,393,4,8,39,43,74,78,109,113,144,148,179,183,319,323,354,358,389,393,5,9,40,44,75,79,110,114,145,149,180,184,320,324,355,359,390,394,6,9,10,41,44,45,76,79,80,111,114,115,146,149,150,181,184,185,321,324,325,356,359,360,391,394,395,7,10,11,42,45,46,77,80,81,112,115,116,147,150,151,182,185,186,322,325,326,357,360,361,392,395,396,8,11,43,46,78,81,113,116,148,151,183,186,323,326,358,361,393,396,9,12,44,47,79,82,114,117,149,152,184,187,324,327,359,362,394,397,10,12,13,45,47,48,80,82,83,115,117,118,150,152,153,185,187,188,325,327,328,360,362,363,395,397,398,11,13,46,48,81,83,116,118,151,153,186,188,326,328,361,363,396,398,12,14,47,49,82,84,117,119,152,154,187,189,327,329,362,364,397,399,13,14,48,49,83,84,118,119,153,154,188,189,328,329,363,364,398,399,14,49,84,119,154,189,329,364,399,4,18,39,53,74,88,109,123,144,158,179,193,319,333,354,368,389,403,8,18,21,43,53,56,78,88,91,113,123,126,148,158,161,183,193,196,323,333,336,358,368,371,393,403,406,11,21,23,46,56,58,81,91,93,116,126,128,151,161,163,186,196,198,326,336,338,361,371,373,396,406,408,13,23,24,48,58,59,83,93,94,118,128,129,153,163,164,188,198,199,328,338,339,363,373,374,398,408,409,14,24,49,59,84,94,119,129,154,164,189,199,329,339,364,374,399,409,0,15,35,50,70,85,105,120,140,155,175,190,315,330,350,365,385,400,1,15,16,36,50,51,71,85,86,106,120,121,141,155,156,176,190,191,316,330,331,351,365,366,386,400,401,5,15,19,40,50,54,75,85,89,110,120,124,145,155,159,180,190,194,320,330,334,355,365,369,390,400,404,2,16,17,37,51,52,72,86,87,107,121,122,142,156,157,177,191,192,317,331,332,352,366,367,387,401,402,6,16,19,20,41,51,54,55,76,86,89,90,111,121,124,125,146,156,159,160,181,191,194,195,321,331,334,335,356,366,369,370,391,401,404,405,9,19,22,44,54,57,79,89,92,114,124,127,149,159,162,184,194,197,324,334,337,359,369,372,394,404,407,3,17,18,38,52,53,73,87,88,108,122,123,143,157,158,178,192,193,318,332,333,353,367,368,388,402,403,7,17,20,21,42,52,55,56,77,87,90,91,112,122,125,126,147,157,160,161,182,192,195,196,322,332,335,336,357,367,370,371,392,402,405,406,10,20,22,23,45,55,57,58,80,90,92,93,115,125,127,128,150,160,162,163,185,195,197,198,325,335,337,338,360,370,372,373,395,405,407,408,12,22,24,47,57,59,82,92,94,117,127,129,152,162,164,187,197,199,327,337,339,362,372,374,397,407,409 };
	static const int coeffs1_ind[] = { 34,69,104,139,174,209,349,384,419,31,34,66,69,101,104,136,139,171,174,206,209,346,349,381,384,416,419,25,31,60,66,95,101,130,136,165,171,200,206,340,346,375,381,410,416,15,25,50,60,85,95,120,130,155,165,190,200,330,340,365,375,400,410,16,25,26,51,60,61,86,95,96,121,130,131,156,165,166,191,200,201,331,340,341,366,375,376,401,410,411,19,25,28,54,60,63,89,95,98,124,130,133,159,165,168,194,200,203,334,340,343,369,375,378,404,410,413,26,31,32,61,66,67,96,101,102,131,136,137,166,171,172,201,206,207,341,346,347,376,381,382,411,416,417,17,26,27,52,61,62,87,96,97,122,131,132,157,166,167,192,201,202,332,341,342,367,376,377,402,411,412,20,26,28,29,55,61,63,64,90,96,98,99,125,131,133,134,160,166,168,169,195,201,203,204,335,341,343,344,370,376,378,379,405,411,413,414,28,31,33,63,66,68,98,101,103,133,136,138,168,171,173,203,206,208,343,346,348,378,381,383,413,416,418,22,28,30,57,63,65,92,98,100,127,133,135,162,168,170,197,203,205,337,343,345,372,378,380,407,413,415,32,34,67,69,102,104,137,139,172,174,207,209,347,349,382,384,417,419,27,32,62,67,97,102,132,137,167,172,202,207,342,347,377,382,412,417,18,27,53,62,88,97,123,132,158,167,193,202,333,342,368,377,403,412,21,27,29,56,62,64,91,97,99,126,132,134,161,167,169,196,202,204,336,342,344,371,377,379,406,412,414,29,32,33,64,67,68,99,102,103,134,137,138,169,172,173,204,207,208,344,347,348,379,382,383,414,417,418,23,29,30,58,64,65,93,99,100,128,134,135,163,169,170,198,204,205,338,344,345,373,379,380,408,414,415,33,34,68,69,103,104,138,139,173,174,208,209,348,349,383,384,418,419,30,33,65,68,100,103,135,138,170,173,205,208,345,348,380,383,415,418,24,30,59,65,94,100,129,135,164,170,199,205,339,345,374,380,409,415 };
	static const int C0_ind[] = { 3,7,11,15,19,23,27,31,35,38,39,42,43,46,47,50,51,54,55,58,59,62,63,66,67,70,71,74,75,78,79,82,83,86,87,90,91,94,95,98,99,102,103,106,107,110,111,114,115,118,119,122,123,126,127,130,131,134,135,138,139,142,143,146,147,150,151,154,155,158,159,162,163,166,167,170,171,174,175,178,179,182,186,190,194,198,202,206,210,214,217,219,221,223,225,227,229,231,233,235,237,239,241,243,245,247,249,251,253,254,255,257,258,259,261,262,263,265,266,267,269,270,271,273,274,275,277,278,279,281,282,283,285,286,287,289,290,291,293,294,295,297,298,299,301,302,303,305,306,307,309,310,311,313,314,315,317,318,319,321,322,323,325,326,327,329,330,331,333,334,335,337,338,339,341,342,343,345,346,347,349,350,351,353,354,355,357,358,359,361,362,365,366,369,370,373,374,377,378,381,382,385,386,389,390,393,394,397,399,401,403,405,407,409,411,413,415,417,419,421,423,425,427,429,431,433,434,435,437,438,439,441,442,443,445,446,447,449,450,451,453,454,455,457,458,459,461,462,463,465,466,467,469,470,471,473,474,475,477,478,479,481,482,483,485,486,487,489,490,491,493,494,495,497,498,499,501,502,503,505,506,509,510,513,514,517,518,521,522,525,526,529,530,533,534,537,538,541,543,545,547,549,551,553,555,557,559,561,563,565,567,569,571,573,575,577,578,579,581,582,583,585,586,587,589,590,591,593,594,595,597,598,599,601,602,603,605,606,607,609,610,611,613,614,617,618,621,622,625,626,629,630,633,634,637,638,641,642,645,646,649,651,653,655,657,659,661,663,665,667,669,671,673,675,677,679,681,683,685,686,689,690,693,694,697,698,701,702,705,706,709,710,713,714,717,718,721,725,729,733,737,741,745,749,753,756,758,760,762,764,766,768,770,772,774,776,778,780,782,784,786,788,790,792,793,794,796,797,798,800,801,802,804,805,806,808,809,810,812,813,814,816,817,818,820,821,822,824,825,826,828,829,830,832,833,834,836,837,838,840,841,842,844,845,846,848,849,850,852,853,854,856,857,858,860,861,862,864,865,866,868,869,870,872,873,874,876,877,878,880,881,882,884,885,886,888,889,890,892,893,894,896,897,898,900,901,904,905,908,909,912,913,916,917,920,921,924,925,928,929,932,933,936,939,940,943,944,947,948,951,952,955,956,959,960,963,964,967,968,971,972,974,975,976,978,979,980,982,983,984,986,987,988,990,991,992,994,995,996,998,999,1000,1002,1003,1004,1006,1007,1008,1009,1011,1012,1013,1015,1016,1017,1019,1020,1021,1023,1024,1025,1027,1028,1029,1031,1032,1033,1035,1036,1037,1039,1040,1041,1043,1044,1046,1047,1048,1050,1051,1052,1054,1055,1056,1058,1059,1060,1062,1063,1064,1066,1067,1068,1070,1071,1072,1074,1075,1076,1078,1079,1080,1081,1082,1083,1084,1085,1086,1087,1088,1089,1090,1091,1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1111,1112,1113,1114,1115,1116,1117,1119,1120,1121,1123,1124,1125,1127,1128,1129,1131,1132,1133,1135,1136,1137,1139,1140,1141,1143,1144,1145,1147,1148,1149,1151,1152,1154,1155,1156,1158,1159,1160,1162,1163,1164,1166,1167,1168,1170,1171,1172,1174,1175,1176,1178,1179,1180,1182,1183,1184,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1213,1214,1215,1216,1217,1218,1219,1220,1221,1222,1223,1224,1225,1226,1227,1228,1229,1230,1231,1232,1233,1234,1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1263,1264,1265,1267,1268,1269,1271,1272,1273,1275,1276,1277,1279,1280,1281,1283,1284,1285,1287,1288,1289,1291,1292,1293,1295 } ;
	static const int C1_ind[] = { 0,4,8,12,16,20,24,28,32,36,39,40,43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,87,88,91,92,95,96,99,100,103,104,107,108,111,112,115,116,119,120,123,124,127,128,131,132,135,136,139,140,143,144,146,147,148,150,151,152,154,155,156,158,159,160,162,163,164,166,167,168,170,171,172,174,175,176,178,179,180,181,183,184,185,187,188,189,191,192,193,195,196,197,199,200,201,203,204,205,207,208,209,211,212,213,215,216,218,219,220,222,223,224,226,227,228,230,231,232,234,235,236,238,239,240,242,243,244,246,247,248,250,251,252,254,255,256,258,259,260,262,263,264,266,267,268,270,271,272,274,275,276,278,279,280,282,283,284,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,327,328,329,331,332,333,335,336,337,339,340,341,343,344,345,347,348,349,351,352,353,355,356,357,359,360,361,363,364,365,367,368,369,371,372,373,375,376,377,379,380,381,383,384,385,387,388,389,391,392,393,395,396,398,400,402,404,406,408,410,412,414,416,418,420,422,424,426,428,430,432,434,436,438,440,442,444,446,448,450,452,454,456,458,460,462,464,466,468,470,472,474,476,478,480,482,484,486,488,490,492,494,496,498,500,502,504,505,506,508,509,510,512,513,514,516,517,518,520,521,522,524,525,526,528,529,530,532,533,534,536,537,538,540,541,542,544,545,546,548,549,550,552,553,554,556,557,558,560,561,562,564,565,566,568,569,570,572,573,574,576,577,578,580,581,582,584,585,586,588,589,590,592,593,594,596,597,598,600,601,602,604,605,606,608,609,610,612,613,616,617,620,621,624,625,628,629,632,633,636,637,640,641,644,645,648,649,652,653,656,657,660,661,664,665,668,669,672,673,676,677,680,681,684,685,688,689,692,693,696,697,700,701,704,705,708,709,712,713,716,717 };

	Matrix<double,36,36> C0; C0.setZero();
	Matrix<double,36,20> C1; C1.setZero();
	for (int i = 0; i < 810; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 450; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

	Matrix<double,36,20> C12 = C0.partialPivLu().solve(C1);

	// Setup action matrix
	Matrix<double,30, 20> RR;
	RR << -C12.bottomRows(10), Matrix<double,20,20>::Identity(20, 20);

	static const int AM_ind[] = { 11,12,13,0,1,2,14,3,4,15,5,16,17,6,7,18,8,19,20,9 };
	Matrix<double, 20, 20> AM;
	for (int i = 0; i < 20; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 3, 20> sols;
	sols.setZero();

	// Solve eigenvalue problem
	EigenSolver<Matrix<double, 20, 20> > es(AM);
	ArrayXcd D = es.eigenvalues();
	ArrayXXcd V = es.eigenvectors();
    V = (V / V.row(0).array().replicate(20, 1)).eval();

    sols.row(0) = D.transpose().array();
    sols.row(1) = V.row(11).array();
    sols.row(2) = V.row(17).array();

	return sols;
}

