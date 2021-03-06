function [result] = HelicalDriveStrength(T1, z1, u12, KA, Lh1, Lh2, ZL1, ZL2, ZX, Zw, ZNmax, YNmax, NstH, NBH, NstF, NBF, SHmin, SFmin, YR1, YR2, Hi1, Hi2, niu1, niu2, E1, E2, MatTT1, MatTT2, HB1, HB2, sigma021, sigma022, sigma_Hlim1, sigma_Hlim2, sigma_Flim1, sigma_Flim2, Raf1, Raf2, Rar1, Rar2, PozAngr, CPrec, mn, aw, alfa_t, alfa_wt, beta, beta_b, eps_alfa, eps_beta, eps_alfa_n, b1, b2, zn1, zn2, xn1, xn2, d1, d2, n1, n2, HDS)
%         // T1			 - Momentul de torsiune corespunzator arborelui de intrare, Nmm
% 		// z1			 - Numerele de dinti ai pinionului
% 		// u12			 - Raportul deal de angrenare corespunzator treptei I
% 		// KA			 - Factorul regimului de functionare
% 		// Lh1,Lh2		 - Durata minima de functionare a pinionului/rotii dintate, ore
% 		// ZL1,ZL2		 - Factorul de ungere corespunzator pinionului/rotii dintate
% 		// Zx			 - Factorul de marime
% 		// ZW			 - Factorul raportului duritatilor flancurilor dintilor
% 		// ZNmax		 - Factorul durabilitatii pentru solicitarea de contact
% 		// YNmax		 - Factorul durabilitatii pentru solicitarea de incovoiere
% 		// NstH,NstF	 - 
% 		// NBH,NBF		 -
% 		// SHmin		 - Coeficientul de siguranta minim pentru solicitarea de contact
% 		// SFmin		 - Coeficientul de siguranta minim pentru solicitarea de incovoiere
% 		// YR1,YR2		 - Factorul rugozitatii flancului dintelui pinionului/rotii dintate pentru solicitarea de incovoiere
% 		// Hi1,Hi2		 - Numarul de roti cu care vine in contact pinionul/roata dintata
% 		// niu1,niu2	 - Coeficientii lui Poisson
% 		// E1,E2		 - Modulul de elasticitate longitudinal, MPa
% 		// MatTT1,MatTT2 -
% 		// HB1,HB2		 - Duritatea materialului pinionului/rotii dintate, MPa
% 		// sigma021,2	 - Limita de curegere a materialului pinionului/rotii dintate, MPa
% 		// sigma_Hlim1,2 - Tensiunea hertziana limita pentru materialul pinionului/rotii dintate, MPa
% 		// sigma_Flim1,2 - Tensiunea de incovoiere limita pentru materialul pinionului/rotii dintate, MPa
% 		// Raf1,Raf2	 - Rugozitatea flancului dintelui pinionului/rotii dintate, micrometri
% 		// Rar1,Rar2	 - Rugozitatea razei de racordare a dintelui pinionului/rotii dintate, micrometri
% 		// PozAngr		 - Pozitia angrenajului
% 		// CPrec		 - Clasa de precizie a angrenajului
% 		// mn			 - Modulul normal al angrenajului, mm
% 		// aw			 - Distanta axiala standardizata/impusa, mm
% 		// alfa_t		 - Unghiul de angrenare de referinta in plan frontal, rad
% 		// alfa_wt		 - Unghiul real de angrenare in plan frontal, rad
% 		// beta			 - Unghiul de inclinare al danturii pe cilindrul de divizare, grade
% 		// beta_b		 - Unghiul de inclinare al danturii pe cilindrul de baza, rad
% 		// eps_alfa		 - Gradul de acoperire al angrenajului in plan frontal
% 		// eps_beta		 - Gradul de acoperire al angrenajului in plan axial
% 		// eps_alfa_n	 - Gradul de acoperire al angrenajului echivalent
% 		// b1,b2		 - Latimea pinionului/rotii dintate, mm
% 		// zn1,zn2		 - Numarele de dinti ale angrenajului echivalent
% 		// xn1,xn2		 - Coeficiemtii deplasarilor de profil in plan normal corespunzatori pinionului/rotii dintate
% 		// d1,d2		 - Diametrele cercurilor de divizare ale pinionului/rotii dintate, mm
% 		// n1,n2		 - Turatia arborelui 1 si 2, rot/min 

    PI = vpa(pi, 100);
% Factorul rugozitatii flancurilor Zr pentru solicitarea de contact
    GZR = [1.689,1.918,2.106,2.315,2.607,2.815,3.149,3.462,3.796,4.172,4.526,4.881,5.299,5.695,6.112,6.446,6.864,7.323,7.782,8.262,8.742,9.242,9.743,10.223,10.766,11.266,11.788,12.268,12.685,13.042,13.520,14.000, 1.100,1.079,1.060,1.040,1.022,1.007,0.993,0.977,0.962,0.948,0.937,0.927,0.915,0.905,0.896,0.887,0.878,0.868,0.859,0.850,0.843,0.835,0.829,0.822,0.817,0.813,0.809,0.805,0.803,0.801,0.800,0.799];
    PZR(1) = [850];
    
% Factorul de viteza Zv
    GZv = [1.000,1.165,1.357,1.605,1.870,2.193,2.585,3.047,3.561,4.102,4.866,5.595,6.487,7.520,8.863,9.943,11.247,12.722,14.522,16.412,18.871,21.385,24.189,27.763,32.187,37.317,43.263,50.340,58.576,69.098,83.467,100.000,
    0.901,0.905,0.909,0.914,0.919,0.925,0.931,0.937,0.944,0.951,0.959,0.966,0.974,0.982,0.993,1.000,1.008,1.015,1.023,1.031,1.039,1.047,1.055,1.064,1.074,1.083,1.092,1.100,1.108,1.115,1.121,1.124];
    PZv(1) = [850];
    
% Factorul de forma al dintelui YFa pentru solicitarea de incovoiere 
    GYFa = [15.423,15.955,16.416,16.931,17.542,18.134,18.876,19.935,21.002,22.001,23.051,24.525,26.301,27.663,28.990,30.819,32.505,34.620,36.701,39.306,41.723,46.297,50.832,56.764,65.714,75.439,89.558,112.077,147.021,208.569,274.782,400,
     1.854,1.856,1.859,1.861,1.864,1.866,1.869,1.874,1.878,1.883,1.887,1.891,1.897,1.901,1.906,1.910,1.914,1.920,1.924,1.931,1.936,1.945,1.951,1.959,1.970,1.980,1.990,2.002,2.013,2.027,2.035,2.043,
     13.442,14.060,14.526,15.062,15.681,16.233,16.821,17.708,18.685,19.824,20.937,22.089,23.169,24.658,26.402,28.204,30.337,32.721,35.202,38.295,41.920,45.487,49.703,56.211,64.457,75.392,86.796,102.781,123.693,158.903,202.822,400,
     1.900,1.901,1.902,1.903,1.905,1.906,1.907,1.908,1.911,1.913,1.917,1.919,1.922,1.926,1.929,1.933,1.938,1.941,1.947,1.952,1.959,1.963,1.970,1.976,1.984,1.993,2.000,2.008,2.015,2.024,2.031,2.048,
     11.652,12.065,12.468,12.926,13.532,14.178,14.915,15.754,16.718,17.480,18.184,19.010,19.979,21.127,22.282,23.619,25.003,27.369,30.231,32.473,35.953,40.293,45.010,50.484,59.705,68.911,79.915,94.883,113.951,147.553,227.989,400,
     1.980,1.977,1.975,1.973,1.971,1.969,1.968,1.966,1.965,1.965,1.965,1.964,1.965,1.964,1.965,1.966,1.967,1.970,1.973,1.976,1.979,1.983,1.987,1.993,1.999,2.006,2.011,2.017,2.023,2.031,2.040,2.049,
     10.033,10.364,10.754,11.147,11.677,12.263,12.873,13.473,14.126,15.071,16.054,16.994,18.114,19.109,20.669,21.934,23.594,25.457,27.535,30.063,32.110,35.318,37.877,41.288,45.286,50.742,59.556,73.435,92.445,125.273,198.877,400,
     2.106,2.099,2.091,2.085,2.076,2.068,2.061,2.054,2.048,2.044,2.037,2.033,2.029,2.026,2.023,2.021,2.018,2.017,2.015,2.015,2.015,2.015,2.016,2.016,2.016,2.019,2.022,2.026,2.032,2.037,2.046,2.050,
     10.042,10.486,10.872,11.319,11.636,12.039,12.456,13.040,13.738,14.504,15.494,16.400,17.523,18.707,19.944,21.136,22.488,24.396,26.107,28.748,31.238,34.134,37.019,41.924,49.846,57.942,70.564,89.904,115.967,151.435,212.942,400,
     2.240,2.224,2.212,2.200,2.192,2.182,2.173,2.162,2.150,2.140,2.128,2.119,2.109,2.101,2.094,2.088,2.082,2.075,2.071,2.066,2.063,2.058,2.056,2.054,2.049,2.049,2.050,2.049,2.050,2.051,2.053,2.056,
     10.036,10.288,10.690,11.103,11.578,12.109,12.540,13.161,13.578,14.128,14.755,15.493,16.335,16.910,17.584,18.617,20.194,21.688,23.244,25.267,27.521,30.043,33.874,38.641,46.351,57.010,73.127,89.580,113.932,147.074,208.324,400,
     2.399,2.387,2.368,2.350,2.331,2.312,2.298,2.280,2.269,2.257,2.244,2.231,2.218,2.209,2.200,2.187,2.172,2.160,2.148,2.137,2.128,2.119,2.108,2.098,2.088,2.080,2.074,2.070,2.068,2.066,2.064,2.063,
     10.349,10.727,11.132,11.554,12.020,12.639,13.205,13.790,14.498,15.300,16.149,17.219,18.330,19.438,20.941,22.493,24.179,26.110,28.476,30.386,32.415,35.337,39.201,44.690,50.362,58.596,71.669,89.108,111.527,138.737,208.061,400,
     2.572,2.545,2.520,2.496,2.472,2.445,2.422,2.402,2.379,2.357,2.336,2.313,2.292,2.275,2.255,2.237,2.221,2.206,2.191,2.179,2.171,2.160,2.148,2.135,2.124,2.114,2.102,2.093,2.085,2.079,2.073,2.066,
     12.041,12.477,12.984,13.473,14.041,14.646,15.344,16.308,17.063,18.052,19.029,20.112,21.342,22.986,24.739,26.499,28.833,31.493,33.856,37.234,40.678,45.309,51.179,60.370,70.608,79.736,92.222,110.265,143.858,202.528,265.057,400,
     2.659,2.631,2.602,2.576,2.551,2.525,2.498,2.465,2.443,2.417,2.395,2.372,2.350,2.325,2.302,2.284,2.263,2.243,2.230,2.211,2.197,2.182,2.166,2.148,2.135,2.125,2.116,2.106,2.095,2.085,2.078,2.073,
     13.768,14.105,14.503,14.961,15.509,15.917,16.425,16.887,17.583,18.259,18.748,19.334,20.128,21.432,22.678,23.899,25.484,27.099,28.825,30.457,33.110,35.116,38.476,42.829,48.870,56.132,68.319,80.656,105.411,148.204,218.751,400,
     2.752,2.731,2.708,2.682,2.657,2.637,2.616,2.596,2.572,2.548,2.532,2.515,2.494,2.464,2.438,2.415,2.390,2.367,2.348,2.330,2.307,2.293,2.269,2.246,2.221,2.199,2.174,2.154,2.132,2.110,2.095,2.080,
     15.465,15.848,16.311,16.743,17.104,17.653,18.261,18.953,19.918,20.787,21.790,22.950,24.298,25.590,27.284,29.077,31.395,34.061,37.493,40.480,44.209,49.224,55.665,64.086,73.217,84.182,99.224,111.339,139.624,190.518,259.934,400,
     2.851,2.828,2.801,2.777,2.758,2.732,2.706,2.678,2.642,2.613,2.584,2.553,2.522,2.496,2.467,2.439,2.409,2.380,2.348,2.326,2.300,2.274,2.249,2.223,2.202,2.182,2.163,2.151,2.133,2.114,2.098,2.085,
     17.221,17.660,18.159,18.727,19.296,20.027,20.799,21.698,22.651,23.597,24.936,26.210,27.551,29.247,31.159,33.341,36.166,38.700,42.221,45.466,49.937,54.589,59.855,67.699,76.628,90.699,109.059,133.679,165.646,221.848,282.297,400,
     2.951,2.923,2.892,2.863,2.836,2.800,2.768,2.734,2.700,2.672,2.635,2.604,2.575,2.542,2.510,2.480,2.445,2.418,2.386,2.360,2.330,2.307,2.285,2.258,2.234,2.206,2.180,2.159,2.139,2.120,2.106,2.093,
     18.950,19.570,20.192,20.929,21.742,22.810,23.859,24.830,25.962,27.109,28.141,29.576,30.900,32.718,34.349,36.674,38.948,41.547,45.445,49.399,54.149,59.666,66.136,74.109,83.914,96.341,111.680,133.580,165.066,200.730,258.259,400,
     3.055,3.013,2.979,2.940,2.903,2.857,2.816,2.781,2.747,2.713,2.685,2.653,2.626,2.593,2.565,2.532,2.502,2.473,2.431,2.401,2.370,2.340,2.311,2.282,2.256,2.230,2.206,2.182,2.159,2.141,2.124,2.101,
     20.604,21.158,21.776,22.678,23.449,24.481,25.741,26.729,28.051,29.279,30.365,31.962,33.650,35.531,37.922,41.149,44.809,48.298,51.547,55.593,59.852,65.086,72.044,78.815,88.559,99.499,111.901,131.231,166.180,195.383,284.662,400,
     3.167,3.128,3.093,3.044,3.006,2.957,2.909,2.874,2.832,2.795,2.764,2.729,2.690,2.655,2.615,2.565,2.521,2.485,2.456,2.425,2.400,2.369,2.340,2.313,2.287,2.260,2.235,2.210,2.179,2.160,2.130,2.111,
     22.351,23.018,23.654,24.403,25.142,26.035,26.926,27.982,29.209,30.312,31.748,32.983,34.385,36.295,38.006,39.849,42.027,44.073,47.437,50.863,56.192,62.863,69.143,77.013,90.394,107.812,127.317,154.703,192.180,231.985,286.242,400,
     3.278,3.235,3.195,3.154,3.112,3.076,3.036,2.992,2.948,2.908,2.867,2.837,2.800,2.759,2.722,2.691,2.653,2.620,2.575,2.538,2.492,2.444,2.408,2.371,2.320,2.277,2.243,2.211,2.183,2.159,2.141,2.122,
     24.080,25.056,26.001,27.003,28.081,29.266,30.692,31.560,32.536,33.627,34.891,36.451,37.863,39.336,41.431,44.451,48.098,51.757,55.631,59.244,64.471,70.188,76.387,83.818,92.281,104.174,119.311,137.936,168.401,204.888,250.783,400,
     3.392,3.330,3.279,3.227,3.168,3.114,3.059,3.027,2.995,2.961,2.927,2.887,2.848,2.815,2.773,2.715,2.660,2.615,2.575,2.541,2.500,2.465,2.425,2.390,2.358,2.323,2.287,2.258,2.222,2.192,2.166,2.133,
     25.864,26.755,27.770,28.729,29.763,30.836,31.879,32.698,33.749,34.905,36.250,37.798,39.652,41.009,42.643,44.284,46.237,49.396,51.375,53.999,57.521,62.289,67.729,74.591,83.763,95.208,110.371,134.681,161.356,211.304,250.164,400,
     3.517,3.453,3.392,3.339,3.287,3.234,3.190,3.157,3.120,3.078,3.036,2.987,2.939,2.904,2.868,2.832,2.790,2.741,2.712,2.675,2.640,2.589,2.543,2.498,2.445,2.397,2.346,2.296,2.253,2.210,2.182,2.139];
     PYFa = [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5];

    % Factorul de corectie al tensiunilor de incovoiere la baza dintelui YSa 
    GYSa = [25.675,28.003,30.046,32.245,34.567,37.389,40.161,43.782,47.010,50.338,54.475,58.397,62.728,68.178,73.388,80.022,87.826,94.433,102.989,112.076,119.936,129.049,140.550,155.145,170.557,191.672,210.689,239.346,274.808,309.771,371.515,400,
     1.382,1.400,1.415,1.428,1.442,1.457,1.472,1.489,1.504,1.518,1.534,1.547,1.561,1.577,1.592,1.609,1.627,1.641,1.657,1.673,1.685,1.698,1.712,1.729,1.745,1.761,1.776,1.794,1.811,1.829,1.847,1.853,
     23.923,25.253,26.649,28.141,29.567,31.303,33.455,35.789,38.316,40.572,43.257,45.921,48.851,52.042,54.958,59.099,63.828,70.486,76.407,82.538,91.269,98.843,110.258,121.052,133.836,148.950,168.383,197.204,219.421,253.331,297.763,400,
     1.409,1.420,1.430,1.441,1.451,1.462,1.474,1.487,1.501,1.512,1.524,1.536,1.548,1.560,1.571,1.584,1.599,1.618,1.633,1.648,1.666,1.681,1.699,1.715,1.732,1.750,1.767,1.787,1.803,1.820,1.843,1.864,
     22.204,23.597,25.192,26.544,27.868,29.329,30.897,32.692,34.889,37.179,39.566,42.226,44.992,48.480,52.623,56.639,61.914,67.556,74.031,80.757,89.354,100.883,110.544,121.755,136.503,154.933,177.846,203.976,234.861,268.826,314.365,400,
     1.436,1.448,1.461,1.471,1.481,1.491,1.500,1.511,1.522,1.534,1.546,1.558,1.571,1.585,1.600,1.614,1.629,1.646,1.661,1.677,1.694,1.715,1.730,1.744,1.761,1.780,1.797,1.816,1.831,1.845,1.862,1.883,
     20.535,21.571,22.664,23.701,24.852,26.280,27.645,28.972,30.742,32.674,34.633,36.733,39.101,41.774,44.894,49.322,54.590,59.724,66.064,74.470,83.458,94.470,106.026,117.638,131.135,150.340,172.169,198.584,234.610,279.152,335.925,400,
     1.465,1.475,1.485,1.493,1.502,1.513,1.522,1.531,1.541,1.552,1.562,1.573,1.585,1.597,1.610,1.627,1.645,1.660,1.677,1.696,1.715,1.735,1.752,1.767,1.782,1.800,1.817,1.834,1.850,1.865,1.883,1.897,
     18.763,20.165,21.661,23.150,24.339,25.491,26.906,28.394,30.018,31.651,33.370,35.199,37.607,40.710,44.088,48.422,52.631,57.745,63.510,69.839,76.779,85.370,94.489,105.517,118.114,131.751,148.449,171.101,199.355,237.241,294.504,400,
     1.495,1.509,1.522,1.535,1.544,1.552,1.562,1.571,1.581,1.590,1.599,1.608,1.620,1.633,1.647,1.663,1.677,1.691,1.706,1.721,1.735,1.751,1.766,1.781,1.795,1.809,1.823,1.838,1.853,1.869,1.886,1.908,
     16.993,17.801,18.642,19.699,20.635,21.900,23.565,25.528,27.212,29.146,31.072,33.448,36.180,39.249,42.771,45.625,49.091,53.136,57.652,64.523,72.903,82.608,93.321,106.338,123.609,148.988,170.765,195.365,230.612,263.927,311.154,400,
     1.524,1.533,1.542,1.553,1.561,1.572,1.585,1.599,1.610,1.621,1.632,1.644,1.656,1.669,1.682,1.693,1.704,1.716,1.729,1.744,1.761,1.778,1.794,1.810,1.827,1.847,1.860,1.871,1.884,1.895,1.906,1.918,
     15.393,16.386,17.399,18.405,19.634,21.047,22.357,23.886,25.294,26.942,28.946,31.233,33.045,35.182,37.131,39.438,42.325,45.848,48.832,52.334,55.847,61.341,68.468,77.231,89.726,103.306,120.856,138.801,168.565,214.881,275.814,400,
     1.556,1.568,1.579,1.590,1.602,1.614,1.624,1.635,1.645,1.655,1.666,1.677,1.686,1.695,1.703,1.712,1.722,1.734,1.743,1.753,1.762,1.773,1.787,1.802,1.819,1.834,1.850,1.863,1.879,1.896,1.912,1.929,
     13.736,14.215,14.697,15.241,16.021,16.910,17.721,19.028,20.193,21.498,22.837,24.372,26.240,27.977,29.969,32.611,35.708,39.332,42.978,47.526,53.409,60.441,67.884,78.481,91.259,106.318,123.544,144.350,176.857,225.844,288.043,400,
     1.586,1.593,1.600,1.607,1.616,1.626,1.635,1.647,1.656,1.667,1.677,1.687,1.698,1.707,1.717,1.729,1.741,1.754,1.766,1.779,1.793,1.807,1.820,1.834,1.850,1.864,1.876,1.887,1.901,1.915,1.927,1.939,
     12.030,12.517,13.021,13.484,14.027,14.689,15.513,16.418,17.531,18.894,20.319,21.613,23.008,24.575,26.747,28.987,31.569,33.823,36.669,40.134,44.431,49.106,54.295,59.541,69.172,80.694,93.712,114.791,137.890,177.612,245.979,400,
     1.615,1.623,1.630,1.637,1.645,1.653,1.663,1.673,1.684,1.697,1.708,1.717,1.726,1.736,1.747,1.758,1.769,1.777,1.786,1.797,1.809,1.820,1.830,1.839,1.853,1.867,1.879,1.893,1.905,1.918,1.932,1.947,
     10.314,10.672,10.953,11.270,11.673,12.157,12.639,13.330,14.128,15.037,15.968,16.965,18.058,19.251,20.476,22.102,23.832,25.820,27.925,30.642,33.758,37.617,42.292,48.084,55.082,64.978,77.953,95.182,119.369,167.223,231.894,400,
     1.637,1.645,1.650,1.656,1.663,1.671,1.678,1.688,1.699,1.710,1.720,1.730,1.740,1.749,1.758,1.768,1.778,1.788,1.797,1.808,1.818,1.829,1.841,1.853,1.865,1.878,1.891,1.903,1.916,1.931,1.942,1.956,
     10.026,10.501,11.010,11.630,12.269,12.853,13.414,14.022,14.678,15.349,16.179,17.190,18.134,18.891,20.507,22.005,23.718,25.846,28.026,31.243,35.132,39.899,45.484,51.387,60.141,69.894,83.047,100.440,127.357,173.717,256.675,400,
     1.684,1.694,1.704,1.714,1.724,1.732,1.740,1.747,1.755,1.763,1.771,1.780,1.788,1.794,1.804,1.813,1.821,1.831,1.839,1.850,1.860,1.871,1.881,1.890,1.900,1.909,1.918,1.927,1.936,1.946,1.954,1.962,
     10.026,10.379,10.693,11.089,11.475,11.981,12.400,12.971,13.728,14.429,15.173,16.053,17.045,18.167,19.300,20.571,22.122,23.631,25.338,27.378,30.015,32.665,35.675,39.653,45.088,52.135,60.279,73.578,96.455,134.612,217.316,400,
     1.732,1.739,1.746,1.753,1.759,1.767,1.773,1.781,1.790,1.798,1.806,1.814,1.823,1.831,1.839,1.846,1.854,1.861,1.867,1.874,1.882,1.888,1.895,1.902,1.910,1.918,1.925,1.933,1.943,1.951,1.960,1.967,
     10.026,10.478,10.883,11.379,11.913,12.411,12.871,13.506,14.119,14.944,15.882,16.977,18.096,19.386,20.935,22.398,24.111,25.658,27.670,29.867,31.877,34.720,38.102,42.898,49.041,57.206,69.174,84.068,106.561,144.825,229.234,400,
     1.773,1.782,1.790,1.798,1.806,1.813,1.819,1.827,1.833,1.842,1.850,1.859,1.867,1.875,1.883,1.889,1.895,1.900,1.906,1.911,1.915,1.920,1.925,1.931,1.937,1.943,1.949,1.953,1.958,1.963,1.967,1.970,
     11.646,12.125,12.458,12.838,13.237,13.760,14.267,14.746,15.227,15.959,16.727,17.602,18.626,19.962,21.062,22.327,23.697,25.394,27.112,28.990,31.602,34.523,37.840,41.687,47.398,55.081,64.679,78.424,97.194,129.376,182.978,400,
     1.832,1.838,1.843,1.847,1.852,1.858,1.863,1.868,1.872,1.878,1.884,1.890,1.896,1.903,1.907,1.912,1.917,1.922,1.926,1.930,1.935,1.939,1.943,1.947,1.951,1.955,1.959,1.962,1.965,1.969,1.971,1.972,
     14.175,14.786,15.335,16.035,16.797,17.679,18.495,19.221,19.947,21.477,23.124,24.958,26.394,27.927,29.393,31.202,33.526,36.218,39.459,43.011,47.080,52.325,57.551,64.439,74.842,86.671,101.224,120.066,145.429,185.306,252.549,400,
     1.874,1.879,1.883,1.888,1.893,1.899,1.904,1.907,1.910,1.917,1.922,1.927,1.931,1.934,1.937,1.940,1.943,1.946,1.950,1.953,1.955,1.958,1.960,1.962,1.965,1.967,1.968,1.969,1.970,1.971,1.972,1.972,
     10.026,10.306,10.617,10.899,11.315,11.803,12.360,12.921,13.751,14.615,15.455,16.414,17.563,18.763,19.931,21.339,22.949,24.520,26.117,28.028,30.093,32.525,35.376,38.899,43.697,49.549,56.906,72.970,97.878,130.866,207.623,400,
     1.773,1.785,1.798,1.808,1.821,1.836,1.851,1.861,1.872,1.880,1.887,1.895,1.902,1.909,1.915,1.920,1.926,1.930,1.934,1.938,1.942,1.945,1.949,1.952,1.956,1.960,1.962,1.966,1.969,1.971,1.973,1.973];
    PYSa = [-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];

   % Factorul de durabilitate al materialului la concentratorul de tensiuni de la baza dintelui Yd
    GYdelta = [1.371,2.000,
     0.953,1.030,
     1.369,2.000,
     0.961,1.026,
     1.366,2.000,
     0.977,1.016];
    PYdelta = [500,600,800];

   % Factorul relativ de sensibilitate al materialului la concentratorul de tensiune de la baza dintelui la solicitarea statica Ydelta_st
    GX = [1.400,2.000,
     1.483,2.000,
     1.400,2.000,
     1.527,2.109,
     1.400,2.000,
     1.572,2.247,
     1.400,2.000,
     1.604,2.382,
     1.400,2.000,
     1.637,2.528,
     1.400,1.958,
     1.678,2.598];
    PX = [1,1.2,1.4,1.6,1.8,2];	

    GYdelta_st = [1.483,2.600,
     0.762,1.298,
     1.483,2.600,
     0.767,1.287,
     1.483,2.600,
     0.776,1.278];
    PYdelta_st = [500,600,800];	

    % Factorul de marime Yx
    GYx = [0.000,5.018,29.935,44.935,
     1.000,1.000,0.852,0.852];
    PYx(1) = [1];

    % Factorul dinamic Kvalfa
    GKvalfa = [0.,6.018,
     1.,1.800];
    PKvalfa(1) = [8];

    % Factorul dinamic Kvbeta
    GKvbeta = [0.000,6.643,
     1.000,1.446];
    PKvbeta(1) = [8];

   % Factorul de repartizare a sarcinii pe latimea danturii KHbeta pentru solicitarea de contact
    GKHbeta = [0.000,0.055,0.119,0.184,0.243,0.286,0.352,0.417,0.481,0.557,0.621,0.728,0.785,0.843,0.907,0.972,1.055,1.138,1.217,1.299,1.366,1.454,1.527,1.593,1.658,1.717,1.773,1.829,1.877,1.926,1.967,2.000,
     1.000,1.005,1.007,1.012,1.017,1.021,1.024,1.030,1.034,1.041,1.046,1.055,1.062,1.066,1.073,1.080,1.088,1.097,1.106,1.117,1.126,1.138,1.149,1.161,1.171,1.182,1.193,1.203,1.215,1.226,1.235,1.244];
    PKHbeta(1) = [5];

    % Factorul de repartizare a sarcinii pe latimea danturii KFbeta pentru solicitarea de incovoiere
    GKFbeta = [0.000,0.105,0.183,0.280,0.345,0.407,0.469,0.542,0.621,0.692,0.767,0.857,0.964,1.063,1.155,1.244,1.310,1.392,1.458,1.528,1.561,1.595,1.655,1.688,1.751,1.786,1.822,1.852,1.895,1.929,1.965,2.001,
     1.000,1.015,1.025,1.034,1.042,1.050,1.060,1.070,1.081,1.094,1.108,1.132,1.156,1.179,1.204,1.231,1.255,1.283,1.306,1.340,1.354,1.370,1.396,1.413,1.443,1.459,1.480,1.499,1.522,1.540,1.559,1.579];
    PKFbeta(1) = [5];


    % Calculul de rezistenta

    validHDS = 1;
    while(validHDS)
        Rzf1=exp(log(4.4)+0.97*log(Raf1)); % Rugozitatea Rz a flancului dintelui pinionului, micrometri
        Rzf2=exp(log(4.4)+0.97*log(Raf2)); % Rugozitatea Rz a flancului dintelui rotii, micrometri
        Rz100=(Rzf1+Rzf2)*sqrt(100/aw)/2; % Rugozitatea ehivalenta, micrometri

        % Factorii rugozitatii flancurilor pentru solicitarea de contact
        ZR1=grafparam(sigma_Hlim1, Rz100);
        if (ZR1==0) validHDS=0; break;
        end
        ZR2=grafparam(sigma_Hlim2, Rz100);
        if(ZR2==0) validHDS=0; break;
        end
        v1=PI*d1*n1/60000.;
        if(v1<1) v1=1; % Viteza periferica pe cercul de divizare al pinionului, m/s
        end
        v2=PI*d2*n2/60000.;
        if(v2<1) v2=1; % Viteza periferica pe cercul de divizare al pinionului, m/s
        end
        % Factorii de viteza pentru solicitarea de contact
        Zv1=grafparam(sigma_Hlim1, v1);
        if (Zv1==0) validHDS=0; break;
        end
        Zv2=grafparam(sigma_Hlim2, v2);
        if (Zv2==0) validHDS=0; break;
        end
        % Factorii de forma pentru solicitarea de incovoiere
        YFa1=grafparam(xn1, zn1);
        if ((YFa1<1.7)||(YFa1>3.7)) validHDS=0; break;
        end
        YFa2=grafparam(xn2, zn2);
        if ((YFa2<1.7)||(YFa2>3.7)) validHDS=0; break;
        end
        % Factorii de corectie a tensiunilor de incovoiere la baza dintelui
        YSa1=grafparam(xn1, zn1);
        if ((YSa1<1.2)||(YSa1>2)) validHDS=0; break;
        end
        YSa2=grafparam(xn2, zn2);
        if ((YSa2<1.2)||(YSa2>2)) validHDS=0; break;
        end
        % Factorii relativi de sensibilitate a materialului la concentratorul de tensiuni de la baza dintelui, la durabilitate nelimitata
        Ydelta1=grafparam(sigma021, YSa1);
        if (Ydelta1==0) validHDS=0; break;
        end
        Ydelta2=grafparam(sigma022, YSa2);
        if (Ydelta2==0) validHDS=0; break;
        end
        % Factorii relativi de sensibilitate a materialului la concentratorul de tensiuni de la baza dintelui, la solicitarea statica
        XX1=grafparam(eps_alfa_n, YSa1);
        if (XX1==0) validHDS=0; break;
        end
        XX2=grafparam(eps_alfa_n, YSa2);
        if (XX2==0) validHDS=0; break;
        end
        Ydelta_st1=grafparam(sigma021, XX1);
        if (Ydelta_st1==0) validHDS=0; break;
        end
        Ydelta_st2=grafparam(sigma022, XX2); 
        if (Ydelta_st2==0) validHDS=0; break;
        end
        % Factorii de marime pentru solicitarea de contact si de incovoiere
        Yx1=grafparam(MatTT1, mn);
        if (Yx1==0) validHDS=0; break;
        end
        Yx2=grafparam(MatTT2, mn);
        if (Yx2==0) validHDS=0; break;
        end
        % Factorii durabilitatii pentru solicitarea de contact
        NL1=60*n1*Lh1*Hi1; % Numarul de cicluri de solicitare al pinionului
        mH1=log(NBH/NstH)/log(ZNmax/ZL1/ZR1/Zv1);
        mF1=log(NBF/NstF)/log(YNmax*Ydelta_st1/Ydelta1/YR1/Yx1); 
        ZN1 =1;
        if (NL1<=NstH) ZN1 = ZNmax;
        end
        if ((NstH<NL1)&&(NL1<NBH)) ZN1=exp(log(NBH/NL1)/mH1);
        end
        NL2=60*n2*Lh2*Hi2; % Numarul de cicluri de solicitare al rotii dintate  
        mH2=log(NBH/NstH)/log(ZNmax/ZL2/ZR2/Zv2); 
        mF2=log(NBF/NstF)/log(YNmax*Ydelta_st2/Ydelta2/YR2/Yx2); 
        ZN2 =1; 
        if (NL2<=NstH) ZN2 =ZNmax;
        end
        if ((NstH<NL2)&&(NL2<NBH)) ZN2=exp(log(NBH/NL2)/mH2);
        end
        % Factorii durabilitatii pentru solicitarea de incovoiere
        YN1 =1; 
        if (NL1<=NstF) YN1 =YNmax;
        end
        if ((NstF<NL1)&&(NL1<NBF)) YN1=exp(log(NBF/NL1)/mF1);
        end
        YN2 =1; 
        if (NL2<=NstF) YN2 =YNmax;
        end
        if ((NstF<NL2)&&(NL2<NBF)) YN2=exp(log(NBF/NL2)/mF2);
        end

        sigma_HP1= sigma_Hlim1*ZN1*Zw*ZL1*ZR1*Zv1*ZX/SHmin; % Tensiunea hertziana admisibila pt materialul pinionului, MPa
        sigma_HP2= sigma_Hlim2*ZN2*Zw*ZL2*ZR2*Zv2*ZX/SHmin; % Tensiunea hertziana admisibila pt materialul rotii, MPa

        sigma_HP=sigma_HP1; % Tensiunea hertziana admisibila, MPa
        if (sigma_HP1>sigma_HP2) sigma_HP=sigma_HP2;
        end
        sigma_FP1=sigma_Flim1*YN1*Ydelta1*YR1*Yx1/SFmin; % Tensiunea de incovoiere admisibila pt materialul pinionului, MPa
        sigma_FP2=sigma_Flim2*YN2*Ydelta2*YR2*Yx2/SFmin; % Tensiunea de incovoiere admisibila pt materialul rotii, MPa

        % Factorul gradului de acoperire pentru solicitarea de contact
        Zeps=sqrt((4-eps_beta)*(1-eps_beta)/3+eps_beta/eps_alfa); 
        if (eps_beta>=1) Zeps=sqrt(1/eps_alfa);
        end
        % Factorul gradului de acoperire pentru solicitarea de incovoiere
        Yeps=0.25+0.75/eps_alfa_n;  

        % Factorul zonei de contact pentru solicitarea de contact
        ZH=sqrt(2*cos(beta_b)/sin(alfa_wt)/cos(alfa_wt)); 

        % Factorul inclinarii dintilor pentru solicitarea de contact
        Zbeta=sqrt(cos(beta)); 

        % Factorii dinamici
        v1z1=v1*zn1*cos(beta)*cos(beta)*cos(beta)/100.; % In loc de z1 am scris zn1*cos(beta)*cos(beta)*cos(beta)
        Kvalfa=grafparam(GKvalfa, 2, 2, PKvalfa, CPrec, v1z1);
        if (Kvalfa==0) validHDS=0; break;
        end
        Kvbeta=grafparam(GKvbeta, 2, 2, PKvbeta, CPrec, v1z1);
        if (Kvbeta==0) validHDS=0; break;
        end
        Kv=Kvbeta; 
        if (eps_beta<1) Kv=Kvbeta-eps_beta*(Kvbeta-Kvalfa);
        end
        % Factorii de repartizare a sarcini pe latimea danturii pentru solicitarea de contact, respectiv de incovoiere
        psi_d=b2/d1;
        KHbeta=grafparam(GKHbeta, 2, 32, PKHbeta, PozAngr, psi_d);
        if (KHbeta==0) validHDS=0; break;
        end
        KFbeta=grafparam(GKFbeta, 2, 32, PKFbeta, PozAngr, psi_d);
        if (KFbeta==0) validHDS=0; break;
        end
        % Abaterea efectiva a pasului de baza
        fpbr=0;
        if ((CPrec==7)&&((1<=mn)&&(mn<3.5))&&(d2<=125)) fpbr=13;
        end
        if ((CPrec==7)&&((1<=mn)&&(mn<3.5))&&((d2>125)&&(d2<=400))) fpbr=15;
        end
        if ((CPrec==7)&&((1<=mn)&&(mn<3.5))&&((d2>400)&&(d2<=800))) fpbr=17;
        end

        if ((CPrec==7)&&((3.5<=mn)&&(mn<6.3))&&(d2<=125)) fpbr=17;
        end
        if ((CPrec==7)&&((3.5<=mn)&&(mn<6.3))&&((d2>125)&&(d2<=400))) fpbr=19;
        end
        if ((CPrec==7)&&((3.5<=mn)&&(mn<6.3))&&((d2>400)&&(d2<=800))) fpbr=19;
        end

        if ((CPrec==7)&&((6.3<=mn)&&(mn<=10))&&(d2<=125)) fpbr=19;
        end
        if ((CPrec==7)&&((6.3<=mn)&&(mn<=10))&&((d2>125)&&(d2<=400))) fpbr=21;
        end
        if ((CPrec==7)&&((6.3<=mn)&&(mn<=10))&&((d2>400)&&(d2<=800))) fpbr=24;
        end

        if ((CPrec==8)&&((1<=mn)&&(mn<3.5))&&(d2<=125)) fpbr=19;
        end
        if ((CPrec==8)&&((1<=mn)&&(mn<3.5))&&((d2>125)&&(d2<=400))) fpbr=21;
        end
        if ((CPrec==8)&&((1<=mn)&&(mn<3.5))&&((d2>400)&&(d2<=800))) fpbr=24;
        end        

        if ((CPrec==8)&&((3.5<=mn)&&(mn<6.3))&&(d2<=125)) fpbr=24; 
        end
        if ((CPrec==8)&&((3.5<=mn)&&(mn<6.3))&&((d2>125)&&(d2<=400))) fpbr=26;
        end
        if ((CPrec==8)&&((3.5<=mn)&&(mn<6.3))&&((d2>400)&&(d2<=800))) fpbr=26;
        end

        if ((CPrec==8)&&((6.3<=mn)&&(mn<=10))&&(d2<=125)) fpbr=26;  
        end
        if ((CPrec==8)&&((6.3<=mn)&&(mn<=10))&&((d2>125)&&(d2<=400))) fpbr=30;  
        end
        if ((CPrec==8)&&((6.3<=mn)&&(mn<=10))&&((d2>400)&&(d2<=800))) fpbr=34; 
        end

        if ((CPrec==9)&&((1<=mn)&&(mn<3.5))&&(d2<=125)) fpbr=26;    
        end
        if ((CPrec==9)&&((1<=mn)&&(mn<3.5))&&((d2>125)&&(d2<=400))) fpbr=30;    
        end
        if ((CPrec==9)&&((1<=mn)&&(mn<3.5))&&((d2>400)&&(d2<=800))) fpbr=34;   
        end

        if ((CPrec==9)&&((3.5<=mn)&&(mn<6.3))&&(d2<=125)) fpbr=34;
        end
        if ((CPrec==9)&&((3.5<=mn)&&(mn<6.3))&&((d2>125)&&(d2<=400))) fpbr=38;
        end
        if ((CPrec==9)&&((3.5<=mn)&&(mn<6.3))&&((d2>400)&&(d2<=800))) fpbr=38;
        end

        if ((CPrec==9)&&((6.3<=mn)&&(mn<=10))&&(d2<=125)) fpbr=32;
        end
        if ((CPrec==9)&&((6.3<=mn)&&(mn<=10))&&((d2>125)&&(d2<=400))) fpbr=42;
        end
        if ((CPrec==9)&&((6.3<=mn)&&(mn<=10))&&((d2>400)&&(d2<=800))) fpbr=48;
        end

        Ft1=2*T1/d1; % Forta tangentiala corespunzatoare diametrului de divizare, N
        qalfa=1; % Factorul auxiliar
        if (4*(0.1+b2*(fpbr-4)/Ft1)<=0.5) qalfa=0.5;
        end

        % Factorul de repartizare a sarcinii in plan frontal pe perechile de dintii aflate simultan in angrenare, pentru solicitarea de contact
        KHalfa=1+2*(qalfa-0.5)*(1/Zeps/Zeps-1); 

        % Factorul de repartizare a sarcinii in plan frontal pe perechile de dintii aflate simultan in angrenare, pentru solicitarea de incovoiere
        KFalfa=qalfa*eps_alfa; 

        % Factorul minim al inclinarii dintilor pentu solicitarea de incovoiere
        Ybeta_min=0.75; 
        if (eps_beta<=1) Ybeta_min=1-0.25*eps_beta;
        end

        % Factorul inclinarii dintilor pentu solicitarea de incovoiere
        Ybeta=1-eps_beta*(beta*180/PI)/120;
        if(Ybeta<Ybeta_min)Ybeta=Ybeta_min; 
        end

        % Factorul de elasticitate al materialului
        ZE=sqrt( 1/( (1-niu1*niu1)/E1+(1-niu2*niu2)/E2 )/PI );

        sigma_H=(u12+1)*ZE*Zeps*ZH*Zbeta*cos(alfa_t)*sqrt(T1*KA*Kv*KHbeta*KHalfa*(u12+1)/b2/u12/2)/aw/cos(alfa_wt); % Tensiunea hertziana

        sigma_F1=T1*z1*(u12+1)*(u12+1)*KA*Kv*KFbeta*KFalfa*Yeps*Ybeta*YFa1*YSa1*cos(alfa_t)*cos(alfa_t)/cos(alfa_wt)/cos(alfa_wt)/cos(beta)/aw/aw/b1/2; % Tensiunea de incovoiere la baza dintelui pinionului, MPa
        sigma_F2=sigma_F1*b1*YFa2*YSa2/b2/YFa1/YSa1; % Tensiunea de incovoiere la baza dintelui rotii, MPa

        break;
    end

    if(validHDS)
        HDS(1)=1;
        HDS(2)=sigma_H;
        HDS(3)=sigma_HP;
        HDS(4)=sigma_F1;
        HDS(5)=sigma_F2;
        HDS(6)=sigma_FP1;
        HDS(7)=sigma_FP2;
    else HDS(1)=0;
        result = HDS;
end

