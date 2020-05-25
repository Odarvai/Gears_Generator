function [result] = ControlHelicalDrive(mn, aw, z1, z2, xn1, xn2, beta, beta_b, beta_w, alfa_n, alfa_t, alfa_wt, alfa_at1, alfa_at2, d1, da1, d2, db1, db2, da2, CHD)
% Returneaza elementele de control ale angrenajului cu r.d.c.d.i
		% mn     - Modulul normal, mm (in cazul r.d.c.c.d.d. modulul normal este modulul m)
		% aw	  - Distanta axiala standardizata/impusa, mm
		% z1,z2  - Numarele de dinti
		% beta   - Unghiul de inclinare al danturii pe cilindrul de divizare, grade
		% beta_b - Unghiul de inclinare al danturii pe cilindrul de baza, rad
		% beta_w - Unghiul de inclinare al danturii pe cilindrul de rostogolire, rad
		% alfa_n - Unghiul de presiune de referinta in plan normal, grade
		% alfa_wt- Unghiul real de angrenare in plan frontal, rad
		% d1     - Diametrul cercului de divizare corespunzator pinionului, mm
		% da1	  - Diametrul cercului de cap corespunzator pinionului, mm
		% d2     - Diametrul cercului de divizare corespunzator rotii dintate, mm
		% da2	  - Diametrul cercului de cap corespunzator rotii dintate, mm
		
		digits(100);
        PI = vpa(pi);

		validCHD = 1;
		
		while(validCHD == 1)
            temp = abs(z1 * cos(alfa_t) / (z1 + 2 * xn1 * cos(beta)));
            if (temp > 1) validCHD = 0; break;
            end
            alfa_Nt1 = acos(z1 * cos(alfa_t) / (z1 + 2 * xn1 * cos(beta))); % rad
            temp = abs(z2 * cos(alfa_t) / (z2 + 2 * xn2 * cos(beta)));
            if (temp > 1) validCHD = 0; break;
            end
            alfa_Nt2 = acos(z2 * cos(alfa_t) / (z2 + 2 * xn2 * cos(beta))); % rad
            Ncalc1 = 0.5 + z1 * (tan(alfa_Nt1) / cos(beta) / cos(beta) - 2 * xn1 * tan(alfa_n) / z1 - involut(alfa_t)) / PI;
            Ncalc2 = 0.5 + z2 * (tan(alfa_Nt2) / cos(beta) / cos(beta) - 2 * xn2 * tan(alfa_n) / z2 - involut(alfa_t)) / PI;
            N1 = floor(Ncalc1); % Numarul de dinti pentru masurarea cotei peste dinti, pentru pinion
            if (Ncalc1 - floor(Ncalc1) >= 0.5) N1 = N1 + 1;
            end
                N2 = floor(Ncalc2); % Numarul de dinti pentru masurarea cotei peste dinti, pentru roata
            if (Ncalc2 - floor(Ncalc2) >= 0.5) N2 = N2 + 1;
            end
            WNn1 = 2 * xn1 * mn * sin(alfa_n) + mn * cos(alfa_n) * ((N1-0.5) * PI + z1 * involut(alfa_t)); % Cota peste N1 dinti in plan normal pt angrenaje fara joc intre flancuri, pentru pinion, mm
            WNn2 = 2 * xn2 * mn * sin(alfa_n) + mn * cos(alfa_n) * ((N2-0.5) * PI + z2 * involut(alfa_t)); % Cota peste N2 dinti in plan normal pt angrenaje fara joc intre flancuri, pentru roata, mm
            WNt1 = WNn1 / cos(beta_b); % Cota peste N1 dinti in plan frontal pt angrenaje fara joc intre flancuri, pentru pinion, mm
            WNt2 = WNn2 / cos(beta_b); % Cota peste N2 dinti in plan frontal pt angrenaje fara joc intre flancuri, pentru roata, mm
            roNt1 = 0.5 * WNt1; % Raza de curbura a profilului in punctele simetrice de masurare a lungimii peste N1 dinti in planul frontal, pentru pinion, mm
            roNt2 = 0.5 * WNt2; % Raza de curbura a profilului in punctele simetrice de masurare a lungimii peste N2 dinti in planul frontal, pentru roata, mm
            roAt1 = aw * sin(alfa_wt) - 0.5 * db2 * tan(alfa_at2); % Raza de curbura a profilului dintelui pinionului in punctul de intrare din angrenare, mm
            if (roAt1 <= 0) validCHD = 0; break;
            end
            roEt2 = aw * sin(alfa_wt) - 0.5 * db1 * tan(alfa_at1); % Raza de curbura a profilului dintelui rotii in punctul de iesire din angrenare, mm
            if (roEt2 <= 0) validCHD = 0; break;
            end
            roat1 = 0.5 * da1 * sin(alfa_at1); % Raza de curbura a profilului la capul dintelui pinionului, mm
            roat2 = 0.5 * da2 * sin(alfa_at2); % Raza de curbura a profilului la capul dintelui rotii, mm
            break;
        end

		if(validCHD)
            CHD(1) = 1;
            CHD(2) = alfa_Nt1;      CHD(3) = alfa_Nt2;      CHD(4) = Ncalc1;    CHD(5) = Ncalc2; 
            CHD(6) = N1;            CHD(7) = N2;            CHD(8) = WNn1;      CHD(9) = WNn2; 
            CHD(10) = WNt1;         CHD(11) = WNt2;         CHD(12) = roNt1;    CHD(13) = roNt2; 
            CHD(14) = roAt1;        CHD(15) = roEt2;        CHD(16) = roat1;    CHD(17) = roat2;

		else CHD(1) = 0;
        end
        result = CHD;
end

