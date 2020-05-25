function [result] = HelicalDriveGeometry(mn, aw, z1, z2, xn1, beta, b2, Delta_b, HDG)
% Returneaza in HDG geomatria unui angrenaj cu r.d.c.d.i.
		% mn      - Modulul normal, mm (in cazul r.d.c.c.d.d. modulul normal este modulul m)
		% aw      - Distanta axiala standardizata/impusa, mm
		% z1,z2   - Numarele de dinti
		% xn1     - Coeficientul deplasarii de profil in plan normal, pentru pinion (in cazul r.d.c.c.d.d. xn1 este x1)
		% beta    - Unghiul de inclinare al danturii pe cilindrul de divizare, grade
		% b2	   - Latimea danturii rotii, mm
		% Delta_b - Diferenta de latime dintre latimea pinionului si a rotii, mm (Delta_b=5)
	
		% Cremaliera de referinta ISO53 (STAS 821)
        
        alfa_n_g = 20; % Unghiul de presiune de referinta in plan normal, grade
		ha_n = 1;		% Coeficientul inaltimii capului dintelui
		cs_n = 0.25;	% Coeficientul jocului la capul dintelui de referinta
        
        digits(100);
        PI = vpa(pi); % crestem numarul de zecimale al variabilei pi
        
        alfa_n = alfa_n_g * PI / 180;
		alfa_t = atan(tan(alfa_n) / cos(beta)); % Unghiul de angrenare de referinta in plan frontal, rad
		validHDG = 1;
		b1 = b2 + Delta_b; % Latimea pinionului, mm
		while (validHDG == 1)
        
            a = mn * (z1 + z2) / cos(beta) / 2; % Distanta axiala elementara, mm
            if (((aw - a) * cos(beta) / mn < -0.4) || ((aw - a) * cos(beta) / mn > 2.5))    validHDG = 0; break;
                temp = abs(a * cos(alfa_t) / aw);
            end
            if (temp > 1) validHDG = 0; break;
            end
            alfa_wt = acos(a * cos(alfa_t) / aw); % Unghiul real  de  angrenare  in  plan  frontal, rad

            xsn = (involut(alfa_wt) - involut(alfa_t)) * (z1 + z2) / tan(alfa_n) / 2; % Suma coeficientilor  deplasarilor  de  profil  in  plan  normal
            xn2 = xsn - xn1; % Coeficientul deplasarii de profil al dintelui rotii, in plan normal

            Delta_yn = xsn - (aw - a) / mn;

            HelicalGear(alfa_n, ha_n, cs_n, mn,z1, xn1, beta, alfa_t, Delta_yn, HG1);
            if(HG1(1) == 0)    validHDG=0; break;
            end
            xt1 = HG1(2); % Coeficientul deplasarii de profil al dintelui pinionului, in plan frontal
            d1 = HG1(3); % Diametrul cercului de divizare pentru pinion, mm
            db1 = HG1(4); % Diametrul cercului de baza pentru pinion, mm
            df1 = HG1(5); % Diametrul cercului de picior pentru pinion, mm
            da1 = HG1(6); % Diametrul cercului de cap pentru pinion, mm
            alfa_at1 = HG1(7); % Unghiul de presiune de referinta pe cercul de cap al pinionului, rad
            st1 = HG1(8); % Lungimea arcului dintelui pe cercul de divizare al pinionului, in plan frontal, mm
            sat1 = HG1(9); % Lungimea arcul dintelui pe cercul de cap al pinionului in plan frontal, mm
            beta_a1 = HG1(10); % Unghiul de inclinare a danturii pe cilindrul de cap al pinionului, rad
            beta_b = HG1(11); % Unghiul de inclinare a danturii pe cilindrul de baza, rad
            sn1 = HG1(12); % Lungimea arcului dintelui pe cercul de divizare al pinionului, in plan normal, mm
            san1 = HG1(13); % Lungimea arcul dintelui pe cercul de cap al pinionului in plan normal, mm
            ADinte1 = HG1(14); % Aria suprafetei frontale a dintelui pinionului, mm^2

            HelicalGear(alfa_n, ha_n, cs_n, mn, z2, xn2, beta, alfa_t, Delta_yn, HG2);
            if(HG2(1) == 0)   validHDG=0; break;
            end
            xt2 = HG2(2); % Coeficientul deplasarii de profil al dintelui rotii, in plan frontal
            d2 = HG2(3); % Diametrul cercului de divizare pentru roata, mm
            db2 = HG2(4); % Diametrul cercului de baza pentru roata, mm
            df2 = HG2(5); % Diametrul cercului de picior pentru roata, mm
            da2 = HG2(6); % Diametrul cercului de cap pentru roata, mm
            alfa_at2 = HG2(7); % Unghiul de presiune de referinta pe cercul de cap al rotii, rad
            st2 = HG2(8); % Lungimea arcului dintelui pe cercul de divizare al rotii, in plan frontal, mm
            sat2 = HG2(9); % Lungimea arcul dintelui pe cercul de cap al rotii in plan frontal, mm
            beta_a2 = HG2(10); % Unghiul de inclinare a danturii pe cilindrul de cap al rotii, rad
            beta_b = HG2(11); % Unghiul de inclinare a danturii pe cilindrul de baza, rad
            sn2 = HG2(12); % Lungimea arcului dintelui pe cercul de divizare al rotii, in plan normal, mm
            san2 = HG2(13); % Lungimea arcul dintelui pe cercul de cap al rotii in plan normal, mm
            ADinte2 = HG2(14); % Aria suprafetei frontale a dintelui rotii, mm^2

            VolDantura1 = ADinte1 * z1 * b1;
            VolDantura2 = ADinte2 * z2 * b2;

            dw1 = d1 * cos(alfa_t) / cos(alfa_wt); % Diametrul cercului de rostogolire pentru pinion, mm
            dw2 = d2 * cos(alfa_t) / cos(alfa_wt); % Diametrul cercului de rostogolire pentru roata, mm

            beta_w = atan(dw1 * tan(beta) / d1); % Unghiul de inclinare a danturii pe cilindrul de rostogolire, rad

            eps_alfa = (sqrt(da1 * da1 - db1 * db1) + sqrt(da2 * da2 - db2 * db2) - 2 * aw * sin(alfa_wt)) * cos(beta) / cos(alfa_t) / mn / PI / 2; % Gradul de acoperie frontal
            if (eps_alfa <= 0) validHDG = 0; break;
            end
            eps_beta = b2 * sin(beta) / mn / PI; % Gradul de acoperire axial
            eps_gama = eps_alfa + eps_beta; % Gradul de acoperire total
            break;
        end
       
		if(validHDG)
		HDG(1) = 1;
		HDG(2) = alfa_t;        HDG(3) = alfa_at1;      HDG(4) = beta_a1;       HDG(5) = beta_b; 
		HDG(6) = d1;            HDG(7) = db1;           HDG(8) = df1;           HDG(9) = da1;
		HDG(10) = st1;          HDG(11) = sat1;         HDG(12) = sn1;          HDG(13) = san1;  
		HDG(14) = xt1;          HDG(15) = ADinte1;      HDG(16) = a;            HDG(17) = alfa_wt;
		HDG(18) = xsn;          HDG(19) = xn2;          HDG(20) = alfa_at2;     HDG(21) = beta_a2;
		HDG(22) = d2;           HDG(23) = db2;          HDG(24) = df2;          HDG(25) = da2;
		HDG(26) = st2;          HDG(27) = sat2;         HDG(28) = sn2;          HDG(29) = san2;
		HDG(30) = xt2;          HDG(31) = ADinte2;      HDG(32) = dw1;          HDG(33) = dw2;
		HDG(34) = beta_w;       HDG(35) = eps_alfa;     HDG(36) = eps_beta;     HDG(37) = eps_gama;
		HDG(38) = b1; 
		HDG(39) = cs_n;
		HDG(40) = alfa_n;
		HDG(41) = VolDantura1;
		HDG(42) = VolDantura2;
		
		else HDG(1) = 0;
        end
        result = HDG;
end

