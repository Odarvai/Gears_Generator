function [result] = EquivalentHelicalDrive(mn, a, z1, z2, beta, beta_b, beta_w, alfa_n, alfa_wt, d1, da1, d2, da2, EHD)
% Returneaza elementele geometrice ale angrenajului echivalent cu r.d.c.c.d. al unui angrenaj cu r.d.c.d.i. echivalarea in plan normal
		% mn     - Modulul normal, mm (in cazul r.d.c.c.d.d. modulul normal este modulul m)
		% a	  - Distanta axiala standardizata/impusa, mm
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

		alfa_n_g = 20; % Unghiul de presiune de referinta in plan normal, grade
        
		alfa_n = alfa_n_g * PI / 180;
		zn1 = z1 / cos(beta) / cos(beta) / cos(beta); % Numarul de dinti ai pinionului echivalent
		zn2 = z2 / cos(beta) / cos(beta) / cos(beta); % Numarul de dinti ai rotii echivalente
		xn1min = (14 - zn1) / 17; % Valoarea minima a coeficientului deplasarii de profil a dintelui pinionului, in plan normal
		xn2min = (14 - zn2) / 17; % Valoarea minima a coeficientului deplasarii de profil a dintelui rotii, in plan normal
		dn1 = mn * zn1; % Diametrul cercului de divizare al pinionului echivalent, mm
		dn2 = mn * zn2; % Diametrul cercului de divizare al rotii echivalente, mm
		dbn1 = dn1 * cos(alfa_n); % Diametrul cercului de baza all pinionului echivalent, mm
		dbn2 = dn2 * cos(alfa_n); % Diametrul cercului de baza a rotii echivalente, mm
		dan1 = dn1 + da1 - d1; % Diametrul cercului de cap al pinionului echivalent, mm
		dan2 = dn2 + da2 - d2; % Diametrul cercului de cap al rotii echivalente, mm
		alfa_wn = acos(cos(alfa_wt) * cos(beta_b) / cos(beta_w)); % Unghiul de presiune al angrenajului echivalent, rad
		awn = a * cos(alfa_n) / cos(beta_b) / cos(beta_b) / cos(alfa_wn); % Distanta dintre axe a angrenajului echivalent, mm
		eps_alfa_n = (sqrt(dan1 * dan1 - dbn1 * dbn1) + sqrt(dan2 * dan2 - dbn2 * dbn2) - 2 * awn * sin(alfa_wn)) / cos(alfa_n) / mn / PI / 2; % Gradul de acoperie frontal al angrenajului echivalent
		
		EHD(1) = 1;
		EHD(2) = zn1;           EHD(3)=zn2;     EHD(4)=xn1min;          EHD(5) = xn2min;      EHD(6) = dn1; 
		EHD(7) = dn2;           EHD(8)=dbn1;    EHD(9)=dbn2;            EHD(10) = dan1;       EHD(11) = dan2; 
		EHD(12) = alfa_wn;      EHD(13)=awn;    EHD(14)=eps_alfa_n;
        result = EHD;
end

