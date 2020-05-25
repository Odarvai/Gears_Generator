function [HG] = HelicalGear(alfa_n, ha_n, cs_n, mn, z, xn, beta, alfa_t, Delta_yn, HG)
% Returneaza in HG geomatria unei r.d.c.d.i.
		% ha_n		- Coeficientul inaltimii capului dintelui
		% cs_n		- Coeficientul jocului la capul dintelui de referinta
		% mn		- Modulul normal, mm (in cazul r.d.c.c.d.d. modulul normal este modulul m al rotii)
		% z		- Numarul de dinti
		% xn		- Coeficientul deplasarii de profil in plan normal, (in cazul r.d.c.c.d.d. xn este x-ul rotii)
		% beta		- Unghiul de inclinare al danturii pe cilindrul de divizare, grade
	    % Delta_yn - Coeficientul de scurtare al dintilor in plan normal
    validHG=1;
    digits(100);
    p = vpa(pi);  % am crescut numarul de zecimale al numarului pi
    
    xt = xn * cos(beta); % Coeficientul deplasarii de profil al dintelui pinionului, in plan frontal
    d = mn * z / cos(beta); % Diametrul cercului de divizare, mm
    db = d * cos(alfa_t); % Diametrul cercului de baza, mm

    while validHG == 1
 
        df = mn * (z / cos(beta) - 2 * (ha_n + cs_n - xn)); % Diametrul cercului de picior, mm
        if (df <= 0) validHG = 0; break;
        end
        da = mn * (z / cos(beta) + 2 * (ha_n + xn - Delta_yn)); % Diametrul cercului de cap, mm
        if (da <= 0) validHG = 0; break;
        end
        temp = abs(d * cos(alfa_t) / da);
        if (temp > 1) validHG = 0; break;
        end
        alfa_at = alfa_yt(da, d, alfa_t); % Unghiul de presiune de referinta pe cercul de cap, rad

        st = syt(mn, z, alfa_t, alfa_t, beta, xn); % Lungimea arcului dintelui pe cercul de divizare, in plan frontal, mm
        if (st <= 0) validHG = 0; break;
        end

        sat = syt(mn, z, alfa_t, alfa_at, beta, xn); %Lungimea arcul dintelui pe cercul de cap in plan frontal, mm 
        if (sat <= 0) validHG = 0; break;
        end
        beta_a = beta_y(da, d, beta); % Unghiul de inclinare a danturii pe cilindrul de cap, rad

        beta_b = beta_y(db, d, beta); % Unghiul de inclinare a danturii pe cilindrul de baza, rad

        sn = syn(st, beta); %Lungimea arcului dintelui pe cercul de divizare, in plan normal, mm
        if (sn <= 0) validHG = 0; break;
        end
        
        san = syn(sat, beta_a); % Lungimea arcul dintelui pe cercul de cap in plan normal, mm
        if (san <= 0) validHG = 0; break;
        end
        break;
    end
    
    if(validHG)

        HG(1) = 1;          HG(2) = xt;          HG(3) = d;      HG(4) = db;     HG(5) = df;
        HG(6) = da;         HG(7) = alfa_at;     HG(8) = st;	 HG(9) = sat;    HG(10) = beta_a;
        HG(11) = beta_b;    HG(12) = sn;         HG(13) = san;	 HG(14) = ArieDinte(da, db, df, xn, z, sat, alfa_n);
    
    else HG(1) = 0;
    end
end

