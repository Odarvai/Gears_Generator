function [arie] = ArieDinte(da, db, df, xn, z, sat, alfa_n)
% Returneaza aria unui dinte cu profil evolventice
		% A-inceputul evolventei (deasupra racordarii)
		% B-sfarsitul evolventei (pe capul dintelui)
    digits(100);
    p = vpa(pi);    %am crescut numarul de zecimale a variabilei pi din matlab
   
    alfa_a = atan(sqrt(da * da - db * db) / db); %Unghiul de sfarsit al evolventei, rad
    alfa_u = atan((z * sin(alfa_n) * sin(alfa_n) - 2 * (1 - xn)) / z / cos(alfa_n) / sin(alfa_n)); %Unghiul punctului de inceput al evolventei, rad
    xA = absEv(db, alfa_u);
    yA = ordEv(db, alfa_u);
    xB = absEv(db, alfa_a);
    yB = ordEv(db, alfa_a);
    du = 2 * sqrt(xA * xA + yA * yA); %Diametrul cercului de inceput al evolventei, mm
    psi = 2 * sat / da;
    alfa_AB = acos((xA * xB + yA * yB) / sqrt(xA * xA + yA * yA) / sqrt(xB * xB + yB * yB));
    AOCA = xA * yA / 2;
    AODB = xB * yB / 2;
    AOAM = ASect(du, alfa_AB);
    AOMN = ASect(du, psi);
    AOBE = ASect(da, psi);
    AOAF = ASect(du, 2 * alfa_AB + psi);
    AOHK = ASect(df, 2 * alfa_AB + psi);
    ABMNE = AOBE - AOMN;
    AAHKF = AOAF - AOHK;
    ror = (du * du - df * df) / df / 4; % Raza aprox. a racordarii, mm
    AOAo = ror * du / 4;
    fi = asin(2 * ror / (df + 2 * ror));
    AOTH = ASect(df, fi);
    AoTA = ASect(2 * ror, p / 2 - fi);
    AR = AOAo - AOTH - AoTA; %Aria aprox a zonei de racord a dintelui
    Ae = db * db * (primEv(tan(alfa_a))-primEv(tan(alfa_u))) / 48; % Aria marginita de evolventa
	Aev = Ae + AOCA;
	AAMB = Aev - AODB - AOAM;
    
    arie = 2 * AAMB + ABMNE + AAHKF + 2 * AR;
end

