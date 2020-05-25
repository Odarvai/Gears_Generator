function [TC] = Tangenta_Comuna(R1, R2, C, TC)
% Returneaza lungimea tangentei comune a doua cercuri si unghiul facut de aceasta cu linia centrelor, mm
% R1, R2 - Razele cercurilor, mm
% C	  - Distanta dintre centrele cercurilor, mm

temp = C * C - (R2 - R1) * (R2 - R1);

if ( temp >= 0 )
    tangenta_comuna = sqrt(temp);
    unghi = atan((R2-R1) / tangenta_comuna);
    TC(1) = 1;
    TC(2) = tangenta_comuna;
    TC(3) = unghi;
else
    TC(1) = 0;
end
end

