function result = Volum_calota_sferica(R, h)
%Returneaza volumul unui calote sferice sectionate la distanta h de centru, mm^3
%R  - Raza sferei, mm
%h  - distanta de la centrul sferei la planul de sectiune, mm
digits(100);
p = vpa(pi);
result = p * (R - h) * (R - h) * (2 * R + h) / 3;
end

