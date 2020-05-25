function volum = Volum_tr_con_gaurit(D, d, di, h)
% Returneaza volumul unui trunchi de con - din care se scade 
% volumul unui cilindru de diametru di si inaltime h, mm^2
% D  - Diametrul bazei mari al trunchiului de con, mm
% d  - Diametrul bazei mici al trunchiului de con, mm
% di - Diametrul interior al cilindrului, mm
digits(100);
p = vpa(pi);
volum = p * ((D * D + d * d + D * d) / 12 - di * di / 4) * h; 
end

