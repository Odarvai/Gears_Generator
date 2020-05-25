function volum = Volum_saiba(De, Di, l)
%Returneaza volumul unei saibe, mm^3
%De - Diametrul exterior, mm
%Di - Diametrul interior, mm
%l  - Grosimea saibei, mm
digits(100);
p = vpa(pi);
volum = (De * De - Di * Di) * p * l / 4;
end

