function result = Volum_cutie_paralelip(tip_cutie, Li, li, hi, gL, gl, gc, gf)
%Returneaza volumul unuei cutii peralelipipedice, mm^3
% tip_cutie = 1 - cutie fara capac si fara fund (gc=0, gf=0)
% tip_cutie = 2 - cutie fara capac, dar cu fund (gc=0)
% tip_cutie = 3 - cutie cu capac, dar fara fund (gf=0)
% tip_cutie = 4 - cutie cu capac si cu fund
% Li - Lungimea interioara a sectiunii, mm
% li - Latimea interioara a sectiunii, mm
% hi - Inaltimea interioara a sectiunii, mm
% gL - Grosimea peretelui cu dimensiunea interioara Li, mm
% gl - Grosimea peretelui cu dimensiunea interioara li, mm
% gc - Grosimea peretelui capacului, mm
% gf - Grosimea peretelui fundului, mm

%Asp, Vol_pereti, Vol_capac, Vol_fund;
Asp = 2 * (Li * gL + li * gl + 2 * gL * gl); %Aria sectiunii peretilor, mm^2
Vol_pereti = Asp * hi;
Vol_capac = (Asp + Li * li) * gc;
Vol_fund = (Asp + Li * li) * gf;
if (tip_cutie == 1) result = Vol_pereti;
end
if (tip_cutie == 2) result = Vol_pereti + Vol_fund;
end
if (tip_cutie == 3) result = Vol_pereti + Vol_capac;
end
if (tip_cutie == 4) result = Vol_pereti + Vol_fund + Vol_fund;
end
end

