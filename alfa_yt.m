function [unghi_de_presiune] = alfa_yt(dy, d, alfa_t)
% Returneaza valoarea unghiului de presiune al profilului pe un cerc de diametru oarecare in plan frontal, rad
		% dy - cercul de diametrul oarecare, mm
		% d  - diametrul de divizare, mm
    unghi_de_presiune = acos(d * cos(alfa_t) / dy);
end

