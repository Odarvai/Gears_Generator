function [unghi_de_inclinare] = beta_y(dy, d, beta)
		% Returneaza valoarea unghiului de inclinare al danturii pe un cilindru oarecare
		% dy - cilindrul de diametrul oarecare, mm
		% d  - diametrul de divizare, mm
		% beta - unghiul de inclinare al danturii pe cilindrul de divizare, rad
    unghi_de_inclinare = atan(dy * tan(beta) / d);
end

