function [Varc] = syt(mn, z, alfa_t, alfa_yt, beta, xn)
% Returneaza valoarea arcului dintelui pe un cerc oarecare in plan frontal
		% mn - modulul normal, mm
		% z  - numarul de dinti
		% 
		% st - lungimea arcului dintelui pe cercul de divizare in plan frontal
		% alfa_yt - unghiul de presiune al profilului pe un cerc de diametru oarecare in plan frontal, rad
		% beta - unghiul de inclinare al danturii pe cilindrul de divizare, rad
    digits(100)
    p = vpa(pi); % am crescut numarul de zecimale al variabilei pi din matlab
    st = (0.5 * p / cos(beta) + 2 * xn * tan(alfa_t)) * mn;
		
	Varc = ((involut(alfa_t) - involut(alfa_yt)) * mn * z / cos(beta) + st) * cos(alfa_t) / cos(alfa_yt);
end

