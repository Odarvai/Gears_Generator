function [arie] = Arie_supracoarda(R, h)
% Returneaza aria supracorzii, mm^2
	% r - Raza sectorului de cerc, mm
	% h - distanta de la centru la coarda, mm

    if(h >= R) arie = 0;
    end
    alfa = acos(h / R);
    Asector = alfa * R * R; % Aria sectorului de cerc, rad
    Atriunghi = h * R * sin(alfa); % Aria triunghiului, mm^2
    arie = Asector - Atriunghi;
end

