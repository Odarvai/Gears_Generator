function [abscisa] = absEv(db, alfa)
%Returneaza abscisa unui punct de pe evolventa cu diametrul cercului
	% de baza db [mm] si determinat de unghiul alfa [rad]
    abscisa = db * (sin(tan(alfa)) - tan(alfa) * cos(tan(alfa))) / 2;
end

