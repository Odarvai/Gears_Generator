function [ordonata] = ordEv(db, alfa)
% Returneaza ordonata unui punct de pe evolventa cu diametrulcercului
		% de baza db [mm] si determinat de unghiul alfa [rad]
    ordonata = db * (cos(tan(alfa)) + tan(alfa) * sin(tan(alfa))) / 2;
end

