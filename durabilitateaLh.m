function Lh = durabilitateaLh( HA, VA, HB, VB, Fa, Y, e, C, n, ff)
%durabilitateaLh - returenaza durabilitatea efectiva a rulmentilor unui arbore

%
%


FRA = sqrt(HA * HA + VA * VA);
FRB = sqrt(HA * HA + VB * VB);


FaAprim = 0.5 * FRA / Y; % Se introduce cu semn
FaBprim =-0.5 * FRB / Y; % Se introduce cu semn

Farb = FaAprim + Fa + FaBprim; % Forta rezultanta din arbore

if (Farb >= 0)
    FaA = FaAprim;
    FaB = Farb - FaBprim;
else
    FaA = -Farb + FaAprim;
    FaB = -FaBprim;
end

if (FaA / FRA <= e)
    PeA = FRA;
else
    PeA = 0.4 * FRA + Y * FaA;
end

if (FaB / FRB <= e)
    PeB = FRB;
else
    PeB = 0.4 * FRB + Y * FaB;
end

if (PeA > PeB)
    Pe = PeA;
else
    Pe = PeB;
end

Lh = ((exp(3.33*log(1000*C/ff/Pe)) * 10e05) / 60 / n);

end



