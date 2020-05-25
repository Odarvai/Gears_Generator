function [ Wz ] = modululderezistentaWz( ttr, d, b, t1 )
%modululderezistentaWz - returneaza Wz al unei sectiuni circulare cu sau fara pana

%ttr = 0 - sectiune circulara (fara canal de pana) de diametru d
%ttr = 1 - sectiune circulara de diametru d cu canal de pana bxt1 

if (ttr == 0) Wz = pi * d * d * d / 32;
else Wz=pi * d * d * d / 32 - (b * t1 * (d - t1) * (d - t1) / 2 / d);
end

end

