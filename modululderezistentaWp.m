function [ Wp ] = modululderezistentaWp( ttr, d, b, t1 )
%modululderezistentaWp - returneaza Wp al unei sectiuni circulare cu sau fara pana

%ttr = 0 - sectiune circulara (fara canal de pana) de diametru d
%ttr = 1 - sectiune circulara de diametru d cu canal de pana bxt1 

if (ttr == 0) Wp = 2 * modululderezistentaWz(ttr, d, b, t1);
else Wp = pi * d * d * d / 16 - (b * t1 * (d - t1) * (d - t1) / 2 / d);    
end

end