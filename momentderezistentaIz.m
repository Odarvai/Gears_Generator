function [ Iz ] = momentderezistentaIz( ttr, d, b, t1 )
%momentderezistentaIz - returneaza Iz al unei sectiuni circulare cu sau fara pana

%ttr = 0 - sectiune circulara (fara canal de pana) de diametru d
%ttr = 1 - sectiune circulara de diametru d cu canal de pana bxt1 

if (ttr == 0) Iz = pi * d * d * d * d / 64;
else Iz = pi * d * d * d * d / 64 - (b * t1 * (d - t1) * (d - t1) / 4);    
end

end