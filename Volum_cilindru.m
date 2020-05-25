function result = Volum_cilindru(D, l)
%Returneaza volumul unui cilindru, mm^3
%D  - Diametrul, mm
%l  - Lungimea, mm
digits(100);
p = vpa(pi);
result = D*D*p*l/4;
end

