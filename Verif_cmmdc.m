function conditie = Verif_cmmdc(a, b)
%Daca numerele a si b sunt prime intre ele functia returneaza -1
% In caz contrar returneaza +1

if(gcd(a, b) == 1)  conditie = -1;
else
    conditie = 1;
end
end
    
