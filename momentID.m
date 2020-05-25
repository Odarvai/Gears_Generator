function MOMD = momentID( x, n, u, F, Mic )
%moment_i_d - returneaza momentul incovoietor la distanta x de nodul 0, consoderat in dreapta lui

% n   - numarul de noduri ale arborelui
% u   - sirul ordonatelor nodurilor
% F   - sirul REACTUALIZAT al sarcinilor din noduri (pozitive daca sunt de jos in sus) !!! S-au introdus reactiunile in nodul (indexul) corespunzator
% Mic - sirul momentelor concentrate din noduri (pozitive daca sunt in sens orar) 

Md = 0;

for i = 1 : n
    if(u(i) <= x)
        Md = Md + F(i) * (x - u(i)) + Mic(i);
    else break;
    end
end

MOMD = Md;

end

