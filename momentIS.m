function MOMS = momentIS( x, n, u, F, Mic )
%momentIS - returneaza momentul incovoietor la distanta x de nodul 0, consoderat in stanga lui

% n   - numarul de noduri ale arborelui
% u   - sirul ordonatelor nodurilor
% F   - sirul REACTUALIZAT al sarcinilor din noduri (pozitive daca sunt de jos in sus) !!! S-au introdus reactiunile in nodul (indexul) corespunzator
% Mic - sirul momentelor concentrate din noduri (pozitive daca sunt in sens orar) 

Ms = 0;

for i = 1 : n
    if(u(i) < x)
        Ms = Ms + F(i) * (x - u(i)) + Mic(i);
    else break;
    end
end

MOMS = Ms;

end

