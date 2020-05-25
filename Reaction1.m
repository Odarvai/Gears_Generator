function [result] = Reaction1 (n, r1, r2, u, F, Mic)
%  Returneaza reactiunea din nodul r1	
%  n   - Numarul de noduri ale arborelui
%  r1  - Numarul de ordine al nodului in care se considera reazemul 1
%  r2  - Numarul de ordine al nodului in care se considera reazemul 2, r2>r1
%  u   - Sirul ordonatelor nodurilor
%  F   - Sirul sarcinilor din noduri (pozitive daca sunt de jos in sus)
%  Mic - Sirul momentelor concentrate din noduri (pozitive daca sunt in sens orar)
 
    M=0;
    
    for i = 1:n  
        M = M + F(i) * (u(r2) - u(i)) + Mic(i);
    end
    result = M / (u(r1) - u(r2));
end

