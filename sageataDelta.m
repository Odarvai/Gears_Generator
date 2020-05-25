function [ DELTA] = sageataDelta( x, n, nr1, nr2, u, Iz_s, Iz_d, Mi_s, Mi_d, E)
%sageataDelta - returneaza valoarea sagetii in sectiunea x
%   Detailed explanation goes here

j = n - 1;

for i = 0 : n - 1
    if ( (u(i) <= x) && ( x < u(i+1) ) )
        j = i;
        break;
    else continue;
    end
end

if (j ~= n - 1)
    for i = 0 : n + 1
        if (i <= j)
            v(i) = u(i);
            vIz_s(i) = Iz_s(i);
            vMi_s(i) = Mi_s(i);
            vMi_d(i) = Mi_d(i);
        end
        if (i == j + 1)
            v(i) = x;
            vIz_s(i) = Iz_d(i - 1);
            vMi_s(i) = ( u(j + 1) - u(j) ~= 0 ? ((Mi_s[j+1]-Mi_s[j])*x+Mi_s[j]*u[j+1]-Mi_s[j+1]*u[j])/(u[j+1]-u[j]) : Mi_s[j]);
            vMi_d(i) = ( u(j + 1) - u(j) ~= 0 ? ((Mi_d[j+1]-Mi_d[j])*x+Mi_d[j]*u[j+1]-Mi_d[j+1]*u[j])/(u[j+1]-u[j]) : Mi_d[j]);
        end
        if (i > j + 1)
            v(i) = u(i - 1);
            vIz_s(i) = Iz_s(i - 1);
            vMi_s(i) = Mi_s(i - 1);
            vMi_d(i) = Mi_d(i - 1);
        end
    end
else
    for i = 0 : n
        v(i) = u(i);
        vIz_s(i) = Iz_s(i);
        vMi_s(i) = Mi_s(i);
        vMi_d(i) = Mi_d(i);
    end
    v(n) = u(n - 1);
    vIz_s(n) = Iz_s(n - 1);
	vMi_s(n) = Mi_s(n - 1);
	vMi_d(n) = Mi_d(n - 1);
end

if (j ~= n - 1 )
    if (j < nr1)
        vnr1 = nr1 + 1;
        vnr2 = nr2 + 1;
    end
    if ((j >= nr1) && (j < nr2))
        vnr1 = nr1;
        vnr2 = nr2 + 1;
    end
    if (j >= nr2)
        vnr1 = nr1;
        vnr2 = nr2;
    end
else
    if (nr2 ~= n - 1)
        vnr1 = nr1;
        vnr2 = nr2;
    else
        vnr1 = nr1;
        vnr2 = nr2;
    end
end

if(j ~= n - 1)
    ix = j + 1;
else
    ix = n;
end

for i = 0 : n + 1
    F(i) = 0;
    Mic(i) = 0;
end

F(ix) = 1;

R1 = reactiuneRA(n + 1, vnr1, vnr2, v, F, Mic);
R2 = reactiuneRB(n + 1, vnr1, vnr2, v, F, Mic);
F(vnr1) = R1;
F(vnr2) = R2;

for i = 0 : n + 1
    Mi_unit_s(i) = momentIS(v(i), n + 1, v, F, Mic);
    Mi_unit_d(i) = momentID(v(i), n + 1, v, F, Mic);
end

for i = 0 : n
    DELTA = DELTA + (v(i + 1) - v(i)) * (Mi_unit_d(i) * (2 * vMi_d(i) + vMi_s(i + 1)) + Mi_unit_s(i+1) * (vMi_d(i) + 2 * vMi_s(i + 1))) / vIz_s(i + 1); 
    DELTA = DELTA / E;
end
end

