function  RA = reactiuneRA(n, r1, r2, u, F, Mic)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


RA = 0;

for i = 1 : n 
    RA = RA + F(i) * (u(r2) - u(i)) + Mic(i);
end

RA = RA / (u(r1) - u(r2));

end

