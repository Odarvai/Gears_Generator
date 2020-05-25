function  RB = reactiuneRB(n, r1, r2, u, F, Mic)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


RB = 0;

for i = 1 : n
    RB = RB + F(i) * ( u(r1)-u(i) ) + Mic(i);
    
      
    
    
end

RB = RB / ( u(r2)-u(r1) );
end

