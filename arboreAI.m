function [ HA, HB, VA, VB, MiH_s, MiH_d ] = arboreAI( n, r1, r2, u, F, Mic )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%Calculul reactiunilor
%--------------------------------------------------------------------------
HA = reactiuneRA(n, r1, r2, u, F, Mic);
HB = reactiuneRB(n, r1, r2, u, F, Mic);

VA = reactiuneRA(n, r1, r2, u, F, Mic);
VB = reactiuneRB(n, r1, r2, u, F, Mic);
%--------------------------------------------------------------------------

%FV=VA;	FV=VB;

%Momente incovoietoare
%--------------------------------------------------------------------------
MiH_s = zeros(1,n);
MiH_d = zeros(1,n);
for i = 1 : n
    MiH_s(i) = momentIS( u(i), n, u, F, Mic );
    MiH_d(i) = momentID( u(i), n, u, F, Mic );
end


    
    
   
    
end

