function [ COEFSO ] = coeficientsigursantalaobsealaC( tipconc, tiptrat, Rm,  )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

% BetaKSigma
if (tipconc == 0)
    BetaKSigma = 1;
elseif (tipconc == 1)
    BetaKSigma = 5.31067 + Rm * (-0.0365703 + Rm * (0.000129762 + Rm * (-2.34808e-007 + Rm * (2.35451e-010 + Rm *(-1.24372e-013 + 2.70484e-017 * Rm))))));
else
    BetaKSigma = fbetaKSigma(Rm, D, d, rr);
end

%BetaKTau
if (tipconc == 0) 
    BetaKTau = 1;
elseif (tipconc == 1)
    BetaKTau = (11.1196 + Rm * (-0.0882135 + Rm * (0.000310523 + Rm * (-5.55491e-007 + Rm * (5.42524e-010 + Rm * (-2.75733e-013 + 5.72111e-017 * Rm))))));
else
    BetaKTau = fbetaKTau(Rm, D, d, rr);
end
%-------------------------------------------------------------------------------

%Coeficient dimensional EpsilonSigma
if (Rm <= 800)
    EpsilonSigma = 1.00985 + d * (-0.0111193 + d * (0.000194058 + d * (-2.55797e-006 + d * (1.94378e-008 + d * (-7.50335e-011 + 1.14459e-013 * d)))));
else
    EpsilonSigma = 0.887348 + d * (-0.00961303 + d * (8.49069e-005 + d * (-3.96427e-007 + d * (1.61692e-009 + d * (-7.74794e-012 + 1.7666e-014 * d)))));
end
%-------------------------------------------------------------------------------

%Coeficientul dimensional EpsilonTau
EpsilonTau = (0.922898 + d * (-0.00545707 + d * (-7.62714e-005 + d * (2.28387e-006 + d * (-1.97118e-008 + d * (-7.29049e-011-9.91386e-014 * d)))))); 
%-------------------------------------------------------------------------------

if (Ra <= 0.8)
    Beta1 = 1;
elseif (Ra == 1.6)
    Beta1 = 0.981253 - 5.10827e-005 * Rm;
elseif (Ra == 3.2)
    Beta1 = 0.961508 - 7.07299e-005 * Rm;
else
    Beta1 = 1.01072 - 0.000227908 * Rm;
end
%-------------------------------------------------------------------------------

if (tiptrat == 0)
    Beta2 = 1;
elseif ((tiptrat == 1) && ((1 <= BetaKSigma) && (BetaKSigma < 1.5)))
    Beta2 = 1.4;
elseif ((tiptrat == 1) && ((1.5 <= BetaKSigma) && (BetaKSigma < 1.8)))
    Beta2 = 1.6;
elseif ((tiptrat == 1) && (BetaKSigma >= 1.8)) 
    Beta2 = 2.2;
end
%-------------------------------------------------------------------------------

if(Sigmai ~= 0)
    cSigma = Sigmaminus1 * EpsilonSigma * Beta1 * Beta2 / BetaKSigma / Sigmai;
else cSigma = 0;
end

if(Taut ~= 0)
    cTau = 2 * Tauminus1 / (BetaKTau / EpsilonTau / Beta1 / Beta2 + PsiTau) / Taut;
else cTau = 0;

%-------------------------------------------------------------------------------
if ((cSigma == 0) && (cTau == 0))  
    COEFSO = 1e10;
elseif ((cSigma == 0) && (cTau ~= 0))
    COEFSO = cTau;
elseif ((cSigma ~= 0)&&(cTau == 0))
    COEFSO = cSigma;
elseif ((cSigma ~= 0) && (cTau ~= 0))  
    COEFSO = cSigma * cTau / sqrt(cSigma * cSigma + cTau * cTau);
end
%-------------------------------------------------------------------------------









end

