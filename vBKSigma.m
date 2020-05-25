function [ BKS ] = vBKSigma( Rm, D, d, rr )
%vBKSigma - returneaza valoarea lui Beta_k_Sigma

%j - numar de ordine
%D - diametrul mare
%d - diametrul mic
%rr - raza de racordare

if ( 500 >= Rm ) 
    BKS = vdiagBetaKSigma( 1, D, d, rr );
elseif ( (500 < Rm ) && ( Rm <= 600) )
    BKS = 
    
else BKS = 1;
end



if ((500<sigma_r)&&(sigma_r<=600)) return (interp(500,600,val_diag_beta_k_sigma(0,D,d,rr),val_diag_beta_k_sigma(1,D,d,rr),sigma_r));
		if ((600<sigma_r)&&(sigma_r<=700)) return (interp(600,700,val_diag_beta_k_sigma(1,D,d,rr),val_diag_beta_k_sigma(2,D,d,rr),sigma_r));
		if ((700<sigma_r)&&(sigma_r<=800)) return (interp(700,800,val_diag_beta_k_sigma(2,D,d,rr),val_diag_beta_k_sigma(3,D,d,rr),sigma_r));
		if ((800<sigma_r)&&(sigma_r<1200)) return (interp(800,1200,val_diag_beta_k_sigma(3,D,d,rr),val_diag_beta_k_sigma(4,D,d,rr),sigma_r));
		if ((sigma_r>=1200)) return (val_diag_beta_k_sigma(4,D,d,rr));




end

