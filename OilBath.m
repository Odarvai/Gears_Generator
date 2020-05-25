function [result] = OilBath(n2, n3, dw2, dw4, df2, df4, da2, da4, DeltaH)
        %Returneaza diferenta minima dintre nivelul maxim si minim al baii de ulei al unui reductor cu r.d.c.c.d.i. cu 2 trepte
		digits(100);
		PI = vpa(pi);

		vw2=PI*dw2*n2/60000; % Viteza rotii conduse 2 pe cilindrul de rostogolire
		vw4=PI*dw4*n3/60000; % Viteza rotii conduse 4 pe cilindrul de rostogolire

		k2=6;
		if(vw2<=2) k2=3;
        end
		k4=6;
		if(vw4<=2) k4=3;
        end
		H2max=0.95*df2/2; % Corespunde nivelului cel mai coborat al uleiului din baie
		H2min=(k2-2)*da2/k2/2; 
		
		H4max=0.95*df4/2; % Corespunde nivelului cel mai coborat al uleiului din baie
		H4min=(k4-2)*da4/k4/2;

		Hmax=min(H2max,H4max);
		Hmin=max(H2min,H4min);

		DeltaH=0.1;
		if ( Hmax>=Hmin) DeltaH=Hmax-Hmin;
        end
		result = DeltaH;
end

