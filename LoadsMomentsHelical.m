function [result] = LoadsMomentsHelical(T1, dw1, dw2, beta_w, alfa_wt, LMH)
%  Returneaza valoarea fortelor tangentiale,axiale si radiale precum si valoarea momentului concentrat
% 		 dw1		- Diametrul de rostogolire al pinionului, mm
% 		 dw2		- Diametrul de rostogolire al rotii, mm
% 		 beta_w	- Unghiul de inclinare al danturii pe cilindrul de rostogolire, rad
% 		 alfa_wt	- Unghiul real de angrenare in plan frontal, rad
        Ft1=2*T1/dw1;
		Fa1=Ft1*tan(beta_w);
		Fr1=Ft1*tan(alfa_wt);

		Mic1=Fa1*dw1/2;
		
		Ft2=Ft1;
		Fa2=Fa1;
		Fr2=Fr1;

		Mic2=Fa2*dw2/2;

		LMH(1)=1; % Se va folosi pentru control
		LMH(2)=Ft1; LMH(3)=Fa1; LMH(4)=Fr1; LMH(5)=Mic1;
		LMH(6)=Ft2; LMH(7)=Fa2; LMH(8)=Fr2; LMH(9)=Mic2;
        result = LMH;
end

