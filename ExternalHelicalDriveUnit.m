function [result] = ExternalHelicalDriveUnit(mn, aw, z1, z2, xn1, beta, b2, Delta_b, ...
							   T1, u12, KA, Lh1, Lh2, ZL1, ZL2, ZX, ...
							   Zw, ZNmax, YNmax, NstH, NBH, NstF, NBF, SHmin, ...
							   SFmin, YR1, YR2, Hi1, Hi2, niu1, niu2, E1, E2, ...
							   MatTT1, MatTT2, HB1, HB2, sigma021, sigma022, sigma_Hlim1, ...
							   sigma_Hlim2, sigma_Flim1, sigma_Flim2, Raf1, Raf2, Rar1, Rar2, ...
							   PozAngr, CPrec, n1, n2, EHDU)
                               
%         // mn			 - Modulul normal al angrenajului, mm
% 		// aw			 - Distanta axiala standardizata/impusa, mm
% 		// z1,z2		 - Numerele de dinti ale pinionului si ale rotii dintate
% 		// xn1,xn2		 - Coeficiemtii deplasarilor de profil in plan normal corespunzatori pinionului
% 		// beta			 - Unghiul de inclinare al danturii pe cilindrul de divizare, grade
% 		// b2			 - Latimea rotii dintate, mm
% 		// Delta_b		 - Diferenta de latime dintre pinion si roata dintata, mm 
% 		// T1			 - Momentul de torsiune corespunzator arborelui de intrare, Nmm
% 		// u12			 - Raportul deal de angrenare corespunzator treptei I
% 		// KA			 - Factorul regimului de functionare
% 		// Lh1,Lh2		 - Durata minima de functionare a pinionului/rotii dintate, ore
% 		// ZL1,ZL2		 - Factorul de ungere corespunzator pinionului/rotii dintate
% 		// Zx			 - Factorul de marime
% 		// ZW			 - Factorul raportului duritatilor flancurilor dintilor
% 		// ZNmax		 - Factorul durabilitatii pentru solicitarea de contact
% 		// YNmax		 - Factorul durabilitatii pentru solicitarea de incovoiere
% 		// NstH,NstF	 - 
% 		// NBH,NBF		 -
% 		// SHmin		 - Coeficientul de siguranta minim pentru solicitarea de contact
% 		// SFmin		 - Coeficientul de siguranta minim pentru solicitarea de incovoiere
% 		// YR1,YR2		 - Factorul rugozitatii flancului dintelui pinionului/rotii dintate pentru solicitarea de incovoiere
% 		// Hi1,Hi2		 - Numarul de roti cu care vine in contact pinionul/roata dintata
% 		// niu1,niu2	 - Coeficientii lui Poisson
% 		// E1,E2		 - Modulul de elasticitate longitudinal, MPa
% 		// MatTT1,MatTT2 -
% 		// HB1,HB2		 - Duritatea materialului pinionului/rotii dintate, MPa
% 		// sigma021,2	 - Limita de curegere a materialului pinionului/rotii dintate, MPa
% 		// sigma_Hlim1,2 - Tensiunea hertziana limita pentru materialul pinionului/rotii dintate, MPa
% 		// sigma_Flim1,2 - Tensiunea de incovoiere limita pentru materialul pinionului/rotii dintate, MPa
% 		// Raf1,Raf2	 - Rugozitatea flancului dintelui pinionului/rotii dintate, micrometri
% 		// Rar1,Rar2	 - Rugozitatea razei de racordare a dintelui pinionului/rotii dintate, micrometri
% 		// PozAngr		 - Pozitia angrenajului
% 		// CPrec		 - Clasa de precizie a angrenajului
% 		// n1,n2		 - Turatia arborelui 1 si 2, rot/min 

        HDG = zeros(1, 42, 'double');
        EHD = zeros(1, 14, 'double');
        CHD = zeros(1, 17, 'double');
        HDS = zeros(1, 7, 'double');
        LMH = zeros(1, 9, 'double');
        
        validEHDU=1;
		while (validEHDU==1)

            % Elementele geometrice ale angrenajului
            HelicalDriveGeometry(mn,aw,z1,z2,xn1,beta,b2,Delta_b,HDG);
            if(HDG(1)==0) validEHDU=0; break;
            end
            alfa_t=HDG(2); 
            alfa_at1=HDG(3);  
            beta_a1=HDG(4);   
            beta_b=HDG(5); 
            d1=HDG(6);      
            db1=HDG(7);       
            df1=HDG(8);       
            da1=HDG(9);
            st1=HDG(10);     
            sat1=HDG(11);     
            sn1=HDG(12);      
            san1=HDG(13);  
            xt1=HDG(14);    
            ADinte1=HDG(15);  
            a=HDG(16);        
            alfa_wt=HDG(17);
            xsn=HDG(18);    
            xn2=HDG(19);      
            alfa_at2=HDG(20); 
            beta_a2=HDG(21);
            d2=HDG(22);     
            db2=HDG(23);      
            df2=HDG(24);      
            da2=HDG(25);
            st2=HDG(26);    
            sat2=HDG(27);	  
            sn2=HDG(28);	   
            san2=HDG(29);
            xt2=HDG(30);    
            ADinte2=HDG(31); 
            dw1=HDG(32);	   
            dw2=HDG(33);
            beta_w=HDG(34); 
            eps_alfa=HDG(35);
            eps_beta=HDG(36);
            eps_gama=HDG(37);
            b1=HDG(38);
            cs_n=HDG(39);
            alfa_n=HDG(40);
            VolDantura1=HDG(41);
            VolDantura2=HDG(42);

            % Angrenajul echivalent
            EquivalentHelicalDrive(mn,a,z1,z2,beta,beta_b,beta_w,alfa_n,alfa_wt,d1,da1,d2,da2,EHD);
            % Elementele angrenajului echivalent
            zn1=EHD(1);
            zn2=EHD(2);
            xn1min=EHD(3);
            xn2min=EHD(4);
            dn1=EHD(5);
            dn2=EHD(6);
            dbn1=EHD(7);
            dbn2=EHD(8);
            dan1=EHD(9);
            dan2=EHD(10);
            alfa_wn=EHD(11);
            awn=EHD(12);
            eps_alfa_n=EHD(13);

            % Elementele de control ale angrenajului
            ControlHelicalDrive(mn,aw,z1,z2,xn1,xn2,beta,beta_b,beta_w,alfa_n,alfa_t,alfa_wt,alfa_at1,alfa_at2,d1,da1,d2,db1,db2,da2,CHD);	
            alfa_Nt1=CHD(1); 
            alfa_Nt2=CHD(2); 
            Ncalc1=CHD(3); 
            Ncalc2=CHD(4); 
            N1=CHD(5);	     
            N2=CHD(6);       
            WNn1=CHD(7);   
            WNn2=CHD(8); 
            WNt1=CHD(9);	 
            WNt2=CHD(10);	  
            roNt1=CHD(11); 
            roNt2=CHD(12); 
            roAt1=CHD(13);   
            roEt2=CHD(14); 
            roat1=CHD(15);
            roat2=CHD(16);

            % Calculul fortelor si momentelor de incovoiere concentrate (datorate fortelor axiale)
            LoadsMomentsHelical (T1,dw1,dw2,beta_w,alfa_wt,LMH);
            Ft1=LMH(1); 
            Fa1=LMH(2);
            Fr1=LMH(3); 
            Mic1=LMH(4);

            Ft2=LMH(5); 
            Fa2=LMH(6); 
            Fr2=LMH(7); 
            Mic2=LMH(8);

            % Calculul de rezistenta
            HelicalDriveStrength(T1,z1,u12,KA,Lh1,Lh2,ZL1,ZL2,ZX,Zw,ZNmax,YNmax,NstH,NBH,NstF,NBF,SHmin,SFmin,YR1,YR2,Hi1,Hi2,niu1,niu2,E1,E2,MatTT1,MatTT2,HB1,HB2,sigma021,sigma022,sigma_Hlim1,sigma_Hlim2,sigma_Flim1,sigma_Flim2,Raf1,Raf2,Rar1,Rar2,PozAngr,CPrec,mn,aw,alfa_t,alfa_wt,beta,beta_b,eps_alfa,eps_beta,eps_alfa_n,b1,b2,zn1,zn2,xn1,xn2,d1,d2,n1,n2,HDS);	
            if(HDS(1)==0) validEHDU=0; break;
            end
            sigma_H=HDS(2);
            sigma_HP=HDS(3);
            sigma_F1=HDS(4);
            sigma_F2=HDS(5);
            sigma_FP1=HDS(6);
            sigma_FP2=HDS(7);
            break;
        end
		if (validEHDU)
			EHDU(1)=1;

            % Elementele geometrice ale angrenajului
            EHDU(2)=alfa_t; 
            EHDU(3)=alfa_at1;  
            EHDU(4)=beta_a1;   
            EHDU(5)=beta_b; 
            EHDU(6)=d1;      
            EHDU(7)=db1;      
            EHDU(8)=df1;       
            EHDU(9)=da1;
            EHDU(10)=st1;     
            EHDU(11)=sat1;     
            EHDU(12)=sn1;      
            EHDU(13)=san1;  
            EHDU(14)=xt1;    
            EHDU(15)=ADinte1;  
            EHDU(16)=a;        
            EHDU(17)=alfa_wt;
            EHDU(18)=xsn;    
            EHDU(19)=xn2;      
            EHDU(20)=alfa_at2; 
            EHDU(21)=beta_a2;
            EHDU(22)=d2;     
            EHDU(23)=db2;      
            EHDU(24)=df2;      
            EHDU(25)=da2;
            EHDU(26)=st2;    
            EHDU(27)=sat2;	  
            EHDU(28)=sn2;	   
            EHDU(29)=san2;
            EHDU(30)=xt2;    
            EHDU(31)=ADinte2; 
            EHDU(32)=dw1;	   
            EHDU(33)=dw2;
            EHDU(34)=beta_w; 
            EHDU(35)=eps_alfa;
            EHDU(36)=eps_beta;
            EHDU(37)=eps_gama;
            EHDU(38)=b1;
            EHDU(39)=cs_n;
            EHDU(40)=alfa_n;

            % Elementele geometrice ale angrenajului echivalent
            EHDU(41)=zn1;
            EHDU(42)=zn2;
            EHDU(43)=xn1min;
            EHDU(44)=xn2min;
            EHDU(45)=dn1;
            EHDU(46)=dn2;
            EHDU(47)=dbn1;
            EHDU(48)=dbn2;
            EHDU(49)=dan1;
            EHDU(50)=dan2;
            EHDU(51)=alfa_wn;
            EHDU(52)=awn;
            EHDU(53)=eps_alfa_n;

            % Elementele de control ale angrenajului
            EHDU(54)=alfa_Nt1; 
            EHDU(55)=alfa_Nt2; 
            EHDU(56)=Ncalc1; 
            EHDU(57)=Ncalc2; 
            EHDU(58)=N1;	     
            EHDU(59)=N2;       
            EHDU(60)=WNn1;   
            EHDU(61)=WNn2; 
            EHDU(62)=WNt1;	 
            EHDU(63)=WNt2;	  
            EHDU(64)=roNt1; 
            EHDU(65)=roNt2;
            EHDU(66)=roAt1;   
            EHDU(67)=roEt2; 
            EHDU(68)=roat1;
            EHDU(69)=roat2;

            % Tensiuni si tensiuni admisibile
            EHDU(70)=sigma_H;
            EHDU(71)=sigma_HP;
            EHDU(72)=sigma_F1;
            EHDU(73)=sigma_F2;
            EHDU(74)=sigma_FP1;
            EHDU(75)=sigma_FP2;

            EHDU(76)=VolDantura1;
            EHDU(77)=VolDantura2;

            EHDU(78)=Ft1; 
            EHDU(79)=Fa1; 
            EHDU(80)=Fr1; 
            EHDU(81)=Mic1;
            EHDU(82)=Ft2; 
            EHDU(83)=Fa2; 
            EHDU(84)=Fr2; 
            EHDU(85)=Mic2;
        
		else EHDU(1)=0;
        end
      result = EHDU;
end

