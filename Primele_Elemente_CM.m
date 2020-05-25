function [result] = Primele_Elemente_CM(Aliniere_Rulmenti, aw_1, aw_2, CCr_1, CCr_3, Dr_1, Dr_3, Tr_1, Tr_3, Br_1, Br_3, c, lr, cb, bm_1, bm_3, ls_min, lplm, PECM)
    ls2=5; % distanta de siguranta dintre un corp in rotatie si un cap de surub, mm

	delta_corp=round(0.025*max(aw_1, aw_2)+5,0); % Grosimea peretelui corpului reductorului, mm
	delta_capac=round(0.825*delta_corp,0); % Grosimea peretelui capacului reductorului, mm
	
	d_sf=round(1.65*delta_corp,0); % Diametrul nominal al suruburilor de fixare a reductorului pe fundatie, mm
	if (d_sf<9)					d_sf= 8;a_t=18; b_t=15; end % a_t, b_t - distantele de la centrul gaurii de fixare a carcasei reductorului de fundatie la peretele carcasei respectiv la marginea aceteia
	if ((d_sf>=9)&&(d_sf<11))	d_sf=10;a_t=20; b_t=20; end
	if ((d_sf>=11)&&(d_sf<15))	d_sf=12;a_t=25; b_t=22; end
	if ((d_sf>=15)&&(d_sf<17))	d_sf=16;a_t=30; b_t=26; end
	if ((d_sf>=17)&&(d_sf<21))	d_sf=20;a_t=40; b_t=32; end
	if ((d_sf>=21)&&(d_sf<23))	d_sf=22;a_t=43; b_t=33; end
	if ((d_sf>=23)&&(d_sf<26))	d_sf=24;a_t=45; b_t=33; end
	if ((d_sf>=26)&&(d_sf<29))	d_sf=27;a_t=50; b_t=35; end
	if (d_sf>=29)				d_sf=30;a_t=55; b_t=40; end

	t=round(2.25*delta_corp,0); % Inaltimea talpii corpului reductorului - varianta fara bosaje,mm

	if ((1.25<=t/delta_corp)&&(t/delta_corp<1.8)) x_t=t-delta_corp; end % Zona de trecere, mm - vezi Atlas reductoare, pag 47, fig 4.13, tabel 4.12
	if ((1.8 <=t/delta_corp)&&(t/delta_corp<2.5)) x_t=0.8*(t-delta_corp);   end
	if  (2.5 <=t/delta_corp) x_t=0.7*(t-delta_corp);    end

	l_t=a_t+b_t+x_t; % Latimea talpii corpului reductorului, mm

	d_s1=round(0.75*d_sf,0); % Diametrul nominal al suruburilor de la lagare (M8 - M24), mm
	if (d_s1<9)					d_s1= 8; a_f=18; b_f=15; Scheie1=13; ep=11; mp=13; Pp=28;   end % a_f, b_f - distantele de la centrul gaurii de fixare a carcasei reductorului de fundatie la peretele carcasei respectiv la marginea aceteia
	if ((d_s1>= 9)&&(d_s1<11))	d_s1=10; a_f=20; b_f=20; Scheie1=17; ep=14; mp=18; Pp=36;   end % ep,mp- distante de la axa surubului la perete, Pp-distanta minima intre 2 suruburi , mm - vezi SATS 776/89 
	if ((d_s1>=11)&&(d_s1<15))	d_s1=12; a_f=25; b_f=22; Scheie1=19; ep=15; mp=20; Pp=40;   end
	if ((d_s1>=15)&&(d_s1<17))	d_s1=16; a_f=30; b_f=26; Scheie1=24; ep=17; mp=23; Pp=49;   end
	if ((d_s1>=17)&&(d_sf<21))	d_s1=20; a_f=40; b_f=32; Scheie1=30; ep=21; mp=29; Pp=60;   end
	if ((d_s1>=21)&&(d_sf<23))	d_s1=22; a_f=43; b_f=33; Scheie1=32; ep=23; mp=30; Pp=64;   end
	if (d_s1>=23)				d_s1=24; a_f=45; b_f=33; Scheie1=36; ep=25; mp=33; Pp=72;   end

	delta_as=d_s1; % Distanta de la peretele alezajului la axa surubului de la lagare, mm
	
	d_s2=round(0.5*d_sf,0); % Diametrul nominal al suruburilor de imbinare a capacului rulmentului cu corpul, mm
	if (d_s2<=7)				d_s2= 6; ks2=4;   gws2=1.6; end
	if ((d_s2>7) && (d_s2<=9))	d_s2= 8; ks2=5.5; gws2=2;   end
	if ((d_s2>9) && (d_s2<=11))	d_s2=10; ks2=7;   gws2=2.5; end
	if (d_s2>11)				d_s2=12; ks2=8;   gws2=3;   end

	hs_1=round(bm_1+1.2,0);
	hs_3=round(bm_3+1.2,0);

	ecr = floor (1.2*d_s2); % Grosimea capacului rulmentului in zona suruburilor de fixare, mm (preliminar)
	
	K_f=mp+ep+delta_corp; % Latimea flansei, mm (preliminar)
	K_fb=K_f+cb; % Latimea flansei in zona bosajelor, mm (preliminar) 
    
    if (0.5*sqrt(Dr_1)>5)
        lc_1_nec = 0.5*sqrt(Dr_1);      % Lungimea de centrare neceasara a capacului 1, mm
    else
        lc_1_nec = 5;
    end
    
    if (0.5*sqrt(Dr_3)>5)
        lc_3_nec = 0.5*sqrt(Dr_3);      % Lungimea de centrare neceasara a capacului 3, mm
    else
        lc_3_nec = 5;
        
	if (Aliniere_Rulmenti==0)
        
		if (K_fb-CCr_1-c < lc_1_nec) K_fb=ceil(CCr_1+c+lc_1_nec); end
		if (K_fb-CCr_3-c < lc_3_nec) K_fb=ceil(CCr_3+c+lc_3_nec); end % Latimea flansei in zona bosajelor, mm (definitiv)

		lc_1=ceil(K_fb-CCr_1-c); % lungimea de centrare, mm (final)
		lc_3=ceil(K_fb-CCr_3-c); % lungimea de centrare, mm (final)	
    end

	if (Aliniere_Rulmenti==1)
	
		if (K_fb-c-Tr_1 < lc_1_nec) K_fb=ceil(Tr_1+c+lc_1_nec); end
		if (K_fb-c-Tr_3 < lc_3_nec) K_fb=ceil(Tr_3+c+lc_3_nec); end %Latimea flansei in zona bosajelor, mm (definitiv)

		lc_1=ceil(K_fb-Tr_1-c); % lungimea de centrare, mm (final)
		lc_3=ceil(K_fb-Tr_3-c); % lungimea de centrare, mm (final)
    end
	
	K_f=K_fb-cb; % Latimea flansei, mm (definitiv)

	if (ecr+lr+lc_1+c < ls_min+hs_1+lplm) ecr=ls_min+hs_1+lplm-lr-lc_1-c; end
	if (ecr+lr+lc_3+c < ls_min+hs_3+lplm) ecr=ls_min+hs_1+lplm-lr-lc_1-c; end %Grosimea capacului rulmentului in zona suruburilor de fixare, mm (definitiv)

	if (ecr+lr+lc_1+c > ls_min+hs_1+lplm)
        ls_1=ecr+lr+lc_1+c-hs_1-lplm; 
    end
	if (ecr+lr+lc_3+c > ls_min+hs_3+lplm)
        ls_3=ecr+lr+lc_3+c-hs_3-lplm;  %distanta dintre manseta si rulment, mm (definitiv)
    end
     
	le_1= ceil( ls2+ks2+gws2+lplm+hs_1+ls_1+Tr_1-Br_1 );
	le_3=ceil( ls2+ks2+gws2+lplm+hs_3+ls_3+Tr_3-Br_3 );


	PECM(1)=1; % pentru control
	PECM(2)=t; % Inaltimea talpii corpului reductorului - varianta fara bosaje,mm
	PECM(3)=l_t; % Latimea talpii corpului reductorului, mm
	PECM(4)=K_f; % Latimea flansei, mm
	PECM(5)=d_sf; % Diametrul nominal al suruburilor de fixare a reductorului pe fundatie, mm
	PECM(6)=d_s1; % Diametrul nominal al suruburilor de la lagare (M8 - M24), mm
	PECM(7)=Pp; % Distanta minima intre 2 suruburi, mm - vezi SATS 776/89
	PECM(8)=d_s2; % Diametrul nominal al suruburilor de imbinare a capacului rulmentului cu corpul, mm

	PECM(9)=delta_corp; % Grosimea peretelui corpului, mm
	PECM(10)=delta_capac; % Grosimea peretelui capacului, mm
	PECM(11)=delta_as; % Distanta de la peretele alezajului la axa surubului de la lagare, mm


	PECM(12)=ecr; % Grosimea capacului rulmentului in zona suruburilor de fixare, mm
	PECM(13)=hs_1; % latimea locasului mansetei, mm
	PECM(14)=hs_3; % latimea locasului mansetei, mm
	PECM(15)=ls_1; % distanta dintre manseta si rulment, mm
	PECM(16)=ls_3; % distanta dintre manseta si rulment, mm
	PECM(17)=lc_1; % lungimea de centrare a capacului, mm
	PECM(18)=lc_3; % lungimea de centrare a capacului, mm
	PECM(19)=le_1; % lungimea tronsonului de etansare a arborelui 1, mm
	PECM(20)=le_3; % lungimea tronsonului de etansare a arborelui 3, mm
    result = PECM;
end

