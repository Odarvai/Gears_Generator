function [result] = RoataDintataCuPana(b, df, darb, mn, beta, leds, RDCP)
% Returneaza lungimea butucului, voulmul rotii (fara dantura) si dimensiunile penei alese
% b	- Latimea danturii, mm
% df	- Diametrul de picior al danturii, mm
% da	- Diametrul arborelui pe care se monteaza roata, mm
% mn	- Modulul normal al danturii, mm (pentru r,d,c.c.d.d se va introduce m)
% beta	- Unghiul de inclinare al danturii pe clindrul de divizare, rad (pentru r,d,c.c.d.d se va introduce 0)

    digits(100);
    PI = vpa(pi);

    lh=ceil(max(1.1*darb, b)); % Lungimea butucului, mm
    
    AP = double.empty();
    
    alegepana(darb, lh-leds, AP);
        bp=AP(1);
        hp=AP(2);
        lp=AP(3);
        t1p=AP(4);
        t2p=AP(5);

    dh=ceil(Maximum(2*(5+sqrt((darb/2+t2p)*(darb/2+t2p)+bp*bp/4)), 1.5*darb)); %Diametrul butucului rotii dintate, mm

    S=ceil(2.2*mn/cos(beta)+0.05*b); % Grosimea coroanei, mm

    if(df-2*S-dh>=20) C=ceil(Maximum(0.25*(2*S+dh-darb), 0.25*b)); 
    else S=0; C=b;
    end

    VolRoata=b*(df*df-(df-2*S)*(df-2*S)) ... % Volumul coroanei, mm^3 (fara PI/4)
            + lh*(dh*dh-darb*darb) ... % Volumul butucului, mm^3 (fara PI/4)
            + C*((df-2*S)*(df-2*S)-dh*dh); % Volumul discului, mm^3 (fara PI/4)
    VolRoata = ValRoata * PI/4;

    RDCP(1)=1; % Se va folosi pentru control
    RDCP(2)=lh;
    RDCP(3)=VolRoata;
    RDCP(4)=bp;
    RDCP(5)=hp;
    RDCP(6)=lp;
    RDCP(7)=t1p;
    RDCP(8)=t2p;
    RDCP(9)=dh;
    result = RDCP;
end

