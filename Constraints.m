%Constraint function
function [c, c_eq] = Constraints(x)

    x = load('temporary.mat');

    PosibilAngrenaj = x.PosibilAngrenaj;
    if (~PosibilAngrenaj)
        return;
    u12 = x.u12;
    i12STAS = x.i12STAS;
    sigma_H = x.sigma_H;
    sigma_HP = x.sigma_HP;
    sigma_F1 = x.sigma_F1;
    sigma_FP1 = x.sigma_FP1;
    sigma_F2 = x.sigma_F2;
    sigma_FP2 = x.sigma_FP2;
    xn1 = x.xn1;
    zn1 = x.zn1;
    xn2 = x.xn2;
    zn2 = x.zn2;
    cs_n = x.cs_n;
    mn = x.mn;
    san1 = x.san1;
    san2 = x.san2;
    eps_alfa_lim = x.eps_alfa_lim;
    eps_alfa = x.eps_alfa;
    WNn1 = x.WNn1;
    beta_b = x.beta_b;
    b1 = x.b1;
    %WNn2 = x.WNn2;
    b2 = x.b2;
    roAt1 = x.roAt1;
    roNt1 = x.roNt1;
    roat1 = x.roat1;
    roEt2 = x.roEt2;
    roNt2 = x.roNt2;
    roat2 = x.roat2;

    digits(100);
    PI = vpa(pi);

    tolStrOfMat 	=		  0.002;     % toleranta cu care se doreste verificarea: sigma <= (1 + tolStrOfMat)*sigma_all


    if (PosibilAngrenaj)
        %R1 Eroarea relativa a raportului de transmitere trebuie sa fie in intervalul [-2,5%...+2,5%], daca u12<4 sau in [-3%...+3%] in caz contrar.
        if (u12 < 4) 
            g(1) = abs(i12STAS - u12) * 40 / i12STAS - 1;
        else
            g(1) = abs(i12STAS - u12) * 100 / i12STAS / 3 - 1;
        end
        % R2 Verificarea la presiunea de contact.
        g(2) = sigma_H / sigma_HP / (1 + tolStrOfMat) - 1;
        % R3 Verificarea la incovoiere a dintelui pinionului.
        g(3) = sigma_F1 / sigma_FP1 / (1 + tolStrOfMat) - 1;
        % R4 Verificarea la incovoiere a dintelui rotii.
        g(4) = sigma_F2 / sigma_FP2 / (1 + tolStrOfMat) - 1;
        % R5 Verificarea danturii pinionului la subtaiere.
        if (xn1 > 0) 
            g(5) = (14 - zn1) / 17 / xn1 - 1;
        elseif (xn1 < 0) 
                g(5) = 1 - (14 - zn1) / 17 / xn1;
        else
            g(5) = 14 / zn1 - 1;            
        end
        %R6 Verificarea danturii rotii la subtaiere.
        if (xn2 > 0) 
            g(6) = (14 - zn2) / 17 / xn2 - 1;
        elseif (xn2 < 0)
                g(6) = 1 - (14 - zn2) / 17 / xn2;  
        else
            g(6) = 14 / zn2 - 1;
        end
        % R7 Verificarea danturii pinionului la ascutire.
        g(7) = cs_n * mn / san1 - 1;
        % R8 Verificarea danturii rotii la ascutire.
        g(8) = cs_n * mn / san2 - 1;
        % R9 Gradul de acoperire frontal trebuie sa fie mai mare decat o valoare minima impusa (in general in functie de viteza anfrenajului).
        g(9) = eps_alfa_lim / eps_alfa - 1;
        % R10 Se verifica daca xn2 este in intervalul [-0.5, 1].
        g(10) = abs(xn2 - 0.25) / 0.75 - 1;
        % R11-16 Pentru masurarea cotei peste dinti trebuie indeplinite conditiile.
        g(11) = (WNn1 * sin(beta_b) + 5) / b1 - 1;
        g(12) = (WNn2 * sin(beta_b) + 5) / b2 - 1;
        g(13) = roAt1 / roNt1 - 1;
        g(14) = roNt1 / roat1 - 1;
        g(15) = roEt2 / roNt2 - 1;
        g(16) = roNt2 / roat2 - 1;
        % R17 z1 si z2 trebuie sa fie prime intre ele.
        g(17) = Verif_cmmdc(z1,z2);
        c = g;
        %c_eg = 0;
    else
        return;
    end  
end