function gearsMass = gearFitness(x)

%marim precizia variabilei PI la 100 de zecimale
digits(100);
PI = vpa(pi); 
EHDU1 = zeros(85);

Lista_mn = [1.0, 1.125, 1.25, 1.375, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 25.0, 28.0, 32.0, 36.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0];

Lista_aw = [56.0, 63.0, 71.0, 80.0, 90.0, 100.0, 112.0, 125.0, 140.0, 160.0, 180.0, 200.0, 225.0, 250.0, 280.0, 315.0];

Lista_iSTAS = [1.12, 1.25, 1.40, 1.60, 1.80, 2.00, 2.24, 2.50, 2.80, 3.15, 3.55, 4.00, 4.50, 5.00, 5.60, 6.30, 7.10, 8.00, 9.00, 10.00, 11.2, 12.50, 14.00, 16.00, 18.00, 20.00, 22.40, 25.00, 28.00, 31.50, 35.50, 40.00];

Lista_z1 = [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56];

P = 2.9;
i12STAS = 3.55;
n1 = 925;

%z1 = Lista_z1(x(1));
%iaw = x(2);
%xn1 = x(3);
%psi_a = x(4);
%beta_g = x(5);
z1 = 31;
iaw = 6;
xn1 = 0.9764;
psi_a = 0.2142;
beta_g = 15;


ZL1             =		  1.05;		% Factor de ungere
ZL2             =	      1.05;		% Factor de ungere
Lh1             =         9000;		% Durata minima de functionare a pinionului, ore
Lh2             =         9000;		% Durata minima de functionare a rotii, ore
Hi1             =		  1;		% Numarul de roti cu care vine in contact pinionul
Hi2             =		  1;		% Numarul de roti cu care vine in contact roata

Delta_b         =		  4;		% Diferenta de latime a rotilor, mm

PozAngr         =		  5;		% Pozitia rotilor angrenajului in raport cu lagarele arborilor transmisiei
CPrec           =		  8; 		% Clasa de precizie

eps_alfa_lim    =		  1.2;		% Valoarea minima admisibila a gradului de acoperie frontal

% Caracteristicile materialelor pentru trepta I
%ro_otel         =         7.85e-6;	% Densitatea otelului, kg/mm^3
MatTT1          =         1;		% Pinion: 41MoCr11 imbunatatit
HB1             =         3000;		% Duritatea flancului dintelui pinionului, MPa
sigma021		=	      690;		% Limita de curgere a materialului pinionului, MPa
Raf1			=		  0.8;		% Rugozitatea flancului dintelui pinionului, micrometri
Rar1			=		  1.6;		% Rugozitatea zonei de racordare a dintelui pinionului, micrometri
sigma_Hlim1     =         760;		% Tensiunea hertziana limita a materialului pinionului, MPa
sigma_Flim1     =		  580;		% Tensiunea de incovoiere limita a materialului pinionului, MPa
YR1             =		  1.02;		% Factorul rugozitatii flancului dintelui pinionului pentru solicitarea de incovoiere
niu1			=		  0.3;		% Coeficientul lui Poisson pentru materialul pinionului
E1              =		  2.05e5;   % Modulul de elasticitate longitudinal pentru materialul pinionului

MatTT2          =      	  1;    	% Roata: 40Cr10 imbunatatit
HB2             =         2700;     % Duritatea flancului dintelui rotii, MPa
sigma022		=	      560;		% Limita de curgere a materialului rotii, MPa
Raf2			=		  0.8;		% Rugozitatea flancului dintelui rotii, micrometri
Rar2			=		  1.6;		% Rugozitatea zonei de racordare a dintelui rotii, micrometri
sigma_Hlim2     =   	  720;		% Tensiunea hertziana limita a materialului rotii, MPa
sigma_Flim2     =		  560;		% Tensiunea de incovoiere limita a materialului rotii, MPa
YR2             =		  1.02;		% Factorul rugozitatii flancului dintelui rotii pentru solicitarea de incovoiere
niu2			=		  0.3;		% Coeficientul lui Poisson pentru materialul rotii
E2              =		  2.05e5;   % Modulul de elasticitate longitudinal pentru materialul rotii


ZX              =		  1;		% Factorul de forma
Zw              =	      1;		% Factorul raportului duritatilor flancurilor dintilor

NstH			=	      1.e5;		%
NBH             =		  5.e7;		%
ZNmax           =		  1.6;		%

NstF			=	      1.e4;		%
NBF             =		  3.e6;		%
YNmax           =		  2.5;		%

SHmin           =		  1.15;		% Coeficientul de siguranta minim pentru solicitarea de contact
SFmin           =		  1.25;		% Coeficientul de siguranta minim pentru solicitarea de incovoiere

eta_rul         =	      0.99;		% Randamentul unei perechi de rulmenti
eta_angr		=		  0.97;		% Randamentul amgrenajului cilindric cu dinti inclinati
eta_ungere      =         0.995;	% randamentul datorita antrenarii uleiului (pt o trepata)

KA              =		  1.25;		% Factorul regimului de functionare

romat           =	      7.85e-6;  % Densitatea materialului capacelor, kg/mm^3

valid_angrenaj = 1;

%%%%%%%%%%%%%%%

while (valid_angrenaj == 1)
    aw = Lista_aw(iaw);
    beta = beta_g * PI / 180; % Unghiul de inclinare al danturii pe cilindrul de divizare, treapta I, rad

    %Cred ca trebuie aproximat prin scadere
    z2 = int32(round(z1 * i12STAS, 0)); % Numarul de dinti ai rotii z2
    u12 = z2 * 1 / z1; % Rapotrul de real de angrenare al primei trepte. Atentie la tipul variabilelor. am pus z2*1./z1 ca sa-l fac real

    mn = GetList(2 * aw * cos(beta) / z1 / (1 + u12), 41, Lista_mn);    %%Aici am facut modificare

    n2 = n1 / u12; % Turatia rotii z2


    T1 = 3.e7 * P * eta_rul / n1 / PI;				% Momentul de torsiune de pe arborele 1, N*mm
    T2 = T1 * eta_rul * eta_angr * eta_ungere * u12; % Momentul de torsiune de pe arborele 2, N*mm

    b2 = floor(psi_a * aw); % Latimea rotii, mm
    ExternalHelicalDriveUnit(mn,aw,z1,z2,xn1,beta,b2,Delta_b,T1,u12,KA,Lh1,Lh2,ZL1,ZL2,ZX,Zw,ZNmax,YNmax,NstH,NBH,NstF,NBF,SHmin, SFmin, ...
                             YR1,YR2,Hi1,Hi2,niu1,niu2,E1,E2,MatTT1,MatTT2,HB1,HB2,sigma021,sigma022,sigma_Hlim1,sigma_Hlim2,sigma_Flim1, ...
                             sigma_Flim2,Raf1,Raf2,Rar1,Rar2,PozAngr,CPrec,n1,n2,EHDU1);
    
    PosibilAngrenaj = int16(EHDU1(1));
   
    if (PosibilAngrenaj == 0) 
        valid_angrenaj = 0;
        break
    
    %%%% Aici am modificat algoritmul 
        %return;
    end

% Elementele geometrice ale angrenajului de pe treapta I
    alfa_t      = EHDU1(2);
    alfa_at1    = EHDU1(3);
    beta_a1     = EHDU1(4);
    beta_b      = EHDU1(5);
    d1          = EHDU1(6);
    db1         = EHDU1(7);
    df1         = EHDU1(8);
    da1         = EHDU1(9);
    st1         = EHDU1(10);
    sat1        = EHDU1(11);
    sn1         = EHDU1(12);
    san1        = EHDU1(13);
    xt1         = EHDU1(14);
    ADinte1     = EHDU1(15);
    a           = EHDU1(16);
    alfa_wt     = EHDU1(17);
    xsn         = EHDU1(18);
    xn2         = EHDU1(19);
    alfa_at2    = EHDU1(20);
    beta_a2     = EHDU1(21);
    d2          = EHDU1(22);
    db2         = EHDU1(23);
    df2         = EHDU1(24);
    da2         = EHDU1(25);
    st2         = EHDU1(26);
    sat2        = EHDU1(27);
    sn2         = EHDU1(28);
    san2        = EHDU1(29);
    xt2         = EHDU1(30);
    ADinte2     = EHDU1(31);
    dw1         = EHDU1(32);
    dw2         = EHDU1(33);
    beta_w      = EHDU1(34);
    eps_alfa    = EHDU1(35);
    eps_beta    = EHDU1(36);
    eps_gama    = EHDU1(37);
    b1          = EHDU1(38);
    cs_n        = EHDU1(39);
    alfa_n      = EHDU1(40);

    % Elementele geometrice ale angrenajului echivalent pentru treapta I
    zn1         = EHDU1(41);
    zn2         = EHDU1(42);
    xn1min      = EHDU1(43);
    xn2min      = EHDU1(44);
    dn1         = EHDU1(45);
    dn2         = EHDU1(46);
    dbn1        = EHDU1(47);
    dbn2        = EHDU1(48);
    dan1        = EHDU1(49);
    dan2        = EHDU1(50);
    alfa_wn     = EHDU1(51);
    awn         = EHDU1(52);
    eps_alfa_n  = EHDU1(53);

    % Elementele de control ale angrenajului treptei I
    alfa_Nt1    = EHDU1(54);
    alfa_Nt2    = EHDU1(55);
    Ncalc1      = EHDU1(56);
    Ncalc2      = EHDU1(57);
    N1          = EHDU1(58);
    N2          = EHDU1(59);
    WNn1        = EHDU1(60);
    WNn2        = EHDU1(61);
    WNt1        = EHDU1(62);
    WNt2        = EHDU1(63);
    roNt1       = EHDU1(64);
    roNt2       = EHDU1(65);
    roAt1       = EHDU1(66);
    roEt2       = EHDU1(67);
    roat1       = EHDU1(68);
    roat2       = EHDU1(69);

    % Tensiuni si tensiuni admisibile pentru materialele treptei I
    sigma_H     = EHDU1(70);
    sigma_HP    = EHDU1(71);
    sigma_F1    = EHDU1(72);
    sigma_F2    = EHDU1(73);
    sigma_FP1   = EHDU1(74);
    sigma_FP2   = EHDU1(75);

    % Volumul danturiilor, mm^3
    VolDantura1 = EHDU1(76);
    VolDantura2 = EHDU1(77);
   
    % Forte si momente
    Ft1         = EHDU1(78);
    Fa1         = EHDU1(79);
    Fr1         = EHDU1(80);
    Mic1        = EHDU1(81); MH7=Mic1;

    Ft2         = EHDU1(82);
    Fa2         = EHDU1(83);
    Fr2         = EHDU1(84);
    Mic2        = EHDU1(85); MH3=Mic2;
    break
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Masaz1 = (VolDantura1 + Volum_cilindru(df1,b1)) * romat;
Masaz2 = (VolDantura2 + Volum_cilindru(df2,b2)) * romat;

gearsMass = Masaz1 + Masaz2;

save('temporary.mat', 'PosibilAngrenaj', 'u12', 'i12STAS', 'sigma_H', 'sigma_HP', 'sigma_F1', ...
'sigma_FP1', 'sigma_F2', 'sigma_FP2', 'xn1', 'zn1', 'xn2', ...
'zn2', 'cs_n', 'mn', 'san1', 'san2', 'eps_alfa_lim', 'eps_alfa', ...
'WNn1', 'beta_b', 'b1', 'WNn2', 'b2', 'roAt1', 'roNt1', 'roat1', ...
'roEt2', 'roNt2', 'roat2');  

end


