// Angrenaj cu roti dintate cilindrice cu dinti inclinati - 15.07.2011



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <array>
#include <vector>
using namespace std;

double RandomNumber()
// Functia returneaza un numar aleator intre 0 si 1
{
	return 1.0*rand()/RAND_MAX;

}

double Maximum(double x,double y)
{
	if (x<=y) return y;
	else return x;
}

double Minimum (double x,double y)
{
	if (x<=y) return x;
	else return y;
}

double sign(double x)
{
	if (x>0) return 1.;
	else if (x<0) return -1.;
	else return 0.;
}

double round(double x,int prec)
// Functia returneaza rotunjirea numarului real "x" (in sus sau in jos, dupa caz) cu "prec" zecimale
// Pentru prec=0, "x" se rotunjeste la cel mai apropiat intreg
{
	int m,i;
	double a,b;

	m=1;
	if (prec!=0)
		{
			for (i=0;i<prec;i++) m*=10;
		}
	a=floor(fabs(x)*m);
	b=floor(fabs(x)*m*10);
	if (b-10*a<5) return sign(x)*a/m;
	else return sign(x)*(a+1)/m;
}

int Verif_cmmdc(int a,int b)
// Daca numerele a si b sunt prime intre ele functia returneaza -1 - in caz contrar returneaza +1
{
	int c,conditie; int cmmdc=1;
    int lim;

    lim=(a<b)?a:b;

	for (c=2; c<=lim; c++)
	{
		if ((a%c==0)&&(b%c==0))
			{
			cmmdc=c;
			break;
			}
	}
	conditie=1;
	if (cmmdc ==1) conditie=-1;
	return(conditie);
}

double GetList(double x, int NrElemList,double List [])
// Gaseste elementul cel mai apropiat de x din sirul List (care are NrElemList elemente)
{
	double dif_min,dif_i,result;
	int i;
    cout << __func__ << '\n';
	dif_min=1.e10;
	for(i=0;i<NrElemList;i++)
	{
		dif_i=fabs(x-List[i]);
		if (dif_i<dif_min)
		{
			result=List[i];
			dif_min=dif_i;
		}
	}
	return result;
}


double Lungime_arc_cerc (double R,double u)
// Functia returneaza lungimea arcului de cerc (mm) cu raza R (mm) si unghiul la centru u (rad)
{
	return R*u;
}

double ASect(double d,double alfa)
// Returneaza aria unui sector de cerc cu diametrul d si unghiul la centru alfa [rad]
{
		return (d*d*alfa/8.);
}

double ord(double x1,double y1,double x2,double y2,double x)
// Returneaza ordonata punctului de abscisa x de pe dreapta determinata de punctele M1(x1,y1) si M2(x2,y2)
{
	return(y1+(x-x1)*(y2-y1)/(x2-x1));
}

double Volum_cilindru(double D,double l)
// Returneaza volumul unui cilindru, mm^3
// D  - Diametrul, mm
// l  - Lungimea, mm
{
	double PI=3.1415926535897932384626433832795;
	return(D*D*PI*l/4);
}

double inv(double x)
// Involuta de unghiul x (in radiani)
{
	return(tan(x)-x);
}

double absEv(double db,double alfa)
// Returneaza abscisa unui punct de pe evolventa cu diametrulcercului de baza db [mm] si determinat de unghiul alfa [rad]
{
	return (db*(sin(tan(alfa))-tan(alfa)*cos(tan(alfa)))/2.);
}

double ordEv(double db,double alfa)
// Returneaza ordonata unui punct de pe evolventa cu diametrulcercului de baza db [mm] si determinat de unghiul alfa [rad]
{
	return (db*(cos(tan(alfa))+tan(alfa)*sin(tan(alfa)))/2.);
}

double primEv(double x)
// Primitiva f evolventice
{
	return(2*x*x*x+3*(1-x*x)*sin(2*x)-6*x*cos(2*x));
}

double ArieDinte(double da,double db,double df,double xn,int z,double sat,double alfa_n)
// Returneaza aria unui dinte cu profil evolventice
// A-inceputul evolventei (deasupra racordarii)
// B-sfarsitul evolventei (pe capul dintelui)
{
	double PI;
	double alfa_a,alfa_u,alfa_AB,xA,yA,du,xB,yB,psi;
	double AOCA,AODB,AOAM,AOMN,AOBE,AOAF,AOHK,ABMNE,AAHKF;
	double ror,AOAo,fi,AOTH,AoTA,AR,Ae,Aev,AAMB;

	PI=3.1415926535897932384626433832795;
    cout << __func__ << '\n';
	alfa_a=atan(sqrt(da*da-db*db)/db); // Unghiul de sfarsit al evolventei, rad
	alfa_u=atan((z*sin(alfa_n)*sin(alfa_n)-2*(1-xn))/z/cos(alfa_n)/sin(alfa_n)); // Unghiul punctului de inceput al evolventei, rad
	xA=absEv(db,alfa_u);
	yA=ordEv(db,alfa_u);
	xB=absEv(db,alfa_a);
	yB=ordEv(db,alfa_a);
	du=2*sqrt(xA*xA+yA*yA); // Diametrul cercului de inceput al evolventei, mm
	psi=2.*sat/da;
	alfa_AB=acos((xA*xB+yA*yB)/sqrt(xA*xA+yA*yA)/sqrt(xB*xB+yB*yB));
	AOCA=xA*yA/2.;
	AODB=xB*yB/2.;
	AOAM=ASect(du,alfa_AB);
	AOMN=ASect(du,psi);
	AOBE=ASect(da,psi);
	AOAF=ASect(du,2*alfa_AB+psi);
	AOHK=ASect(df,2*alfa_AB+psi);
	ABMNE=AOBE-AOMN;
	AAHKF=AOAF-AOHK;

	ror=(du*du-df*df)/df/4.; // Raza aprox. a racordarii, mm
	AOAo=ror*du/4.;
	fi=asin(2*ror/(df+2*ror));
	AOTH=ASect(df,fi);
	AoTA=ASect(2*ror,PI/2.-fi);
	AR=AOAo-AOTH-AoTA; // Aria aprox a zonei de racord a dintelui

	Ae=db*db*(primEv(tan(alfa_a))-primEv(tan(alfa_u)))/48.;	// Aria marginita de evolventa
	Aev=Ae+AOCA;

	AAMB=Aev-AODB-AOAM;

	return(2.*AAMB+ABMNE+AAHKF+2.*AR);
}

double grafparam(double Gu[],int lin,int col,double P[],double parameter,double x)
{
    cout << lin << ' ' << col << ' ' << parameter << ' ' << x << ' ';


	int i, j;
	double result, Y1, Y2;
	int valid, found;
	found=0;
	result=0;
	valid=1;
    cout << __func__ << ' ';

	for (i=0; i<lin/2; i++)
	{
		if (P[i]==parameter)
		{

			for (j=0; j<col; j++)
			{
				if (Gu[2*i*col+j]==x)
				{
				    cout << "if-ul 1" << "Gu[2*i*col+j]" << ' ' << "Gu[" << 2*i*col+j << "] ";
					result=Gu[(2*i+1)*col+j];
					found=1;
					break;
				}
				if ((j!=col-1) && (Gu[2*i*col+j]<x) && (Gu[2*i*col+j+1]>x))
				{
				    cout << "if-ul 2 " << Gu[2*i*col+j] << ' ' <<
                    i << ' ' << col << ' ' << j << ' ';
					result=ord(Gu[2*i*col+j],Gu[(2*i+1)*col+j],Gu[2*i*col+j+1],Gu[(2*i+1)*col+j+1],x);
					found=1;
					break;
				}
			}
		}
		if ((found==1)||(i==lin/2-1)) break;

		if (((P[i]<parameter)&&(P[i+1]>parameter))||((P[i+1]<parameter)&&(P[i]>parameter)))
		{
		    cout << 2*i*col << ' ' << (2*i+2)*col << " max " << Maximum(Gu[2*i*col],Gu[(2*i+2)*col]) << " min " << Minimum(Gu[(2*i+1)*col-1],Gu[(2*i+3)*col-1]) << ' ';
			if((x<Maximum(Gu[2*i*col],Gu[(2*i+2)*col]))||(x>Minimum(Gu[(2*i+1)*col-1],Gu[(2*i+3)*col-1])))
			{
				valid=0;
				break;
			}
			else
			{
				for (j=0; j<col-1; j++)
				{
					if ((Gu[2*i*col+j]<=x) && (Gu[2*i*col+j+1]>=x))
					{
					    cout << ' ' <<i << " else 1 for 1 ";
					    cout << Gu[2*i*col+j] << ' ' << Gu[2*i*col+j+1] << ' ' << i << ' ' << j << ' ';
						Y1=ord(Gu[2*i*col+j],Gu[(2*i+1)*col+j],Gu[2*i*col+j+1],Gu[(2*i+1)*col+j+1],x);
						break;
					}
				}
				for (j=0; j<col-1; j++)
				{
					if ((Gu[(2*i+2)*col+j]<=x) && (Gu[(2*i+2)*col+j+1]>=x))
					{
					    cout << "else 1 for 2 ";
					    cout << Gu[(2*i+2)*col+j] << ' ' << Gu[(2*i+2)*col+j+1] << ' ' << i << ' ' << j << ' ';
						Y2=ord(Gu[(2*i+2)*col+j],Gu[(2*i+3)*col+j],Gu[(2*i+2)*col+j+1],Gu[(2*i+3)*col+j+1],x);
						break;
					}
				}
				result=ord(P[i],Y1,P[i+1],Y2,parameter);
				found=1;
				break;
			}
		}
	}

	if ((valid==1)&&(found==0))
	{
		if (lin!=2)
		{
			if((x<Maximum(Gu[0],Gu[(lin-1)*col]))||(x>Minimum(Gu[col-1],Gu[(lin-1)*col-1]))) valid=0;
			else
			{
				for (j=0; j<col-1; j++)
				{
					if ((Gu[0*col+j]<=x) && (Gu[0*col+j+1]>=x))
					{
					    cout << "else 2 if 1 ";

						Y1=ord(Gu[0*col+j],Gu[1*col+j],Gu[0*col+j+1],Gu[1*col+j+1],x);
						break;
					}
				}
				for (j=0; j<col-1; j++)
				{
					if ((Gu[(lin-2)*col+j]<=x) && (Gu[(lin-2)*col+j+1]>=x))
					{
					    cout << "else 2 if 1 2";
						Y2=ord(Gu[(lin-2)*col+j],Gu[(lin-1)*col+j],Gu[(lin-2)*col+j+1],Gu[(lin-1)*col+j+1],x);
						break;
					}
				}
				result=ord(P[0],Y1,P[lin/2-1],Y2,parameter);
				found=1;
			}
		}
		else // else 3
		{
		    cout << "intrare " << col - 1 << ' ' << Gu[col-1] << ' ';
			if ((x<Gu[0])||(x>Gu[col-1])) valid=0;

			else
			{
				for (j=0; j<col-1; j++)
				{
					if ((Gu[j]<=x) && (Gu[j+1]>=x))
					{
					    cout << Gu[col-1] << ' ';
					    cout << " indecsi " << j << ' ' << col + j << ' ' << j + 1 << ' ' << col + j + 1 << '\n';
					    cout << "valori " << "else 4 if 1 " << Gu[j] << ' ' << Gu[col+j] << ' ' << Gu[j+1] << ' ' << Gu[col+j+1] << ' ';
						result=ord(Gu[j],Gu[col+j],Gu[j+1],Gu[col+j+1],x);
						found=1;
						break;
					}
				}
			}
		}
	}
	cout << result << '\n';
	return result;
}

double alfa_yt(double dy,double d,double alfa_t)
// Returneaza valoarea unghiului de presiune al profilului pe un cerc de diametru oarecare in plan frontal, rad
// dy - cercul de diametrul oarecare, mm
// d  - diametrul de divizare, mm
{
	return(acos(d*cos(alfa_t)/dy));
}

double beta_y(double dy,double d,double beta)
// Returneaza valoarea unghiului de inclinare al danturii pe un cilindru oarecare
// dy - cilindrul de diametrul oarecare, mm
// d  - diametrul de divizare, mm
// beta - unghiul de inclinare al danturii pe cilindrul de divizare, rad
{
	return(atan(dy*tan(beta)/d));
}

double syt(double mn,double z,double alfa_t,double alfa_yt,double beta,double xn)
// Returneaza valoarea arcului dintelui pe un cerc oarecare in plan frontal
// mn - modulul normal, mm
// z  - numarul de dinti
//
// st - lungimea arcului dintelui pe cercul de divizare in plan frontal
// alfa_yt - unghiul de presiune al profilului pe un cerc de diametru oarecare in plan frontal, rad
// beta - unghiul de inclinare al danturii pe cilindrul de divizare, rad
{
	double PI,st;
	PI=3.1415926535897932384626433832795;
	st=(0.5*PI/cos(beta)+2*xn*tan(alfa_t))*mn;

	return(((inv(alfa_t)-inv(alfa_yt))*mn*z/cos(beta)+st)*cos(alfa_t)/cos(alfa_yt));
}

double syn(double syt,double beta_y)
// Returneaza valoarea arcului dintelui pe un cerc oarecare in plan normal
{
	return(syt*cos(beta_y));
}

void HelicalGear(double alfa_n,double ha_n,double cs_n,double mn,int z,double xn,double beta,double alfa_t,double Delta_yn,double *HG)
// Returneaza in HG geomatria unei r.d.c.d.i.
// ha_n		- Coeficientul inaltimii capului dintelui
// cs_n		- Coeficientul jocului la capul dintelui de referinta
// mn		- Modulul normal, mm (in cazul r.d.c.c.d.d. modulul normal este modulul m al rotii)
// z		- Numarul de dinti
// xn		- Coeficientul deplasarii de profil in plan normal, (in cazul r.d.c.c.d.d. xn este x-ul rotii)
// beta		- Unghiul de inclinare al danturii pe cilindrul de divizare, grade
// Delta_yn - Coeficientul de scurtare al dintilor in plan normal
{
	double PI,alfa_at,beta_a,beta_b,xt,d,db,df,da,sn,san,st,sat,temp;
	int validHG;
	validHG=1;

	PI=3.1415926535897932384626433832795;

	xt=xn*cos(beta); // Coeficientul deplasarii de profil al dintelui pinionului, in plan frontal
	d=mn*z/cos(beta); // Diametrul cercului de divizare, mm
	db=d*cos(alfa_t); // Diametrul cercului de baza, mm

	while(validHG==1){
            cout << __func__ << '\n';
	df=mn*(z/cos(beta)-2*(ha_n+cs_n-xn)); // Diametrul cercului de picior, mm
	if (df<=0.) {validHG=0;break;}

	da=mn*(z/cos(beta)+2*(ha_n+xn-Delta_yn)); // Diametrul cercului de cap, mm
	if (da<=0.) {validHG=0;break;}

	temp=fabs(d*cos(alfa_t)/da);
	if (temp>1.) {validHG=0;break;}
	alfa_at=alfa_yt(da,d,alfa_t); // Unghiul de presiune de referinta pe cercul de cap, rad

	st=syt(mn,z,alfa_t,alfa_t,beta,xn); // Lungimea arcului dintelui pe cercul de divizare, in plan frontal, mm
	if (st<=0.) {validHG=0;break;}

	sat=syt(mn,z,alfa_t,alfa_at,beta,xn); // Lungimea arcul dintelui pe cercul de cap in plan frontal, mm
	if (sat<=0.) {validHG=0;break;}

	beta_a=beta_y(da,d,beta); // Unghiul de inclinare a danturii pe cilindrul de cap, rad

	beta_b=beta_y(db,d,beta); // Unghiul de inclinare a danturii pe cilindrul de baza, rad

	sn=syn(st,beta); // Lungimea arcului dintelui pe cercul de divizare, in plan normal, mm
	if (sn<=0.) {validHG=0;break;}

	san=syn(sat,beta_a); // Lungimea arcul dintelui pe cercul de cap in plan normal, mm
	if (san<=0.) {validHG=0;break;}

	break;
	}

	if(validHG)
	{

		HG[0]=1;
		HG[1]=xt;
		HG[2]=d;
		HG[3]=db;
		HG[4]=df;
		HG[5]=da;
		HG[6]=alfa_at;
		HG[7]=st;
		HG[8]=sat;
		HG[9]=beta_a;
		HG[10]=beta_b;
		HG[11]=sn;
		HG[12]=san;
		HG[13]=ArieDinte(da,db,df,xn,z,sat,alfa_n);

	}
	else HG[0]=0;
}

void HelicalDriveGeometry(double mn,double aw,int z1,int z2,double xn1,double beta,double b2,double Delta_b,double *HDG)
// Returneaza in HDG geomatria unui angrenaj cu r.d.c.d.i.
// mn      - Modulul normal, mm (in cazul r.d.c.c.d.d. modulul normal este modulul m)
// aw      - Distanta axiala standardizata/impusa, mm
// z1,z2   - Numarele de dinti
// xn1     - Coeficientul deplasarii de profil in plan normal, pentru pinion (in cazul r.d.c.c.d.d. xn1 este x1)
// beta    - Unghiul de inclinare al danturii pe cilindrul de divizare, grade
// b2	   - Latimea danturii rotii, mm
// Delta_b - Diferenta de latime dintre latimea pinionului si a rotii, mm (Delta_b=5)
{
	double alfa_n_g,ha_n,cs_n;
	double PI;
	double temp;
	double HG1[14],HG2[14];
	double alfa_t,alfa_n,beta_b,beta_w,a,alfa_wt,xsn,eps_alfa,eps_beta,eps_gama,Delta_yn;
	double xt1,alfa_at1,d1,db1,df1,da1,dw1,st1,sat1,beta_a1,sn1,san1,ADinte1,b1,VolDantura1;
	double xn2,xt2,alfa_at2,d2,db2,df2,da2,dw2,st2,sat2,beta_a2,sn2,san2,ADinte2,VolDantura2;
	int validHDG;

	// Cremaliera de referinta ISO53 (STAS 821)

	alfa_n_g=20; // Unghiul de presiune de referinta in plan normal, grade
	ha_n=1;		// Coeficientul inaltimii capului dintelui
	cs_n=0.25;	// Coeficientul jocului la capul dintelui de referinta



	PI=3.1415926535897932384626433832795;
	alfa_n=alfa_n_g*PI/180.;
	alfa_t=atan(tan(alfa_n)/cos(beta)); // Unghiul de angrenare de referinta in plan frontal, rad
	validHDG=1;
	b1=b2+Delta_b; // Latimea pinionului, mm
	while(validHDG==1){
            cout << __func__ << '\n';
	a=mn*(z1+z2)/cos(beta)/2.; // Distanta axiala elementara, mm
	if (((aw-a)*cos(beta)/mn<-0.4)||((aw-a)*cos(beta)/mn>2.5))  {validHDG=0;break;}
	temp = fabs(a*cos(alfa_t)/aw);
	if (temp > 1.) {validHDG=0;break;}
	alfa_wt=acos(a*cos(alfa_t)/aw); // Unghiul real  de  angrenare  in  plan  frontal, rad

	xsn=(inv(alfa_wt)-inv(alfa_t))*(z1+z2)/tan(alfa_n)/2.; // Suma coeficientilor  deplasarilor  de  profil  in  plan  normal
	xn2=xsn-xn1; // Coeficientul deplasarii de profil al dintelui rotii, in plan normal

	Delta_yn=xsn-(aw-a)/mn;

	HelicalGear(alfa_n,ha_n,cs_n,mn,z1,xn1,beta,alfa_t,Delta_yn,HG1);
	if(HG1[0]==0){validHDG=0;break;}

	xt1=HG1[1]; // Coeficientul deplasarii de profil al dintelui pinionului, in plan frontal
	d1=HG1[2]; // Diametrul cercului de divizare pentru pinion, mm
	db1=HG1[3]; // Diametrul cercului de baza pentru pinion, mm
	df1=HG1[4]; // Diametrul cercului de picior pentru pinion, mm
	da1=HG1[5]; // Diametrul cercului de cap pentru pinion, mm
	alfa_at1=HG1[6]; // Unghiul de presiune de referinta pe cercul de cap al pinionului, rad
	st1=HG1[7]; // Lungimea arcului dintelui pe cercul de divizare al pinionului, in plan frontal, mm
	sat1=HG1[8]; // Lungimea arcul dintelui pe cercul de cap al pinionului in plan frontal, mm
	beta_a1=HG1[9]; // Unghiul de inclinare a danturii pe cilindrul de cap al pinionului, rad
	beta_b=HG1[10]; // Unghiul de inclinare a danturii pe cilindrul de baza, rad
	sn1=HG1[11]; // Lungimea arcului dintelui pe cercul de divizare al pinionului, in plan normal, mm
	san1=HG1[12]; // Lungimea arcul dintelui pe cercul de cap al pinionului in plan normal, mm
	ADinte1=HG1[13]; // Aria suprafetei frontale a dintelui pinionului, mm^2

	HelicalGear(alfa_n,ha_n,cs_n,mn,z2,xn2,beta,alfa_t,Delta_yn,HG2);
	if(HG2[0]==0){validHDG=0;break;}
	xt2=HG2[1]; // Coeficientul deplasarii de profil al dintelui rotii, in plan frontal
	d2=HG2[2]; // Diametrul cercului de divizare pentru roata, mm
	db2=HG2[3]; // Diametrul cercului de baza pentru roata, mm
	df2=HG2[4]; // Diametrul cercului de picior pentru roata, mm
	da2=HG2[5]; // Diametrul cercului de cap pentru roata, mm
	alfa_at2=HG2[6]; // Unghiul de presiune de referinta pe cercul de cap al rotii, rad
	st2=HG2[7]; // Lungimea arcului dintelui pe cercul de divizare al rotii, in plan frontal, mm
	sat2=HG2[8]; // Lungimea arcul dintelui pe cercul de cap al rotii in plan frontal, mm
	beta_a2=HG2[9]; // Unghiul de inclinare a danturii pe cilindrul de cap al rotii, rad
	beta_b=HG2[10]; // Unghiul de inclinare a danturii pe cilindrul de baza, rad
	sn2=HG2[11]; // Lungimea arcului dintelui pe cercul de divizare al rotii, in plan normal, mm
	san2=HG2[12]; // Lungimea arcul dintelui pe cercul de cap al rotii in plan normal, mm
	ADinte2=HG2[13]; // Aria suprafetei frontale a dintelui rotii, mm^2

	VolDantura1=ADinte1*z1*b1;
	VolDantura2=ADinte2*z2*b2;

	dw1=d1*cos(alfa_t)/cos(alfa_wt); // Diametrul cercului de rostogolire pentru pinion, mm
	dw2=d2*cos(alfa_t)/cos(alfa_wt); // Diametrul cercului de rostogolire pentru roata, mm

	beta_w=atan(dw1*tan(beta)/d1); // Unghiul de inclinare a danturii pe cilindrul de rostogolire, rad

	eps_alfa=(sqrt(da1*da1-db1*db1)+sqrt(da2*da2-db2*db2)-2*aw*sin(alfa_wt))*cos(beta)/cos(alfa_t)/mn/PI/2.; // Gradul de acoperie frontal
	if (eps_alfa<=0) {validHDG=0;break;}

	eps_beta=b2*sin(beta)/mn/PI; // Gradul de acoperire axial
	eps_gama=eps_alfa+eps_beta; // Gradul de acoperire total
	break;
	}

	if(validHDG)
	{

		HDG[0]=1;
		HDG[1]=alfa_t;
		HDG[2]=alfa_at1;
		HDG[3]=beta_a1;
		HDG[4]=beta_b;
		HDG[5]=d1;
		HDG[6]=db1;
		HDG[7]=df1;
		HDG[8]=da1;
		HDG[9]=st1;
		HDG[10]=sat1;
		HDG[11]=sn1;
		HDG[12]=san1;
		HDG[13]=xt1;
		HDG[14]=ADinte1;
		HDG[15]=a;
		HDG[16]=alfa_wt;
		HDG[17]=xsn;
		HDG[18]=xn2;
		HDG[19]=alfa_at2;
		HDG[20]=beta_a2;
		HDG[21]=d2;
		HDG[22]=db2;
		HDG[23]=df2;
		HDG[24]=da2;
		HDG[25]=st2;
		HDG[26]=sat2;
		HDG[27]=sn2;
		HDG[28]=san2;
		HDG[29]=xt2;
		HDG[30]=ADinte2;
		HDG[31]=dw1;
		HDG[32]=dw2;
		HDG[33]=beta_w;
		HDG[34]=eps_alfa;
		HDG[35]=eps_beta;
		HDG[36]=eps_gama;
		HDG[37]=b1;
		HDG[38]=cs_n;
		HDG[39]=alfa_n;
		HDG[40]=VolDantura1;
		HDG[41]=VolDantura2;

	}
	else HDG[0]=0;
}

void EquivalentHelicalDrive(double mn,double a,int z1,int z2,double beta,double beta_b,double beta_w,double alfa_n,double alfa_wt,double d1,double da1,double d2,double da2,double *EHD)
// Returneaza elementele geometrice ale angrenajului echivalent cu r.d.c.c.d. al unui angrenaj cu r.d.c.d.i. echivalarea in plan normal
// mn     - Modulul normal, mm (in cazul r.d.c.c.d.d. modulul normal este modulul m)
// a	  - Distanta axiala standardizata/impusa, mm
// z1,z2  - Numarele de dinti
// beta   - Unghiul de inclinare al danturii pe cilindrul de divizare, grade
// beta_b - Unghiul de inclinare al danturii pe cilindrul de baza, rad
// beta_w - Unghiul de inclinare al danturii pe cilindrul de rostogolire, rad
// alfa_n - Unghiul de presiune de referinta in plan normal, grade
// alfa_wt- Unghiul real de angrenare in plan frontal, rad
// d1     - Diametrul cercului de divizare corespunzator pinionului, mm
// da1	  - Diametrul cercului de cap corespunzator pinionului, mm
// d2     - Diametrul cercului de divizare corespunzator rotii dintate, mm
// da2	  - Diametrul cercului de cap corespunzator rotii dintate, mm
{
	double PI;
	double zn1,zn2,xn1min,xn2min,dn1,dn2,dbn1,dbn2,dan1,dan2;
	double alfa_wn,awn,eps_alfa_n;
	double alfa_n_g;

	alfa_n_g=20; // Unghiul de presiune de referinta in plan normal, grade
	PI=3.1415926535897932384626433832795;
	alfa_n=alfa_n_g*PI/180.;
	zn1=z1/cos(beta)/cos(beta)/cos(beta); // Numarul de dinti ai pinionului echivalent
	zn2=z2/cos(beta)/cos(beta)/cos(beta); // Numarul de dinti ai rotii echivalente
	xn1min=(14-zn1)/17.; // Valoarea minima a coeficientului deplasarii de profil a dintelui pinionului, in plan normal
	xn2min=(14-zn2)/17.; // Valoarea minima a coeficientului deplasarii de profil a dintelui rotii, in plan normal
	dn1=mn*zn1; // Diametrul cercului de divizare al pinionului echivalent, mm
	dn2=mn*zn2; // Diametrul cercului de divizare al rotii echivalente, mm
	dbn1=dn1*cos(alfa_n); // Diametrul cercului de baza all pinionului echivalent, mm
	dbn2=dn2*cos(alfa_n); // Diametrul cercului de baza a rotii echivalente, mm
	dan1=dn1+da1-d1; // Diametrul cercului de cap al pinionului echivalent, mm
	dan2=dn2+da2-d2; // Diametrul cercului de cap al rotii echivalente, mm
	alfa_wn=acos(cos(alfa_wt)*cos(beta_b)/cos(beta_w)); // Unghiul de presiune al angrenajului echivalent, rad
	awn=a*cos(alfa_n)/cos(beta_b)/cos(beta_b)/cos(alfa_wn); // Distanta dintre axe a angrenajului echivalent, mm
	eps_alfa_n=(sqrt(dan1*dan1-dbn1*dbn1)+sqrt(dan2*dan2-dbn2*dbn2)-2*awn*sin(alfa_wn))/cos(alfa_n)/mn/PI/2.; // Gradul de acoperie frontal al angrenajului echivalent

	EHD[0]=1;
	EHD[1]=zn1;
	EHD[2]=zn2;
	EHD[3]=xn1min;
	EHD[4]=xn2min;
	EHD[5]=dn1;
	EHD[6]=dn2;
	EHD[7]=dbn1;
	EHD[8]=dbn2;
	EHD[9]=dan1;
	EHD[10]=dan2;
	EHD[11]=alfa_wn;
	EHD[12]=awn;
	EHD[13]=eps_alfa_n;
}

void ControlHelicalDrive(double mn,double aw,int z1,int z2,double xn1,double xn2,double beta,double beta_b,double beta_w,double alfa_n,double alfa_t,double alfa_wt,double alfa_at1,double alfa_at2,double d1,double da1,double d2,double db1,double db2,double da2,double *CHD)
// Returneaza elementele de control ale angrenajului cu r.d.c.d.i
// mn     - Modulul normal, mm (in cazul r.d.c.c.d.d. modulul normal este modulul m)
// aw	  - Distanta axiala standardizata/impusa, mm
// z1,z2  - Numarele de dinti
// beta   - Unghiul de inclinare al danturii pe cilindrul de divizare, grade
// beta_b - Unghiul de inclinare al danturii pe cilindrul de baza, rad
// beta_w - Unghiul de inclinare al danturii pe cilindrul de rostogolire, rad
// alfa_n - Unghiul de presiune de referinta in plan normal, grade
// alfa_wt- Unghiul real de angrenare in plan frontal, rad
// d1     - Diametrul cercului de divizare corespunzator pinionului, mm
// da1	  - Diametrul cercului de cap corespunzator pinionului, mm
// d2     - Diametrul cercului de divizare corespunzator rotii dintate, mm
// da2	  - Diametrul cercului de cap corespunzator rotii dintate, mm
{
	double PI,temp;
	double alfa_Nt1,Ncalc1,N1,WNn1,WNt1,roat1,roNt1,roAt1;
	double alfa_Nt2,Ncalc2,N2,WNn2,WNt2,roat2,roNt2,roEt2;
	int validCHD;

	validCHD=1;
	PI=3.1415926535897932384626433832795;

	while(validCHD==1){
            cout << __func__ << '\n';
	temp = fabs(z1*cos(alfa_t)/(z1+2.*xn1*cos(beta)));
	if (temp>1.) {validCHD=0; break;}
	alfa_Nt1=acos(z1*cos(alfa_t)/(z1+2.*xn1*cos(beta))); // rad
	temp = fabs(z2*cos(alfa_t)/(z2+2.*xn2*cos(beta)));
	if (temp>1.) {validCHD=0; break;}
	alfa_Nt2=acos(z2*cos(alfa_t)/(z2+2.*xn2*cos(beta))); // rad
	Ncalc1=0.5+z1*(tan(alfa_Nt1)/cos(beta)/cos(beta)-2*xn1*tan(alfa_n)/z1-inv(alfa_t))/PI;
	Ncalc2=0.5+z2*(tan(alfa_Nt2)/cos(beta)/cos(beta)-2*xn2*tan(alfa_n)/z2-inv(alfa_t))/PI;
	N1=floor(Ncalc1); // Numarul de dinti pentru masurarea cotei peste dinti, pentru pinion
	if (Ncalc1-floor(Ncalc1)>=0.5) N1+=1;
	N2=floor(Ncalc2); // Numarul de dinti pentru masurarea cotei peste dinti, pentru roata
	if (Ncalc2-floor(Ncalc2)>=0.5) N2+=1;
	WNn1=2*xn1*mn*sin(alfa_n)+mn*cos(alfa_n)*((N1-0.5)*PI+z1*inv(alfa_t)); // Cota peste N1 dinti in plan normal pt angrenaje fara joc intre flancuri, pentru pinion, mm
	WNn2=2*xn2*mn*sin(alfa_n)+mn*cos(alfa_n)*((N2-0.5)*PI+z2*inv(alfa_t)); // Cota peste N2 dinti in plan normal pt angrenaje fara joc intre flancuri, pentru roata, mm
	WNt1=WNn1/cos(beta_b); // Cota peste N1 dinti in plan frontal pt angrenaje fara joc intre flancuri, pentru pinion, mm
	WNt2=WNn2/cos(beta_b); // Cota peste N2 dinti in plan frontal pt angrenaje fara joc intre flancuri, pentru roata, mm
	roNt1=0.5*WNt1; // Raza de curbura a profilului in punctele simetrice de masurare a lungimii peste N1 dinti in planul frontal, pentru pinion, mm
	roNt2=0.5*WNt2; // Raza de curbura a profilului in punctele simetrice de masurare a lungimii peste N2 dinti in planul frontal, pentru roata, mm
	roAt1=aw*sin(alfa_wt)-0.5*db2*tan(alfa_at2); // Raza de curbura a profilului dintelui pinionului in punctul de intrare din angrenare, mm
	if (roAt1<=0.) {validCHD=0; break;}
	roEt2=aw*sin(alfa_wt)-0.5*db1*tan(alfa_at1); // Raza de curbura a profilului dintelui rotii in punctul de iesire din angrenare, mm
	if (roEt2<=0.) {validCHD=0; break;}
	roat1=0.5*da1*sin(alfa_at1); // Raza de curbura a profilului la capul dintelui pinionului, mm
	roat2=0.5*da2*sin(alfa_at2); // Raza de curbura a profilului la capul dintelui rotii, mm
	break;
	}

	if(validCHD)
	{

		CHD[0]=1;
		CHD[1]=alfa_Nt1;
		CHD[2]=alfa_Nt2;
		CHD[3]=Ncalc1;
		CHD[4]=Ncalc2;
		CHD[5]=N1;
		CHD[6]=N2;
		CHD[7]=WNn1;
		CHD[8]=WNn2;
		CHD[9]=WNt1;
		CHD[10]=WNt2;
		CHD[11]=roNt1;
		CHD[12]=roNt2;
		CHD[13]=roAt1;
		CHD[14]=roEt2;
		CHD[15]=roat1;
		CHD[16]=roat2;

	}
	else CHD[0]=0;
}

void HelicalDriveStrength(double T1,double z1,double u12,double KA,double Lh1,double Lh2,double ZL1,double ZL2,double ZX,double Zw,double ZNmax,double YNmax,double NstH,double NBH,double NstF,double NBF,double SHmin,double SFmin,double YR1,double YR2,double Hi1,double Hi2,double niu1,double niu2,double E1,double E2,double MatTT1,double MatTT2,double HB1,double HB2,double sigma021,double sigma022,double sigma_Hlim1,double sigma_Hlim2,double sigma_Flim1,double sigma_Flim2,double Raf1,double Raf2,double Rar1,double Rar2,double PozAngr,double CPrec,double mn,double aw,double alfa_t,double alfa_wt,double beta,double beta_b,double eps_alfa,double eps_beta,double eps_alfa_n,double b1,double b2,double zn1,double zn2,double xn1,double xn2,double d1,double d2,double n1,double n2,double *HDS)
// T1			 - Momentul de torsiune corespunzator arborelui de intrare, Nmm
// z1			 - Numerele de dinti ai pinionului
// u12			 - Raportul deal de angrenare corespunzator treptei I
// KA			 - Factorul regimului de functionare
// Lh1,Lh2		 - Durata minima de functionare a pinionului/rotii dintate, ore
// ZL1,ZL2		 - Factorul de ungere corespunzator pinionului/rotii dintate
// Zx			 - Factorul de marime
// ZW			 - Factorul raportului duritatilor flancurilor dintilor
// ZNmax		 - Factorul durabilitatii pentru solicitarea de contact
// YNmax		 - Factorul durabilitatii pentru solicitarea de incovoiere
// NstH,NstF	 -
// NBH,NBF		 -
// SHmin		 - Coeficientul de siguranta minim pentru solicitarea de contact
// SFmin		 - Coeficientul de siguranta minim pentru solicitarea de incovoiere
// YR1,YR2		 - Factorul rugozitatii flancului dintelui pinionului/rotii dintate pentru solicitarea de incovoiere
// Hi1,Hi2		 - Numarul de roti cu care vine in contact pinionul/roata dintata
// niu1,niu2	 - Coeficientii lui Poisson
// E1,E2		 - Modulul de elasticitate longitudinal, MPa
// MatTT1,MatTT2 -
// HB1,HB2		 - Duritatea materialului pinionului/rotii dintate, MPa
// sigma021,2	 - Limita de curegere a materialului pinionului/rotii dintate, MPa
// sigma_Hlim1,2 - Tensiunea hertziana limita pentru materialul pinionului/rotii dintate, MPa
// sigma_Flim1,2 - Tensiunea de incovoiere limita pentru materialul pinionului/rotii dintate, MPa
// Raf1,Raf2	 - Rugozitatea flancului dintelui pinionului/rotii dintate, micrometri
// Rar1,Rar2	 - Rugozitatea razei de racordare a dintelui pinionului/rotii dintate, micrometri
// PozAngr		 - Pozitia angrenajului
// CPrec		 - Clasa de precizie a angrenajului
// mn			 - Modulul normal al angrenajului, mm
// aw			 - Distanta axiala standardizata/impusa, mm
// alfa_t		 - Unghiul de angrenare de referinta in plan frontal, rad
// alfa_wt		 - Unghiul real de angrenare in plan frontal, rad
// beta			 - Unghiul de inclinare al danturii pe cilindrul de divizare, grade
// beta_b		 - Unghiul de inclinare al danturii pe cilindrul de baza, rad
// eps_alfa		 - Gradul de acoperire al angrenajului in plan frontal
// eps_beta		 - Gradul de acoperire al angrenajului in plan axial
// eps_alfa_n	 - Gradul de acoperire al angrenajului echivalent
// b1,b2		 - Latimea pinionului/rotii dintate, mm
// zn1,zn2		 - Numarele de dinti ale angrenajului echivalent
// xn1,xn2		 - Coeficiemtii deplasarilor de profil in plan normal corespunzatori pinionului/rotii dintate
// d1,d2		 - Diametrele cercurilor de divizare ale pinionului/rotii dintate, mm
// n1,n2		 - Turatia arborelui 1 si 2, rot/min
{
	double PI;
	double Rzf1,Rzf2,Rz100,ZR1,ZR2,v1,v2,Zv1,Zv2,YFa1,YFa2,YSa1,YSa2;
	double Ydelta1,Ydelta2,XX1,XX2,Ydelta_st1,Ydelta_st2,Yx1,Yx2,NL1,NL2;
	double mH1,mH2,mF1,mF2,ZN1,ZN2,YN1,YN2,sigma_HP1,sigma_HP2,sigma_HP,sigma_H,sigma_FP1,sigma_FP2,sigma_F1,sigma_F2;
	double v1z1,ZE,Zeps,Yeps,ZH,Zbeta,Kvalfa,Kvbeta,Kv,psi_d,KHbeta,KFbeta,Ft1,qalfa;
	double KHalfa,KFalfa,Ybeta_min,Ybeta,fpbr;
	int validHDS;

	// Factorul rugozitatii flancurilor Zr pentru solicitarea de contact
	double GZR[64]=
	{1.689,1.918,2.106,2.315,2.607,2.815,3.149,3.462,3.796,4.172,4.526,4.881,5.299,5.695,6.112,6.446,6.864,7.323,7.782,8.262,8.742,9.242,9.743,10.223,10.766,11.266,11.788,12.268,12.685,13.042,13.520,14.000,
	1.100,1.079,1.060,1.040,1.022,1.007,0.993,0.977,0.962,0.948,0.937,0.927,0.915,0.905,0.896,0.887,0.878,0.868,0.859,0.850,0.843,0.835,0.829,0.822,0.817,0.813,0.809,0.805,0.803,0.801,0.800,0.799};
	double PZR[1]={850};

	// Factorul de viteza Zv
	double GZv[64]=
	{1.000,1.165,1.357,1.605,1.870,2.193,2.585,3.047,3.561,4.102,4.866,5.595,6.487,7.520,8.863,9.943,11.247,12.722,14.522,16.412,18.871,21.385,24.189,27.763,32.187,37.317,43.263,50.340,58.576,69.098,83.467,100.000,
	0.901,0.905,0.909,0.914,0.919,0.925,0.931,0.937,0.944,0.951,0.959,0.966,0.974,0.982,0.993,1.000,1.008,1.015,1.023,1.031,1.039,1.047,1.055,1.064,1.074,1.083,1.092,1.100,1.108,1.115,1.121,1.124};
	double PZv[1]={850};

	// Factorul de forma al dintelui YFa pentru solicitarea de incovoiere
	double GYFa[1024]=
	{15.423,15.955,16.416,16.931,17.542,18.134,18.876,19.935,21.002,22.001,23.051,24.525,26.301,27.663,28.990,30.819,32.505,34.620,36.701,39.306,41.723,46.297,50.832,56.764,65.714,75.439,89.558,112.077,147.021,208.569,274.782,400,
	 1.854,1.856,1.859,1.861,1.864,1.866,1.869,1.874,1.878,1.883,1.887,1.891,1.897,1.901,1.906,1.910,1.914,1.920,1.924,1.931,1.936,1.945,1.951,1.959,1.970,1.980,1.990,2.002,2.013,2.027,2.035,2.043,
	 13.442,14.060,14.526,15.062,15.681,16.233,16.821,17.708,18.685,19.824,20.937,22.089,23.169,24.658,26.402,28.204,30.337,32.721,35.202,38.295,41.920,45.487,49.703,56.211,64.457,75.392,86.796,102.781,123.693,158.903,202.822,400,
	 1.900,1.901,1.902,1.903,1.905,1.906,1.907,1.908,1.911,1.913,1.917,1.919,1.922,1.926,1.929,1.933,1.938,1.941,1.947,1.952,1.959,1.963,1.970,1.976,1.984,1.993,2.000,2.008,2.015,2.024,2.031,2.048,
	 11.652,12.065,12.468,12.926,13.532,14.178,14.915,15.754,16.718,17.480,18.184,19.010,19.979,21.127,22.282,23.619,25.003,27.369,30.231,32.473,35.953,40.293,45.010,50.484,59.705,68.911,79.915,94.883,113.951,147.553,227.989,400,
	 1.980,1.977,1.975,1.973,1.971,1.969,1.968,1.966,1.965,1.965,1.965,1.964,1.965,1.964,1.965,1.966,1.967,1.970,1.973,1.976,1.979,1.983,1.987,1.993,1.999,2.006,2.011,2.017,2.023,2.031,2.040,2.049,
	 10.033,10.364,10.754,11.147,11.677,12.263,12.873,13.473,14.126,15.071,16.054,16.994,18.114,19.109,20.669,21.934,23.594,25.457,27.535,30.063,32.110,35.318,37.877,41.288,45.286,50.742,59.556,73.435,92.445,125.273,198.877,400,
	 2.106,2.099,2.091,2.085,2.076,2.068,2.061,2.054,2.048,2.044,2.037,2.033,2.029,2.026,2.023,2.021,2.018,2.017,2.015,2.015,2.015,2.015,2.016,2.016,2.016,2.019,2.022,2.026,2.032,2.037,2.046,2.050,
	 10.042,10.486,10.872,11.319,11.636,12.039,12.456,13.040,13.738,14.504,15.494,16.400,17.523,18.707,19.944,21.136,22.488,24.396,26.107,28.748,31.238,34.134,37.019,41.924,49.846,57.942,70.564,89.904,115.967,151.435,212.942,400,
	 2.240,2.224,2.212,2.200,2.192,2.182,2.173,2.162,2.150,2.140,2.128,2.119,2.109,2.101,2.094,2.088,2.082,2.075,2.071,2.066,2.063,2.058,2.056,2.054,2.049,2.049,2.050,2.049,2.050,2.051,2.053,2.056,
	 10.036,10.288,10.690,11.103,11.578,12.109,12.540,13.161,13.578,14.128,14.755,15.493,16.335,16.910,17.584,18.617,20.194,21.688,23.244,25.267,27.521,30.043,33.874,38.641,46.351,57.010,73.127,89.580,113.932,147.074,208.324,400,
	 2.399,2.387,2.368,2.350,2.331,2.312,2.298,2.280,2.269,2.257,2.244,2.231,2.218,2.209,2.200,2.187,2.172,2.160,2.148,2.137,2.128,2.119,2.108,2.098,2.088,2.080,2.074,2.070,2.068,2.066,2.064,2.063,
	 10.349,10.727,11.132,11.554,12.020,12.639,13.205,13.790,14.498,15.300,16.149,17.219,18.330,19.438,20.941,22.493,24.179,26.110,28.476,30.386,32.415,35.337,39.201,44.690,50.362,58.596,71.669,89.108,111.527,138.737,208.061,400,
	 2.572,2.545,2.520,2.496,2.472,2.445,2.422,2.402,2.379,2.357,2.336,2.313,2.292,2.275,2.255,2.237,2.221,2.206,2.191,2.179,2.171,2.160,2.148,2.135,2.124,2.114,2.102,2.093,2.085,2.079,2.073,2.066,
	 12.041,12.477,12.984,13.473,14.041,14.646,15.344,16.308,17.063,18.052,19.029,20.112,21.342,22.986,24.739,26.499,28.833,31.493,33.856,37.234,40.678,45.309,51.179,60.370,70.608,79.736,92.222,110.265,143.858,202.528,265.057,400,
	 2.659,2.631,2.602,2.576,2.551,2.525,2.498,2.465,2.443,2.417,2.395,2.372,2.350,2.325,2.302,2.284,2.263,2.243,2.230,2.211,2.197,2.182,2.166,2.148,2.135,2.125,2.116,2.106,2.095,2.085,2.078,2.073,
	 13.768,14.105,14.503,14.961,15.509,15.917,16.425,16.887,17.583,18.259,18.748,19.334,20.128,21.432,22.678,23.899,25.484,27.099,28.825,30.457,33.110,35.116,38.476,42.829,48.870,56.132,68.319,80.656,105.411,148.204,218.751,400,
	 2.752,2.731,2.708,2.682,2.657,2.637,2.616,2.596,2.572,2.548,2.532,2.515,2.494,2.464,2.438,2.415,2.390,2.367,2.348,2.330,2.307,2.293,2.269,2.246,2.221,2.199,2.174,2.154,2.132,2.110,2.095,2.080,
	 15.465,15.848,16.311,16.743,17.104,17.653,18.261,18.953,19.918,20.787,21.790,22.950,24.298,25.590,27.284,29.077,31.395,34.061,37.493,40.480,44.209,49.224,55.665,64.086,73.217,84.182,99.224,111.339,139.624,190.518,259.934,400,
	 2.851,2.828,2.801,2.777,2.758,2.732,2.706,2.678,2.642,2.613,2.584,2.553,2.522,2.496,2.467,2.439,2.409,2.380,2.348,2.326,2.300,2.274,2.249,2.223,2.202,2.182,2.163,2.151,2.133,2.114,2.098,2.085,
	 17.221,17.660,18.159,18.727,19.296,20.027,20.799,21.698,22.651,23.597,24.936,26.210,27.551,29.247,31.159,33.341,36.166,38.700,42.221,45.466,49.937,54.589,59.855,67.699,76.628,90.699,109.059,133.679,165.646,221.848,282.297,400,
	 2.951,2.923,2.892,2.863,2.836,2.800,2.768,2.734,2.700,2.672,2.635,2.604,2.575,2.542,2.510,2.480,2.445,2.418,2.386,2.360,2.330,2.307,2.285,2.258,2.234,2.206,2.180,2.159,2.139,2.120,2.106,2.093,
	 18.950,19.570,20.192,20.929,21.742,22.810,23.859,24.830,25.962,27.109,28.141,29.576,30.900,32.718,34.349,36.674,38.948,41.547,45.445,49.399,54.149,59.666,66.136,74.109,83.914,96.341,111.680,133.580,165.066,200.730,258.259,400,
	 3.055,3.013,2.979,2.940,2.903,2.857,2.816,2.781,2.747,2.713,2.685,2.653,2.626,2.593,2.565,2.532,2.502,2.473,2.431,2.401,2.370,2.340,2.311,2.282,2.256,2.230,2.206,2.182,2.159,2.141,2.124,2.101,
	 20.604,21.158,21.776,22.678,23.449,24.481,25.741,26.729,28.051,29.279,30.365,31.962,33.650,35.531,37.922,41.149,44.809,48.298,51.547,55.593,59.852,65.086,72.044,78.815,88.559,99.499,111.901,131.231,166.180,195.383,284.662,400,
	 3.167,3.128,3.093,3.044,3.006,2.957,2.909,2.874,2.832,2.795,2.764,2.729,2.690,2.655,2.615,2.565,2.521,2.485,2.456,2.425,2.400,2.369,2.340,2.313,2.287,2.260,2.235,2.210,2.179,2.160,2.130,2.111,
	 22.351,23.018,23.654,24.403,25.142,26.035,26.926,27.982,29.209,30.312,31.748,32.983,34.385,36.295,38.006,39.849,42.027,44.073,47.437,50.863,56.192,62.863,69.143,77.013,90.394,107.812,127.317,154.703,192.180,231.985,286.242,400,
	 3.278,3.235,3.195,3.154,3.112,3.076,3.036,2.992,2.948,2.908,2.867,2.837,2.800,2.759,2.722,2.691,2.653,2.620,2.575,2.538,2.492,2.444,2.408,2.371,2.320,2.277,2.243,2.211,2.183,2.159,2.141,2.122,
	 24.080,25.056,26.001,27.003,28.081,29.266,30.692,31.560,32.536,33.627,34.891,36.451,37.863,39.336,41.431,44.451,48.098,51.757,55.631,59.244,64.471,70.188,76.387,83.818,92.281,104.174,119.311,137.936,168.401,204.888,250.783,400,
	 3.392,3.330,3.279,3.227,3.168,3.114,3.059,3.027,2.995,2.961,2.927,2.887,2.848,2.815,2.773,2.715,2.660,2.615,2.575,2.541,2.500,2.465,2.425,2.390,2.358,2.323,2.287,2.258,2.222,2.192,2.166,2.133,
	 25.864,26.755,27.770,28.729,29.763,30.836,31.879,32.698,33.749,34.905,36.250,37.798,39.652,41.009,42.643,44.284,46.237,49.396,51.375,53.999,57.521,62.289,67.729,74.591,83.763,95.208,110.371,134.681,161.356,211.304,250.164,400,
	 3.517,3.453,3.392,3.339,3.287,3.234,3.190,3.157,3.120,3.078,3.036,2.987,2.939,2.904,2.868,2.832,2.790,2.741,2.712,2.675,2.640,2.589,2.543,2.498,2.445,2.397,2.346,2.296,2.253,2.210,2.182,2.139};
	double PYFa[16]={1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5};

	// Factorul de corectie al tensiunilor de incovoiere la baza dintelui YSa
	double GYSa[1024]=
	{25.675,28.003,30.046,32.245,34.567,37.389,40.161,43.782,47.010,50.338,54.475,58.397,62.728,68.178,73.388,80.022,87.826,94.433,102.989,112.076,119.936,129.049,140.550,155.145,170.557,191.672,210.689,239.346,274.808,309.771,371.515,400,
	 1.382,1.400,1.415,1.428,1.442,1.457,1.472,1.489,1.504,1.518,1.534,1.547,1.561,1.577,1.592,1.609,1.627,1.641,1.657,1.673,1.685,1.698,1.712,1.729,1.745,1.761,1.776,1.794,1.811,1.829,1.847,1.853,
	 23.923,25.253,26.649,28.141,29.567,31.303,33.455,35.789,38.316,40.572,43.257,45.921,48.851,52.042,54.958,59.099,63.828,70.486,76.407,82.538,91.269,98.843,110.258,121.052,133.836,148.950,168.383,197.204,219.421,253.331,297.763,400,
	 1.409,1.420,1.430,1.441,1.451,1.462,1.474,1.487,1.501,1.512,1.524,1.536,1.548,1.560,1.571,1.584,1.599,1.618,1.633,1.648,1.666,1.681,1.699,1.715,1.732,1.750,1.767,1.787,1.803,1.820,1.843,1.864,
	 22.204,23.597,25.192,26.544,27.868,29.329,30.897,32.692,34.889,37.179,39.566,42.226,44.992,48.480,52.623,56.639,61.914,67.556,74.031,80.757,89.354,100.883,110.544,121.755,136.503,154.933,177.846,203.976,234.861,268.826,314.365,400,
	 1.436,1.448,1.461,1.471,1.481,1.491,1.500,1.511,1.522,1.534,1.546,1.558,1.571,1.585,1.600,1.614,1.629,1.646,1.661,1.677,1.694,1.715,1.730,1.744,1.761,1.780,1.797,1.816,1.831,1.845,1.862,1.883,
	 20.535,21.571,22.664,23.701,24.852,26.280,27.645,28.972,30.742,32.674,34.633,36.733,39.101,41.774,44.894,49.322,54.590,59.724,66.064,74.470,83.458,94.470,106.026,117.638,131.135,150.340,172.169,198.584,234.610,279.152,335.925,400,
	 1.465,1.475,1.485,1.493,1.502,1.513,1.522,1.531,1.541,1.552,1.562,1.573,1.585,1.597,1.610,1.627,1.645,1.660,1.677,1.696,1.715,1.735,1.752,1.767,1.782,1.800,1.817,1.834,1.850,1.865,1.883,1.897,
	 18.763,20.165,21.661,23.150,24.339,25.491,26.906,28.394,30.018,31.651,33.370,35.199,37.607,40.710,44.088,48.422,52.631,57.745,63.510,69.839,76.779,85.370,94.489,105.517,118.114,131.751,148.449,171.101,199.355,237.241,294.504,400,
	 1.495,1.509,1.522,1.535,1.544,1.552,1.562,1.571,1.581,1.590,1.599,1.608,1.620,1.633,1.647,1.663,1.677,1.691,1.706,1.721,1.735,1.751,1.766,1.781,1.795,1.809,1.823,1.838,1.853,1.869,1.886,1.908,
	 16.993,17.801,18.642,19.699,20.635,21.900,23.565,25.528,27.212,29.146,31.072,33.448,36.180,39.249,42.771,45.625,49.091,53.136,57.652,64.523,72.903,82.608,93.321,106.338,123.609,148.988,170.765,195.365,230.612,263.927,311.154,400,
	 1.524,1.533,1.542,1.553,1.561,1.572,1.585,1.599,1.610,1.621,1.632,1.644,1.656,1.669,1.682,1.693,1.704,1.716,1.729,1.744,1.761,1.778,1.794,1.810,1.827,1.847,1.860,1.871,1.884,1.895,1.906,1.918,
	 15.393,16.386,17.399,18.405,19.634,21.047,22.357,23.886,25.294,26.942,28.946,31.233,33.045,35.182,37.131,39.438,42.325,45.848,48.832,52.334,55.847,61.341,68.468,77.231,89.726,103.306,120.856,138.801,168.565,214.881,275.814,400,
	 1.556,1.568,1.579,1.590,1.602,1.614,1.624,1.635,1.645,1.655,1.666,1.677,1.686,1.695,1.703,1.712,1.722,1.734,1.743,1.753,1.762,1.773,1.787,1.802,1.819,1.834,1.850,1.863,1.879,1.896,1.912,1.929,
	 13.736,14.215,14.697,15.241,16.021,16.910,17.721,19.028,20.193,21.498,22.837,24.372,26.240,27.977,29.969,32.611,35.708,39.332,42.978,47.526,53.409,60.441,67.884,78.481,91.259,106.318,123.544,144.350,176.857,225.844,288.043,400,
	 1.586,1.593,1.600,1.607,1.616,1.626,1.635,1.647,1.656,1.667,1.677,1.687,1.698,1.707,1.717,1.729,1.741,1.754,1.766,1.779,1.793,1.807,1.820,1.834,1.850,1.864,1.876,1.887,1.901,1.915,1.927,1.939,
	 12.030,12.517,13.021,13.484,14.027,14.689,15.513,16.418,17.531,18.894,20.319,21.613,23.008,24.575,26.747,28.987,31.569,33.823,36.669,40.134,44.431,49.106,54.295,59.541,69.172,80.694,93.712,114.791,137.890,177.612,245.979,400,
	 1.615,1.623,1.630,1.637,1.645,1.653,1.663,1.673,1.684,1.697,1.708,1.717,1.726,1.736,1.747,1.758,1.769,1.777,1.786,1.797,1.809,1.820,1.830,1.839,1.853,1.867,1.879,1.893,1.905,1.918,1.932,1.947,
	 10.314,10.672,10.953,11.270,11.673,12.157,12.639,13.330,14.128,15.037,15.968,16.965,18.058,19.251,20.476,22.102,23.832,25.820,27.925,30.642,33.758,37.617,42.292,48.084,55.082,64.978,77.953,95.182,119.369,167.223,231.894,400,
	 1.637,1.645,1.650,1.656,1.663,1.671,1.678,1.688,1.699,1.710,1.720,1.730,1.740,1.749,1.758,1.768,1.778,1.788,1.797,1.808,1.818,1.829,1.841,1.853,1.865,1.878,1.891,1.903,1.916,1.931,1.942,1.956,
	 10.026,10.501,11.010,11.630,12.269,12.853,13.414,14.022,14.678,15.349,16.179,17.190,18.134,18.891,20.507,22.005,23.718,25.846,28.026,31.243,35.132,39.899,45.484,51.387,60.141,69.894,83.047,100.440,127.357,173.717,256.675,400,
	 1.684,1.694,1.704,1.714,1.724,1.732,1.740,1.747,1.755,1.763,1.771,1.780,1.788,1.794,1.804,1.813,1.821,1.831,1.839,1.850,1.860,1.871,1.881,1.890,1.900,1.909,1.918,1.927,1.936,1.946,1.954,1.962,
	 10.026,10.379,10.693,11.089,11.475,11.981,12.400,12.971,13.728,14.429,15.173,16.053,17.045,18.167,19.300,20.571,22.122,23.631,25.338,27.378,30.015,32.665,35.675,39.653,45.088,52.135,60.279,73.578,96.455,134.612,217.316,400,
	 1.732,1.739,1.746,1.753,1.759,1.767,1.773,1.781,1.790,1.798,1.806,1.814,1.823,1.831,1.839,1.846,1.854,1.861,1.867,1.874,1.882,1.888,1.895,1.902,1.910,1.918,1.925,1.933,1.943,1.951,1.960,1.967,
	 10.026,10.478,10.883,11.379,11.913,12.411,12.871,13.506,14.119,14.944,15.882,16.977,18.096,19.386,20.935,22.398,24.111,25.658,27.670,29.867,31.877,34.720,38.102,42.898,49.041,57.206,69.174,84.068,106.561,144.825,229.234,400,
	 1.773,1.782,1.790,1.798,1.806,1.813,1.819,1.827,1.833,1.842,1.850,1.859,1.867,1.875,1.883,1.889,1.895,1.900,1.906,1.911,1.915,1.920,1.925,1.931,1.937,1.943,1.949,1.953,1.958,1.963,1.967,1.970,
	 11.646,12.125,12.458,12.838,13.237,13.760,14.267,14.746,15.227,15.959,16.727,17.602,18.626,19.962,21.062,22.327,23.697,25.394,27.112,28.990,31.602,34.523,37.840,41.687,47.398,55.081,64.679,78.424,97.194,129.376,182.978,400,
	 1.832,1.838,1.843,1.847,1.852,1.858,1.863,1.868,1.872,1.878,1.884,1.890,1.896,1.903,1.907,1.912,1.917,1.922,1.926,1.930,1.935,1.939,1.943,1.947,1.951,1.955,1.959,1.962,1.965,1.969,1.971,1.972,
	 14.175,14.786,15.335,16.035,16.797,17.679,18.495,19.221,19.947,21.477,23.124,24.958,26.394,27.927,29.393,31.202,33.526,36.218,39.459,43.011,47.080,52.325,57.551,64.439,74.842,86.671,101.224,120.066,145.429,185.306,252.549,400,
	 1.874,1.879,1.883,1.888,1.893,1.899,1.904,1.907,1.910,1.917,1.922,1.927,1.931,1.934,1.937,1.940,1.943,1.946,1.950,1.953,1.955,1.958,1.960,1.962,1.965,1.967,1.968,1.969,1.970,1.971,1.972,1.972,
	 10.026,10.306,10.617,10.899,11.315,11.803,12.360,12.921,13.751,14.615,15.455,16.414,17.563,18.763,19.931,21.339,22.949,24.520,26.117,28.028,30.093,32.525,35.376,38.899,43.697,49.549,56.906,72.970,97.878,130.866,207.623,400,
	 1.773,1.785,1.798,1.808,1.821,1.836,1.851,1.861,1.872,1.880,1.887,1.895,1.902,1.909,1.915,1.920,1.926,1.930,1.934,1.938,1.942,1.945,1.949,1.952,1.956,1.960,1.962,1.966,1.969,1.971,1.973,1.973};
	double PYSa[16]={-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};

	// Factorul de durabilitate al materialului la concentratorul de tensiuni de la baza dintelui Yd
	double GYdelta[12]=
	{1.371,2.000,
	 0.953,1.030,
	 1.369,2.000,
	 0.961,1.026,
	 1.366,2.000,
	 0.977,1.016};
	double PYdelta[3]={500,600,800};

	// Factorul relativ de sensibilitate al materialului la concentratorul de tensiune de la baza dintelui la solicitarea statica Ydelta_st
	double GX[24]=
	{1.400,2.000,
	 1.483,2.000,
	 1.400,2.000,
	 1.527,2.109,
	 1.400,2.000,
	 1.572,2.247,
	 1.400,2.000,
	 1.604,2.382,
	 1.400,2.000,
	 1.637,2.528,
	 1.400,1.958,
	 1.678,2.598};
	double PX[6]={1,1.2,1.4,1.6,1.8,2};

	double GYdelta_st[12]=
	{1.483,2.600,
	 0.762,1.298,
	 1.483,2.600,
	 0.767,1.287,
	 1.483,2.600,
	 0.776,1.278};
	double PYdelta_st[3]={500,600,800};

	// Factorul de marime Yx
	double GYx[8]=
	{0.000,5.018,29.935,44.935,
	 1.000,1.000,0.852,0.852};
	double PYx[1]={1};

	// Factorul dinamic Kvalfa
	double GKvalfa[4]=
	{0.,6.018,
	 1.,1.800};
	double PKvalfa[1]={8};

	// Factorul dinamic Kvbeta
	double GKvbeta[4]=
	{0.000,6.643,
	 1.000,1.446};
	double PKvbeta[1]={8};

	// Factorul de repartizare a sarcinii pe latimea danturii KHbeta pentru solicitarea de contact
	double GKHbeta[64]=
	{0.000,0.055,0.119,0.184,0.243,0.286,0.352,0.417,0.481,0.557,0.621,0.728,0.785,0.843,0.907,0.972,1.055,1.138,1.217,1.299,1.366,1.454,1.527,1.593,1.658,1.717,1.773,1.829,1.877,1.926,1.967,2.000,
	 1.000,1.005,1.007,1.012,1.017,1.021,1.024,1.030,1.034,1.041,1.046,1.055,1.062,1.066,1.073,1.080,1.088,1.097,1.106,1.117,1.126,1.138,1.149,1.161,1.171,1.182,1.193,1.203,1.215,1.226,1.235,1.244};
	double PKHbeta[1]={5};

	// Factorul de repartizare a sarcinii pe latimea danturii KFbeta pentru solicitarea de incovoiere
	double GKFbeta[64]=
	{0.000,0.105,0.183,0.280,0.345,0.407,0.469,0.542,0.621,0.692,0.767,0.857,0.964,1.063,1.155,1.244,1.310,1.392,1.458,1.528,1.561,1.595,1.655,1.688,1.751,1.786,1.822,1.852,1.895,1.929,1.965,2.001,
	 1.000,1.015,1.025,1.034,1.042,1.050,1.060,1.070,1.081,1.094,1.108,1.132,1.156,1.179,1.204,1.231,1.255,1.283,1.306,1.340,1.354,1.370,1.396,1.413,1.443,1.459,1.480,1.499,1.522,1.540,1.559,1.579};
	double PKFbeta[1]={5};


	PI=3.1415926535897932384626433832795;

	validHDS=1;
	while(validHDS){
            cout << __func__ << '\n';
	Rzf1=exp(log(4.4)+0.97*log(Raf1)); // Rugozitatea Rz a flancului dintelui pinionului, micrometri
	Rzf2=exp(log(4.4)+0.97*log(Raf2)); // Rugozitatea Rz a flancului dintelui rotii, micrometri
	Rz100=(Rzf1+Rzf2)*sqrt(100./aw)/2.; // Rugozitatea ehivalenta, micrometri

	// Factorii rugozitatii flancurilor pentru solicitarea de contact
	ZR1=grafparam((double *)GZR, 2, 32, PZR, sigma_Hlim1, Rz100);
	if (ZR1==0.) {validHDS=0; break;}
	ZR2=grafparam((double *)GZR, 2, 32, PZR, sigma_Hlim2, Rz100);
	if(ZR2==0.) {validHDS=0; break;}

	v1=PI*d1*n1/60000.;
	if(v1<1) v1=1.; // Viteza periferica pe cercul de divizare al pinionului, m/s
	v2=PI*d2*n2/60000.;
	if(v2<1) v2=1.; // Viteza periferica pe cercul de divizare al pinionului, m/s

	// Factorii de viteza pentru solicitarea de contact
	Zv1=grafparam((double *)GZv, 2, 32, PZv, sigma_Hlim1, v1);
	if (Zv1==0.) {validHDS=0; break;}
	Zv2=grafparam((double *)GZv, 2, 32, PZv, sigma_Hlim2, v2);
	if (Zv2==0.) {validHDS=0; break;}

	xn1=0.2;zn1=27.5;
	// Factorii de forma pentru solicitarea de incovoiere
	YFa1=grafparam((double *)GYFa, 32, 32, PYFa, xn1, zn1);
	if ((YFa1<1.7)||(YFa1>3.7)) {validHDS=0; break;}
	YFa2=grafparam((double *)GYFa, 32, 32, PYFa, xn2, zn2);
	if ((YFa2<1.7)||(YFa2>3.7)) {validHDS=0; break;}

	// Factorii de corectie a tensiunilor de incovoiere la baza dintelui
	YSa1=grafparam((double *)GYSa, 32, 32, PYSa, xn1, zn1);
	if ((YSa1<1.2)||(YSa1>2)) {validHDS=0; break;}
	YSa2=grafparam((double *)GYSa, 32, 32, PYSa, xn2, zn2);
	if ((YSa2<1.2)||(YSa2>2)) {validHDS=0; break;}

	// Factorii relativi de sensibilitate a materialului la concentratorul de tensiuni de la baza dintelui, la durabilitate nelimitata
	Ydelta1=grafparam((double *)GYdelta, 6, 2, PYdelta, sigma021, YSa1);
	if (Ydelta1==0.) {validHDS=0; break;}
	Ydelta2=grafparam((double *)GYdelta, 6, 2, PYdelta, sigma022, YSa2);
	if (Ydelta2==0.) {validHDS=0; break;}

	// Factorii relativi de sensibilitate a materialului la concentratorul de tensiuni de la baza dintelui, la solicitarea statica
	XX1=grafparam((double *)GX, 12, 2, PX, eps_alfa_n, YSa1);
	if (XX1==0.) {validHDS=0; break;}
	XX2=grafparam((double *)GX, 12, 2, PX, eps_alfa_n, YSa2);
	if (XX2==0.) {validHDS=0; break;}

	Ydelta_st1=grafparam((double *)GYdelta_st, 6, 2, PYdelta_st, sigma021, XX1);
	if (Ydelta_st1==0.) {validHDS=0; break;}
	Ydelta_st2=grafparam((double *)GYdelta_st, 6, 2, PYdelta_st, sigma022, XX2);
	if (Ydelta_st2==0.) {validHDS=0; break;}

	// Factorii de marime pentru solicitarea de contact si de incovoiere
	Yx1=grafparam((double *)GYx, 2, 4, PYx, MatTT1, mn);
	if (Yx1==0.) {validHDS=0; break;}
	Yx2=grafparam((double *)GYx, 2, 4, PYx, MatTT2, mn);
	if (Yx2==0.) {validHDS=0; break;}

	// Factorii durabilitatii pentru solicitarea de contact
	NL1=60.*n1*Lh1*Hi1; // Numarul de cicluri de solicitare al pinionului
	mH1=log(NBH/NstH)/log(ZNmax/ZL1/ZR1/Zv1);
	mF1=log(NBF/NstF)/log(YNmax*Ydelta_st1/Ydelta1/YR1/Yx1);
	ZN1 =1.;
	if (NL1<=NstH) ZN1 =ZNmax;
	if ((NstH<NL1)&&(NL1<NBH)) ZN1=exp(log(NBH/NL1)/mH1);

	NL2=60.*n2*Lh2*Hi2; // Numarul de cicluri de solicitare al rotii dintate
	mH2=log(NBH/NstH)/log(ZNmax/ZL2/ZR2/Zv2); //
	mF2=log(NBF/NstF)/log(YNmax*Ydelta_st2/Ydelta2/YR2/Yx2); //
	ZN2 =1.;
	if (NL2<=NstH) ZN2 =ZNmax;
	if ((NstH<NL2)&&(NL2<NBH)) ZN2=exp(log(NBH/NL2)/mH2);

	// Factorii durabilitatii pentru solicitarea de incovoiere
	YN1 =1.;
	if (NL1<=NstF) YN1 =YNmax;
	if ((NstF<NL1)&&(NL1<NBF)) YN1=exp(log(NBF/NL1)/mF1);
	YN2 =1.;
	if (NL2<=NstF) YN2 =YNmax;
	if ((NstF<NL2)&&(NL2<NBF)) YN2=exp(log(NBF/NL2)/mF2);

	sigma_HP1= sigma_Hlim1*ZN1*Zw*ZL1*ZR1*Zv1*ZX/SHmin; // Tensiunea hertziana admisibila pt materialul pinionului, MPa
	sigma_HP2= sigma_Hlim2*ZN2*Zw*ZL2*ZR2*Zv2*ZX/SHmin; // Tensiunea hertziana admisibila pt materialul rotii, MPa

	sigma_HP=sigma_HP1; // Tensiunea hertziana admisibila, MPa
	if (sigma_HP1>sigma_HP2) sigma_HP=sigma_HP2;

	sigma_FP1=sigma_Flim1*YN1*Ydelta1*YR1*Yx1/SFmin; // Tensiunea de incovoiere admisibila pt materialul pinionului, MPa
	sigma_FP2=sigma_Flim2*YN2*Ydelta2*YR2*Yx2/SFmin; // Tensiunea de incovoiere admisibila pt materialul rotii, MPa

	// Factorul gradului de acoperire pentru solicitarea de contact
	Zeps=sqrt((4.-eps_beta)*(1.-eps_beta)/3.+eps_beta/eps_alfa);
	if (eps_beta>=1.) Zeps=sqrt(1/eps_alfa);

	// Factorul gradului de acoperire pentru solicitarea de incovoiere
	Yeps=0.25+0.75/eps_alfa_n;

	// Factorul zonei de contact pentru solicitarea de contact
	ZH=sqrt(2*cos(beta_b)/sin(alfa_wt)/cos(alfa_wt));

	// Factorul inclinarii dintilor pentru solicitarea de contact
	Zbeta=sqrt(cos(beta));

	// Factorii dinamici
	v1z1=v1*zn1*cos(beta)*cos(beta)*cos(beta)/100.; // In loc de z1 am scris zn1*cos(beta)*cos(beta)*cos(beta)
	Kvalfa=grafparam((double *)GKvalfa, 2, 2, PKvalfa, CPrec, v1z1);
	if (Kvalfa==0.) {validHDS=0; break;}
	Kvbeta=grafparam((double *)GKvbeta, 2, 2, PKvbeta, CPrec, v1z1);
	if (Kvbeta==0.) {validHDS=0; break;}

	Kv=Kvbeta;
	if (eps_beta<1.) Kv=Kvbeta-eps_beta*(Kvbeta-Kvalfa);

	// Factorii de repartizare a sarcini pe latimea danturii pentru solicitarea de contact, respectiv de incovoiere
	psi_d=b2/d1;
	KHbeta=grafparam((double *)GKHbeta, 2, 32, PKHbeta, PozAngr, psi_d);
	if (KHbeta==0.) {validHDS=0; break;}
	KFbeta=grafparam((double *)GKFbeta, 2, 32, PKFbeta, PozAngr, psi_d);
	if (KFbeta==0.) {validHDS=0; break;}

	// Abaterea efectiva a pasului de baza
	fpbr=0;
	if ((CPrec==7)&&((1<=mn)&&(mn<3.5))&&(d2<=125)) fpbr=13;
	if ((CPrec==7)&&((1<=mn)&&(mn<3.5))&&((d2>125)&&(d2<=400))) fpbr=15;
	if ((CPrec==7)&&((1<=mn)&&(mn<3.5))&&((d2>400)&&(d2<=800))) fpbr=17;

	if ((CPrec==7)&&((3.5<=mn)&&(mn<6.3))&&(d2<=125)) fpbr=17;
	if ((CPrec==7)&&((3.5<=mn)&&(mn<6.3))&&((d2>125)&&(d2<=400))) fpbr=19;
	if ((CPrec==7)&&((3.5<=mn)&&(mn<6.3))&&((d2>400)&&(d2<=800))) fpbr=19;

	if ((CPrec==7)&&((6.3<=mn)&&(mn<=10))&&(d2<=125)) fpbr=19;
	if ((CPrec==7)&&((6.3<=mn)&&(mn<=10))&&((d2>125)&&(d2<=400))) fpbr=21;
	if ((CPrec==7)&&((6.3<=mn)&&(mn<=10))&&((d2>400)&&(d2<=800))) fpbr=24;

	if ((CPrec==8)&&((1<=mn)&&(mn<3.5))&&(d2<=125)) fpbr=19;
	if ((CPrec==8)&&((1<=mn)&&(mn<3.5))&&((d2>125)&&(d2<=400))) fpbr=21;
	if ((CPrec==8)&&((1<=mn)&&(mn<3.5))&&((d2>400)&&(d2<=800))) fpbr=24;

	if ((CPrec==8)&&((3.5<=mn)&&(mn<6.3))&&(d2<=125)) fpbr=24;
	if ((CPrec==8)&&((3.5<=mn)&&(mn<6.3))&&((d2>125)&&(d2<=400))) fpbr=26;
	if ((CPrec==8)&&((3.5<=mn)&&(mn<6.3))&&((d2>400)&&(d2<=800))) fpbr=26;

	if ((CPrec==8)&&((6.3<=mn)&&(mn<=10))&&(d2<=125)) fpbr=26;
	if ((CPrec==8)&&((6.3<=mn)&&(mn<=10))&&((d2>125)&&(d2<=400))) fpbr=30;
	if ((CPrec==8)&&((6.3<=mn)&&(mn<=10))&&((d2>400)&&(d2<=800))) fpbr=34;

	if ((CPrec==9)&&((1<=mn)&&(mn<3.5))&&(d2<=125)) fpbr=26;
	if ((CPrec==9)&&((1<=mn)&&(mn<3.5))&&((d2>125)&&(d2<=400))) fpbr=30;
	if ((CPrec==9)&&((1<=mn)&&(mn<3.5))&&((d2>400)&&(d2<=800))) fpbr=34;

	if ((CPrec==9)&&((3.5<=mn)&&(mn<6.3))&&(d2<=125)) fpbr=34;
	if ((CPrec==9)&&((3.5<=mn)&&(mn<6.3))&&((d2>125)&&(d2<=400))) fpbr=38;
	if ((CPrec==9)&&((3.5<=mn)&&(mn<6.3))&&((d2>400)&&(d2<=800))) fpbr=38;

	if ((CPrec==9)&&((6.3<=mn)&&(mn<=10))&&(d2<=125)) fpbr=32;
	if ((CPrec==9)&&((6.3<=mn)&&(mn<=10))&&((d2>125)&&(d2<=400))) fpbr=42;
	if ((CPrec==9)&&((6.3<=mn)&&(mn<=10))&&((d2>400)&&(d2<=800))) fpbr=48;

	Ft1=2*T1/d1; // Forta tangentiala corespunzatoare diametrului de divizare, N
	qalfa=1.; // Factorul auxiliar
	if (4*(0.1+b2*(fpbr-4.)/Ft1)<=0.5) qalfa=0.5;

	// Factorul de repartizare a sarcinii in plan frontal pe perechile de dintii aflate simultan in angrenare, pentru solicitarea de contact
	KHalfa=1.+2.*(qalfa-0.5)*(1./Zeps/Zeps-1.);

	// Factorul de repartizare a sarcinii in plan frontal pe perechile de dintii aflate simultan in angrenare, pentru solicitarea de incovoiere
	KFalfa=qalfa*eps_alfa;

	// Factorul minim al inclinarii dintilor pentu solicitarea de incovoiere
	Ybeta_min=0.75;
	if (eps_beta<=1.) Ybeta_min=1.-0.25*eps_beta;

	// Factorul inclinarii dintilor pentu solicitarea de incovoiere
	Ybeta=1.-eps_beta*(beta*180/PI)/120.;
	if(Ybeta<Ybeta_min)Ybeta=Ybeta_min;

	// Factorul de elasticitate al materialului
	ZE=sqrt( 1/( (1-niu1*niu1)/E1+(1-niu2*niu2)/E2 )/PI );

	sigma_H=(u12+1)*ZE*Zeps*ZH*Zbeta*cos(alfa_t)*sqrt(T1*KA*Kv*KHbeta*KHalfa*(u12+1)/b2/u12/2.)/aw/cos(alfa_wt); // Tensiunea hertziana

	sigma_F1=T1*z1*(u12+1)*(u12+1)*KA*Kv*KFbeta*KFalfa*Yeps*Ybeta*YFa1*YSa1*cos(alfa_t)*cos(alfa_t)/cos(alfa_wt)/cos(alfa_wt)/cos(beta)/aw/aw/b1/2.; // Tensiunea de incovoiere la baza dintelui pinionului, MPa
	sigma_F2=sigma_F1*b1*YFa2*YSa2/b2/YFa1/YSa1; // Tensiunea de incovoiere la baza dintelui rotii, MPa

	break;
	}

	if(validHDS)
	{

		HDS[0]=1;
		HDS[1]=sigma_H;
		HDS[2]=sigma_HP;
		HDS[3]=sigma_F1;
		HDS[4]=sigma_F2;
		HDS[5]=sigma_FP1;
		HDS[6]=sigma_FP2;

	}
	else HDS[0]=0;
}

void LoadsMomentsHelical(double T1,double dw1,double dw2,double beta_w,double alfa_wt,double *LMH)
// Returneaza valoarea fortelor tangentiale,axiale si radiale precum si valoarea momentului concentrat
// dw1		- Diametrul de rostogolire al pinionului, mm
// dw2		- Diametrul de rostogolire al rotii, mm
// beta_w	- Unghiul de inclinare al danturii pe cilindrul de rostogolire, rad
// alfa_wt	- Unghiul real de angrenare in plan frontal, rad
{
	double Ft1,Fa1,Fr1,Mic1;
	double Ft2,Fa2,Fr2,Mic2;

	Ft1=2*T1/dw1;
	Fa1=Ft1*tan(beta_w);
	Fr1=Ft1*tan(alfa_wt);

	Mic1=Fa1*dw1/2;

	Ft2=Ft1;
	Fa2=Fa1;
	Fr2=Fr1;

	Mic2=Fa2*dw2/2;

	LMH[0]=1; // Se va folosi pentru control
	LMH[1]=Ft1;
	LMH[2]=Fa1;
	LMH[3]=Fr1;
	LMH[4]=Mic1;
	LMH[5]=Ft2;
	LMH[6]=Fa2;
	LMH[7]=Fr2;
	LMH[8]=Mic2;
}

void ExternalHelicalDriveUnit (double mn,double aw,int z1,int z2,double xn1,double beta,double b2,double Delta_b,double T1,double u12,double KA,double Lh1,double Lh2, double ZL1, double ZL2,double ZX,double Zw,double ZNmax,double YNmax, double NstH,double NBH,double NstF,double NBF,double SHmin,double SFmin,double YR1,double YR2,double Hi1,double Hi2,double niu1,double niu2,double E1,double E2,double MatTT1,double MatTT2,double HB1,double HB2,double sigma021,double sigma022,double sigma_Hlim1,double sigma_Hlim2,double sigma_Flim1,double sigma_Flim2,double Raf1,double Raf2,double Rar1,double Rar2,double PozAngr,double CPrec,double n1,double n2,double *EHDU)
// mn			 - Modulul normal al angrenajului, mm
// aw			 - Distanta axiala standardizata/impusa, mm
// z1,z2		 - Numerele de dinti ale pinionului si ale rotii dintate
// xn1,xn2		 - Coeficiemtii deplasarilor de profil in plan normal corespunzatori pinionului
// beta			 - Unghiul de inclinare al danturii pe cilindrul de divizare, grade
// b2			 - Latimea rotii dintate, mm
// Delta_b		 - Diferenta de latime dintre pinion si roata dintata, mm
// T1			 - Momentul de torsiune corespunzator arborelui de intrare, Nmm
// u12			 - Raportul deal de angrenare corespunzator treptei I
// KA			 - Factorul regimului de functionare
// Lh1,Lh2		 - Durata minima de functionare a pinionului/rotii dintate, ore
// ZL1,ZL2		 - Factorul de ungere corespunzator pinionului/rotii dintate
// Zx			 - Factorul de marime
// ZW			 - Factorul raportului duritatilor flancurilor dintilor
// ZNmax		 - Factorul durabilitatii pentru solicitarea de contact
// YNmax		 - Factorul durabilitatii pentru solicitarea de incovoiere
// NstH,NstF	 -
// NBH,NBF		 -
// SHmin		 - Coeficientul de siguranta minim pentru solicitarea de contact
// SFmin		 - Coeficientul de siguranta minim pentru solicitarea de incovoiere
// YR1,YR2		 - Factorul rugozitatii flancului dintelui pinionului/rotii dintate pentru solicitarea de incovoiere
// Hi1,Hi2		 - Numarul de roti cu care vine in contact pinionul/roata dintata
// niu1,niu2	 - Coeficientii lui Poisson
// E1,E2		 - Modulul de elasticitate longitudinal, MPa
// MatTT1,MatTT2 -
// HB1,HB2		 - Duritatea materialului pinionului/rotii dintate, MPa
// sigma021,2	 - Limita de curegere a materialului pinionului/rotii dintate, MPa
// sigma_Hlim1,2 - Tensiunea hertziana limita pentru materialul pinionului/rotii dintate, MPa
// sigma_Flim1,2 - Tensiunea de incovoiere limita pentru materialul pinionului/rotii dintate, MPa
// Raf1,Raf2	 - Rugozitatea flancului dintelui pinionului/rotii dintate, micrometri
// Rar1,Rar2	 - Rugozitatea razei de racordare a dintelui pinionului/rotii dintate, micrometri
// PozAngr		 - Pozitia angrenajului
// CPrec		 - Clasa de precizie a angrenajului
// n1,n2		 - Turatia arborelui 1 si 2, rot/min
{
	double HDG[42],EHD[14],CHD[17],HDS[7],LMH[9];
	double cs_n;
	double alfa_t,alfa_n,beta_b,beta_w,a,alfa_wt,xsn,eps_alfa,eps_beta,eps_gama;
	double xt1,alfa_at1,d1,db1,df1,da1,dw1,st1,sat1,beta_a1,sn1,san1,ADinte1,b1,VolDantura1;
	double xn2,xt2,alfa_at2,d2,db2,df2,da2,dw2,st2,sat2,beta_a2,sn2,san2,ADinte2,VolDantura2;
	double zn1,zn2,xn1min,xn2min,dn1,dn2,dbn1,dbn2,dan1,dan2,alfa_wn,awn,eps_alfa_n;
	double alfa_Nt1,Ncalc1,N1,WNn1,WNt1,roat1,roNt1,roAt1,alfa_Nt2,Ncalc2,N2,WNn2,WNt2,roat2,roNt2,roEt2;
	double sigma_HP,sigma_H,sigma_FP1,sigma_FP2,sigma_F1,sigma_F2;
	double Ft1,Fa1,Fr1,Mic1;
	double Ft2,Fa2,Fr2,Mic2;

	int validEHDU;
	validEHDU=1;
	while (validEHDU==1){
    cout << __func__ << '\n';
	// Elementele geometrice ale angrenajului
	HelicalDriveGeometry(mn,aw,z1,z2,xn1,beta,b2,Delta_b,HDG);
	if(HDG[0]==0) {validEHDU=0; break;}
	alfa_t=HDG[1];
	alfa_at1=HDG[2];
	beta_a1=HDG[3];
	beta_b=HDG[4];
	d1=HDG[5];
	db1=HDG[6];
	df1=HDG[7];
	da1=HDG[8];
	st1=HDG[9];
	sat1=HDG[10];
	sn1=HDG[11];
	san1=HDG[12];
	xt1=HDG[13];
	ADinte1=HDG[14];
	a=HDG[15];
	alfa_wt=HDG[16];
	xsn=HDG[17];
	xn2=HDG[18];
	alfa_at2=HDG[19];
	beta_a2=HDG[20];
	d2=HDG[21];
	db2=HDG[22];
	df2=HDG[23];
	da2=HDG[24];
	st2=HDG[25];
	sat2=HDG[26];
	sn2=HDG[27];
	san2=HDG[28];
	xt2=HDG[29];
	ADinte2=HDG[30];
	dw1=HDG[31];
	dw2=HDG[32];
	beta_w=HDG[33];
	eps_alfa=HDG[34];
	eps_beta=HDG[35];
	eps_gama=HDG[36];
	b1=HDG[37];
	cs_n=HDG[38];
	alfa_n=HDG[39];
	VolDantura1=HDG[40];
	VolDantura2=HDG[41];

	// Angrenajul echivalent
	EquivalentHelicalDrive(mn,a,z1,z2,beta,beta_b,beta_w,alfa_n,alfa_wt,d1,da1,d2,da2,EHD);
	// Elementele angrenajului echivalent
	zn1=EHD[1];
	zn2=EHD[2];
	xn1min=EHD[3];
	xn2min=EHD[4];
	dn1=EHD[5];
	dn2=EHD[6];
	dbn1=EHD[7];
	dbn2=EHD[8];
	dan1=EHD[9];
	dan2=EHD[10];
	alfa_wn=EHD[11];
	awn=EHD[12];
	eps_alfa_n=EHD[13];

	// Elementele de control ale angrenajului
	ControlHelicalDrive(mn,aw,z1,z2,xn1,xn2,beta,beta_b,beta_w,alfa_n,alfa_t,alfa_wt,alfa_at1,alfa_at2,d1,da1,d2,db1,db2,da2,CHD);
	alfa_Nt1=CHD[1];
	alfa_Nt2=CHD[2];
	Ncalc1=CHD[3];
	Ncalc2=CHD[4];
	N1=CHD[5];
	N2=CHD[6];
	WNn1=CHD[7];
	WNn2=CHD[8];
	WNt1=CHD[9];
	WNt2=CHD[10];
	roNt1=CHD[11];
	roNt2=CHD[12];
	roAt1=CHD[13];
	roEt2=CHD[14];
	roat1=CHD[15];
	roat2=CHD[16];

	// Calculul fortelor si momentelor de incovoiere concentrate (datorate fortelor axiale)
	LoadsMomentsHelical (T1,dw1,dw2,beta_w,alfa_wt,LMH);
	Ft1=LMH[1];
	Fa1=LMH[2];
	Fr1=LMH[3];
	Mic1=LMH[4];

	Ft2=LMH[5];
	Fa2=LMH[6];
	Fr2=LMH[7];
	Mic2=LMH[8];

	// Calculul de rezistenta
	HelicalDriveStrength(T1,z1,u12,KA,Lh1,Lh2,ZL1,ZL2,ZX,Zw,ZNmax,YNmax,NstH,NBH,NstF,NBF,SHmin,SFmin,YR1,YR2,Hi1,Hi2,niu1,niu2,E1,E2,MatTT1,MatTT2,HB1,HB2,sigma021,sigma022,sigma_Hlim1,sigma_Hlim2,sigma_Flim1,sigma_Flim2,Raf1,Raf2,Rar1,Rar2,PozAngr,CPrec,mn,aw,alfa_t,alfa_wt,beta,beta_b,eps_alfa,eps_beta,eps_alfa_n,b1,b2,zn1,zn2,xn1,xn2,d1,d2,n1,n2,HDS);
	if(HDS[0]==0) {validEHDU=0; break;}
	sigma_H=HDS[1];
	sigma_HP=HDS[2];
	sigma_F1=HDS[3];
	sigma_F2=HDS[4];
	sigma_FP1=HDS[5];
	sigma_FP2=HDS[6];
	break;
	}

	if (validEHDU)
	{

		EHDU[0]=1;
		// Elementele geometrice ale angrenajului
		EHDU[1]=alfa_t;
		EHDU[2]=alfa_at1;
		EHDU[3]=beta_a1;
		EHDU[4]=beta_b;
		EHDU[5]=d1;
		EHDU[6]=db1;
		EHDU[7]=df1;
		EHDU[8]=da1;
		EHDU[9]=st1;
		EHDU[10]=sat1;
		EHDU[11]=sn1;
		EHDU[12]=san1;
		EHDU[13]=xt1;
		EHDU[14]=ADinte1;
		EHDU[15]=a;
		EHDU[16]=alfa_wt;
		EHDU[17]=xsn;
		EHDU[18]=xn2;
		EHDU[19]=alfa_at2;
		EHDU[20]=beta_a2;
		EHDU[21]=d2;
		EHDU[22]=db2;
		EHDU[23]=df2;
		EHDU[24]=da2;
		EHDU[25]=st2;
		EHDU[26]=sat2;
		EHDU[27]=sn2;
		EHDU[28]=san2;
		EHDU[29]=xt2;
		EHDU[30]=ADinte2;
		EHDU[31]=dw1;
		EHDU[32]=dw2;
		EHDU[33]=beta_w;
		EHDU[34]=eps_alfa;
		EHDU[35]=eps_beta;
		EHDU[36]=eps_gama;
		EHDU[37]=b1;
		EHDU[38]=cs_n;
		EHDU[39]=alfa_n;

		// Elementele geometrice ale angrenajului echivalent
		EHDU[40]=zn1;
		EHDU[41]=zn2;
		EHDU[42]=xn1min;
		EHDU[43]=xn2min;
		EHDU[44]=dn1;
		EHDU[45]=dn2;
		EHDU[46]=dbn1;
		EHDU[47]=dbn2;
		EHDU[48]=dan1;
		EHDU[49]=dan2;
		EHDU[50]=alfa_wn;
		EHDU[51]=awn;
		EHDU[52]=eps_alfa_n;

		// Elementele de control ale angrenajului
		EHDU[53]=alfa_Nt1;
		EHDU[54]=alfa_Nt2;
		EHDU[55]=Ncalc1;
		EHDU[56]=Ncalc2;
		EHDU[57]=N1;
		EHDU[58]=N2;
		EHDU[59]=WNn1;
		EHDU[60]=WNn2;
		EHDU[61]=WNt1;
		EHDU[62]=WNt2;
		EHDU[63]=roNt1;
		EHDU[64]=roNt2;
		EHDU[65]=roAt1;
		EHDU[66]=roEt2;
		EHDU[67]=roat1;
		EHDU[68]=roat2;

		// Tensiuni si tensiuni admisibile
		EHDU[69]=sigma_H;
		EHDU[70]=sigma_HP;
		EHDU[71]=sigma_F1;
		EHDU[72]=sigma_F2;
		EHDU[73]=sigma_FP1;
		EHDU[74]=sigma_FP2;

		EHDU[75]=VolDantura1;
		EHDU[76]=VolDantura2;

		EHDU[77]=Ft1;
		EHDU[78]=Fa1;
		EHDU[79]=Fr1;
		EHDU[80]=Mic1;
		EHDU[81]=Ft2;
		EHDU[82]=Fa2;
		EHDU[83]=Fr2;
		EHDU[84]=Mic2;

	}
	else EHDU[0]=0;
}

double afiseaza()
{

		int i,z2;

		double T1,T2;
		double cs_n;
		double u12,n2;
		double alfa_n,alfa_t,alfa_wt,alfa_wn,beta,beta_b,beta_w;
		double alfa_at1,alfa_Nt1,beta_a1;
		double alfa_at2,alfa_Nt2,beta_a2;
		double mn,a,aw,awn,b1,b2;
		double xsn,xn2,xt1,xt2,xn1min,xn2min;
		double zn1,zn2;
		double d1,db1,dw1,df1,da1;
		double d2,db2,dw2,df2,da2;
		double dn1,dbn1,dan1;
		double dn2,dbn2,dan2;
		double sn1,st1,sat1,san1;
		double sn2,st2,sat2,san2;
		double eps_alfa,eps_beta,eps_gama,eps_alfa_n;
		double sigma_HP,sigma_H;
		double sigma_FP1,sigma_F1;
		double sigma_FP2,sigma_F2;
		double Ncalc1,N1;
		double Ncalc2,N2;
		double WNn1,WNt1,roNt1,roAt1,roat1;
		double WNn2,WNt2,roNt2,roEt2,roat2;
		double Ft1,Ft2,Fa1,Fa2,Fr1,Fr2,Mic1,Mic2,MH3,MH7;

		double ADinte1,VolDantura1;
		double ADinte2,VolDantura2;
		double Masaz1,Masaz2,MasaAngrenaj;


		double EHDU1[85];
		double ValIesire;


		double PI			=		  3.1415926535897932384626433832795;
		double CoefPen		=	      1e20;		// Coeficientul de penalizare

		double tolStrOfMat	=		  0.002; // toleranta cu care se doreste verificarea: sigma <= (1 + tolStrOfMat)*sigma_all





		//=====================================================================================================================
		// DATE DE INTRARE
		//=====================================================================================================================





		double P			=         2.9;		// Puterea de intrare, kW
		double i12STAS		=	      3.55;		// Raportul de transmitere total
		double n1			=		  925;		// Turatia arborelui 1, rot/min

		// Tip lubrifiant TIN 125 EP STAS 10588-76 avand vascozitatea cinematica 125-140 [mm2/s] la 500 C
		double ZL1			=		  1.05;		// Factor de ungere
		double ZL2			=	      1.05;		// Factor de ungere
		double Lh1			=	   9000;		// Durata minima de functionare a pinionului, ore
		double Lh2			=      9000;		// Durata minima de functionare a rotii, ore
		double Hi1			=		  1;		// Numarul de roti cu care vine in contact pinionul
		double Hi2			=		  1;		// Numarul de roti cu care vine in contact roata

		double Delta_b      =		  4;		// Diferenta de latime a rotilor, mm

		double PozAngr      =		  5;		// Pozitia rotilor angrenajului in raport cu lagarele arborilor transmisiei
		double CPrec        =		  8; 		// Clasa de precizie

		double eps_alfa_lim =		  1.2;		// Valoarea minima admisibila a gradului de acoperie frontal

		// Caracteristicile materialelor pentru trepta I
		double ro_otel		=		  7.85e-6;	// Densitatea otelului, kg/mm^3
		double MatTT1		=		  1;		// Pinion: 41MoCr11 imbunatatit
		double HB1			=	   3000;		// Duritatea flancului dintelui pinionului, MPa
		double sigma021		=	    690;		// Limita de curgere a materialului pinionului, MPa
		double Raf1			=		  0.8;		// Rugozitatea flancului dintelui pinionului, micrometri
		double Rar1			=		  1.6;		// Rugozitatea zonei de racordare a dintelui pinionului, micrometri
		double sigma_Hlim1	=		760;		// Tensiunea hertziana limita a materialului pinionului, MPa
		double sigma_Flim1	=		580;		// Tensiunea de incovoiere limita a materialului pinionului, MPa
		double YR1			=		  1.02;		// Factorul rugozitatii flancului dintelui pinionului pentru solicitarea de incovoiere
		double niu1			=		  0.3;		// Coeficientul lui Poisson pentru materialul pinionului
		double E1			=		  2.05e5;    // Modulul de elasticitate longitudinal pentru materialul pinionului

		double MatTT2		=   	  1;    	// Roata: 40Cr10 imbunatatit
		double HB2			=	   2700;        // Duritatea flancului dintelui rotii, MPa
		double sigma022		=	    560;		// Limita de curgere a materialului rotii, MPa
		double Raf2			=		  0.8;		// Rugozitatea flancului dintelui rotii, micrometri
		double Rar2			=		  1.6;		// Rugozitatea zonei de racordare a dintelui rotii, micrometri
		double sigma_Hlim2	=		720;		// Tensiunea hertziana limita a materialului rotii, MPa
		double sigma_Flim2	=		560;		// Tensiunea de incovoiere limita a materialului rotii, MPa
		double YR2			=		  1.02;		// Factorul rugozitatii flancului dintelui rotii pentru solicitarea de incovoiere
		double niu2			=		  0.3;		// Coeficientul lui Poisson pentru materialul rotii
		double E2			=		  2.05e5;    // Modulul de elasticitate longitudinal pentru materialul rotii


		double ZX			=		  1.;		// Factorul de forma
		double Zw			=	      1.;		// Factorul raportului duritatilor flancurilor dintilor

		double NstH			=	      1.e5;		//
		double NBH			=		  5.e7;		//
		double ZNmax		=		  1.6;		//

		double NstF			=	      1.e4;		//
		double NBF			=		  3.e6;		//
		double YNmax		=		  2.5;		//

		double SHmin		=		  1.15;		// Coeficientul de siguranta minim pentru solicitarea de contact
		double SFmin		=		  1.25;		// Coeficientul de siguranta minim pentru solicitarea de incovoiere

		double eta_rul		=	      0.99;		// Randamentul unei perechi de rulmenti
		double eta_angr		=		  0.97;		// Randamentul amgrenajului cilindric cu dinti inclinati
		double eta_ungere	=			0.995;	// randamentul datorita antrenarii uleiului (pt o trepata)

		double KA			=		  1.25;		// Factorul regimului de functionare

		double romat		=	      7.85e-6;  // Densitatea materialului capacelor, kg/mm^3

		int PosibilAngrenaj;

		int nr_restrictii;
		double pen, obj;
		double g[100];


		int valid_angrenaj=1;

		double Lista_mn[41]=
		{1.0,1.125,1.25,1.375,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,8.0,9.0,10.0,
		11.0,12.0,14.0,16.0,18.0,20.0,22.0,25.0,28.0,32.0,36.0,40.0,45.0,50.0,55.0,60.0,70.0,80.0,90.0,100.0};

		double Lista_aw[16]=
		{56.0,63.0,71.0,80.0,90.0,100.0,112.0,125.0,140.0,160.0,180.0,200.0,225.0,250.0,280.0,315.0};

		double Lista_iSTAS[32]=
		{1.12,1.25,1.40,1.60,1.80,2.00,2.24,2.50,2.80,3.15,3.55,4.00,4.50,5.00,5.60,6.30,7.10,8.00,9.00,10.00,
		11.2,12.50,14.00,16.00,18.00,20.00,22.40,25.00,28.00,31.50,35.50,40.00};






		//=====================================================================================================================
		// GENELE
		//=====================================================================================================================


/*
		//----------------------------------------------------------------------------------------------------------------------
		int		z1=(int) vect[0];		// Numarul de dinti ai pinionului (intreg, 25 ... 56) -> 32 valori
		int		iaw=(int) vect[1];		// Indicele distantei dintre axe, mm (intreg, 0...15) -> 16 valori
		double	xn1=round(vect[2],2);	// Coeficientul deplasarii de profil a profilului dintelui pinionului, in plan normal (real, -0.5 ... +1.0) -> 128 valori
		double	psi_a=round(vect[3],2);	// Coeficientul de latime a pinionului (b/aw) (real, 0.2 ... +0.8) -> 128 valori
		double	beta_g=vect[4];			// Unghiul de inclinare al danturii pe cilindrul de divizare, grade (real, 7.25 ... +15.00) -> 32 valori
		//----------------------------------------------------------------------------------------------------------------------

//31   6  0.9882  0.2142 15.00  4.6675e+000
//31   6  0.9764  0.2142 15.00  4.6684e+000
*/
		//Varianata optima
		//----------------------------------------------------------------------------------------------------------------------
		int		z1=31;		// Numarul de dinti ai pinionului (intreg, 25 ... 56) -> 32 valori
		int		iaw=6;		// Indicele distantei dintre axe, mm (intreg, 0...15) -> 16 valori
		double	xn1=0.9764;		// Coeficientul deplasarii de profil a profilului dintelui pinionului, in plan normal (real, -0.5 ... +1.0) -> 128 valori
		double	psi_a=0.2142;		// Coeficientul de latime a pinionului (b/aw) (real, 0.2 ... +0.8) -> 128 valori
		double	beta_g=15;	// Unghiul de inclinare al danturii pe cilindrul de divizare, grade (real, 7.25 ... +15.00) -> 32 valori
		//----------------------------------------------------------------------------------------------------------------------
/*
		//Varianta clasica
		//----------------------------------------------------------------------------------------------------------------------
		int		z1=25;		// Numarul de dinti ai pinionului (intreg, 25 ... 56) -> 32 valori
		int		iaw=5;		// Indicele distantei dintre axe, mm (intreg, 0...15) -> 16 valori
		double	xn1=0.079;		// Coeficientul deplasarii de profil a profilului dintelui pinionului, in plan normal (real, -0.5 ... +1.0) -> 128 valori
		double	psi_a=0.25;		// Coeficientul de latime a pinionului (b/aw) (real, 0.2 ... +0.8) -> 128 valori
		double	beta_g=8;	// Unghiul de inclinare al danturii pe cilindrul de divizare, grade (real, 7.25 ... +15.00) -> 32 valori
		//----------------------------------------------------------------------------------------------------------------------

*/


		ValIesire=1000.*CoefPen*(1.+RandomNumber());

		while (valid_angrenaj==1)
        {



            aw=Lista_aw[iaw];
            beta=beta_g*PI/180.; // Unghiul de inclinare al danturii pe cilindrul de divizare, treapta I, rad

            z2=(int)round(z1*i12STAS, 0); // Numarul de dinti ai rotii z2
            u12=z2*1./z1; // Rapotrul de real de angrenare al primei trepte. Atentie la tipul variabilelor. am pus z2*1./z1 ca sa-l fac real

            mn=GetList(2.*aw*cos(beta)/z1/(1+u12), 67, Lista_mn);

            n2=n1/u12; // Turatia rotii z2


            T1=3.e7*P*eta_rul/n1/PI;				// Momentul de torsiune de pe arborele 1, N*mm
            T2=T1*eta_rul*eta_angr*eta_ungere*u12; // Momentul de torsiune de pe arborele 2, N*mm

            b2=floor(psi_a*aw); // Latimea rotii, mm
            cout << mn << ' ' << aw  << ' ' << z1 << ' ' << z2 << ' ' << xn1 << ' ' << beta << ' ' << b2 << ' ' << Delta_b << ' ' << T1 << ' ' <<u12
            << ' ' << KA << ' ' <<Lh1 << ' ' <<Lh2 << ' ' <<ZL1 << ' ' <<ZL2 << ' ' <<ZX << ' ' <<Zw << ' ' <<ZNmax << ' ' <<YNmax << ' ' <<NstH << ' '
            <<NBH << ' ' <<NstF << ' ' <<NBF << ' ' <<SHmin << ' ' << SFmin << ' ' <<
            YR1 << ' ' <<YR2 << ' ' <<Hi1 << ' ' <<Hi2 << ' ' <<niu1 << ' ' <<niu2 << ' ' <<E1 << ' ' <<E2 << ' ' <<MatTT1 << ' ' <<MatTT2 << ' ' <<HB1
            << ' ' <<HB2 << ' ' <<sigma021 << ' ' <<sigma022 << ' ' <<sigma_Hlim1 << ' ' <<sigma_Hlim2 << ' ' <<sigma_Flim1 << ' ' <<
            sigma_Flim2 << ' ' <<Raf1 << ' ' <<Raf2 << ' ' <<Rar1 << ' ' <<Rar2 << ' ' <<PozAngr << ' ' <<CPrec << ' ' <<n1 << ' ' <<n2;

            ExternalHelicalDriveUnit(mn,aw,z1,z2,xn1,beta,b2,Delta_b,T1,u12,KA,Lh1,Lh2,ZL1,ZL2,ZX,Zw,ZNmax,YNmax,NstH,NBH,NstF,NBF,SHmin, SFmin,
                                     YR1,YR2,Hi1,Hi2,niu1,niu2,E1,E2,MatTT1,MatTT2,HB1,HB2,sigma021,sigma022,sigma_Hlim1,sigma_Hlim2,sigma_Flim1,
                                     sigma_Flim2,Raf1,Raf2,Rar1,Rar2,PozAngr,CPrec,n1,n2,EHDU1);
            cout << "Ne oprim aici\n";

            cout << "Rezultatul " << EHDU1[0];
            PosibilAngrenaj=(int)EHDU1[0];
            if (PosibilAngrenaj==0) {valid_angrenaj=0; break;}

            // Elementele geometrice ale angrenajului de pe treapta I
            alfa_t=EHDU1[1];
            alfa_at1=EHDU1[2];
            beta_a1=EHDU1[3];
            beta_b=EHDU1[4];
            d1=EHDU1[5];
            db1=EHDU1[6];
            df1=EHDU1[7];
            da1=EHDU1[8];
            st1=EHDU1[9];
            sat1=EHDU1[10];
            sn1=EHDU1[11];
            san1=EHDU1[12];
            xt1=EHDU1[13];
            ADinte1=EHDU1[14];
            a=EHDU1[15];
            alfa_wt=EHDU1[16];
            xsn=EHDU1[17];
            xn2=EHDU1[18];
            alfa_at2=EHDU1[19];
            beta_a2=EHDU1[20];
            d2=EHDU1[21];
            db2=EHDU1[22];
            df2=EHDU1[23];
            da2=EHDU1[24];
            st2=EHDU1[25];
            sat2=EHDU1[26];
            sn2=EHDU1[27];
            san2=EHDU1[28];
            xt2=EHDU1[29];
            ADinte2=EHDU1[30];
            dw1=EHDU1[31];
            dw2=EHDU1[32];
            beta_w=EHDU1[33];
            eps_alfa=EHDU1[34];
            eps_beta=EHDU1[35];
            eps_gama=EHDU1[36];
            b1=EHDU1[37];
            cs_n=EHDU1[38];
            alfa_n=EHDU1[39];

            // Elementele geometrice ale angrenajului echivalent pentru treapta I
            zn1=EHDU1[40];
            zn2=EHDU1[41];
            xn1min=EHDU1[42];
            xn2min=EHDU1[43];
            dn1=EHDU1[44];
            dn2=EHDU1[45];
            dbn1=EHDU1[46];
            dbn2=EHDU1[47];
            dan1=EHDU1[48];
            dan2=EHDU1[49];
            alfa_wn=EHDU1[50];
            awn=EHDU1[51];
            eps_alfa_n=EHDU1[52];

            // Elementele de control ale angrenajului treptei I
            alfa_Nt1=EHDU1[53];
            alfa_Nt2=EHDU1[54];
            Ncalc1=EHDU1[55];
            Ncalc2=EHDU1[56];
            N1=EHDU1[57];
            N2=EHDU1[58];
            WNn1=EHDU1[59];
            WNn2=EHDU1[60];
            WNt1=EHDU1[61];
            WNt2=EHDU1[62];
            roNt1=EHDU1[63];
            roNt2=EHDU1[64];
            roAt1=EHDU1[65];
            roEt2=EHDU1[66];
            roat1=EHDU1[67];
            roat2=EHDU1[68];

            // Tensiuni si tensiuni admisibile pentru materialele treptei I
            sigma_H=EHDU1[69];
            sigma_HP=EHDU1[70];
            sigma_F1=EHDU1[71];
            sigma_F2=EHDU1[72];
            sigma_FP1=EHDU1[73];
            sigma_FP2=EHDU1[74];

            // Volumul danturiilor, mm^3
            VolDantura1=EHDU1[75];
            VolDantura2=EHDU1[76];

            // Forte si momente
            Ft1=EHDU1[77];
            Fa1=EHDU1[78];
            Fr1=EHDU1[79];
            Mic1=EHDU1[80]; MH7=Mic1;

            Ft2=EHDU1[81];
            Fa2=EHDU1[82];
            Fr2=EHDU1[83];
            Mic2=EHDU1[84]; MH3=Mic2;
            break;
        } // aici se termina while-ul


		if (valid_angrenaj==1){


            Masaz1=(VolDantura1+Volum_cilindru(df1,b1))*romat;
            Masaz2=(VolDantura2+Volum_cilindru(df2,b2))*romat;



            MasaAngrenaj=Masaz1+Masaz2;





            //=====================================================================================================================
            // FUNCTIA OBIECTIV
            //=====================================================================================================================





            obj=MasaAngrenaj;
            //obj=PI*romat*(dw1*dw1*b1+dw2*dw2*b2)/4.;





            //=====================================================================================================================
            // RESTRICTIILE PROBLEMEI
            //=====================================================================================================================





            // R1 Eroarea relativa a raportului de transmitere trebuie sa fie in intervalul [-2,5%...+2,5%], daca u12<4 sau in [-3%...+3%] in caz contrar.
            if (u12<4.) g[1]=fabs(i12STAS-u12)*40./i12STAS-1.;
            else g[1]=fabs(i12STAS-u12)*100./i12STAS/3.-1.;
            // R2 Verificarea la presiunea de contact.
            g[2]=sigma_H/sigma_HP/(1+tolStrOfMat)-1.;
            // R3 Verificarea la incovoiere a dintelui pinionului.
            g[3]=sigma_F1/sigma_FP1/(1+tolStrOfMat)-1.;
            // R4 Verificarea la incovoiere a dintelui rotii.
            g[4]=sigma_F2/sigma_FP2/(1+tolStrOfMat)-1.;
            // R5 Verificarea danturii pinionului la subtaiere.
            if (xn1>0) g[5]=(14-zn1)/17/xn1-1.;
            else if (xn1<0) g[5]=1.-(14-zn1)/17./xn1;
            else g[5]=14./zn1-1.;
            // R6 Verificarea danturii rotii la subtaiere.
            if (xn2>0) g[6]=(14-zn2)/17./xn2-1.;
            else if (xn2<0) g[6]=1.-(14-zn2)/17./xn2;
            else g[6]=14./zn2-1.;
            // R7 Verificarea danturii pinionului la ascutire.
            g[7]=cs_n*mn/san1-1.;
            // R8 Verificarea danturii rotii la ascutire.
            g[8]=cs_n*mn/san2-1.;
            // R9 Gradul de acoperire frontal trebuie sa fie mai mare decat o valoare minima impusa (in general in functie de viteza anfrenajului).
            g[9]=eps_alfa_lim/eps_alfa-1.;
            // R10 Se verifica daca xn2 este in intervalul [-0.5, 1].
            g[10]=fabs(xn2-0.25)/0.75-1.;
            // R11-16 Pentru masurarea cotei peste dinti trebuie indeplinite conditiile.
            g[11]=(WNn1*sin(beta_b)+5.)/b1-1.;
            g[12]=(WNn2*sin(beta_b)+5.)/b2-1.;
            g[13]=roAt1/roNt1-1.;
            g[14]=roNt1/roat1-1.;
            g[15]=roEt2/roNt2-1.;
            g[16]=roNt2/roat2-1.;
            // R17 z1 si z2 trebuie sa fie prime intre ele.
            g[17]=Verif_cmmdc(z1,z2);

            nr_restrictii=18;
            pen=0;
            for (i=1;i<=nr_restrictii; i++) {if (g[i]>0) pen +=g[i];}
            obj += pen*CoefPen;
            cout << "Merge" << '\n';
            return obj;
		}
		else
        {
            cout << 2 << '\n';
            return ValIesire;
        }
}
int main()
{
    double result = afiseaza();
    cout << result;
    return 0;
}
// ************************* End of file **********************
