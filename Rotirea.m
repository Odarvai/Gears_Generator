function [result] = Rotirea (x, n, nr1, nr2, u[20], Iz_s[20], Iz_d[20], Mi_s[20], Mi_d[20], EE)
%Returneaza valoarea rotirii in punctul de abscisa x pentru un arbore discretizat cu u[n], Iz_s[], Iz_d[] si incarcat cu Mi_s[],Mi_d[],

        j=n-1; % indicele nodului a.i. u[j]<=x<u[j+1]. Se initializeaza ca fiind ultimul nod
		double fi=0;
		
		% Determinarea lui j
		for i = 1:n-1
			if ((u(i)<=x)&&(x<u(i+1)))
				j=i;break;
            end
			else continue;
        end

		% Se insereaza noul nod in siruri (nodul de abscisa x)
		if (j~=n-1)
            for i = 1:n+1
                if (i<=j)
                    v(i)=u(i);
                    vIz_s(i)=Iz_s(i);
                    vMi_s(i)=Mi_s(i);
                    vMi_d(i)=Mi_d(i);
                end
                if (i==j+1) 
                    v(i)=x;
                    vIz_s(i)=Iz_d(i-1);
                    if (u[j+1]-u[j] ~= 0)
                        vMi_s(i) = (Mi_s[j+1]-Mi_s[j])*x+Mi_s[j]*u[j+1]-Mi_s[j+1]*u[j])/(u[j+1]-u[j]);
                    else
                        vMi_s(i) = Mi_s[j];
                    end
                    if (u[j+1]-u[j] ~=0)
                        vMi_d(i) = (Mi_d[j+1]-Mi_d[j])*x+Mi_d[j]*u[j+1]-Mi_d[j+1]*u[j])/(u[j+1]-u[j]);
                    else
                        vMi_d(i) = Mi_d[j];
                    end
                end
                if(i>j+1)
                    v(i)=u(i-1);
                    vIz_s(i)=Iz_s(i-1);
                    vMi_s(i)=Mi_s(i-1);
                    vMi_d(i)=Mi_d(i-1);
                end
            end
        else
			for i = 1:n
				v(i)=u(i);
				vIz_s(i)=Iz_s(i);
				vMi_s(i)=Mi_s(i);
				vMi_d(i)=Mi_d(i);
            end

			v(n)=u(n-1);
			vIz_s(n)=Iz_s(n-1);
			vMi_s(n)=Mi_s(n-1);
			vMi_d(n)=Mi_d(n-1);
        end
		% Determiarea indecsilor reazemelor in noul sir
		if (j~=n-1)
			if(j<nr1)
				vnr1=nr1+1;
				vnr2=nr2+1;
            end
			if((j>=nr1)&&(j<nr2))
				vnr1=nr1;
				vnr2=nr2+1;
            end
			if(j>=nr2)
				vnr1=nr1;
				vnr2=nr2;
            end
        else
			if(nr2~=n-1)
				vnr1=nr1;
				vnr2=nr2;
            else
				vnr1=nr1;
				vnr2=nr2+1;
            end
        end

		% Determinarea indicelui punctului de abscisa x in noile siruri
		if (j~=n-1) ix=j+1;
		else ix=n;
        end
		%Generarea matricei de incarcare cu sarcina unitara
		for i= 1:n+1
			F(i)=0;
			Mic(i)=0;
        end
		Mic(ix)=1;
		% Calculul momentelor in fiecare nod, date de momentul unitar
		R1=Reaction1(n+1,vnr1,vnr2,v,F,Mic);
		R2=Reaction2(n+1,vnr1,vnr2,v,F,Mic);
		F(vnr1)=R1;
		F(vnr2)=R2;

		for i = 1:n+1  
			Mi_unit_s(i)=Moment_i_s (v(i),n+1,v,F,Mic);
			Mi_unit_d(i)=Moment_i_d (v(i),n+1,v,F,Mic);
        end
        
		for i = 1:n 
            fi = fi + (v(i+1)-v(i))*(Mi_unit_d(i)*(2*vMi_d(i)+vMi_s(i+1))+Mi_unit_s(i+1)*(vMi_d(i)+2*vMi_s(i+1)))/vIz_s(i+1);
		fi = fi / (6*EE);
		
		result = fi;

end

