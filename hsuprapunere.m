function [result] = hsuprapunere(R1, R2, C, HS)
	% R1 - Raza bosajului din partea stanga, mm
	% R2 - Raza dosajului din partea dreapta, mm
	% C  - Distanta dintre axe, mm
    if(C - R1 - R2 >= 0) HS(1) = 0; HS(2) = R1; HS(3) = R2;
    else
        temp = abs((R1 * R1 - R2 * R2) / C / 2);
        if ( R1 >= R2)  h1 = 0.5 * C + temp;
        else h1 = 0.5 * C - temp;
        end        
        h2 = C - h1;
        HS(1) = 1;
        HS(2) = h1;
        HS(3) = h2;
        result = HS;
    end
end

