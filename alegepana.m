function [ APh, APb, APl, APt1, APt2 ] = alegepana( dtr, ltr )
%alegepana - retureanza dimensiunile geometrice ale penei (b x h x l, t1, t2) corespunzatoare diametrului de arbore dtr

%dtr - diametrul tronsonului arborelui pe care se monteaza pana
%ltr - lungimea tronsonului de arbore 


A = [
        6,8,2,2,6,20,1.2,1;
		8,10,3,3,6,36,1.8,1.4;
		10,12,4,4,8,45,2.5,1.8;
		12,17,5,5,10,56,3,2.3;
		17,22,6,6,14,70,3.5,2.8;
		22,30,8,7,18,90,4,3.3;
		30,38,10,8,22,110,5,3.3;
		38,44,12,8,28,140,5,3.3;
		44,50,14,9,36,110,5.5,3.8;
		50,58,16,10,45,180,6,4.3;
		58,65,18,11,50,200,7,4.4;
		65,75,20,12,56,220,7.5,4.9;
		75,85,22,14,63,250,9,5.4;
		85,95,25,14,70,280,9,5.4;
		95,110,28,16,80,320,10,6.4;
		110,130,32,18,90,360,11,7.4
     ];
 
[row , ~] = size(A);
 
L = [6,8,10,12,14,16,18,20,22,25,28,32,36,40,45,50,56,63,70,80,90,100,110,125,140,160,180,200,220,250,280,320,360,400];

for i = 1 : row
    if ((A(i,1) < dtr) && ( dtr <= A(i,2)))
        b = A(i,3);
        h = A(i,4);
        lmin = A(i,5);
        lmax = A(i,6);
        t1 = A(i,7);
        t2 = A(i,8);
        break;
    end
end

l = lmin;

for i = 1 : 34
    if (lmin == L(i))
        for j = i : 34
            if ((L(j) <= lmax) && (L(j) < ltr))
                l = L(j);
            end
        end
    end
end

APb = b; APh = h; APl = l; APt1 = t1; APt2 = t2;

end



