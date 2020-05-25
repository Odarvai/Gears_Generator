function [abscisa] = ord(x1, y1, x2, y2, x)
%Returneaza ordonata punctului de abscisa x de pe dreapta determinata
%de punctele M1(x1,y1) si M2(x2,y2)7
   abscisa = y1 + (x - x1) * (y2 - y1) / (x2 - x1);
end

