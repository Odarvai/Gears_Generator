function [result] = GetList(x, NrElemList, List)
%Gaseste elementul cel mai apropiat de x din sirul List (care are NrElemList elemente)
%dif_min, dif_i, result;

dif_min = 1.e10;
    for i = 1:NrElemList
        dif_i = abs( x - List(i) );
        if ( dif_i < dif_min )
            result = List(i);
            dif_min = dif_i;
        end
    end

end

