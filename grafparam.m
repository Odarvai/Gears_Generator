function result = grafparam(Gu, lin, col, P, parameter, x)
    found = 0;
    result = 0;
    valid = 1;
    for i = 1 : lin/2
        if P(i) == parameter
            for j = 1 : col
                disp(i + "  " + j);
                if Gu(2 * (i - 1) * col + j ) == x
                    disp("if 1 ");
                    result = Gu((2 * i + 1) * col + j);
                    found = 1;
                    break;
                end
                
                if ((j ~= col) && (Gu(2 * (i-1) * col + j) < x) && (Gu(2 * (i-1) * col + j + 1) > x)) %Aici e gresit
                    disp("if 2 " + Gu(2 * (i-1) * col + j) + " " + Gu((2 * (i-1) + 1) * col + j) + " " + Gu(2 * (i-1) * col + j + 1) + " " + Gu((2 * (i-1) + 1) * col + j + 1) + " " + i + " " + col + " " + j + " " + (2*(i-1)*col+j) );
                    
                    result = ord(Gu(2 * (i-1) * col + j), Gu((2 * (i-1) + 1) * col + j), Gu(2 * (i-1) * col + j + 1), Gu((2 * (i-1) + 1) * col + j + 1), x); %Si aici e gresit
                    found = 1;
                    break;
                end
            end
            
        end
        if ((found == 1) || (i == lin / 2))
            break;
        end
        if ( ((P(i) < parameter) && (P(i+1) > parameter)) || ((P(i+1) < parameter) && (P(i) > parameter)) )
            try
                if((x < max(Gu(2*(i-1)*col+1), Gu((2*(i-1)+2)*col+1))) || (x > min(Gu((2*(i-1)+1)*col), Gu((2*(i-1)+3)*col)))) %scadem 1 din fiecare i
                disp("o puscat in primul if " + i + ' ' + max(Gu(2*(i-1)*col+1), Gu((2*(i-1)+2)*col+1)));
                    valid = 0;
                    break;
                end
            catch
                warning("Index exceeds the number of array elements (12).");
                %continue;
            end
                for j = 1 : col
                    if ((Gu(2 * (i-1) * col + j ) <= x) && (Gu(2 * (i-1) * col + j + 1) >= x ) ) %scadem 1 din i; 
                        disp("else 1 for 1 " + i + " " + j);
                        Y1 = ord(Gu(2 * (i-1) * col + j ), Gu((2 * (i-1) + 1) * col + j ), Gu(2 * (i-1) * col + j +1), Gu((2 * (i-1) + 1) * col + j +1), x); %scadem 1 din i; 
                        break;
                    end
                end
                for j = 1 : col
                    if ( (Gu((2 * (i-1) + 2) * col + j ) <= x) && (Gu((2 * (i-1) + 2) * col + j +1) >= x)) %scadem 1 din i; 
                        disp("else 1 for 2 " + i + " " + j);
                        Y2 = ord(Gu((2 * (i-1) + 2) * col + j ), Gu((2 * (i-1) + 3) * col + j ), Gu((2 * (i-1) + 2) * col + j + 1), Gu((2 * (i-1) + 3) * col + j +1), x); %scadem 1 din i; 
                        break;
                    end
                end
                disp("o mers in al doilea if mare");
                result = ord(P(i), Y1, P(i + 1), Y2, parameter);
                found = 1;
                break;
            %end
        end
    end
    if ((valid == 1) && (found == 0))
        if (lin ~= 2)
            if (( x < max(Gu(1), Gu(lin * col))) || (x > min(Gu(col), Gu(lin * col-1)))) %am sters -1
                disp("o puscat in al doilea if");
                valid = 0;
              
            else
                for j = 1 : col
					if ((Gu(0 * col + j) <= x) && (Gu(0 * col + j + 1) >= x))
                        disp("else 2 if 1 ");
						Y1 = ord(Gu(0 * col + j), Gu(1 * col + j), Gu(0 * col + j + 1), Gu(1 * col + j + 1), x);
						break;
                    end
                end
                for j = 1 : col
                    if ((Gu((lin - 1) * col + j) <= x) && (Gu((lin - 1) * col + j + 1) >= x)) %Aici am pus -1 in loc de -2 cand am scazut din lin
                        disp("else 2 if ");
                        Y2 = ord(Gu((lin - 1) * col + j), Gu((lin) * col + j), Gu((lin - 1) * col + j + 1), Gu((lin) * col + j + 1), x); %Si aici am adunat peste tot cu 1 unde era -1 sau -2
                        break;
                    end
                end
                result = ord(P(0), Y1, P(lin/2), Y2, parameter); %Am inlocuit P(lin/2-1)
                found = 1;
            end
        else % else 3
            if ((x < Gu(1)) || (x > Gu(col)))  %Am inlocuit Gu(col-1) si Gu(0) cu Gu(1)
                valid = 0;
                disp("o puscat in al treilea if "+ col + " " +Gu(col));
                if x < Gu(1)
                    disp(x + "<" + Gu(1));
                end
                if x > Gu(col)
                    disp(col);
                    digits(10);
                    disp(vpa(Gu(col), 5));
                    disp(Gu(1)+ " si " + Gu(col));
                end
              
            else
                for j = 1 : col
					if ((Gu(j) <= x) && (Gu(j+1) >= x))
                        V = [j, col+j, j+1, col+j+1];
                        X = [Gu(j), Gu(col + j), Gu(j + 1), Gu(col + j + 1)];
                        disp(V);
                        disp(X);
                        disp("rezultatu iese din ultimu else");
						result = ord(Gu(j), Gu(col + j), Gu(j + 1), Gu(col + j + 1), x);  %Aici poate trebe peste tot adaugat 1 la j dar nu-s sigur
						found = 1;
						break;
                    end
                end
            end
        end
    end   
    
end