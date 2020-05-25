function [primitiva] = primEv(x)
%Primitiva f evolventice
    primitiva = 2 * x * x * x + 3 * (1 - x * x) * sin(2 * x) - 6 * x * cos(2 * x);
end

