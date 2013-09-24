function [Z_p] = parallell(Z1, Z2)
%PARALLELL Compute impedance of parallel coupling
%   PARALLELL(Z1, Z2) computes the equivalent impedance of the parallel
%   coupling of impedances Z1 and Z2. If one of Z1 and Z2 is inf, the one
%   which is not inf is returned. If both are inf, then inf is returned. If
%   one or both are zero, zero is returned. The calculation is done
%   elementwise to allow for vector inputs.

Z_p = zeros(length(Z1), 1);
for i = 1:length(Z1)
    if (Z1(i) == inf)
        Z_p(i) = Z2(i);
    elseif (Z2(i) == inf)
        Z_p(i) = Z1(i);
    elseif (Z1(i) == 0 || Z2(i) == 0)
        Z_p(i) = 0;
    else
        Z_p(i) = Z1(i)*Z2(i)/(Z1(i) + Z2(i));
    end
end
        

end

