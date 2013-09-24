function [V_out] = voltage_division(V_in, Z_in, Z_out)
%EO.VOLTAGE_DIVISION Compute output voltage from potential divider
%   eo.voltage_division(V_in, Z_in, Z_out) computes the output voltage from a
%   system with input impedance Z_in, output impedance Z_out and fed by an
%   input voltage V_in. The calculation is done elementwise to allow for
%   vector inputs.

V_out = zeros(length(V_in), 1);
for i = 1:length(V_in)
    if (Z_in(i) == inf || Z_out(i) == 0)
        V_out(i) = 0;
    elseif (Z_out(i) == inf || Z_in(i) == 0)
        V_out(i) = V_in(i);
    else
        V_out(i) = V_in(i)*Z_out(i)/(Z_in(i) + Z_out(i));
    end

end

