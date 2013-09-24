function [f_p] = plasmafreq(n_e)
%PLASMAFREQ Compute the plasma frequency (in Hertz) from n_e.

%--------------------------------------------------------------------------
% Physical constants

e = 1.602*10^(-19);
e_0 = 8.854*10^(-12);
m_e = 9.109*10^(-31);

%--------------------------------------------------------------------------
% Plasma parameters

n_e = n_e*10^6; % (conversion from cm^-3)

%--------------------------------------------------------------------------
% Plasma frequency

w_p = sqrt(n_e*e^2/(e_0*m_e));
f_p = w_p/(2*pi);

end

