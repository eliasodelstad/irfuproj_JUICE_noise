function [lambda_D] = debye(n_e, T_e)
%DEBYE Compute the Debye length (in meters) from n_e and T_e.

%--------------------------------------------------------------------------
% Physical constants

e = 1.602*10^(-19);
K = 1.381*10^(-23);
e_0 = 8.854*10^(-12);

%--------------------------------------------------------------------------
% Plasma parameters

n_e = n_e*10^6; % (conversion from cm^-3)
T_e = (e/K)*T_e; % (conversion from eV)

%--------------------------------------------------------------------------
% Debye length

lambda_D = sqrt(e_0*K*T_e./(n_e*e^2));

end

