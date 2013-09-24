function [varargout] = grard(varargin)
%GRARD Computes photoelectron current vs voltage using expression in Grard
%(1973).
%   grard(V, V_p, J_ao, R, A) computes the photoelectron current as
%   predicted by Grard (1973) for a spherical probe at the voltage points 
%   contained in vector V, for photoelectron e-folding energy V_p, saturated
%   current density J_ao (A/m^2 at 1 AU, material property), distance from the sun R (AU) 
%   and sunward projected probe area A (m^2).

% Default parameters
V_B = (-10:0.1:20);
V_p = 1.5;
J_ao = -25*10^(-6);
R = 5.2;
A = pi*0.05^2;

% User-defined parameters
if (nargin >= 1)
    V_B = cell2mat(varargin(1));
end
if (nargin >= 2)
    V_p = cell2mat(varargin(2));
end
if (nargin >= 3)
    J_ao = -cell2mat(varargin(3));
end
if (nargin >= 4)
    R = cell2mat(varargin(4));
end
if (nargin >= 5)
    A = cell2mat(varargin(5));
end

%--------------------------------------------------------------------------
% Photocurrent according to Grard (1973)
% Saturation current
I_ph = A/R^2*J_ao*ones(length(V_B), 1);
% Positive bias
V_B1 = V_B(V_B >= 0);
I_ph(V_B >= 0) = A/R^2*J_ao*(1 + V_B1/V_p).*exp(-V_B1/V_p);

%--------------------------------------------------------------------------
% Give output
if (nargout == 0)
    plot(V_B, I_ph, 'k', 'LineWidth', 1.2)
    set(gca, 'FontSize', 14)
    xlabel('Probe bias voltage [V]', 'FontSize', 14, 'Interpreter', 'latex')
    ylabel('Probe current [A]', 'FontSize', 14, 'Interpreter', 'latex')
else
    varargout = {I_ph};
end
    

end

