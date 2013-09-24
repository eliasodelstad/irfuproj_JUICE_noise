function [varargout] = oml(n_j, T_j, v_Dj, varargin)
%OML Computes current from spherical Langmuir probe due to single particle species under OML approximation.
%   OML(n_j, T_j, v_Dj, {q_j, m_j, V_B, r_p}) computes and plots the particle
%   current from given plasma parameters (from equations (7.16) - (7.21) in
%   H?ymork et al) particle density n_j, temperature T_j, drift velocity
%   v_Dj, charge q_j, mass m_j and probe radius r_p at bias voltage points
%   in vector V_B.


%--------------------------------------------------------------------------
% Physical constants

e = 1.602*10^(-19);
k = 1.381*10^(-23);
m_e = 9.109*10^(-31);
m_H = 1.67*10^(-27);

%--------------------------------------------------------------------------
% Plasma characteristics

% Density (conversion from cm^(-3))
n_j = n_j*10^6;
% Temperature in K (conversion from eV)
T_j = (e/k)*T_j;
% Relative drift velocity (conversion from km/s)
v_Dj = v_Dj*10^3;
% Particle charge
if (nargin >= 4)
    q_j = cell2mat(varargin(1))*e;
else
    q_j = -e;
end
% Particle mass
if (nargin >= 5)
    m_j = cell2mat(varargin(2))*m_H;
else if (q_j == -e)
        m_j = m_e;
    else
        m_j = m_H;
    end
end

%--------------------------------------------------------------------------
% Probe characteristics

% Bias voltage (probe potential)
if (nargin >= 6)
    V_B = cell2mat(varargin(3));
else
    V_B = (-250:0.001:100);
end
% Probe radius
if (nargin >= 7)
    r_p = cell2mat(varargin(4));
else
    r_p = 0.05;
end
% Probe surface area
A = 4*pi*r_p^2;


%--------------------------------------------------------------------------
% Current

% Ratio of relative drift velocity to mean speed of ions
S = v_Dj/sqrt(2*k*T_j/m_j);
% Normalized potential
X_j = q_j*(V_B)/(k*T_j);
% Saturation current
I_j0 = -A*n_j*q_j*sqrt(k*T_j/(2*pi*m_j));
% Current
I_j = zeros(length(V_B), 1);
% Negative bias
X_1 = X_j(q_j*V_B <= 0);
I_j(q_j*V_B <= 0) = I_j0*(1/2)*(sqrt(pi)*(S + 1/(2*S) - X_1/S)*erf(S)...
    + exp(-S^2));
% Positive bias
X_2 = X_j(q_j*V_B >= 0);
I_j(q_j*V_B >= 0) = I_j0*1/(4*S)*(sqrt(pi)*(S^2 + 1/2 - X_2).*(erf(S +...
    sqrt(X_2)) - erf(sqrt(X_2) - S)) + (S - sqrt(X_2)).*exp(-(S + ...
    sqrt(X_2)).^2) + (S + sqrt(X_2)).*exp(-(sqrt(X_2) - S).^2));
    

%--------------------------------------------------------------------------
% Give output

if (nargout == 0)
    % Set color
    if (q_j > 0)
        s = 'g';
    else
        s = 'r';
    end
    figure('Position', [560 450 640 480])
    plot(V_B, I_j, s, 'LineWidth', 1.2)
    set(gca, 'FontSize', 14)
    xlabel('Probe bias voltage [V]', 'FontSize', 14, 'Interpreter', 'latex')
    ylabel('Probe current [A]', 'FontSize', 14, 'Interpreter', 'latex')
else
    varargout = {I_j};
end


end





