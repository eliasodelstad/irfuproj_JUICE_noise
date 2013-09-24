function [varargout] = qtnmod(n_e, T_e, varargin)
%EO.QTNMOD Compute quasi-thermal noise and associated antenna impedance for kappa-distributed electrons.
%   eo.qtnmod(n_e, T_e) computes the quasi-thermal noise in a double-sphere 
%   dipole antenna in a plasma where the electron velocities obey a "Kappa-
%   distribution". n_e is the electron density in cm^-3 and T_e is the
%   electron temperature in eV.

%--------------------------------------------------------------------------
% Physical constants

e = 1.602*10^(-19);
e_0 = 8.854*10^(-12);
K = 1.381*10^(-23);
m_e = 9.109*10^(-31);

%--------------------------------------------------------------------------
% Plasma parameters

n_e = n_e*10^6; % (conversion from cm^-3)
T_e = (e/K)*T_e; % (conversion from eV)
f_p = sqrt(n_e*e^2/(e_0*m_e))/(2*pi); % (plasma frequency)
v_T = sqrt(2*K*T_e/m_e); % (electron thermal velocity)
L_Dm = sqrt(e_0*K*T_e/(n_e*e^2)); % (Maxwellian Debye length)

if (nargin == 4)
    kappa = cell2mat(varargin(2));
else
    kappa = 4;
end

% "Thermal" velocity
v_0 = sqrt((K*T_e)/m_e*(2*kappa-3)/kappa);

%--------------------------------------------------------------------------
% Antenna geometry

a = 0.05; % (probe radius)
L = 6; % (antenna length)

% Antenna geometry factor (current distribution)
F_1 = @(x) (1/4)*(1 - (sin(x))./x);
F_a = @(x) (sin(x)./x).^2;  % For finite radius a
%     F_a = @(x) 1;         % For lim(ka --> 0)
F = @(x) F_1(x).*F_a(x);

%--------------------------------------------------------------------------
% Frequency band

if (nargin >= 3)
    f = cell2mat(varargin(1))*f_p;
else
    f = eo.f_sample(0.1, 10, 1)*f_p;
end
w = 2*pi*f;

%--------------------------------------------------------------------------
% Voltage spectral density

% Initialize
V2 = zeros(length(f), 1);
imp = zeros(length(f), 1);

% Singularity offset
if (nargin == 5)
    epsilon = cell2mat(varargin(3));
else
    epsilon = 0.001;
end

% Warnig logs
F_V = [];
F_Z = [];

% Iterate over frequency range
for i = 1:length(f)
    r = f(i)/f_p;

    % Make z and p same size matrices
    P = @(z,p) p'*ones(1,length(z));
    Z = @(z,p) ones(length(p),1)*z;
    % Implement summation term
    S = @(z,p) factorial(kappa + P(z,p))./(factorial(P(z,p)).*(2*1i)...
        .^(kappa+1+P(z,p)).*(Z(z,p)+1i).^(kappa+1-P(z,p)));
    % Dielectric permittivity
    e_L = @(z) 1 + (z./r).^2.*(2*kappa - 1 + ((-2)^(kappa+1))/...
        prod(1:2:2*kappa-3)*1i*z.*sum(S(z,0:kappa)));
    
    % Debye length in kappa-distributed plasma
    L_D = v_0/(2*pi*f_p)*sqrt(kappa/(2*kappa-1));
    u = L/L_D;

    % Voltage spectral density
    V2(i) = 2^(kappa+3)/(pi^2*e_0)*factorial(kappa)/(prod(1:2:2*kappa-3)...
        *sqrt(kappa))*m_e*v_0/r^2*quadgk(@(z) z.*F(r*u./(z*sqrt(2*kappa-1)))...
        .*((1+z.^2).^kappa.*abs(e_L(z)).^2).^(-1), 0, inf, 'MaxIntervalCount', 1000);
    % Find warning
    [~, s] = lastwarn;
    if (strcmp(s, 'MATLAB:quadgk:MaxIntervalCountReached'))
        F_V = f(i);
        lastwarn('', '');
    end
    
    % Impedance
    Q = quadgk(@(z)F_1(w(i)*L./(sqrt(kappa)*v_0*z)).*F_a(w(i)*a./(sqrt(kappa)...
        *v_0*z))./(z.^2.*e_L(z)), 0+epsilon, inf, 'MaxIntervalCount', 1000);
    imp(i) = 4*1i./(pi^2*e_0*sqrt(kappa)*v_0).*Q;
    % Find warning
    [~, s] = lastwarn;
    if (strcmp(s, 'MATLAB:quadgk:MaxIntervalCountReached'))
        F_Z = f(i);
        lastwarn('', '');
    end
    
end

%--------------------------------------------------------------------------
% Output

if (nargout == 1)
    varargout = {V2};
elseif (nargout == 2)
    varargout = {V2, imp};
elseif (nargout == 4)
    varargout = {V2, imp, F_V, F_Z};
else
    figure('Position', [560 450 640 480])
    loglog(f/f_p, V2, 'k', 'LineWidth', 1.2)
    set(gca, 'FontSize', 14)
    xlabel('Normalized frequency $$f/f_p$$', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
    ylabel('Voltage spectral density [V$$^2/$$Hz]', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
    axis tight
    hold on
    text('Interpreter', 'latex', 'String', ['$$\frac{L}{\lambda_{D_m}} = $$', ...
        num2str(L/L_Dm,3)], 'Units', 'normalized', 'Position', [0.97 0.97], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14)
    hold off

end
















