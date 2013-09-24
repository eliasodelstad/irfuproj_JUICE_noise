function [varargout] = noise(n_e, T_e, V_bias, varargin)
%EO.NOISE Compute and plot the output voltage spectral densities of all the
%noise components in the small signal integrated noise model, together with
%the impedance of the QTN antenna model (Z_A), the IU-curve and the
%associated resistances (R_s, R_ph).
%   eo.noise(n_e, T_e, V_bias) computes and plots as described above with the
%   electron density n_e [cm^-3] and temperature T_e [eV], plus a bias
%   potential of V_bias [V].

%--------------------------------------------------------------------------
% Physical constants
k_B = 1.381*10^(-23);
e = 1.602*10^(-19);

% Settings
V_B = (-250:0.001:100)';
f = eo.f_sample(0.01, 10, 1) + 0.004;     % f/f_p
f_p = eo.plasmafreq(n_e);
L = 6;
L_Dm = eo.debye(n_e, T_e);

%--------------------------------------------------------------------------
% Quasi-thermal noise and antenna impedance
[V2_qtn, Z] = eo.qtnmod(n_e, T_e, f);
V_qtn = sqrt(V2_qtn);

% Probe resistance and currents due to ambient plasma particles
% Electron current
I_E = eo.oml(n_e, T_e, 0.01);
% Ion Current
if (nargin >= 4)
    m_i = cell2mat(varargin(1));
    if (nargin >= 6)
        v_Di = cell2mat(varargin(3));
    else
        v_Di = 0.01;
    end
    I_I = eo.oml(n_e, T_e, v_Di, 1, m_i);
else
    I_I =zeros(length(V_B), 1);
end
% Total ambient probe-plasma current
I_S = I_I + I_E;
% Interpolate
I_i_pp = spline(V_B, I_I);
I_e_pp = spline(V_B, I_E);
I_S_pp = spline(V_B, I_S);
% Differentiate
Cond_S_pp = fnder(I_S_pp);
% Evaluate
I_s = ppval(I_S_pp, V_bias);
I_e = ppval(I_e_pp, V_bias);
I_i = ppval(I_i_pp, V_bias);
R_s = 2*1./(ppval(Cond_S_pp, V_bias))*ones(length(Z), 1); % Serial connection

% Photo-electron contributions
% Photoelectron current
I_PH = eo.grard(V_B);
% Interpolate
I_PH_pp = spline(V_B, I_PH);
% Differentiate
Cond_PH_pp = fnder(I_PH_pp);
% Evaluate
I_ph = ppval(I_PH_pp, V_bias);
R_ph = 2*1./(ppval(Cond_PH_pp, V_bias))*ones(length(Z), 1); % Serial connection
if (isnan(R_ph))
    R_ph = Inf*ones(length(Z), 1);
end

% Johnson-Nyquist noise
T_ph = (e/k_B)*1.5;
V_ph = sqrt(4*k_B*T_ph*R_ph);
T_s = (e/k_B)*T_e;
V_s = sqrt(4*k_B*T_s*R_s);

% Shot noise
I2_shot_ph = 2*e*abs(I_ph);
% I2_shot_s = 2*e*abs(I_s);
I2_shot_s = 2*e*(abs(I_e) + abs(I_i));  % Assume uncorrelated ion and electron contributions
Z_tot = eo.parallell(R_s, eo.parallell(R_ph, Z));

% Output voltages
V_out_qtn = eo.voltage_division(V_qtn, Z, eo.parallell(R_s, R_ph));
V2_out_qtn = V_out_qtn.*conj(V_out_qtn);
V_out_ph = eo.voltage_division(V_ph, R_ph, eo.parallell(R_s, Z));
V2_out_ph = V_out_ph.*conj(V_out_ph);
V_out_s = eo.voltage_division(V_s, R_s, eo.parallell(R_ph, Z));
V2_out_s = V_out_s.*conj(V_out_s);
V2_out_shot_s = I2_shot_s.*Z_tot.*conj(Z_tot);
V2_out_shot_ph = I2_shot_ph.*Z_tot.*conj(Z_tot);

% if (max(V2_out_ph) < 10^(-21))
%     V2_out_ph = zeros(length(V2_out_ph), 1);
% end

if (nargout == 0)
    % Plot noise levels
    figure('Position', [0 450 640 550])
    [l] = loglog(f, V2_out_qtn, f, V2_out_s, f, V2_out_shot_s, f, V2_out_ph, ...
        f, V2_out_shot_ph, 'LineWidth', 1.2);
    set(l(1), 'Color', 'b')
    set(l(2), 'Color', [0 0.5 0])
    set(l(3), 'Color', get(l(2), 'Color'), 'LineStyle', '--')
    set(l(4), 'Color', [0.8 0.8 0])
    set(l(5), 'Color', get(l(4), 'Color'), 'LineStyle', '--')
    YLim0 = get(gca, 'YLim');
    set(gca, 'FontSize', 14, 'YTick', 10.^(ceil(log10(YLim0(1))):floor(log10(YLim0(2)))), ...
        'YMinorTick', 'on')
    xlabel(['Normalized frequency $$f/f_p$$ ($$f_p =$$ ' sprintf('%1.2f', ...
        f_p*10^(-fix(log10(f_p)))) '$$\cdot 10^' num2str(fix(log10(f_p))) '$$ Hz)'], ...
        'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
    ylabel('Voltage spectral density [V$$^2/$$Hz]', 'FontSize', 18, 'interpreter', ...
        'latex', 'unit', 'character')
    if (nargin >= 5)
        if (length(cell2mat(varargin(2))) >= 1)
            axis(cell2mat(varargin(2)));
        else
            axis tight
        end
    else
        axis tight
    end
    leg1 = legend('Quasi-thermal noise', 'Ambient plasma Nyquist noise', ...
        'Ambient plasma shot noise', 'Photoelectron Nyquist noise', ...
        'Photoelectron shot noise', 'location', 'SouthWest');
    set(leg1, 'Interpreter', 'latex', 'box', 'off')
    if (nargin >= 6)
        v_D = [sprintf('\n'), '$$v_{D_i} =$$ ', ...
        num2str(v_Di, '%g'), ' km/s'];
    else
        v_D = '';
    end
    text('Interpreter', 'latex', 'String', ['$$n_e = $$', num2str(n_e, '%g'), ...
        ' cm$$^{-3}$$', sprintf('\n'), '$$T_e =$$ ', num2str(T_e, '%g'), ' eV', ...
        sprintf('\n'), '$$V_{bias} =$$ ', num2str(V_bias, '%g'), ' V' sprintf('\n'), ...
        '$$L/\lambda_{D_m} = $$ ', num2str(L/L_Dm,3), v_D], 'Units', 'normalized', 'Position', [0.97 0.97], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14, 'BackgroundColor', [1 1 1])
    
    % Plot impedances
    figure('Position', [640 450 640 550])
    C = 1./(2*pi*f_p*f.*imag(Z))*10^12; % pF
%     [AX, h1, h2] = plotyy(f, [real(Z) R_s R_ph], f, C, 'loglog', 'semilogx');
    [AX, h1, h2] = plotyy(f-0.004, [real(Z) abs(imag(Z))], f-0.004, C, 'loglog', 'semilogx');
    set(AX, 'FontSize', 14)
    YLim1 = get(AX(1), 'YLim');
    set(AX(1), 'YTick', 10.^(ceil(log10(YLim1(1))):floor(log10(YLim1(2)))), 'YColor', 'k', 'Color', 'none');
    set(h1, 'LineWidth', 1.2)
    set(h1(1), 'Color', 'k')
    set(h1(2), 'Color', 'k', 'LineStyle', '--')
%     set(h1(2), 'Color', [0 0.5 0])
%     set(h1(3), 'Color', [0.8 0.8 0])
    axes(AX(1))
    xlabel(['Normalized frequency $$f/f_p$$ ($$f_p =$$ ' sprintf('%1.2f', ...
        f_p*10^(-fix(log10(f_p)))) '$$\cdot 10^' num2str(fix(log10(f_p))) '$$ Hz)'], ...
        'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
    ylabel('Impedance ($$\Omega$$)', 'FontSize', 18, 'interpreter', ...
        'latex', 'unit', 'character')
    set(h2, 'LineWidth', 1.2, 'Color', 'r')
    set(AX(2), 'YColor', get(h2, 'Color'), 'YTick', [0:10], 'YMinorTick', 'on', 'Color', 'w')
    axes(AX(2))
    ylabel('Capacitance (pF)', 'FontSize', 18, 'interpreter', ...
        'latex', 'unit', 'character')
    text('Interpreter', 'latex', 'String', ['$$n_e = $$', num2str(n_e, '%g'), ...
        ' cm$$^{-3}$$', sprintf('\n'), '$$T_e =$$ ', num2str(T_e, '%g'), ' eV', ...
        sprintf('\n'), '$$V_{bias} =$$ ', num2str(V_bias, '%g'), ' V' sprintf('\n'), ...
        '$$L/\lambda_{D_m} = $$ ', num2str(L/L_Dm,3), v_D], 'Units', 'normalized', 'Position', [0.97 0.97], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14, 'BackgroundColor', [1 1 1])
%     leg2 = legend([h1(1) h2 h1(2) h1(3)], {'Antenna resistance', 'Antenna capacitance', ...
%         'Ambient electron resistance', 'Photoelectron resistance'}, 'location', 'SouthWest');
    leg2 = legend([h1(1) h1(2) h2], {'Antenna resistance (Re$$\{ Z \}$$)', 'Reactance (magnitude, $$|$$Im$$\{ Z \}|$$)', 'Antenna capacitance'});
    set(leg2, 'Interpreter', 'latex', 'box', 'off', 'location', 'NorthWest')
    uistack(AX(1), 'top');
    
    % Plot I-V curve (for oml)
    figure('Position', [0 450 640 550])
    [m] = plot(V_B, [I_I I_E I_PH, I_S+I_PH], 'LineWidth', 1.2);
    set(m(1), 'Color', 'b');
    set(m(2), 'Color', 'r');
    set(m(3), 'Color', [0.8 0.8 0]);
    set(m(4), 'Color', 'k', 'LineStyle', ':');
    set(gca, 'FontSize', 14)
    YLim2 = get(gca, 'Ylim');
    LimX = [max(-5*T_e, min(V_B)) min(5*T_e, max(V_B))];
    axis([LimX(1) LimX(2) (LimX(1)/(-100))*YLim2(1) (LimX(2)/100)*YLim2(2)]);
    xlabel('Probe voltage (V)', 'FontSize', 18, 'interpreter', 'latex')
    ylabel('Probe currents (A)', 'FontSize', 18, 'interpreter', 'latex')
    [~, I] = min(abs(I_S + I_PH));
    V_F = V_B(I);
    if (V_F <= min(V_B))
        V_f = ['$$<$$ ', num2str(min(V_B), '%1.0f'), ' V'];
    elseif (V_F >= max(V_B))
        V_f = ['$$>$$ ', num2str(max(V_B), '%1.0f'), ' V'];
    else
        V_f = ['$$=$$ ', num2str(V_F, '%1.3f'), ' V'];
    end
    text('Interpreter', 'latex', 'String', ['$$n_e = $$', num2str(n_e, '%g'), ...
        ' cm$$^{-3}$$', sprintf('\n'), '$$T_e =$$ ', num2str(T_e, '%g'), 'eV', ...
        v_D, sprintf('\n'), '$$V_{\textrm{float}}$$ ', V_f], 'Units', 'normalized', ...
        'Position', [0.97 0.97], 'HorizontalAlignment', 'right', 'VerticalAlignment', ...
        'top', 'FontSize', 14, 'BackgroundColor', [1 1 1])
    leg3 = legend('Ion current', 'Electron current', 'Photoelectron current', ...
        'Total probe current');
    set(leg3, 'interpreter', 'latex', 'location', 'NorthWest', 'box', 'off')
    
    % Plot resistances
    R_S = 2*1./(ppval(Cond_S_pp, V_B));
    R_PH = 2*1./(ppval(Cond_PH_pp, V_B));
    figure('Position', [640 450 640 550])
    [n] = semilogy(V_B, [R_S R_PH], 'LineWidth', 1.2);
    set(gca, 'FontSize', 14, 'YMinorTick', 'on');
    set(n(1), 'Color', [0 0.5 0]);
    set(n(2), 'Color', [0.8 0.8 0]);
    axis([min(LimX(1), -1) max(LimX(2), 1) min(min([R_S R_PH]))/10 max(max(R_S), 6.2*10^11)]);
    xlabel('Probe voltage (V)', 'FontSize', 18, 'interpreter', 'latex')
    ylabel('Resistance ($$\Omega$$)', 'FontSize', 18, 'interpreter', 'latex')
    text('Interpreter', 'latex', 'String', ['$$n_e = $$', num2str(n_e, '%g'), ...
        ' cm$$^{-3}$$', sprintf('\n'), '$$T_e =$$ ', num2str(T_e, '%g'), 'eV', ...
        v_D, sprintf('\n'), '$$V_{\textrm{float}}$$ ', V_f], 'Units', 'normalized', ...
        'Position', [0.97 0.97], 'HorizontalAlignment', 'right', 'VerticalAlignment', ...
        'top', 'FontSize', 14, 'BackgroundColor', [1 1 1])
    leg4 = legend('Ambient plasma', 'Photoelectrons');
    set(leg4, 'interpreter', 'latex', 'location', 'NorthWest', 'box', 'off')
    
else
    varargout = {real(V2_out_qtn + V2_out_ph + V2_out_s + V2_out_shot_s + V2_out_shot_ph)};
end

end

