% This script plots the quasi-thermal noise (QTN) and associated antenna
% impedance and capacitance for a number of electron densities at
% constant temperature. When the antenna is short, '%1.1f' should be used in
% the legend string formatting. For long antennas this should be replaced
% by '%1.0g'.

f = f_sample(0.1, 10, 1);
n_e = 10.^(4:0.5:5);
f_p = plasmafreq(n_e);
T_e = 1;
V2 = zeros(length(f), length(n_e));
Z = V2;
S = cell(1, length(n_e));
for i = 1:length(n_e);
    [V2(:,i), Z(:,i)] = qtnmod(n_e(i), T_e, f);
    S(i) = {['10$$^{' num2str(log10(n_e(i))) '}$$ eV ($$L/\lambda_{D_m}$$ = ' num2str(6/debye(n_e(i), T_e), '%1.1f') ')']};
end

figure('Position', [0 450 640 480])
loglog(f, V2, 'LineWidth', 1.2)
set(gca, 'FontSize', 14)
xlabel('Normalized frequency $$f/f_p$$', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
ylabel('Voltage spectral density [V$$^2/$$Hz]', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
legend(S, 'interpreter', 'latex', 'location', 'SouthWest', 'box', 'off')

figure('Position', [700 450 640 480])
loglog(f, real(Z), 'LineWidth', 1.2)
set(gca, 'FontSize', 14)
xlabel('Normalized frequency $$f/f_p$$', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
ylabel('Resistance [$$\Omega$$]', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
legend(S, 'interpreter', 'latex', 'location', 'SouthWest', 'box', 'off')

C = 1./(2*pi*f*ones(1, length(T_e))*f_p.*imag(Z));
figure('Position', [0 0 640 480])
semilogx(f, 10^12*C, 'LineWidth', 1.2)
set(gca, 'FontSize', 14)
xlabel('Normalized frequency $$f/f_p$$', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
ylabel('Capacitance [pF]', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character')
legend(S, 'interpreter', 'latex', 'location', 'NorthEast', 'box', 'off')


% Discarded code:
% figure('Position', [0 450 640 480])
% surf(f, n_e, V2', log10(V2'));
% set(gca, 'FontSize', 14, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log', ...
%     'XLim', 10.^[floor(log10(min(f))) ceil(log10(max(f)))], 'YLim', [min(n_e) max(n_e)], 'ZLim', [min(min(V2)) ...
%     max(max(V2))], 'XTick', 10.^(floor(log10(min(f))):1:ceil(log10(max(f)))), ...
%     'YTick', 10.^(ceil(log10(min(n_e))):2:floor(log10(max(n_e)))), 'ZTick', ...
%     10.^(ceil(log10(min(min(V2)))):2:floor(log10(max(max(V2))))));
% view([90 90]);
% X = xlabel('Normalized frequency $$f/f_p$$', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character');
% Y = ylabel('Plasma density $$n_e$$ (cm$$^-3$$)', 'FontSize', 18, 'interpreter', 'latex');
% Z = zlabel('Voltage spectral density [V$$^2/$$Hz]', 'FontSize', 18, 'interpreter', 'latex', 'unit', 'character');
% shading flat
% caxis(log10(get(gca, 'ZLim')));
% set(colorbar, 'Position', [1 2.21 0.03 0.555].*get(gca, 'Position'), 'YTick', []);

% [Y X] = max(V2, [], 1);
% loglog(n_e, Y);
    


