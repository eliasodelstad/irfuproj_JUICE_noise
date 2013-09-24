function [f] = f_sample(f_min, f_max, res)
%EO.F_SAMPLE Generate frequncy samples with constant relaitive resolution
%   eo.f_sample(f_min, f_max, res) returns an Nx1 vector of frequency samples
%   from f_min to f_max in the form of a geometric sequence, assuring that 
%   the relative sampling length dw/w is constant and equal to res(%).

% Number of samples
N = 1 + log10(f_max/f_min)/log10(1 + 0.01*res);

% Frequency samples
f = f_min*(1 + 0.01*res).^(0:N-1)';

end

