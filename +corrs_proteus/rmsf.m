function [alpha, corrTime] = rmsf(u, dt)
%rmsf  Root mean square fluctuation exponent alpha and correlation time of a move-step series sampled
%   every dt s (Viswanathan et al. 1996, Nature 381:413, doi:10.1038/381413a0). The first 3,500 steps are
%   used, or the whole series if shorter. Alpha is the log-log slope from the first lag up to each of 24
%   fitting ranges that grow by 62.5 s (125 steps at 0.5 s); corrTime, in steps, is the longest range
%   fitted with R^2 > 0.98.
    arguments
        u (:,1) double {mustBeReal}
        dt (1,1) double {mustBePositive}
    end
    u = u(1:min(3500, numel(u)));
    Y = cumsum([0; u(:)]);
    dim = numel(Y);
    F = zeros(dim, 1);
    for l = 1:dim
        m = dim - l;
        if m < 1, F(l) = NaN; continue; end
        dY = Y(l+1:dim) - Y(1:m);
        F(l) = sqrt(mean(dY.^2) - mean(dY)^2);
    end
    time = (1:dim)' * dt;
    step = round(62.5 / dt); nWindows = 24;
    junto = zeros(nWindows, 3);
    zz = 0;
    for timeMax = step:step:nWindows*step
        zz = zz + 1;
        idx = find(time(1:timeMax) > 0 & F(1:timeMax) > 0);
        logL = log10(time(idx)); logF = log10(F(idx));
        p = polyfit(logL, logF, 1);
        Rsq = 1 - sum((logF - polyval(p, logL)).^2) / ((numel(logF) - 1) * var(logF));
        junto(zz, :) = [p(1) Rsq timeMax];
    end
    dev = junto(junto(:, 2) > 0.98, :);
    alpha = dev(end, 1);
    corrTime = dev(end, 3);
end
