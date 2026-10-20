function [slope, deltat] = msd(x, y)
%msd  Mean squared displacement exponent beta of a two-dimensional trajectory, lags up to N/4
%   (Einstein 1905; Dieterich et al. 2008, PNAS 105:459, doi:10.1073/pnas.0707603105).
    tau = 0.5;
    x = x(:); y = y(:);
    nLag = floor(numel(x) / 4);
    deltat = (1:nLag)' * tau;
    msdpts = zeros(nLag, 1);
    for dt = 1:nLag
        msdpts(dt) = mean((x(1+dt:end) - x(1:end-dt)).^2 + (y(1+dt:end) - y(1:end-dt)).^2);
    end
    p = polyfit(log10(deltat), log10(msdpts), 1);
    slope = p(1);
end
