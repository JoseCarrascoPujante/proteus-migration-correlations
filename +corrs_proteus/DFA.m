function alpha = DFA(R)
%DFA  Detrended fluctuation analysis exponent gamma (Peng et al. 1994, Phys Rev E 49:1685,
%   doi:10.1103/PhysRevE.49.1685), boxes cut from both ends of the series (Kantelhardt et al. 2001,
%   Physica A 295:441, doi:10.1016/S0378-4371(01)00144-3).
    arguments
        R (:, 1) double {mustBeReal}
    end
    N = numel(R);
    n = round(10 .^ (2:0.1:3));
    assert(floor(N / n(end)) >= 3, 'DFA:short', ...
        'N = %d leaves fewer than 3 boxes per end at %d frames.', N, n(end));
    z = cumsum(R - mean(R));
    D = zeros(numel(n), 1);
    for i = 1:numel(n)
        L = n(i);
        nBox = floor(N / L);
        t = (1:L)';
        X = [ones(L, 1) t];
        Z = [reshape(z(1:nBox * L), L, nBox), reshape(z(N - nBox * L + 1:N), L, nBox)];
        res = Z - X * (X \ Z);
        D(i) = sqrt(mean(res(:) .^ 2));
    end
    p = polyfit(log10(n(:)), log10(D), 1);
    alpha = p(1);
end
