function apen = ApEn(dim, r, s)
%ApEn  Approximate entropy of a series (Pincus 1991, PNAS 88:2297, doi:10.1073/pnas.88.6.2297).
    s = s(:);
    apen = phim(dim, r, s) - phim(dim + 1, r, s);
end

function p = phim(m, r, s)
    N = numel(s); M = N - m + 1;
    v = zeros(M, m);
    for i = 1:m, v(:, i) = s(i:i+M-1); end
    C = zeros(M, 1);
    for i = 1:M
        C(i) = sum(max(abs(v - v(i, :)), [], 2) <= r) / M;
    end
    p = mean(log(C));
end
