function out = roundSigMin3(x)
%roundSigMin3  Round to 3 decimals, keeping at least one significant figure below 1e-3.
    out = x;
    fin = isfinite(x) & x ~= 0;
    big = fin & abs(x) >= 1e-3;
    small = fin & abs(x) < 1e-3;
    out(big) = round(x(big), 3);
    xs = x(small);
    d  = ceil(-log10(abs(xs)));
    for i = 1:numel(xs)
        xs(i) = round(xs(i), d(i));
    end
    out(small) = xs;
end
