function dcorBiasCorrected = distance_correlation(x, y)
%distance_correlation  Bias-corrected squared distance correlation R*_n (Szekely & Rizzo 2013,
%   doi:10.1016/j.jmva.2013.02.012).
    arguments
        x (:,1) double {mustBeReal}
        y (:,1) double {mustBeReal}
    end

    n = numel(x);
    if numel(y) ~= n
        error("distance_correlation:sizeMismatch", ...
            "x and y must have the same number of elements.");
    end

    a = abs(x - x.');
    b = abs(y - y.');

    Au = uCenter(a);
    Bu = uCenter(b);
    innerAB = sum(Au.*Bu, "all");
    innerAA = sum(Au.*Au, "all");
    innerBB = sum(Bu.*Bu, "all");
    denomBC = sqrt(innerAA*innerBB);
    if denomBC > 0
        dcorBiasCorrected = innerAB/denomBC;
    else
        dcorBiasCorrected = 0;
    end
end

function centered = uCenter(d)
    n = size(d, 1);
    centered = d - sum(d, 2)/(n - 2) - sum(d, 1)/(n - 2) + sum(d, "all")/((n - 1)*(n - 2));
    centered(1:n+1:end) = 0;
end
