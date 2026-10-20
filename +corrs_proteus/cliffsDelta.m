function d = cliffsDelta(x, y)
%cliffsDelta  Cliff's delta, P(X>Y) - P(X<Y) (Cliff 1993, Psychol Bull 114:494, doi:10.1037/0033-2909.114.3.494).
    x = x(:); y = y(:);
    d = sum(sign(x - y.'), 'all') / (numel(x) * numel(y));
end
