function T = couplingStrength()
%couplingStrength  R*_n with BCa 95% CI and Spearman rho of every dependent coupling,
%   rounded once, from full precision, to two decimals.
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    derivedDir = fullfile(projectRoot, 'data', 'derived');
    dep = readtable(fullfile(derivedDir, 'pairwiseDependence.csv'), 'TextType', 'string');
    dep = sortrows(dep(ismember(lower(string(dep.Coupled)), ["1", "true"]), :), 'Rank');

    sfx = ["NoStim", "Galv", "Chem"];
    T = table(dep.Pair, 'VariableNames', {'Pair'});
    for c = 1:numel(sfx)
        s = sfx(c);
        T.(char("bcDcor_CI_" + s)) = compose("%.2f [%.2f, %.2f]", dep.(char("bcDcor_" + s)), ...
            dep.(char("bcDcor_lo_" + s)), dep.(char("bcDcor_hi_" + s)));
        T.(char("rho_" + s)) = round(dep.(char("Spearman_" + s)), 2);
    end
    out = fullfile(derivedDir, 'couplingStrength.csv');
    writetable(T, out);
end
