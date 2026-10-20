function Table1 = table1_modeledCouplings()
%table1_modeledCouplings  Table 1: R*_n with BCa 95% CI and Spearman rho of every dependent coupling,
%   rounded once, from full precision, to the two decimals of the paper.
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    derivedDir = fullfile(projectRoot, 'data', 'derived');
    dep = readtable(fullfile(derivedDir, 'pairwiseDependence.csv'), 'TextType', 'string');
    dep = sortrows(dep(ismember(lower(string(dep.Coupled)), ["1", "true"]), :), 'Rank');

    sfx = ["NoStim", "Galv", "Chem"];
    Table1 = table(dep.Pair, 'VariableNames', {'Pair'});
    for c = 1:numel(sfx)
        s = sfx(c);
        Table1.(char("bcDcor_CI_" + s)) = compose("%.2f [%.2f, %.2f]", dep.(char("bcDcor_" + s)), ...
            dep.(char("bcDcor_lo_" + s)), dep.(char("bcDcor_hi_" + s)));
        Table1.(char("rho_" + s)) = round(dep.(char("Spearman_" + s)), 2);
    end
    out = fullfile(derivedDir, 'table1_modeledCouplings.csv');
    writetable(Table1, out);
end
