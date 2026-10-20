function dt = resolvedInterval(coordinates)
%resolvedInterval  Shortest multiple of the sampling interval at which fewer than 1% of the steps are
%   zero (subsampling as in Gorelik & Gautreau 2014, Nat Protoc 9:1931, doi:10.1038/nprot.2014.131).
    arguments
        coordinates (1, 1) struct
    end
    maxZeroShare = 0.01;
    samplingInterval = 0.5;   % positions are recorded every 0.5 s
    conds = fieldnames(coordinates);
    nFrames = size(coordinates.(conds{1}).x, 1);
    for step = 1:floor(nFrames / 2)
        [nZero, nStep] = deal(0);
        for c = 1:numel(conds)
            X = coordinates.(conds{c}).x(1:step:end, :);
            Y = coordinates.(conds{c}).y(1:step:end, :);
            d = hypot(diff(X), diff(Y));
            nZero = nZero + nnz(d == 0);
            nStep = nStep + numel(d);
        end
        if nZero / nStep < maxZeroShare
            dt = step * samplingInterval;
            return
        end
    end
    error('resolvedInterval:none', 'No interval reaches a zero-step share below %g.', maxZeroShare);
end
