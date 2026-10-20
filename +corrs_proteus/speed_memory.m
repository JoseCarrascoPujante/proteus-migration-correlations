function ac1 = speed_memory(x, y, dt)
%speed_memory  Lag-1 speed autocorrelation AC_v(1) of a trajectory sampled every dt s
%   (Shaebani et al. 2022, Biophys J 121:4099, doi:10.1016/j.bpj.2022.09.033).
    arguments
        x (:,1) double {mustBeReal}
        y (:,1) double {mustBeReal}
        dt (1,1) double {mustBePositive}
    end

    step = round(dt/0.5);   % positions are recorded every 0.5 s
    xv = x(1:step:end);
    yv = y(1:step:end);
    speed = hypot(diff(xv), diff(yv));
    speedDev = speed - mean(speed);
    n = numel(speedDev);
    ac1 = mean(speedDev(1:n-1).*speedDev(2:n))/var(speed, 1);
end
