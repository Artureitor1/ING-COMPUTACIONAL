function [maxD] = dampedForcedVibrationMaxDSteady(d0,dampingFactor,MODES,Freq, DOFl,Force,Force_w, t)
d = zeros(length(d0(DOFl)),1);
    for mode = 1:length(Freq)
        f_i = MODES(:,mode)'*Force;
        rho_i = Force_w / Freq(mode);
        q_max = (f_i / (Freq(mode)^2))/((1-rho_i^2)^2 + (2*dampingFactor*rho_i)^2);
        d = d+ MODES(:,mode)*q_max * ...
            ((1-rho_i^2)*sin(Force_w*t) - 2*dampingFactor*rho_i*cos(Force_w*t));
    end
maxD = max(abs(d));
end
