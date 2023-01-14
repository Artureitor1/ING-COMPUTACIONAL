function [d] = dampedForcedVibrationMaxDSteady(d0,dampingFactor,MODES,Freq, DOFl,Force,Force_w, t)
d = zeros(length(d0(DOFl)),1);
    for mode = 1:length(Freq)
        f_i = MODES(:,mode)'*Force;
        rho_i = Force_w / Freq(mode);
        q_max_i = (f_i / (Freq(mode)^2))/((1-rho_i^2)^2 + (2*dampingFactor*rho_i)^2);
        d = d+ MODES(:,mode)*q_max_i * sqrt((1-rho_i^2)^2 +  (2*dampingFactor*rho_i)^2);
    end
end
