function [d] = dampedForcedVibration(d0,dampingFactor,MODES,Freq, DOFl,Force,Force_w, t)
d = zeros(length(d0(DOFl)),1);
    for mode = 1:length(Freq)
        wBar = Freq(mode)*sqrt(1-dampingFactor^2);
        f_i = MODES(:,mode)'*Force;
        rho_i = Force_w / Freq(mode);
        q_max = (f_i / (Freq(mode)^2))/((1-rho_i^2)^2 + (2*dampingFactor*rho_i)^2);
        theta_i = atan(2*dampingFactor*rho_i / (1-rho_i^2));
        a_i = dampingFactor*sin(theta_i) - rho_i*cos(theta_i);
        b_i = sin(theta_i);
        d = d+ MODES(:,mode)*q_max * ...
            (...
            (exp(-dampingFactor*Freq(mode)*t) * (a_i * cos(wBar*t) + b_i*sin(wBar*t)))...
            + ...
            ((1-rho_i^2)*sin(Force_w*t) - 2*dampingFactor*rho_i*cos(Force_w*t))...
            );
    end 
end