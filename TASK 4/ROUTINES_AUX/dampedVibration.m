function [d] = dampedVibration(d0,dampingFactor,MODES,Freq,M, DOFl,t)
d = zeros(length(d0(DOFl)),1);
    for mode = 1:length(Freq)
        wBar = Freq(mode)*sqrt(1-dampingFactor^2);
        q = MODES(:,mode)'*M(DOFl,DOFl)*d0(DOFl);
        d = d+MODES(:,mode)*(exp(-dampingFactor*Freq(mode)*t)*(q*cos(wBar*t)+(q+dampingFactor*q*sin(wBar*t))/(wBar)));
    end 
end