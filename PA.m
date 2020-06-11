function [S2,P_RFI] = PA(Stokes,P,tau)

%Inputs:
%Stokes: Stokes 3'rd or 4'th 
%P: Second order moment, power
%tau:  Threshold determined from MC - simulation
%outputs:
%S2: Stokes in spectral sub-bands after removing RFI
%P_RFI: Power of detected RFI

S2 = Stokes;
%3'rd and 4'th are always close to zero for natural signal
mu = 1e-3;
%Determine upper bound for maximum deviation of Stokes
idx_up = find(Stokes > mu*tau);
%Detmine lower bound and ensure that fractional part is properly determined
if tau > 1
    idx_lo = find(Stokes < mu-(mu*tau-floor(mu*tau)));
else
    idx_lo = find(Stokes < mu-(mu*tau-ceil(mu*tau)));
end

%Power of potentially detected RFI
P_RFIu = sum(P(idx_up));
P_RFIl = sum(P(idx_lo));

%Only consider either lower or upper bound condition as flag for RFI
if (P_RFIl < P_RFIu)
    P_RFI = P_RFIu;
    S2(idx_up) = NaN;
else
    P_RFI = P_RFIl;
    S2(idx_lo) = NaN;
end


end

