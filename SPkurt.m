function [K2,P_RFI] = SPkurt(Kurt,P,tau)

%Inputs:
%Kurt: Spectral Kurtosis H or V calculated from getProducts.m
%P: Second order product from getProducts.m
%tau:  Threshold determined from MC - simulation
%outputs:
%K2: Kurtosis value in spectral sub-bands after removing RFI
%P_RFI: Power of detected RFI

K2 = Kurt;

%Spectral Kurtosis value for Gaussian signal
mu = 2;
%Determine upper bound for maximum deviation of SP-kurtosis
idx_up = find(Kurt > mu*tau);
%Detmine lower bound and ensure that fractional part is properly determined
if tau > 1
    idx_lo = find(Kurt < mu-(mu*tau-floor(mu*tau)));
else
    idx_lo = find(Kurt < mu-(mu*tau-ceil(mu*tau)));
end

%Power of potentially detected RFI
P_RFIu = sum(P(idx_up));
P_RFIl = sum(P(idx_lo));

%Only consider either lower or upper bound condition as flag for RFI
if (P_RFIl < P_RFIu)
    P_RFI = P_RFIu;
    K2(idx_up) = NaN;
else
    P_RFI = P_RFIl;
    K2(idx_lo) = NaN;
end


end

