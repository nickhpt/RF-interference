function [P_acc2,P_detect] = anomalousAmp(Pow,tau)

%Function should be called on output from radioPix.m
%Inputs:
%Pow:  Power values in spectral sub-bands
%tau:  Threshold determined from Monte-Carlo simulation
%output:
%P_acc2: Output data after flagging RFI
%P_detect: Power of flagged RFI

%Determining spectral sub-band with RFI
%NFFT = size(Pow,1)*2; 
P_acc = Pow;
P_acc2 = P_acc;

y_gauss = mean(P_acc)*tau; %Threshold
idx = find(P_acc>y_gauss);
%Output data after flagging RFI
P_acc2(idx) = NaN;
%Power of flagged RFI
P_detect = sum(P_acc(idx));


end

