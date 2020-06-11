function [P_spsb2,P_RFI] = CFPA(dat,tau,plott)
%Delete the output i
%Input: X or H, Stokes3, threshold, plot
%Cross frequency peak algorithm
%Input: Radiometer pixel, outputs: 1)data and RFI, 2)RFI removed

data = dat;
P_spsb2 = data;

cnt = 1;
for i = 1:1e10 %Some high number
    
    mu_spsb = nanmean(P_spsb2);
    idx = find(P_spsb2 == max(P_spsb2));
    
         if P_spsb2(idx) > mu_spsb*tau
           % clf  
            P_spsb2(idx) = NaN;
            RFI_spsb(cnt) = idx;
            P_RFI(cnt) = dat(idx); %Count RFI from second order products
            cnt = cnt+1;    
            
        else
            break; %Break when no more RFI
        end

end

if(exist('P_RFI'))
P_RFI = sum(P_RFI);
else
    P_RFI = 0;
end

%Out:
%P_spsb2: Output data after RFI flagging
%P_RFI: Power of flagged RFI

if plott==1

norm = max(P_spsb2);
P_spsb2 = P_spsb2/norm;
%subplot(2,1,2)
bar(P_spsb2);
hold on
plot((zeros(1,length(P_spsb2))+mu_spsb*tau)/norm,'r')
xlabel('Spectral sub-band');ylabel('Normalized Power');title('Post detection')  

end
end

