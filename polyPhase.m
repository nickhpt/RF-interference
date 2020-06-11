function [H_filt,V_filt] = polyPhase(datH,datV,NFFT,AFFT,WFFT)

%Input: H/V - pol data must be of dim([NFFT]x[AFFT+WFFT-1])
SX = datH;
SY = datV;

M=NFFT*WFFT;
h=zeros(1,M);ham=h;
for i = 1:M   
   h(i) = sinc(0.95*(i-M/2)/(NFFT)); %Multiply with factor to play with mainlobe width
   ham(i) = 0.54-0.46*cos(2*pi*i/M);   
end
h = h.*ham;

%Create NFFT coefficients
N_coeff = h.';

%Serializing data
SX_s = SX(:);
SY_s = SY(:);

%Perform polyphase structure
incr = 0;
for i = 1:AFFT
    N_wx = N_coeff.*SX_s(1+incr*NFFT:(incr+WFFT)*NFFT);
    N_wy = N_coeff.*SY_s(1+incr*NFFT:(incr+WFFT)*NFFT);
    
    N_wx=reshape(N_wx,NFFT,WFFT);
    N_wy=reshape(N_wy,NFFT,WFFT);
    
    N_sumx = sum(N_wx,2);
    N_sumy = sum(N_wy,2); 
    
    FFT_Nx(i,:) = fftshift(fft(N_sumx)/(NFFT));
    FFT_Ny(i,:) = fftshift(fft(N_sumy)/(NFFT));
    incr = incr+1;
  
end

H_filt = FFT_Nx;
V_filt = FFT_Ny;

end

