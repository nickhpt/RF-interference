function [pixH,pixV,pixST3,pixST4,pixKV,pixKH,INR_H,INR_V,noise_x,RFI_x] = radiometerPix(NFFT,WFFT,AFFT,freq,toggleH,toggleV,RFI,stokes,duty,bw)

%Inputs
%Freq:   Vector of relatives frequencies (0-0.5)for arbitary #temporal sub-samples
%toggle: len(toggle)=len(Freq), 1:keeps RFI, 0:Deletes RFI in a sub-sample
%RFI:    String input 'SIN' for sinDuty otherwise QPSK 
%stokes: 2-element vector, with 3 and 4 stokes for RFI and noise
%Outputs
%Radiometer pixel for different products: 
%[powerH,PowerV,stokes3,stokes4,kurtosisV,kurtosisH,noiseH,noiseY,RFIH,RFIV]

%Number of temporal subsamples 
N_band = length(freq);
frequency = freq;
%h = waitbar(0,'Calculating...');
for i = 1:N_band
    type = RFI(i);
   % waitbar(i/N_band,h);
    if type == 1
      [SX, SY] = CimrStokesSinDuty(frequency(i), duty(i), NFFT, WFFT, AFFT, stokes(i,1), stokes(i,2));  
      
    elseif type == 2
       [SX, SY] = CimrStokesQpskBw(bw(i), frequency(i), 1, NFFT, WFFT, AFFT, stokes(i,1), stokes(i,2), 'QPSK sig');
    end
    
    %Add noise to RFI 
    pd = makedist('Normal'); %zero mean, 1 std. 
    [nx, ny] = CimrStokesNoise(pd, NFFT, WFFT, AFFT, 0, 0);
    %Toggle RFI on/off and control intensity
    SX_dat=(SX.*toggleH(i)+nx).'; 
    SY_dat=(SY.*toggleV(i)+ny).'; 
    %Apply polyphase filter to RFI + signal
    [FFT_Nx,FFT_Ny] = polyPhase(SX_dat,SY_dat,NFFT,AFFT,WFFT);
    %Apply polyphase filter to only signal
    [FFT_sigx,FFT_sigy] = polyPhase(nx.',ny.',NFFT,AFFT,WFFT);
    %Apply polyphase filter to only RFI
    [FFT_nox,FFT_noy] = polyPhase((SX.*toggleH(i)).',(SY.*toggleV(i)).',NFFT,AFFT,WFFT);
    
    %Get power products
    [P_V,P_H,P3,P4,K_V,K_H] = getProducts(FFT_Nx,FFT_Ny);
    [P_Vs,P_Hs,~,~,~,~] = getProducts(FFT_sigx,FFT_sigy);
    [P_Vn,P_Hn,~,~,~,~] = getProducts(FFT_nox,FFT_noy);
    
    %Accumulate product, dim = NFFT, and extract only positive    
    %H-pol (,power,kurtosis)
    P1H = sum(P_H); %Gaussian noise + RFI
    K_H = (AFFT.*sum(K_H))./(sum(P_H).^2); %Kurtosis
    PHs= sum(P_Hs); %Gaussian noise
    PHn= sum(P_Hn); %RFI
    
    ST3 = sum(P3);  %Third Stokes
    
    %Get positive frequencies
    ST3 = ST3(NFFT/2+1:end-1);
    P1H = P1H(NFFT/2+1:end-1);
    PHs= PHs(NFFT/2+1:end-1);
    PHn= PHn(NFFT/2+1:end-1);
    K_H = K_H(NFFT/2+2:end);
      
    %V-pol
    P1V = sum(P_V);
    K_V = (AFFT.*sum(K_V))./(sum(P_V).^2);
    PVs = sum(P_Vs);
    PVn = sum(P_Vn);
    
    ST4 = sum(P4); %Fourth Stokes
    
    %Get positive frequencies
    ST4 = ST4(NFFT/2+1:end-1);
    P1V = P1V(NFFT/2+1:end-1);
    PVs = PVs(NFFT/2+1:end-1);
    PVn = PVn(NFFT/2+1:end-1);
    K_V = K_V(NFFT/2+2:end);
    
    %Generate i'th temporal sub-sample
    %Power
    pixH(:,i) = P1H;   
    pixV(:,i) = P1V;    
    %Stokes
    pixST3(:,i) = ST3;
    pixST4(:,i) = ST4;    
    %Kurtosis
    pixKV(:,i) = K_V;
    pixKH(:,i) = K_H;
    
    %Natural signal that went into the system
    noise_x = PHs;
    noise_y = PVs;    
    %RFI that went into the system
    RFI_x = PHn;
    RFI_y = PVn;
    
    %Input INR at detector
    INR_H(i) = 10*log10(sum(RFI_x)/sum(noise_x)) ;
    INR_V(i) = 10*log10(sum(RFI_y)/sum(noise_y)) ;

%     %If Peak INR is needed, then uncomment
%     INR_H(i) = 10*log10(max(RFI_x)/mean(noise_x));
%     INR_V(i) = 10*log10(max(RFI_y)/mean(noise_y));
       
   
end
%close(h)



end

