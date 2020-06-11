close all; clear all; clc

%RFI Characterstics section

NFFT = 512; %Number of inputs to each FFT 
WFFT = 3;   %Number of FFT lengths used for weighting
AFFT = 128; %Number of FFTs used in the simulation

%Frequencies where RFI is seen
frequencies = [0.25 0.15 0.35 0.1];

%Control H/V intensities
toggleV  = [1 1 1 1];            
toggleH = [1 1 1 1];    
%Duty for sin
duty = [0 0.001 0.1 0.5];    
%RFI type (SIN,QPSK) = (1,2)
type = [2 1 1 1];    
%Stokes 3,4
stokes = [1 0; 1 0; -1 -1; 1 1];  
%Bandwidth of QPSK RFI
bw = [0.4 0.1 0.1 0.1];        

%Create radiometer-pixel. Number of temporal sub-samples in this case is 4
[H,V,S3,S4,KV,KH,INRH,INRV] = radiometerPix(NFFT,WFFT,AFFT,frequencies,toggleH,toggleV,type,stokes,duty,bw);
%Define frequency vector:
f = (-NFFT/2:NFFT/2-1)/NFFT;f=f(NFFT/2+1:end-1);

%Plot Power in the first temporal sub-sample
plot(f,10*log10(H(:,1)/max(H(:,1))));
xlabel('Relative Frequency');ylabel('Normalized power [dB]');grid on;

%%
%Generation of filter response to sinusoidal input 
%Update: After code modification natural noise is always turned on
%Noise must be turned off manually inside
%radimeterPix.m to show filter response

NFFT = 32;%Number of inputs to each FFT 
WFFT = 3; %Number of FFT lengths used for weighting
AFFT = 128; %Number of FFTs used in the simulation
loop_frq = 0.001:0.001:0.499;
Temp_P = zeros(length(loop_frq),NFFT/2-1);
incr = 0;
for frq = loop_frq
[X,Y,~,~,~,~,~,~] = radiometerPix(NFFT,WFFT,AFFT,frq,0,0,1,[1 1],1,1);
incr = incr +1;
Temp_P(incr,:) = Y; 
end

for k = 1:NFFT/2-1
    plot(loop_frq,10*log10(Temp_P(:,k)/max(Temp_P(:,k))),'LineWidth',2);   
    hold on
end
grid on
xlabel('Relative frequency');
ylabel('Normalized power [dB]')

%%
close all; clear all; clc
%One example of result generation for figure 19. 
thresholds2 = load('thresholdswgn.mat');
tau2 = thresholds2.tau;
tauval = 1; %Set FAR
%Define frequency vector:
f = (-512/2:512/2-1)/512;f=f(512/2+1:end-1);

for i = 1:100
rng(1) %Set random seed, such that same noise is used
%Create radiometer pixel (1-temporal sub-sample)
%[H,~,S3,~,~,KH,INRH,~,Nx,RFIx] = radiometerPix(512,3,128/4,0.35,3/sqrt(2*i),3/sqrt(2*i),1,[1 0],0.1,0.1);
%[H,~,S3,~,KV,KH,INRH,~,Nx,RFIx] = radiometerPix(512,3,128/4,0.35,1/sqrt(i),1/sqrt(i),1,[1 0],0.001,0.1);
[H,~,S3,~,KV,KH,INRH,~,Nx,RFIx] = radiometerPix(512,3,128/4,0.25,20/sqrt(i),20/sqrt(i),1,[1 0],0.001,0.1);
plot(f,H(:,1))
hold on
grid on
%Run algorithms 
[P_outcf,P_RFIcf] = CFPA(H(:,1),tau2(tauval),0);
[P_outsp,P_RFISP] = SPkurt(KH(:,1),H(:,1),tau2(tauval));
[P_outamp,P_RFIamp] = anomalousAmp(H(:,1),tau2(tauval));
[P_outpa,P_RFIpa] = PA(S3(:,1),H(:,1),tau2(tauval));

Pdetcf(i) = (P_RFIcf); %CFPA
Pdetspk(i) = P_RFISP;%SPK
Pdetamp(i) = P_RFIamp;%A-amp
Pdetpa(i) = P_RFIpa;%SCFA
INR1(i) = INRH; %Before detection SNR at detector 
i
end
%xlabel('Relative frequency');ylabel('Kurtosis H-pol');

figure
plot((Pdetcf/max(Pdetcf)),INR1,'--r','LineWidth',1.4)
hold on
plot((Pdetspk/max(Pdetspk)),INR1,'--g','LineWidth',1.4)
hold on
plot((Pdetamp/max(Pdetamp)),INR1,'--b','LineWidth',1.4)
hold on
plot((Pdetpa/max(Pdetpa)),INR1,'--y','LineWidth',1.4)
xlabel('Normalized detected RFI power');ylabel('INR at detector [dB]');
legend('CFPA','SPK','A-amp','PA');
grid on
