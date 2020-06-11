function [QX, QY] = CimrStokesQpskBw(Bw, FrqQpsk, FrqSamp, NFFT, WFFT, AFFT, S3, S4, CimrFileStr)
%
% QPSK input generating function 
%
% Input parameters:
%    FrqQpsk: Frequency (relative 0.0 - 0.5)
%    Bw: Band width ( relative0.0 - 1.0)
%    NFFT: Number of inputs to each FFT
%    WFFT: Number of FFT lengths used for weighting
%    AFFT: Number of FFTs used in the simulation
%    S3: Third Stokes parameter value
%    S4: Fourth Stokes parameter value
%
% Output parameters
%    SX: X-polarization output matrix of size AFFT+WFFT-1 times NFFT
%    SY: Y-polarization output matrix of size AFFT+WFFT-1 times NFFT
%
% No RFI if zero bandwidth
%
if (Bw == 0.0)
  QX = zeros(AFFT+WFFT-1, NFFT);
  QY = zeros(AFFT+WFFT-1, NFFT);
  return
end
%
% RRC parameters (Root Raised Cosine)
%
RRC_RollOff = 0.2; % 20% bandwidth widening
RRC_Symbols = 10;
RRC_Shape = 'Square root';
PSK_NPhase = 4;
PSK_PhaseOff = pi/4;
%PSK_SymbLng  = round(FrqSamp/Bw/(1+RRC_RollOff));
PSK_SymbLng  = round(FrqSamp/Bw);
Samples = (AFFT+WFFT-1)*NFFT;
%
% PSK modulation parameters
%
SamplesExtra = Samples+2*RRC_Symbols*PSK_SymbLng;
PSK_NSymb = uint32(ceil(SamplesExtra/PSK_SymbLng)+2*RRC_Symbols); % Number of symbols in signal
PSK_NBits = PSK_NSymb*uint32(log2(PSK_NPhase)); % Number of bits to generate N-PSK 
%
% PSK modulation at baseband
%
PSK_BitData = randi([0 1],PSK_NBits,1); % Create binary data for random modulation
hPSKModulator = comm.PSKModulator(PSK_NPhase, 'BitInput', true); % Create a N-PSK modulator System object with bits as inputs and Gray-coded signal constellation
hPSKModulator.PhaseOffset = PSK_PhaseOff; % Change the phase offset
PSK_IQSignal = step(hPSKModulator, PSK_BitData); % Generate complex modulation signal at baseband
%
% RRC (Square Root Raised Cosine) weighting
%
RRC_Gain = sqrt(double(1.0*PSK_SymbLng));
RRC_Output = 1.0*PSK_SymbLng;
hRRCFilter = comm.RaisedCosineTransmitFilter('Shape', RRC_Shape, 'RolloffFactor',RRC_RollOff, ...
  'FilterSpanInSymbols',RRC_Symbols,'OutputSamplesPerSymbol',RRC_Output, 'Gain', RRC_Gain);
%b =coeffs(hRRCFilter);
PSK_Encode = step(hRRCFilter, PSK_IQSignal);

%
% Multiply by carrier
%
F = 1:Samples;
F = exp(1i * 2 * pi * FrqQpsk * F / FrqSamp);
F = F .* PSK_Encode(1+RRC_Symbols*PSK_SymbLng:Samples+RRC_Symbols*PSK_SymbLng)'; % Avoid edge effects
AX = real(F); % Real signal
%
% Scale RFI to power 1.0
%
Pw = dot(AX,AX)/Samples;
AX = AX .* (1/sqrt(Pw));
%
QX = AX;
QY = zeros(size(AX));

if (S3 ~= 0.0)
  if (S3 > 1.0)
    S3 = 0.5;
  elseif (S3 < -1.0)
    S3 = -0.5;
  else
    S3 = S3/2;
  end
  QX = sqrt((1-abs(S3)))*AX;
  QY = sign(S3)*sqrt(abs(S3))*AX;
elseif (S4 ~= 0.0)
  if (S4 > 1.0)
    S4 = 0.5;
  elseif (S4 < -1.0)
    S4 = -0.5;
  else
    S4 = S4/2;
  end
  QX = sqrt((1-abs(S4)))*AX;
  QY = sign(S4)*sqrt(abs(S4))*circshift(AX, +1);
end

%
%CimrMisc;
SavePlot = 0;
%
% Plot QPSK signal
%
clf;
Pdir = 'QpskBwPlot';
frq = (1:Samples/2)*FrqSamp/Samples;
QF = abs(fft(QX));
QS = 20*log10(QF);
QS(1:length(frq)) = QS(1:length(frq)) - max(QS(1:length(frq)));

%added following line 112
%And supressed CimrMisc. line 98 & CimrSave Line 124 + added close all in end
LineWidth = 1; FontSize = 1; 

%Also commented drawnow (line 130)
% plot(frq, QS(1:length(frq)), 'r',  'LineWidth', LineWidth);
% axis([0.0 FrqSamp/2 -100 0]);
% Str = sprintf('QPSK Spectrum');
% title(Str);
% ylabel('Amplitude [dB]', 'FontSize', FontSize);
if (FrqQpsk > 50) % Not L-band
  Str = sprintf('(BW: %3.0f MHz, Frq = %3.0f MHz)         Frequency [MHz]', Bw, FrqQpsk);
else
  Str = sprintf('(BW: %3.1f MHz, Frq = %3.1f MHz)         Frequency [MHz]', Bw, FrqQpsk);
end
xlabel(Str, 'FontSize', FontSize);
set(gca,'FontSize',FontSize);
grid;
SaveStr = sprintf('%s/%s_QPSK_FRQ%1.0f_BW%1.0f_title.png', Pdir, CimrFileStr, FrqQpsk, Bw);
%CimrSave;
%drawnow;
%
% Reshape
%
QX = reshape(QX(1:Samples), NFFT, AFFT+WFFT-1)';
QY = reshape(QY(1:Samples), NFFT, AFFT+WFFT-1)';

close all