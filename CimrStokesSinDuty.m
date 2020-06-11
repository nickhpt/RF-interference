function [SX, SY] = CimrStokesSinDuty(SinFrq, Duty, NFFT, WFFT, AFFT, S3, S4)
%
% Sinusoidal input generating function 
%
% Input parameters:
%    SinFrq: Frequency (relative 0.0 - 0.5)
%    Duty: Duty cycle (0.0 - 1.0)
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
%
% No RFI if zero duty cycle
%
if (Duty == 0.0)
  SX = zeros(AFFT+WFFT-1, NFFT);
  SY = zeros(AFFT+WFFT-1, NFFT);
  return
end
%
% Duty cycle
%
TIdx = 0:(AFFT+WFFT-1)*NFFT-1;
Idx = (TIdx <= length(TIdx)*Duty);
%
% Single polarization
%
AX = cos(2*pi*SinFrq*TIdx);
AX = AX .* Idx;
Pw = sum(dot(AX,AX))/numel(AX);
AX = AX .* (1/sqrt(Pw));
AX = circshift(AX, (WFFT+1)*NFFT);
AX = reshape(AX, NFFT, AFFT+WFFT-1)';
%
% No 3rd or 4th
%
if (S3 == 0.0) && (S4 == 0.0)
  SX = reshape(AX, NFFT, AFFT+WFFT-1)';
  SY = zeros(AFFT+WFFT-1, NFFT);
  return
end
%
% 3rd Stokes
%
if (S3 ~= 0.0)
  if (S3 > 1.0)
    S3 = 0.5;
  elseif (S3 < -1.0)
    S3 = -0.5;
  else
    S3 = S3/2;
  end
  SX = sqrt((1-abs(S3)))*AX; %1-abs
  SY = sign(S3)*sqrt(abs(S3))*AX;
  return;
end
  
if (S4 ~= 0.0)
  if (S4 > 1.0)
    S4 = 0.5;
  elseif (S4 < -1.0)
    S4 = -0.5;
  else
    S4 = S4/2;
  end
  SX = sqrt((1-abs(S4)))*AX;
  SY = sign(S4)*sqrt(abs(S4))*circshift(AX, +1);
  return;
end

SX = zeros(AFFT+WFFT-1, NFFT);
SY = zeros(AFFT+WFFT-1, NFFT);
