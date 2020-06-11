function [NX, NY] = CimrStokesNoise(Noise_pd, NFFT, WFFT, AFFT, S3, S4)
%
% Noise input generating function
%
% Input parameters:
%    Noise_pd: Distribution of the random noise (normal)
%    NFFT: Number of inputs to each FFT
%    WFFT: Number of FFT lengths used for weighting
%    AFFT: Number of FFTs used in the simulation
%    S3: Third Stokes parameter value
%    S4: Fourth Stokes parameter value
%
% Output parameters
%    NX: X-polarization output matrix of size AFFT+WFFT-1 times NFFT
%    NY: Y-polarization output matrix of size AFFT+WFFT-1 times NFFT
%
if (S3 == 0.0) && (S4 == 0.0)
  NX = random(Noise_pd, AFFT+WFFT-1, NFFT);
  NY = random(Noise_pd, AFFT+WFFT-1, NFFT);
else
  QX = random(Noise_pd, NFFT*(AFFT+WFFT-1), 1);
  QY = random(Noise_pd, NFFT*(AFFT+WFFT-1), 1);
  
  if (S3 ~= 0.0)
    if (S3 > 1.0)
      S3 = 0.5;
    elseif (S3 < -1.0)
      S3 = -0.5;
    else
      S3 = S3/2;
    end
    NX = sqrt((1-abs(S3)))*QX + sign(S3)*sqrt(abs(S3))*QY;
    NY = sqrt((1-abs(S3)))*QY + sign(S3)*sqrt(abs(S3))*QX;
  end
  
  if (S4 ~= 0.0)
    if (S3 ~= 0.0)
      disp('ERROR: SimrStokesNoise S3~= 0 and S4~=0 not implemented');
    else
      if (S4 > 1.0)
        S4 = 0.5;
      elseif (S4 < -1.0)
        S4 = -0.5;
      else
        S4 = S4/2;
      end
      NX = sqrt((1-abs(S4)))*QX + sign(S4)*sqrt(abs(S4))*circshift(QY, +1);
      NY = sqrt((1-abs(S4)))*QY - sign(S4)*sqrt(abs(S4))*circshift(QX, +1);
    end
  end
  NX = reshape(NX,  NFFT, AFFT+WFFT-1)';
  NY = reshape(NY,  NFFT, AFFT+WFFT-1)';
end