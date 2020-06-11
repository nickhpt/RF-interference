function [P_V,P_H,P3,P4,K_V,K_H] = getProducts(FFT_Nx,FFT_Ny)

%Inputs
%FFT_Nx: FFT based filtered H-pol
%FFT_Ny: FFT based filtered V-pol
%Outputs
%Power V-H , 3,4 stokes parameter, Kurtosis V-H

%Calculate second and fourth order products
P_V = (real(FFT_Ny).*real(FFT_Ny)+imag(FFT_Ny).*imag(FFT_Ny));
P_H = (real(FFT_Nx).*real(FFT_Nx)+imag(FFT_Nx).*imag(FFT_Nx));
P3 = (real(FFT_Ny).*real(FFT_Nx)+imag(FFT_Ny).*imag(FFT_Nx));
P4 = (real(FFT_Ny).*imag(FFT_Nx)-imag(FFT_Ny).*real(FFT_Nx));
K_V = (P_V.*P_V);
K_H = (P_H.*P_H);


end

