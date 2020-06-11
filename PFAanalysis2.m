close all; clear all; clc
%Threshold setting
%thresh1 = (qfuncinv(PF)./sqrt(NFFT))+ 1; % Theroretical value of threshold
%Input probability of false alarm
%PF = ([0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);
PF = 0.01:0.01:1; 
NFFT = 512; %Number of inputs to each FFT 
WFFT = 3; %Number of FFT lengths used for weighting
AFFT = 128; %Number of FFTs used in the simulation
steps = linspace(0.4,1.4,1e3);
done = false;   
incr = 0;
tic
h = waitbar(0,'Calculating');
for k = 1:length(PF)
    waitbar(k/length(PF),h);
    FAR = PF(k);  
    %Tolerated number of flags
    tol_flag = ceil(NFFT*FAR);
    
    for j = 1:length(steps)
        step = steps(j);
        
        for i = 1:100 
            %Make noise
            pd = makedist('Normal'); %zero mean, 1 std. s
            [nx, ny] = CimrStokesNoise(pd, NFFT, WFFT, AFFT, 0, 0);
            %Filter only signal
            [FFT_sigx,FFT_sigy] = polyPhase(nx.',ny.',NFFT,AFFT,WFFT);
            %Get products
            [P_Vs,~,~,~,~,~] = getProducts(FFT_sigx,FFT_sigy);  
            %Accumulated samples
            P_accs = sum(P_Vs); 
            
            %mu_N = 1+step;
            y_gauss = mean(P_accs)*step;
            idx(i) = numel(find(P_accs >y_gauss));     
        end
        
        cond = find(idx > tol_flag);      
        if isempty(cond)
            done = true;          
            thresh = step;
            break;
        end
            
    end
    
    if done
    tau(k) = thresh;
    done = false;
    end
    k
end
toc
close(h);