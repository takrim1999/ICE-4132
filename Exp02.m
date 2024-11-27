
clc;   clear;   close all;  

%% To generate an AM waveform, we can use the following MATLAB code:

Fs = 256*1e3        % Sampling frequency 40KHz (samples per second) 
dt = 1/Fs;         % seconds per sample 
N = 256;
 
 %% For one cycle get time period
 Fc = 50000;         % Sine wave frequency (hertz) of carrier signal
 Fm = 2000;         % Sine wave frequency (hertz) of modulating signal
 xc = [];           xm = [];    
 wc = 2*pi*Fc;      wm = 2*pi*Fm;   
 % time step for one time period 
 % Carrier Signal
 Ac = 4
 StopTime = 0.005; % seconds 
 tt = 0:dt:StopTime; % seconds
 xc = sin(2*pi*Fc*tt);
 figure(1)
 plot(tt, Ac*xc,'k--o');    grid on;
 title('Carrier Signal')
xlabel(' time s ');
ylabel( ' amplitude ' );

% Modulating Signal
Am = 3
 StopTime = 0.005; % seconds 
 t = 0:dt:StopTime; % seconds
 xm = sin(2*pi*Fm*t);
 figure(2)
 plot(t, Am*xm,'b--o');    grid on;
 title('Modulating Signal')
xlabel(' time s ');
ylabel( ' amplitude ' );


%% DFT using Matlab built-in function: 
% m = length(xc)
% dft_matlab_xc = xc(1,1:N)*dftmtx(N);
% F_analysis = Fs/N  
% 
% dft_matlab_xm = xm(1,1:N)*dftmtx(N);
% F_analysis = Fs/N 
% 
% Xe_mag = [];     Xe_real = [];    Xe_imag = [];
% for ii = 1:N
%    
%     Xe_real(ii,1) = real(dft_matlab_xc(1,ii));
%             if Xe_real(ii,1) > 0 && Xe_real(ii,1) < 1e-10
%                Xe_real(ii,1) = 0;  
%             end  
%             if Xe_real(ii,1) < 0 && Xe_real(ii,1) > -1e-10
%                Xe_real(ii,1) = 0;  
%             end
%         Xe_imag(ii,1) = imag(dft_matlab_xc(1,ii));
%             if Xe_imag(ii,1) > 0 && Xe_imag(ii,1) < 1e-10
%                Xe_imag(ii,1) = 0;  
%             end  
%             if Xe_imag(ii,1) < 0 && Xe_imag(ii,1) > -1e-10
%                Xe_imag(ii,1) = 0;  
%             end
%         Xe_mag_xc(ii,1) = sqrt(Xe_real(ii,1).^2 + Xe_imag(ii,1).^2);
%  
% end
%  
% mf = 0:N-1;
% figure(3); 
% stem(mf,Xe_mag_xc,'LineStyle','--',...
%      'MarkerSize',15,'Marker','s',...
%      'MarkerFaceColor','black',...
%      'MarkerEdgeColor','green')
% grid on;
% title('Magnitude of Xc_exp(m)')
% xlabel('m (KHz)')
% ylabel('Magnitude')
%  
% Xe_mag = [];     Xe_real = [];    Xe_imag = [];
% for ii = 1:N
%    
%     Xe_real(ii,1) = real(dft_matlab_xm(1,ii));
%             if Xe_real(ii,1) > 0 && Xe_real(ii,1) < 1e-10
%                Xe_real(ii,1) = 0;  
%             end  
%             if Xe_real(ii,1) < 0 && Xe_real(ii,1) > -1e-10
%                Xe_real(ii,1) = 0;  
%             end
%         Xe_imag(ii,1) = imag(dft_matlab_xm(1,ii));
%             if Xe_imag(ii,1) > 0 && Xe_imag(ii,1) < 1e-10
%                Xe_imag(ii,1) = 0;  
%             end  
%             if Xe_imag(ii,1) < 0 && Xe_imag(ii,1) > -1e-10
%                Xe_imag(ii,1) = 0;  
%             end
%         Xe_mag_xm(ii,1) = sqrt(Xe_real(ii,1).^2 + Xe_imag(ii,1).^2);
%  
% end
%  
% mf = 0:N-1;
% figure(4); 
% stem(mf,Xe_mag_xm,'LineStyle','--',...
%      'MarkerSize',15,'Marker','s',...
%      'MarkerFaceColor','black',...
%      'MarkerEdgeColor','green')
% grid on;
% title('Magnitude of Xm_exp(m)')
% xlabel('m (KHz)')
% ylabel('Magnitude')

%% AM generation
xc_len = length(xc);
mod_xm = xm(1,1:xc_len);
Modulation_index = Am/Ac

x_AM = Am*mod_xm.*xc + Ac*xc;

figure(5)
plot(tt, x_AM);    grid on;
title('Amplitude Modulated Signal')
xlabel(' time s ');
ylabel( ' amplitude ' ); 
 
%% AM demodulation (envelop detection )

xam_Dmod = Am*Ac*mod_xm + (((Modulation_index*Ac)/2)^2)*sin(2*wm*tt);

figure(6)
plot(tt, xam_Dmod, 'b--o' );    grid on;
xlabel(' time s ');
ylabel( ' amplitude ' );















