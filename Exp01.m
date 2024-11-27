clc;
close all;
clear all;
m=1;
Am=5; 
fa=2000; 
Ta=1/fa;
t=0:Ta/999:6*Ta;
ym=Am*sin(2*pi*fa*t);
figure (1)
subplot (3,1,1);
plot(t,ym)
title ('Modulating signal')

Ac=Am/m;
fc=fa*10;
Tc=1/fc;
yc=Ac*sin(2*pi*fc*t);
subplot (3,1,2)
plot(t, yc)
grid on;
title(' Carrier signal');
%AM Modulation
y=Ac+(1+m*sin(2*pi*fa*t)).*sin(2*pi*fc*t);
subplot (3,1,3)
plot (t,y)
title('Amplitude Modulated Signal')
grid on