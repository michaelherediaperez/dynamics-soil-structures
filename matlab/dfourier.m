function [w,X,Xp,t,i]=dfourier(x,dt,wmax)
%----------------------------------------------------- 
% [w,X,Xp,t,i]=dfourier(x,dt,wmax)
%----------------------------------------------------- 
%
% Calcula la transformada de Fourier de una señal.
%
%
% Por: Jorge E. Hurtado G.
%      Universidad Nacional de Colombia
%
%
% x:     señal
% dt:    intervalo de discretización de la señal
% wmax:  frecuencia angular máxima para la representación gráfica
%
% w:     vector de frecuencias
% X:     Transformada de Fourier
% Xp:    módulo de la transformada
% i:     vector de índices con w <= wmax 
% t:     vector de tiempo 
%
%
%----------------------------------------------------- 
% 
% [w,X,Xp,i,t]=dfourier(elcentro,0.02,100)
% [w,X,Xp,i,t]=dfourier(mexico,0.02,30);
%
X=fft(x);
ws=2*pi/dt;          % sampling frequency
wn=ws/2;             % Nyquist frequency
w=linspace(0,wn,ceil(length(x)/2))';
Xp=abs(X(1:ceil(length(x)/2)));
i=find(w<=wmax);
t=linspace(1,dt*length(x),length(x));
whos
figure
plot(t,x)
xlabel('tiempo')
ylabel('aceleración')
figure
plot(w(i),Xp(i))
xlabel('frecuencia angular')
ylabel('módulo de la transformada')

