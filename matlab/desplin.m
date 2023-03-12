function [T,Sd,Sv,Sa]=desplin(as,Tmin,Tmax,DT,vxi,dt)
%
%----------------------------------------------------- 
% [T,Sd,Sv,Sa]=desplin(as,Tmin,Tmax,DT,vxi,dt)
%----------------------------------------------------- 
%
% Calcula los expectros de respuesta de un sistema 
% sencillo lineal por el metodo de la aceleracion lineal.
%
%
% Por: Jorge E. Hurtado G.
%      Universidad Nacional de Colombia
%
%
% as:    vector (columna) de la aceleracion del suelo
% Tmin:  periodo minimo de calculo
% Tmax:  periodo maximo de calculo
% DT:    incremento del periodo
% vxi:   vector que contiene las fracciones de amortiguamiento 
%        viscoso para las cuales se han de calcular los espectros
% dt:    paso de tiempo del acelerograma
%
% Sd: espectro de desplazamiento
% Sv: espectro de velocidad 
% Sa: espectro de aceleracion 
%
% load vinadelmar.txt;
% [T,Sd,Sv,Sa]=desplin(vinadelmar,0.05,4,0.05,[0.02 0.05],0.005);
% load elcentro.txt;
% [T,Sd,Sv,Sa]=desplin(elcentro,0.05,4,0.05,[0.02 0.05],0.02);
% load mexico.txt
% [T,Sd,Sv,Sa]=desplin(mexico,0.05,4,0.05,[0.02 0.05],0.02);%----------------------------------------------------- 
%
%
% as=input('Acelerograma = ');
% Tmin=input('Periodo minimo = ');
% Tmax=input('Periodo maximo = ');
% DT=input('Incremento de periodo = ');
% vxi=input('Vector de amortiguamientos = ');
% dt=input('Incremento de tiempo = ');
%
l=length(vxi);
m=(Tmax-Tmin)/DT+1;
%
T=linspace(Tmin,Tmax,m)';
W=2*pi./T;
%
for i=1:l
  xi=vxi(i);
  for j=1:m
    w=W(j);
    [d,v,a]=dmaclin1(-as,1,w,xi,dt);
    Sd(i,j)=max(abs(d));
    Sv(i,j)=max(abs(v));
    Sa(i,j)=max(abs(as+a));
   end
end
%
t=linspace(1,dt*length(as),length(as));
figure
plot(t,as)
xlabel('Tiempo')
ylabel('Aceleración del suelo')
figure
plot(T,Sd)
xlabel('Periodo')
ylabel('Espectro de desplazamiento');
legend('\xi=0.02','\xi=0.05')
figure
plot(T,Sv)
xlabel('Periodo')
ylabel('Espectro de velocidad');
legend('\xi=0.02','\xi=0.05')
figure
plot(T,Sa)
xlabel('Periodo')
ylabel('Espectro de aceleracion');
legend('\xi=0.02','\xi=0.05')
%
%
%----------------------------------------------------- fin

