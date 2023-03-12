function [t,d,v,a]=dmaclin(p,m,w,xi,dt)
%
%----------------------------------------------------- 
% [t,d,v,a]=dmaclin(p,m,w,xi,dt)
%----------------------------------------------------- 
%
% Calcula la respuesta de un sistema 
% sencillo lineal por el metodo de la aceleracion lineal.
%
%
% Por: Jorge E. Hurtado G.
%      Universidad Nacional de Colombia
%
%
% p: vector (columna) de carga externa
% m: masa del sistema
% w: frecuencia natural del sistema
% xi: fraccion de amortiguamiento viscoso
% dt: paso de tiempo
%
% t: vector de tiempo
% d: desplazamiento de respuesta
% v: velocidad de respuesta
% a: aceleracion de respuesta
%
% load tokachikenoki.txt
% [t,d,v,a]=dmaclin(-tokachikenoki,1,2*pi,0.05,0.02);
% load elcentro.txt
% [t,d,v,a]=dmaclin(-elcentro,1,2*pi,0.05,0.02);
% load mexico.txt
% [t,d,v,a]=dmaclin(-mexico,1,2*pi,0.05,0.02);
% load vinadelmar.txt
% [t,d,v,a]=dmaclin(-vinadelmar,1,2*pi,0.05,0.005);
%----------------------------------------------------- 
%
%
n=length(p);
tmax=dt*n;
t=linspace(0,tmax,n)';
d0=0;
v0=0;
a0=0;
%
k=m*w^2;
c=2*m*w*xi;
kbar=k+3*c/dt+6*m/(dt^2);
ikbar=1/kbar;
%
for i=1:n
  p1=p(i,:);
  dp=m*(6*d0/dt^2+6*v0/dt+2*a0);
  dp=dp+c*(3*d0/dt+2*v0+dt*a0/2);
  pbar=p1+dp;
  d1=ikbar*pbar;
  v1=3*(d1-d0)/dt-2*v0-dt*a0/2;
  a1=6*(d1-d0)/dt^2 -6*v0/dt-2*a0;
  d(i,1)=d1;
  v(i,1)=v1;
  a(i,1)=a1;
  d0=d1;
  v0=v1;
  a0=a1;
end
%
figure
plot(t,-p/m)
xlabel('Tiempo')
ylabel('Aceleración del suelo');
figure
plot(t,d)
xlabel('Tiempo')
ylabel('Desplazamiento');
figure
plot(t,v)
xlabel('Tiempo')
ylabel('Velocidad');
figure
plot(t,a)
xlabel('Tiempo')
ylabel('Aceleración');
% %
% %
%----------------------------------------------------- fin

