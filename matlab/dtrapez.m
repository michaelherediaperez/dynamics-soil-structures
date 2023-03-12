function [t,d]=dtrapez(p,m,w,xi,dt)
%
%----------------------------------------------------- 
% [t,d]=dtrapez(p,m,w,xi,dt)
%----------------------------------------------------- 
%
% Calcula la integral de Duhamel (respuesta de un sistema 
% sencillo lineal) por la regla de los trapecios.
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
%
%
%----------------------------------------------------- 
%
%
n=length(p);
tmax=dt*n;
t=linspace(0,tmax,n)';
wa=w*sqrt(1-xi^2);
f=p.*cos(wa*t);
g=p.*sin(wa*t);
 f1=[0; f(1:n-1)];
 g1=[0; g(1:n-1)];
pc=f1*exp(-xi*w*dt)+f;
ps=g1*exp(-xi*w*dt)+g;
pc=pc*dt/m/wa/2;
ps=ps*dt/m/wa/2;
for i=1:n
   if i==1
      c(i,1)=0;
      s(i,1)=0;
   else
      c(i,1)=c(i-1,1)*exp(-xi*w*dt)+pc(i,1);
      s(i,1)=s(i-1,1)*exp(-xi*w*dt)+ps(i,1);
   end
end
d=c.*sin(wa*t)-s.*cos(wa*t);
figure
plot(t,d)
xlabel('Tiempo')
ylabel('Desplazamiento')
%
%
%----------------------------------------------------- fin

