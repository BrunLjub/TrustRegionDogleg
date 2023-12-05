clc;
clear;
close all;

f = @(x) (x(1)^2+x(2)^2 - 2)^2 + (exp(x(1)-1) + x(2)^3 - 2)^2;
df = @(x)[2*exp(x(1)-1)*(exp(x(1)-1)+x(2)^3-2)+4*x(1)*(x(1)^2+x(2)^2-2); 6*x(2)^2*(x(2)^3+exp(x(1)-1)-2)+4*x(2)*(x(2)^2+x(1)^2-2)];
d2f = @(x)[2*exp(x(1)-2)*(2*exp(x(1))+exp(1)*(x(2)^3-2))+4*(3*x(1)^2+x(2)^2-2) 6*exp(x(1)-1)*x(2)^2+8*x(1)*x(2);6*exp(x(1)-1)*x(2)^2+8*x(1)*x(2) 18*x(2)^4+12*x(2)*(x(2)^3+exp(x(1)-1)-2)+4*(x(2)^2+x(1)^2-2)+8*x(2)^2];


epsilon = 1e-10;
delta= 1e-10;
dim = 2;
N = 50;

deltac = -1; %poèetni delta je Cauchyjev korak;

%dvije nultoèke su [1;1], [-0.71375;1.22089]


x0 = [0.5;0.4];
[x,k,X,FX,mi,Deltas] = dogleg(f,df,d2f,x0,dim,epsilon,delta,N, deltac)
f(x);
x0 = [-1.2;1.3];
[x,k] = dogleg(f,df,d2f,x0,dim,epsilon,delta,N, deltac)
f(x)
x0 = [-2;3];
[x,k] = dogleg(f,df,d2f,x0,dim,epsilon,delta,N, deltac)
f(x)

%%
%drugi primjer


f = @(x) x(2)^4*sin(x(1))^2+cos(x(2))^2*x(1)^4+(x(1)^2+2*x(2)^2)/2;
df = @(x) [2*x(2)^4*cos(x(1))*sin(x(1))+4*cos(x(2))^2*x(1)^3+x(1);-2*x(1)^4*cos(x(2))*sin(x(2))+4*sin(x(1))^2*x(2)^3+2*x(2)];
d2f = @(x) [-2*x(2)^4*sin(x(1))^2+2*x(2)^4*cos(x(1))^2+12*cos(x(2))^2*x(1)^2+1,8*(x(2)^3*cos(x(1))*sin(x(1))-cos(x(2))*sin(x(2))*x(1)^3);8*(x(2)^3*cos(x(1))*sin(x(1))-cos(x(2))*sin(x(2))*x(1)^3),2*(x(1)^4*sin(x(2))^2-x(1)^4*cos(x(2))^2+6*sin(x(1))^2*x(2)^2+1)];

%Jedini minimum [0,0]
epsilon = 1e-10;
delta= 1e-10;
dim = 2;
N = 50;

deltac = 1; %Poèetni delta = 1
x0 = 20*randn(2,1)
[x,k,X,FX,mi] = dogleg(f,df,d2f,x0,dim,epsilon,delta,N, deltac)

