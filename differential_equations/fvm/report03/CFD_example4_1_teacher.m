% Exemplo 01 do livro: An Introduction to Computational
% Fluid Dynamics_THE FINITE VOLUME METHOD_2ªEd.

clear
clc

L = 0.5;
E = 5; % Número de pontos de divisão
Dx = L/E;
Ta = 100;
Tb = 500;
K = 1000;
A =10^-2;

% Elementos iniciais;

aW(1) = 0;
aE(1)= (K/Dx)*A;
sP(1) = -2*(K/Dx)*A;
aP(1) = aW(1) +aE(1) -sP(1);
sU(1) = 2*Ta*(K/Dx)*A;

% Elementos centrais

i=2:E-1;
aW(i)=(K/Dx)*A;
aE(i)=aW(i);
aP(i) = aW(i) + aE(i);
sP(i) = 0;
sU(i) = 0;

% Elementos finais;

aW(E) = (K/Dx)*A;
aE(E) = 0;
sP(E) = -2*(K/Dx)*A;
aP(E) = aW(E) +aE(E) -sP(E);
sU(E) = 2*Tb*(K/Dx)*A;

% Matriz

h=diag(aP);
for j=1:E-1
    u=j+1;
    h(j,u)=-aE(j);
    h(u,j)=h(j,u);
end
r=inv(h);
Resultado=r*sU';