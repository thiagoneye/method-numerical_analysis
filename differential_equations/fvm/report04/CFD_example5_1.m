% UNIVERSIDADE FEDERAL DA PARAÍBA
% CENTRO DE ENERGIAS ALTERNATIVAS E RENOVÁVEIS
% DEPARTAMENTO DE ENGENHARIA DE ENERGIAS RENOVÁVEIS
% DISCIPLINA DE TÓPICOS ESPECIAIS EM ENGENHARIA DE ENERGIAS RENOVÁVEIS
% ALUNO THIAGO NEY EVARISTO RODRIGUES
%
% Book: An Introduction to Computational Fluid Dynamics - Versteeg, H.K.
% and Malalasekera, W. 
% 
% Chapter 5: The Finite Volume Method for Convection-Diffusion Problems
% 
% Example 5.1.

clear
close all
clc

%% Inputs

n = 20;        % Number of nodes (grid generation)
u = 2.5;      % Velocity [m/s]
L = 1;        % Length [m]
rho = 1;      % Density [kg/m^3]
phi0 = 1;     % Boundary condition (x = 0 m)
phiL = 0;     % Boundary condition (x = L m)
Gamma = 0.1;  % Dynamic viscosity [kg/m.s]

%% Numerical Solution

dx = L/n;
F = rho*u;
D = Gamma/dx;

% Preallocation

aW = zeros(n,1);
aE = zeros(n,1);
SP = zeros(n,1);
Su = zeros(n,1);

% Node 1

aE(1) = D - F/2;
SP(1) = -(2*D + F);
Su(1) = (2*D + F)*phi0;

% Central Nodes

i = 2:(n-1);

aW(i) = D + F/2;
aE(i) = D - F/2;

% Node n

aW(n) = D + F/2;
SP(n) = -(2*D - F);
Su(n) = (2*D - F)*phiL;

% Matrix Equation

aP = aW + aE - SP;
a = diag(aP);

for i = 2:n
    
    j = i-1;
    a(i,j) = -aW(i);
    a(j,i) = -aE(j);
    
end

phi = a\Su;

%% Outputs

x_num = linspace(0,L-dx,n) + dx/2;
x_num = [0; x_num'; L];
f_num = [phi0; phi; phiL];

x_ex = linspace(0,L,100);
f_ex = phi0 + (phiL - phi0)*(exp(rho*u*x_ex/Gamma) - 1)/(exp(rho*u*L/Gamma) - 1);

figure
plot(x_num,f_num, 'sk',  'MarkerFaceColor', 'k')
hold on
plot(x_ex,f_ex, 'k')
legend('Numerical Solution', 'Exact Solution')
legend('Location', 'southwest')
title('u = 2.5 m/s')
xlabel('Distance x [m]')
ylabel('Property \phi')
grid