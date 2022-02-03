% UNIVERSIDADE FEDERAL DA PARAÍBA
% CENTRO DE ENERGIAS ALTERNATIVAS E RENOVÁVEIS
% DEPARTAMENTO DE ENGENHARIA DE ENERGIAS RENOVÁVEIS
% DISCIPLINA DE TÓPICOS ESPECIAIS EM ENGENHARIA DE ENERGIAS RENOVÁVEIS
% ALUNO THIAGO NEY EVARISTO RODRIGUES
%
% Book: An Introduction to Computational Fluid Dynamics - Versteeg, H.K.
% and Malalasekera, W. 
% 
% Chapter 4: The Finite Volume Method for Diffusion Problems
% 
% Example 4.2.

clear
close all
clc

%% Inputs

TA = 100;      % Temperature A [°C]
TB = 200;      % Temperature B [°C]
q = 1000*10^3; % Heat generation [W/m^3]
k = 0.5;       % Thermal conductivity constant [W/m.K]
A = 1;         % Cross-sectional area [m^2]
L = 0.02;      % Length [m]
n = 5;         % Number of nodes (grid generation)

%% Numerical Solution

dx = L/n;

% Preallocation

aW = zeros(n,1);
aE = zeros(n,1);
SP = zeros(n,1);
Su = zeros(n,1);

% Node 1

aE(1) = k*A/dx;
SP(1) = -2*k*A/dx;
Su(1) = (2*k*A/dx)*TA + q*A*dx;

% Central Nodes

i = 2:(n-1);

aW(i) = k*A/dx;
aE(i) = k*A/dx;
Su(i) = q*A*dx;

% Node n

aW(n) = k*A/dx;
SP(n) = -2*k*A/dx;
Su(n) = (2*k*A/dx)*TB + q*A*dx;

% Matrix Equation

aP = aW + aE - SP;
a = diag(aP);

for i = 2:n
    
    j = i-1;
    a(i,j) = -aW(i);
    a(j,i) = -aE(j);
    
end

T = a\Su;

%% Outputs

x_num = linspace(0,L-dx,n) + dx/2;
x_num = [0; x_num'; L];
T_num = [TA; T; TB];

x_ex = linspace(0,L,100);
T_ex = ((TB - TA)/L + q*(L - x_ex)/(2*k)).*x_ex + TA;

figure
plot(x_num,T_num, 'sk',  'MarkerFaceColor', 'k')
hold on
plot(x_ex,T_ex, 'k')
legend({'Numerical Solution', 'Exact Solution'}, 'Location', 'southeast')
xlabel('Distance x [m]')
ylabel('Temperature [ºC]')
grid