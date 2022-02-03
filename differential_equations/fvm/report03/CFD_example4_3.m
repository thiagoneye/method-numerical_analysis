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
% Example 4.3.

clear
close all
clc

%% Inputs

TB = 100;      % Temperature B [°C]
Tinf = 20;     % Ambiente Temperature [ºC]
L = 1;         % Length [m]
h = 25;        % h = n^2 [m^-2]
n = 10;         % Number of nodes (grid generation)

%% Numerical Solution

dx = L/n;

% Preallocation

aW = zeros(n,1);
aE = zeros(n,1);
SP = zeros(n,1);
Su = zeros(n,1);

% Node 1

aE(1) = inv(dx);
SP(1) = -h*dx - 2/dx;
Su(1) = h*dx*Tinf + 2*TB/dx;

% Central Nodes

i = 2:(n-1);

aW(i) = inv(dx);
aE(i) = inv(dx);
SP(i) = -h*dx;
Su(i) = h*dx*Tinf;

% Node n

aW(n) = inv(dx);
SP(n) = -h*dx;
Su(n) = h*dx*Tinf;

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
x_num = [0; x_num'];
T_num = [TB; T];

x_ex = linspace(0,L,100);
T_ex = Tinf + (TB - Tinf)*(cosh(sqrt(h)*(L - x_ex))/cosh(sqrt(h)*L));

figure
plot(x_num,T_num, 'sk',  'MarkerFaceColor', 'k')
hold on
plot(x_ex,T_ex, 'k')
legend('Numerical Solution', 'Exact Solution')
xlabel('Distance x [m]')
ylabel('Temperature [ºC]')
grid