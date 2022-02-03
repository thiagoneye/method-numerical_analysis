% UNIVERSIDADE FEDERAL DA PARAÍBA
% CENTRO DE TECNOLOGIA
% DEPARTAMENTO DE ENGENHARIA MECÂNICA
% DISCIPLINA DE ANÁLISE MATRICIAL E MODELAGEM DE ESTRUTURAS
% ALUNO THIAGO NEY EVARISTO RODRIGUES

% 1ª TAREFA

clear
close all
clc

%% Inputs

F2 = 0; % Force [N]
F3 = 5000; % Force [N]
F4 = 0; % Force [N]
F5 = -2000; % Force [N]
u1 = 0; % Displacement [m]
L1 = 1; % Lenght [m]
L2 = 1; % Lenght [m]
E = 100*10^9; % Young's Modulus [Pa]
A1 = 1*10^(-4); % Area [m^2]
A2 = 2*10^(-4); % Area [m^2]
n = 5; % Degrees of Freedom

%% Characterization

M = [1 -1; -1 1];

% Stiffness

A = [A1; A1; A2; A2];
x = cumsum([0; L1/2; L1/2; L2/2; L2/2]);

k = @(E,x,A) E.*A./(diff(x));
k = k(E,x,A);

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;
K4 = k(4)*M;

% Global Matrix Equation

K = zeros(n);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;
K(4:5,4:5) = K(4:5,4:5) + K4;

% Force Vector

syms F1
F = [F1; F2; F3; F4; F5];

% Displacement Vectors

syms u2 u3 u4 u5
U = [u1; u2; u3; u4; u5];

%% Calculations

AN = solve(F-K*U);
 
F = double([AN.F1; F2; F3; F4; F5]);
U = double([u1; AN.u2; AN.u3; AN.u4; AN.u5]);

f = k.*diff(U);

clear AN A1 A2 F1 F2 F3 F4 F5 L1 L2 M n u1 u2 u3 u4 u5