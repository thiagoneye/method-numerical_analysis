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

F2 = 5000; % Force [N]
F3 = -2000; % Force [N]
u1 = 0; % Displacement [m]
L1 = 1; % Lenght [m]
L2 = 1; % Lenght [m]
E = 100*10^9; % Young's Modulus [Pa]
A1 = 1*10^(-4); % Area [m^2]
A2 = 2*10^(-4); % Area [m^2]
n = 3; % Degrees of Freedom

%% Characterization

M = [1 -1; -1 1];

% Stiffness

A = [A1; A2];
x = cumsum([0; L1; L2]);

k = E.*A./(diff(x));

K1 = k(1)*M;
K2 = k(2)*M;

% Global Matrix Equation

K = zeros(n);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;

% Force Vector

syms F1
F = [F1; F2; F3];

% Displacement Vectors

syms u2 u3
U = [u1; u2; u3];

%% Calculations

AN = solve(F-K*U);
 
F = double([AN.F1; F2; F3]);
U = double([u1; AN.u2; AN.u3]);

f = k.*diff(U);

clear AN A1 A2 F1 F2 F3 L1 L2 M n u1 u2 u3