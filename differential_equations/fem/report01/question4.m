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
F3 = 10000; % Force [N]
F4 = 0; % Force [N]
u1 = 0; % Displacement [m]
u2 = 0; % Displacement [m]
L = 1; % Lenght [m]
E = 100*10^6; % Young's Modulus [Pa]
n = 5; % Degrees of Freedom

%% Characterization

M = [1 -1; -1 1];

% Stiffness

x = linspace(0,L,n)';
r = 0.05 - 0.04.*(x); % Radius [m]
A = @(r) pi*r(1:end-1).*r(2:end); % Area [m^2]
A = A(r);

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

syms F1 F5
F = [F1; F2; F3; F4; F5];

% Displacement Vectors

syms u3 u4 u5
U = [u1; u3; u4; u5; u2];

%% Calculations

AN = solve(F-K*U);
 
F = double([AN.F1; F2; F3; F4; AN.F5]);
U = double([u1; AN.u3; AN.u4; AN.u5; u2]);

f = k.*diff(U);

clear AN F1 F2 F3 F4 F5 K1 K2 K3 K4 M n u1 u2 u3 u4 u5