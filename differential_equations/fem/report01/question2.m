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

F3 = -1000; % Force [N]
F4 = 0; % Force [N]
u1 = 0; % Displacement [mm]
u2 = 0; % Displacement [mm]
u5 = 0; % Displacement [mm]
k1 = 500; % Spring Constant [N/mm]
k2 = 400; % Spring Constant [N/mm]
k3 = 600; % Spring Constant [N/mm]
k4 = 200; % Spring Constant [N/mm]
k5 = 400; % Spring Constant [N/mm]
k6 = 300; % Spring Constant [N/mm]
n = 5; % Degrees of Freedom

%% Characterization

A = [1 -1; -1 1];

% Stiffness

K1 = k1*A;
K2 = k2*A;
K3 = k3*A;
K4 = k4*A;
K5 = k5*A;
K6 = k6*A;

% Global Matrix Equation

K = zeros(n);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:2:4,2:2:4) = K(2:2:4,2:2:4) + K2;
K(2:3,2:3) = K(2:3,2:3) + K3;
K(1:2:3,1:2:3) = K(1:2:3,1:2:3) + K4;
K(3:4,3:4) = K(3:4,3:4) + K5;
K(4:5,4:5) = K(4:5,4:5) + K6;

% Force Vector

syms F1 F2 F5
F = [F1; F2; F3; F4; F5];

% Displacement Vectors

syms u3 u4
U = [u1; u2; u3; u4; u5];

%% Calculations

AN = solve(F-K*U);

F = double([AN.F1; AN.F2; F3; F4; AN.F5]);
U = double([u1; u2; AN.u3; AN.u4; u5]);

clear A AN F1 F2 F3 F4 F5 K1 K2 K3 K4 K5 K6 n u1 u2 u3 u4 u5