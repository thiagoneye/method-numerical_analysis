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
F3 = -1000; % Force [N]
F4 = 0; % Force [N]
u1 = 0; % Displacement [mm]
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

syms F1 F5
F = [F1; F2; F3; F4; F5];

% Displacement Vectors

syms u2 u3 u4
U = [u1; u2; u3; u4; u5];

%% Calculations

AN = solve(F-K*U);

F = double([AN.F1; F2; F3; F4; AN.F5]);
U = double([u1; AN.u2; AN.u3; AN.u4; u5]);

f1 = k1*(U(1) - U(2));
f2 = k2*(U(2) - U(4));
f3 = k3*(U(2) - U(3));
f4 = k4*(U(1) - U(3));
f5 = k5*(U(3) - U(4));
f6 = k6*(U(4) - U(5));

f = [f1; f2; f3; f4; f5; f6];

clear A AN f1 f2 f3 f4 f5 f6 F1 F2 F3 F4 F5 K1 K2 K3 K4 K5 K6 n u1 u2 u3 u4 u5