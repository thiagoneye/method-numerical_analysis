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

syms k1 k2 k3 k4 F1 F2 F3 F4 F5 u1 u2 u3

u4 = 0; % Displacement [mm]
u5 = 0; % Displacement [mm]

%% Characterization

% Stiffness

K = [k1+k3+k4 -k4 -k3 -k1 0; ...
     -k4 k4 0 0 0; ...
     -k3 0 k2+k3 0 -k2; ...
     -k1 0 0 k1 0; ...
     0 0 k2 0 -k2];

% Force Vector

F = [F1; F2; F3; F4; F5];

% Displacement Vectors

U = [u1; u2; u3; u4; u5];

clear F1 F2 F3 F4 F5 u1 u2 u3 u4 u5