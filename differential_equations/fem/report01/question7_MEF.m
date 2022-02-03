% UNIVERSIDADE FEDERAL DA PARAÍBA
% CENTRO DE TECNOLOGIA
% DEPARTAMENTO DE ENGENHARIA MECÂNICA
% DISCIPLINA DE ANÁLISE MATRICIAL E MODELAGEM DE ESTRUTURAS
% ALUNO THIAGO NEY EVARISTO RODRIGUES

% 1ª Tarefa

clear
close all
clc

%% Inputs

u1 = 0; % Boundary Condition
p = 1000; % Load [N/m]
l = 1.5; % Length [m]
r = 0.1; % Ray [m]
E = 207*10^9; % Young's Modulus [Pa]

%% MEF - Método Direto (3 Elementos)

NE = 3; % Number of Elements
ND = NE + 1;
M = [1 -1; -1 1];

% Stiffness

A = pi*r^2;
X = linspace(0,l,ND);
L = diff(X);

k = E*A./L;

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;

% Global Matrix Equation

K = zeros(ND);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;

% Force Vector

syms F1 F2 F3 F4

F = [F1; F2; F3; F4];
aux = p*X;
F(2:4) = aux(2:4);

% Displacement Vectors

syms u2 u3 u4
U = [u1; u2; u3; u4];

% Calculations

AN = solve(F-K*U);
 
F(1) = AN.F1;
U(2:4) = [AN.u2; AN.u3; AN.u4];

F = double(F);
U = double(U);

figure
plot(X,F/A,'-ob','MarkerFaceColor','b')
legend('MEF - Método Direto','Location','southeast')
xlabel('Comprimento [m]')
ylabel('Tensão [Pa]')
grid minor

figure
plot(X,U,'-ob','MarkerFaceColor','b')
legend('MEF - Método Direto','Location','southeast')
xlabel('Comprimento Inicial [m]')
ylabel('Deslocamento [m]')
grid minor

figure
plot(X,F/A,'-ob','MarkerFaceColor','b')

figure
plot(X,U,'-ob','MarkerFaceColor','b')

%% MEF - Método Direto (4 Elementos)

NE = 4; % Number of Elements
ND = NE + 1;
M = [1 -1; -1 1];

% Stiffness

A = pi*r^2;
X = linspace(0,l,ND);
L = diff(X);

k = E*A./L;

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;
K4 = k(4)*M;

% Global Matrix Equation

K = zeros(ND);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;
K(4:5,4:5) = K(4:5,4:5) + K4;

% Force Vector

syms F1 F2 F3 F4 F5

F = [F1; F2; F3; F4; F5];
aux = p*X;
F(2:5) = aux(2:5);

% Displacement Vectors

syms u2 u3 u4 u5
U = [u1; u2; u3; u4; u5];

% Calculations

AN = solve(F-K*U);
 
F(1) = AN.F1;
U(2:5) = [AN.u2; AN.u3; AN.u4; AN.u5];

F = double(F);
U = double(U);

figure(3)
hold on
plot(X,F/A,'-or','MarkerFaceColor','r')

figure(4)
hold on
plot(X,U,'-or','MarkerFaceColor','r')

%% MEF - Método Direto (6 Elementos)

NE = 6; % Number of Elements
ND = NE + 1;
M = [1 -1; -1 1];

% Stiffness

A = pi*r^2;
X = linspace(0,l,ND);
L = diff(X);

k = E*A./L;

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;
K4 = k(4)*M;
K5 = k(5)*M;
K6 = k(6)*M;

% Global Matrix Equation

K = zeros(ND);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;
K(4:5,4:5) = K(4:5,4:5) + K4;
K(5:6,5:6) = K(5:6,5:6) + K5;
K(6:7,6:7) = K(6:7,6:7) + K6;

% Force Vector

syms F1 F2 F3 F4 F5 F6 F7

F = [F1; F2; F3; F4; F5; F6; F7];
aux = p*X;
F(2:7) = aux(2:7);

% Displacement Vectors

syms u2 u3 u4 u5 u6 u7
U = [u1; u2; u3; u4; u5; u6; u7];

% Calculations

AN = solve(F-K*U);
 
F(1) = AN.F1;
U(2:7) = [AN.u2; AN.u3; AN.u4; AN.u5; AN.u6; AN.u7];

F = double(F);
U = double(U);

figure(3)
hold on
plot(X,F/A,'-og','MarkerFaceColor','g')

figure(4)
hold on
plot(X,U,'-og','MarkerFaceColor','g')

%% MEF - Método Direto (8 Elementos)

NE = 8; % Number of Elements
ND = NE + 1;
M = [1 -1; -1 1];

% Stiffness

A = pi*r^2;
X = linspace(0,l,ND);
L = diff(X);

k = E*A./L;

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;
K4 = k(4)*M;
K5 = k(5)*M;
K6 = k(6)*M;
K7 = k(7)*M;
K8 = k(8)*M;

% Global Matrix Equation

K = zeros(ND);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;
K(4:5,4:5) = K(4:5,4:5) + K4;
K(5:6,5:6) = K(5:6,5:6) + K5;
K(6:7,6:7) = K(6:7,6:7) + K6;
K(7:8,7:8) = K(7:8,7:8) + K7;
K(8:9,8:9) = K(8:9,8:9) + K8;

% Force Vector

syms F1 F2 F3 F4 F5 F6 F7 F8 F9

F = [F1; F2; F3; F4; F5; F6; F7; F8; F9];
aux = p*X;
F(2:9) = aux(2:9);

% Displacement Vectors

syms u2 u3 u4 u5 u6 u7 u8 u9
U = [u1; u2; u3; u4; u5; u6; u7; u8; u9];

% Calculations

AN = solve(F-K*U);
 
F(1) = AN.F1;
U(2:9) = [AN.u2; AN.u3; AN.u4; AN.u5; AN.u6; AN.u7; AN.u8; AN.u9];

F = double(F);
U = double(U);

figure(3)
hold on
plot(X,F/A,'-oy','MarkerFaceColor','y')

figure(4)
hold on
plot(X,U,'-oy','MarkerFaceColor','y')

%% MEF - Método Direto (10 Elementos)

NE = 10; % Number of Elements
ND = NE + 1;
M = [1 -1; -1 1];

% Stiffness

A = pi*r^2;
X = linspace(0,l,ND);
L = diff(X);

k = E*A./L;

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;
K4 = k(4)*M;
K5 = k(5)*M;
K6 = k(6)*M;
K7 = k(7)*M;
K8 = k(8)*M;
K9 = k(9)*M;
K10 = k(10)*M;

% Global Matrix Equation

K = zeros(ND);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;
K(4:5,4:5) = K(4:5,4:5) + K4;
K(5:6,5:6) = K(5:6,5:6) + K5;
K(6:7,6:7) = K(6:7,6:7) + K6;
K(7:8,7:8) = K(7:8,7:8) + K7;
K(8:9,8:9) = K(8:9,8:9) + K8;
K(9:10,9:10) = K(9:10,9:10) + K9;
K(10:11,10:11) = K(10:11,10:11) + K10;

% Force Vector

syms F1 F2 F3 F4 F5 F6 F7 F8 F9 F10 F11

F = [F1; F2; F3; F4; F5; F6; F7; F8; F9; F10; F11];
aux = p*X;
F(2:11) = aux(2:11);

% Displacement Vectors

syms u2 u3 u4 u5 u6 u7 u8 u9 u10 u11
U = [u1; u2; u3; u4; u5; u6; u7; u8; u9; u10; u11];

% Calculations

AN = solve(F-K*U);
 
F(1) = AN.F1;
U(2:11) = [AN.u2; AN.u3; AN.u4; AN.u5; AN.u6; AN.u7; AN.u8; AN.u9; AN.u10; AN.u11];

F = double(F);
U = double(U);

figure(3)
hold on
plot(X,F/A,'-oc','MarkerFaceColor','c')
legend('Método Direto (3 Elementos)','Método Direto (4 Elementos)',...
    'Método Direto (6 Elementos)','Método Direto (8 Elementos)',...
    'Método Direto (10 Elementos)','Location','southeast')
xlabel('Comprimento [m]')
ylabel('Tensão [Pa]')
grid minor

figure(4)
hold on
plot(X,U,'-oc','MarkerFaceColor','c')
legend('Método Direto (3 Elementos)','Método Direto (4 Elementos)',...
    'Método Direto (6 Elementos)','Método Direto (8 Elementos)',...
    'Método Direto (10 Elementos)','Location','northwest')
xlabel('Comprimento Inicial [m]')
ylabel('Deslocamento [m]')
grid minor