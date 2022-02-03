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
du15dx = 1; % Boundary Condition
p = 1000; % Load [N/m]
l = 1.5; % Length [m]
r = 0.1; % Ray [m]
E = 207*10^9; % Young's Modulus [Pa]
NE = 3; % Number of Elements

%% FEM - Direct Method

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

% Plots

figure
plot(X,F/A,'-ob','MarkerFaceColor','b')

figure
plot(X,U,'-ob','MarkerFaceColor','b')

%% Garlekin's Method

% Characterization

syms x F1 u2 u3 u4 du0dx
ND = NE + 1;

X = linspace(0,l,ND);
L = diff(X);

for i = 1:ND
    
    if i == 1

        phi(i,1) = (X(i+1) - x)/L(i);
        phi(i,2) = 0;
        phi(i,3) = 0;

    elseif i == ND

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = (x - X(i-1))/L(i-1);

    elseif i == 2
        
        phi(i,1) = (x - X(i-1))/L(i-1);        
        phi(i,2) = (X(i+1) - x)/L(i);
        phi(i,3) = 0;
        
    else

        phi(i,1) = 0;        
        phi(i,2) = (x - X(i-1))/L(i-1);
        phi(i,3) = (X(i+1) - x)/L(i);
        
    end
    
end

% Force

P = p*X;
P = [P(2:4); P(2:4); P(2:4); P(2:4)];
f = P.*phi;
F = int(f(:,1),x,0,0.5) + int(f(:,2),x,0.5,1) + int(f(:,3),x,1,1.5) + du15dx*subs(phi(:,2),x,1.5) - du0dx*subs(phi(:,1),x,0);
F(1) = F1;

% Stiffness

M = [1 -1; -1 1];

A = pi*r^2;
k = E*A./L;

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;

K = zeros(ND);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;

% Displacement

U = [u1; u2; u3; u4];

% Calculations

AN = solve(F - K*U);
F(1) = AN.F1;
U(2:4) = [AN.u2; AN.u3; AN.u4];

% Plot

figure(1)
hold on
plot(X,F/A,'-or','MarkerFaceColor','r')
legend('Método Direto', 'Método de Garlekin','Location','southeast')
xlabel('Comprimento [m]')
ylabel('Tensão [Pa]')
grid minor

figure(2)
hold on
plot(X,U,'-or','MarkerFaceColor','r')
legend('Método Direto', 'Método de Garlekin','Location','southeast')
xlabel('Comprimento Inicial [m]')
ylabel('Deslocamento [m]')
grid minor