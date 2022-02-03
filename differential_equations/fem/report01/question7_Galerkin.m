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

syms x F1 u2 u3 u4 du0dx

u1 = 0; % Boundary Condition
du1dx = 1; % Boundary Condition

p = 1000; % Load [N/m]
l = 1.5; % Length [m]
NE = 3; % Number of Elements

%% Characterization

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

dphidx = diff(phi,x);

% Stiffness

k11 = dphidx(1,:).^2;
K11 = int(k11(1),0,0.5) + int(k11(2),0.5,1) + int(k11(3),1,1.5);

k12 = dphidx(1,:).*dphidx(2,:);
K12 = int(k12(1),0,0.5) + int(k12(2),0.5,1) + int(k12(3),1,1.5);

k13 = dphidx(1,:).*dphidx(3,:);
K13 = int(k13(1),0,0.5) + int(k13(2),0.5,1) + int(k13(3),1,1.5);

k14 = dphidx(1,:).*dphidx(4,:);
K14 = int(k14(1),0,0.5) + int(k14(2),0.5,1) + int(k14(3),1,1.5);

k22 = dphidx(2,:).^2;
K22 = int(k22(1),0,0.5) + int(k22(2),0.5,1) + int(k22(3),1,1.5);

k23 = dphidx(2,:).*dphidx(3,:);
K23 = int(k23(1),0,0.5) + int(k23(2),0.5,1) + int(k23(3),1,1.5);

k24 = dphidx(2,:).*dphidx(4,:);
K24 = int(k24(1),0,0.5) + int(k24(2),0.5,1) + int(k24(3),1,1.5);

k33 = dphidx(3,:).^2;
K33 = int(k33(1),0,0.5) + int(k33(2),0.5,1) + int(k33(3),1,1.5);

k34 = dphidx(3,:).*dphidx(4,:);
K34 = int(k34(1),0,0.5) + int(k34(2),0.5,1) + int(k34(3),1,1.5);

k44 = dphidx(4,:).^2;
K44 = int(k44(1),0,0.5) + int(k44(2),0.5,1) + int(k44(3),1,1.5);

K = double([K11 K12 K13 K14; K12 K22 K23 K24; K13 K23 K33 K34; K14 K24 K34 K44]);

% Force

P = p*X;
P = [P(2:4); P(2:4); P(2:4); P(2:4)];
f = P.*phi;
F = int(f(:,1),x,0,0.5) + int(f(:,2),x,0.5,1) + int(f(:,3),x,1,1.5) + du1dx*subs(phi(:,2),x,1) - du0dx*subs(phi(:,1),x,0);
F(1) = F1;

% Displacement

U = [u1; u2; u3; u4];

%% Calculations

AN = solve(F - K*U);
F(1) = AN.F1;
U(2:4) = [AN.u2; AN.u3; AN.u4];

% clear AN i x F1 u2 u3 u4 du0dx dphidx phi1 phi2 k11 k12 k13 k14 k22 k23 k24 k33 k34 k44 K11 K12 K13 K14 K22 K23 K24 K33 K34 K44 f