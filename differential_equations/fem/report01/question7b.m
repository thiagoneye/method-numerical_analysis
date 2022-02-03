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
du1dx = 1; % Boundary Condition
p = 1000; % Load [N/m]
l = 1.5; % Length [m]
r = 0.1; % Ray [m]
E = 207*10^9; % Young's Modulus [Pa]

%% Garlekin's Method (4 Elements)

% Characterization

syms x F1 u2 u3 u4 u5 du0dx
NE = 4; % Number of Elements
ND = NE + 1;

X = linspace(0,l,ND);
L = diff(X);

for i = 1:ND
    
    if i == 1

        phi(i,1) = (X(i+1) - x)/L(i);
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;

    elseif i == 2
        
        phi(i,1) = (x - X(i-1))/L(i-1);        
        phi(i,2) = (X(i+1) - x)/L(i);
        phi(i,3) = 0;
        phi(i,4) = 0;
        
    elseif i == 3

        phi(i,1) = 0;        
        phi(i,2) = (x - X(i-1))/L(i-1);
        phi(i,3) = (X(i+1) - x)/L(i);
        phi(i,4) = 0;
        
    elseif i == 4

        phi(i,1) = 0;        
        phi(i,2) = 0;
        phi(i,3) = (x - X(i-1))/L(i-1);
        phi(i,4) = (X(i+1) - x)/L(i);

    else

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = (x - X(i-1))/L(i-1);
        
    end
    
end

% Force

P = p*X;
P = [P(2:5); P(2:5); P(2:5); P(2:5); P(2:5)];
f = P.*phi;
F = int(f(:,1),x,X(1),X(2)) + int(f(:,2),x,X(2),X(3)) + int(f(:,3),x,X(3),X(4)) ...
    + int(f(:,4),x,X(4),X(5)) + du1dx*subs(phi(:,2),x,1.5) - du0dx*subs(phi(:,1),x,0);
F(1) = F1;

% Stiffness

M = [1 -1; -1 1];

A = pi*r^2;
k = E*A./L;

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;
K4 = k(4)*M;

K = zeros(ND);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;
K(4:5,4:5) = K(4:5,4:5) + K4;

% Displacement

U = [u1; u2; u3; u4; u5];

% Calculations

AN = solve(F - K*U);
F(1) = AN.F1;
U(2:5) = [AN.u2; AN.u3; AN.u4; AN.u5];

% Plot

figure
plot(X,F/A,'-ob','MarkerFaceColor','b')

figure
plot(X,U,'-ob','MarkerFaceColor','b')

%% Garlekin's Method (6 Elements)

% Characterization

syms x F1 u2 u3 u4 u5 u6 u7 du0dx
NE = 6; % Number of Elements
ND = NE + 1;

X = linspace(0,l,ND);
L = diff(X);

for i = 1:ND
    
    if i == 1

        phi(i,1) = (X(i+1) - x)/L(i);
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;

    elseif i == 2
        
        phi(i,1) = (x - X(i-1))/L(i-1);        
        phi(i,2) = (X(i+1) - x)/L(i);
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;        
        
    elseif i == 3

        phi(i,1) = 0;        
        phi(i,2) = (x - X(i-1))/L(i-1);
        phi(i,3) = (X(i+1) - x)/L(i);
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;        
        
    elseif i == 4

        phi(i,1) = 0;        
        phi(i,2) = 0;
        phi(i,3) = (x - X(i-1))/L(i-1);
        phi(i,4) = (X(i+1) - x)/L(i);
        phi(i,5) = 0;
        phi(i,6) = 0;        

    elseif i == 5

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = (x - X(i-1))/L(i-1);
        phi(i,5) = (X(i+1) - x)/L(i);
        phi(i,6) = 0;
        
    elseif i == 6
        
        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = (x - X(i-1))/L(i-1);
        phi(i,6) = (X(i+1) - x)/L(i);
        
    else
        
        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = (x - X(i-1))/L(i-1);        
        
    end
    
end

% Force

P = p*X;
P = [P(2:7); P(2:7); P(2:7); P(2:7); P(2:7); P(2:7); P(2:7)];
f = P.*phi;
F = int(f(:,1),x,X(1),X(2)) + int(f(:,2),x,X(2),X(3)) + int(f(:,3),x,X(3),X(4)) ...
    + int(f(:,4),x,X(4),X(5)) + int(f(:,5),x,X(5),X(6)) + int(f(:,6),x,X(6),X(7)) ...
    + du1dx*subs(phi(:,2),x,1.5) - du0dx*subs(phi(:,1),x,0);
F(1) = F1;

% Stiffness

M = [1 -1; -1 1];

A = pi*r^2;
k = E*A./L;

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;
K4 = k(4)*M;
K5 = k(5)*M;
K6 = k(6)*M;

K = zeros(ND);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;
K(4:5,4:5) = K(4:5,4:5) + K4;
K(5:6,5:6) = K(5:6,5:6) + K5;
K(6:7,6:7) = K(6:7,6:7) + K6;

% Displacement

U = [u1; u2; u3; u4; u5; u6; u7];

% Calculations

AN = solve(F - K*U);
F(1) = AN.F1;
U(2:7) = [AN.u2; AN.u3; AN.u4; AN.u5; AN.u6; AN.u7];

% Plot

figure(1)
hold on
plot(X,F/A,'-or','MarkerFaceColor','r')

figure(2)
hold on
plot(X,U,'-or','MarkerFaceColor','r')

%% Garlekin's Method (8 Elements)

% Characterization

syms x F1 u2 u3 u4 u5 u6 u7 u8 u9 du0dx
NE = 8; % Number of Elements
ND = NE + 1;

X = linspace(0,l,ND);
L = diff(X);

for i = 1:ND
    
    if i == 1

        phi(i,1) = (X(i+1) - x)/L(i);
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;        

    elseif i == 2
        
        phi(i,1) = (x - X(i-1))/L(i-1);        
        phi(i,2) = (X(i+1) - x)/L(i);
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;        
        
    elseif i == 3

        phi(i,1) = 0;        
        phi(i,2) = (x - X(i-1))/L(i-1);
        phi(i,3) = (X(i+1) - x)/L(i);
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;        
        
    elseif i == 4

        phi(i,1) = 0;        
        phi(i,2) = 0;
        phi(i,3) = (x - X(i-1))/L(i-1);
        phi(i,4) = (X(i+1) - x)/L(i);
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;        

    elseif i == 5

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = (x - X(i-1))/L(i-1);
        phi(i,5) = (X(i+1) - x)/L(i);
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;        
        
    elseif i == 6
        
        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = (x - X(i-1))/L(i-1);
        phi(i,6) = (X(i+1) - x)/L(i);
        phi(i,7) = 0;
        phi(i,8) = 0;        
        
    elseif i == 7
        
        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = (x - X(i-1))/L(i-1);
        phi(i,7) = (X(i+1) - x)/L(i);
        phi(i,8) = 0;        
        
    elseif i == 8

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = (x - X(i-1))/L(i-1);
        phi(i,8) = (X(i+1) - x)/L(i);        
        
    else

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = (x - X(i-1))/L(i-1);        
        
    end
    
end

% Force

P = p*X;
P = [P(2:9); P(2:9); P(2:9); P(2:9); P(2:9); P(2:9); P(2:9); P(2:9); P(2:9)];
f = P.*phi;
F = int(f(:,1),x,X(1),X(2)) + int(f(:,2),x,X(2),X(3)) + int(f(:,3),x,X(3),X(4)) ...
    + int(f(:,4),x,X(4),X(5)) + int(f(:,5),x,X(5),X(6)) + int(f(:,6),x,X(6),X(7)) ...
    + int(f(:,7),x,X(7),X(8)) + int(f(:,8),x,X(8),X(9)) + du1dx*subs(phi(:,2),x,1.5) ...
    - du0dx*subs(phi(:,1),x,0);
F(1) = F1;

% Stiffness

M = [1 -1; -1 1];

A = pi*r^2;
k = E*A./L;

K1 = k(1)*M;
K2 = k(2)*M;
K3 = k(3)*M;
K4 = k(4)*M;
K5 = k(5)*M;
K6 = k(6)*M;
K7 = k(7)*M;
K8 = k(8)*M;

K = zeros(ND);
K(1:2,1:2) = K(1:2,1:2) + K1;
K(2:3,2:3) = K(2:3,2:3) + K2;
K(3:4,3:4) = K(3:4,3:4) + K3;
K(4:5,4:5) = K(4:5,4:5) + K4;
K(5:6,5:6) = K(5:6,5:6) + K5;
K(6:7,6:7) = K(6:7,6:7) + K6;
K(7:8,7:8) = K(7:8,7:8) + K7;
K(8:9,8:9) = K(8:9,8:9) + K8;

% Displacement

U = [u1; u2; u3; u4; u5; u6; u7; u8; u9];

% Calculations

AN = solve(F - K*U);
F(1) = AN.F1;
U(2:9) = [AN.u2; AN.u3; AN.u4; AN.u5; AN.u6; AN.u7; AN.u8; AN.u9];

% Plot

figure(1)
hold on
plot(X,F/A,'-og','MarkerFaceColor','g')

figure(2)
hold on
plot(X,U,'-og','MarkerFaceColor','g')

%% Garlekin's Method (10 Elements)

% Characterization

syms x F1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 du0dx
NE = 10; % Number of Elements
ND = NE + 1;

X = linspace(0,l,ND);
L = diff(X);

for i = 1:ND
    
    if i == 1

        phi(i,1) = (X(i+1) - x)/L(i);
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;
        phi(i,9) = 0;
        phi(i,10) = 0;

    elseif i == 2
        
        phi(i,1) = (x - X(i-1))/L(i-1);        
        phi(i,2) = (X(i+1) - x)/L(i);
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;
        phi(i,9) = 0;
        phi(i,10) = 0;        
        
    elseif i == 3

        phi(i,1) = 0;        
        phi(i,2) = (x - X(i-1))/L(i-1);
        phi(i,3) = (X(i+1) - x)/L(i);
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;
        phi(i,9) = 0;
        phi(i,10) = 0;        
        
    elseif i == 4

        phi(i,1) = 0;        
        phi(i,2) = 0;
        phi(i,3) = (x - X(i-1))/L(i-1);
        phi(i,4) = (X(i+1) - x)/L(i);
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;
        phi(i,9) = 0;
        phi(i,10) = 0;        

    elseif i == 5

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = (x - X(i-1))/L(i-1);
        phi(i,5) = (X(i+1) - x)/L(i);
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;
        phi(i,9) = 0;
        phi(i,10) = 0;        
        
    elseif i == 6
        
        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = (x - X(i-1))/L(i-1);
        phi(i,6) = (X(i+1) - x)/L(i);
        phi(i,7) = 0;
        phi(i,8) = 0;
        phi(i,9) = 0;
        phi(i,10) = 0;        
        
    elseif i == 7
        
        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = (x - X(i-1))/L(i-1);
        phi(i,7) = (X(i+1) - x)/L(i);
        phi(i,8) = 0;
        phi(i,9) = 0;
        phi(i,10) = 0;        
        
    elseif i == 8

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = (x - X(i-1))/L(i-1);
        phi(i,8) = (X(i+1) - x)/L(i);
        phi(i,9) = 0;
        phi(i,10) = 0;        
        
    elseif i == 9

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = (x - X(i-1))/L(i-1);
        phi(i,9) = (X(i+1) - x)/L(i);
        phi(i,10) = 0;
        
    elseif i == 10

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;        
        phi(i,8) = 0;                
        phi(i,9) = (x - X(i-1))/L(i-1);
        phi(i,10) = (X(i+1) - x)/L(i);
        
    else

        phi(i,1) = 0;
        phi(i,2) = 0;
        phi(i,3) = 0;
        phi(i,4) = 0;
        phi(i,5) = 0;
        phi(i,6) = 0;
        phi(i,7) = 0;
        phi(i,8) = 0;        
        phi(i,9) = 0;        
        phi(i,10) = (x - X(i-1))/L(i-1);
        
    end
    
end

% Force

P = p*X;
P = [P(2:11); P(2:11); P(2:11); P(2:11); P(2:11); P(2:11); P(2:11); P(2:11); ...
    P(2:11); P(2:11); P(2:11)];
f = P.*phi;
F = int(f(:,1),x,X(1),X(2)) + int(f(:,2),x,X(2),X(3)) + int(f(:,3),x,X(3),X(4)) ...
    + int(f(:,4),x,X(4),X(5)) + int(f(:,5),x,X(5),X(6)) + int(f(:,6),x,X(6),X(7)) ...
    + int(f(:,7),x,X(7),X(8)) + int(f(:,8),x,X(8),X(9)) + int(f(:,9),x,X(9),X(10)) ...
    + int(f(:,10),x,X(10),X(11)) + du1dx*subs(phi(:,2),x,1.5) - du0dx*subs(phi(:,1),x,0);
F(1) = F1;

% Stiffness

M = [1 -1; -1 1];

A = pi*r^2;
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

% Displacement

U = [u1; u2; u3; u4; u5; u6; u7; u8; u9; u10; u11];

% Calculations

AN = solve(F - K*U);
F(1) = AN.F1;
U(2:11) = [AN.u2; AN.u3; AN.u4; AN.u5; AN.u6; AN.u7; AN.u8; AN.u9; AN.u10; AN.u11];

% Plot

figure(1)
hold on
plot(X,F/A,'-oc','MarkerFaceColor','c')
legend('4 Elementos', '6 Elementos', '8 Elementos', '10 Elementos','Location','southeast')
xlabel('Comprimento [m]')
ylabel('Tensão [Pa]')
grid minor

figure(2)
hold on
plot(X,U,'-oc','MarkerFaceColor','c')
legend('4 Elementos', '6 Elementos', '8 Elementos', '10 Elementos','Location','southeast')
xlabel('Comprimento Inicial [m]')
ylabel('Deslocamento [m]')
grid minor