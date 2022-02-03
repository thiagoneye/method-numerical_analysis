% UNIVERSIDADE FEDERAL DA PARAÍBA
% CENTRO DE TECNOLOGIA
% DEPARTAMENTO DE ENGENHARIA MECÂNICA
% 
% TRABALHO DE CONCLUSÃO DE CURSO
% 
% ANÁLISE NUMÉRICA DE CONDUÇÃO TRANSIENTE COM CONDIÇÃO DE CONTORNO NÃO
% LINEAR PELO MÉTODO DAS LINHAS
% 
% DISCENTE THIAGO NEY EVARISTO RODRIGUES
% ORIENTADOR DR. JACQUES CÉSAR DOS SANTOS

clear
close all
clc

%% Inputs

tab = load('..\paper-cotta2015\tab3\results.txt');

nx = 1001; % Points in spatial grid
nt = 1001; % Points in temporal grid
x0 = 0; % Beginning of the X axis
xl = 1; % End of the X axis (Length L)
t0 = 0; % Start time
tl = 0.5; % End time

%% Calculations

x = linspace(x0,xl,nx); % X axis
t = linspace(t0,tl,nt); % Time
T0 = ones(size(x));       % Initial temperature

[~,T] = ode15s(@d2udx2,t,T0); % Solution

x03 = linspace(x0,xl,11);
[~,idx] = ismember(x,x03);
idx = find(idx);
time03 = find(t == 0.3);

EA = tab(:,3) - T(time03,idx)'; % Absolute error
Er = EA./tab(:,3); % Relative error

%% Plots

[X,Y] = meshgrid(t,x);

figure
colormap('winter')
surf(X,Y,T','EdgeAlpha', 0.05)
xlabel('Time [s]')
ylabel('Distance')
zlabel('Temperature')
% view(-20,40)
view(-25,35)
% view(-15,25)

figure
plot(tab(:,1),tab(:,3), '-.')
hold on
plot(x,T(time03,:), '--')
xlabel('Distance')
ylabel('Temperature')
legend('Cotta', 'Ney')
% legend({'Cotta', 'Ney'}, 'Location', 'southoutside', ...
%     'Orientation', 'horizontal')
% axis equal
grid

%% Function

function ut = d2udx2(t,u)

	% Problem parameters

	x0 = 0;
	xl = 1;
	nx = length(u);
    Bic = 0;
    Bir = 20;
    gamma = 1/3;
	alpha = 1;
    
	% Partial derivative

	dx = (xl - x0)/(nx-1);
    dx2 = dx^2;

	uxx = zeros(nx,1); % Preallocation

    % Partial derivative
    
	for i = 2:nx-1
		
		uxx(i) = (u(i+1) - 2*u(i) + u(i-1))/dx2;
		
	end

	% Boundary condition

	uxx(1) = 2*(u(2) - u(1))/dx2;
	uxx(end) = (u(end-1) + u(end)*(-2*Bic*((u(end))^(1/3))*dx - ...
        2*Bir*dx*((1 + gamma*(u(end)) + ((gamma^2)/2)*((u(end))^2))*(1 ...
        + gamma*(u(end))/2))) - 2*u(end) + u(end-1))/dx2;

    % Partial differential equation

	ut = alpha*uxx;
  
end