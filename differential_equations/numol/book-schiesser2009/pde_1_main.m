% clear
% close all
clc

% Parameters shared with the ODE routine
global ncall ndss

% Initial condition
n=21;
for i=1:n
	u0(i)=sin((pi/2.0)*(i-1)/(n-1));
end

% Independent variable for ODE integration
t0=0.0;
tf=2.5;
tout=linspace(t0,tf,n);
nout=n;
ncall=0;

% ODE integration
mf=1;
reltol=1.0e-04; abstol=1.0e-04;
options=odeset('RelTol',reltol,'AbsTol',abstol);
if(mf==1) % explicit FDs
	[t,u]=ode15s(@pde_1,tout,u0,options);
end
if(mf==2) ndss=4; % ndss = 2, 4, 6, 8 or 10 required
	[t,u]=ode15s(@pde_2,tout,u0,options);
end
if(mf==3) ndss=44; % ndss = 42, 44, 46, 48 or 50 required
	[t,u]=ode15s(@pde_3,tout,u0,options);
end

% Store numerical and analytical solutions, errors at x = 1/2
n2=(n-1)/2.0+1;
sine=sin(pi/2.0*0.5);
for i=1:nout
	u_plot(i)=u(i,n2);
	u_anal(i)=exp(-pi^2/4.0*t(i))*sine;
	err_plot(i)=u_plot(i)-u_anal(i);
end

% Display selected output
fprintf('\n mf = %2d abstol = %8.1e reltol = %8.1e\n', ...
    mf,abstol,reltol);
fprintf('\n t u(0.5,t) u_anal(0.5,t) err u(0.5,t)\n');
for i=1:5:nout
    fprintf('%6.3f%15.6f%15.6f%15.7f\n', ...
        t(i),u_plot(i),u_anal(i),err_plot(i));
end
fprintf('\n ncall = %4d\n',ncall);

% Plot numerical solution and errors at x = 1/2
figure
subplot(1,2,1)
plot(t,u_plot); axis tight
title('u(0.5,t) vs t'); xlabel('t'); ylabel('u(0.5,t)')
subplot(1,2,2)
plot(t,err_plot); axis tight
title('Err u(0.5,t) vs t'); xlabel('t');
ylabel('Err u(0.5,t)');
% print -deps pde.eps; print -dps pde.ps; print -dpng pde.png

% Plot numerical solution in 3D perspective
figure
colormap('Gray');
C=ones(n);
g=linspace(0,1,n); % For distance x
h1=waterfall(t,g,u',C);
axis('tight');
grid off
xlabel('Axis x: t, time')
ylabel('Axis y: x, distance')
zlabel('Axis z: u(x,t)')
s1=sprintf('Diffusion Equation - MOL Solution');
sTmp=sprintf('u(x,0) = sin(\\pi x/2 )');
s2=sprintf('Initial condition: %s',sTmp);
title([{s1}, {s2}],'fontsize',12);
rotate3d on;