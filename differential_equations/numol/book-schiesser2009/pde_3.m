function ut=pde_3(t,u)
%
% Problem parameters
global ncall ndss
xl=0.0;
xu=1.0;
%
% BC at x = 0
u(1)=0.0;
%
% BC at x = 1
n=length(u);
ux(n)=0.0;
%
% Calculate uxx
nl=1; % Dirichlet
nu=2; % Neumann
if (ndss==42) uxx=dss042(xl,xu,n,u,ux,nl,nu);
% second order
elseif(ndss==44) uxx=dss044(xl,xu,n,u,ux,nl,nu);
% fourth order
elseif(ndss==46) uxx=dss046(xl,xu,n,u,ux,nl,nu);
% sixth order
elseif(ndss==48) uxx=dss048(xl,xu,n,u,ux,nl,nu);
% eighth order
elseif(ndss==50) uxx=dss050(xl,xu,n,u,ux,nl,nu);
% tenth order
end
%
% PDE
ut=uxx';
ut(1)=0.0;
%
% Increment calls to pde_3
ncall=ncall+1;