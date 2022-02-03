function ut=pde_2(t,u)
%
% Problem parameters
global ncall ndss
xl=0.0;
xu=1.0;
%
% BC at x = 0 (Dirichlet)
u(1)=0.0;
%
% Calculate ux
n=length(u);
if (ndss== 2) ux=dss002(xl,xu,n,u); % second order
elseif(ndss== 4) ux=dss004(xl,xu,n,u); % fourth order
elseif(ndss== 6) ux=dss006(xl,xu,n,u); % sixth order
elseif(ndss== 8) ux=dss008(xl,xu,n,u); % eighth order
elseif(ndss==10) ux=dss010(xl,xu,n,u); % tenth order
end
%
% BC at x = 1 (Neumann)
ux(n)=0.0;
%
% Calculate uxx
if (ndss== 2) uxx=dss002(xl,xu,n,ux); % second order
elseif(ndss== 4) uxx=dss004(xl,xu,n,ux); % fourth order
elseif(ndss== 6) uxx=dss006(xl,xu,n,ux); % sixth order
elseif(ndss== 8) uxx=dss008(xl,xu,n,ux); % eighth order
elseif(ndss==10) uxx=dss010(xl,xu,n,ux); % tenth order
end
%
% PDE
ut=uxx';
ut(1)=0.0;
%
% Increment calls to pde_2
ncall=ncall+1;