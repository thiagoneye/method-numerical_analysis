function ut=pde_1(t,u)

% Problem parameters
global ncall
xl=0.0;
xu=1.0;

% PDE
n=length(u);
dx2=((xu-xl)/(n-1))^2;
for i=1:n
if(i==1) ut(i)=0.0;
elseif(i==n) ut(i)=2.0*(u(i-1)-u(i))/dx2;
else ut(i)=(u(i+1)-2.0*u(i)+u(i-1))/dx2;
end
end
ut=ut';

% Increment calls to pde_1
ncall=ncall+1;