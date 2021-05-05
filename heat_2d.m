function [ resU ] = heat_2d( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Simulating the 2-D Diffusion equation by the Finite Difference
...Method 
% Numerical scheme used is a first order upwind in time and a second
...order central difference in space (Implicit and Explicit)

%%
%Specifying parameters
a=1;
b=1;
nx=21;                           %Number of steps in space(x)
ny=21;                           %Number of steps in space(y)       
nt=9; % nt==0 mean calculation in T==1*dt; nt==i -> T==(i+1)*dt  %Number of time steps 
dt=0.005;                         %Width of each time step
dx=a/(nx-1);                     %Width of space step(x)
dy=b/(ny-1);                     %Width of space step(y)
x=0:dx:a;                        %Range of x(0,2) and specifying the grid points
y=0:dy:b;                        %Range of y(0,2) and specifying the grid points
u=zeros(nx,ny);                  %Preallocating u
un=zeros(nx,ny);                 %Preallocating un
vis=0.2;                         %Diffusion coefficient/viscocity
UW=0.0;                            %x=0 Dirichlet B.C 
UE=0.0;                            %x=L Dirichlet B.C 
US=0.0;                            %y=0 Dirichlet B.C 
UN=0.0;                            %y=L Dirichlet B.C 
UnW=0.0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0.0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0.0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0.0;                           %y=L Neumann B.C (du/dn=UnN)

problemN = 0;

%%
%Initial Conditions
% for i=1:nx
%     for j=1:ny
%         if ((1<=y(j))&&(y(j)<=1.5)&&(1<=x(i))&&(x(i)<=1.5))
%             u(i,j)=2;
%         else
%             u(i,j)=0;
%         end
%     end
% end

%==============Test0================================================
%  u = repmat(0.2, [nx ny]);
%  u(ceil(nx*0.2):ceil(nx*0.8),ceil(ny*0.2):ceil(ny*0.8)) = 0.95;
%==============Test1================================================
% for i = 1:1:nx
%     for j = 1:1:ny
%         u(i,j) = -0.1*(1+2*pi)*sin(2*pi*(dx*i+dy*j));
%     end
% end
%==============Test2================================================
% for i = 1:1:nx
%     for j = 1:1:ny
%         u(i,j) = exp(-((dx*i - 0.5*2)^2 / (2*0.3*2^2) + (dy*j-0.5*2)^2 / (2*0.1*2^2) ));
%     end
% end
%==============Test3================================================
%u(:,:) = 0.95;
%===================================================================

%==============Test4================================================
% vis=(1.0/8.0);
% u = repmat(0.0, [nx ny]);
% u(:,1:11) = 0.5;
% u(1,:)=UW;
% u(nx,:)=UE;
% u(:,1)=US;
% u(:,ny)=UN;

%==============Test5================================================
% vis=(1.0/5.0);
% u = x'*y; 
% u(1,:)=UW;
% u(nx,:)=UE;
% problemN = 5;
% u(:,1)=US;
% u(:,ny)=UN;
%==============Test6================================================
vis=(1.0/8.0);
for i = 1:1:nx
    for j = 1:1:ny
        u(i,j) = x(i)*(1-x(i))* y(j) * (1-y(j));
    end
end
u(1,:)=UW;
u(nx,:)=UE;
u(:,1)=US;
u(:,ny)=UN;
%==============Test7================================================
% vis=(1.0/3.0);
% u = repmat(0.1, [nx ny]);
% u(1,:)=UW;
% u(nx,:)=UE;
% u(:,1)=US;
% problemN = 7;

%%
%B.C vector
bc=zeros(nx-2,ny-2);
bc(1,:)=UW/dx^2; bc(nx-2,:)=UE/dx^2;  %Dirichlet B.Cs
bc(:,1)=US/dy^2; bc(:,ny-2)=UN/dy^2;  %Dirichlet B.Cs
% bc(1,:)=-UnW/dx; bc(nx-2,:)=UnE/dx;  %Neumann B.Cs
%bc(:,nx-2)=UnN/dy;  bc(:,1)=-UnS/dy;  %Neumann B.Cs
%B.Cs at the corners:
bc(1,1)=UW/dx^2+US/dy^2; bc(nx-2,1)=UE/dx^2+US/dy^2;
bc(1,ny-2)=UW/dx^2+UN/dy^2; bc(nx-2,ny-2)=UE/dx^2+UN/dy^2;
bc=vis*dt*bc;

%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
% Ax(1,1)=-1; Ax(nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
% Ay(ny-2,ny-2)=-1; Ay(1,1)=-1;   %Neumann B.Cs
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
D=speye((nx-2)*(ny-2))-vis*dt*A;

%%
%Calculating the field variable for each time step
i=2:nx-1;
j=2:ny-1;
for it=0:nt
%     pause(1.0);
% if it == nt-1
%     break;
% end

    un=u;
    h=surf(x,y,u','EdgeColor','none');       %plotting the field variable
    shading interp
    axis ([0 2 0 2 0 2])
    title({['2-D Diffusion with {\nu} = ',num2str(vis)];['time (\itt) = ',num2str(it*dt)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('{\leftarrow} Spatial co-ordinate (y)')
    zlabel('Transport property profile (u) \rightarrow')
    drawnow; 
    refreshdata(h)
    %Uncomment as necessary
    %Implicit method:
%     U=un;U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[];
%     U=reshape(U+bc,[],1);
%     U=D\U;
%     U=reshape(U,nx-2,ny-2);
%     u(2:nx-1,2:ny-1)=U;
%     %Boundary conditions
%     %Dirichlet:
%     u(1,:)=UW;
%     u(nx,:)=UE;
%     u(:,1)=US;
%     u(:,ny)=UN;
%     %Neumann:
% %     u(1,:)=u(2,:)-UnW*dx;
% %     u(nx,:)=u(nx-1,:)+UnE*dx;
% %     u(:,1)=u(:,2)-UnS*dy;
% %     u(:,ny)=u(:,ny-1)+UnN*dy;
    %}
    %Explicit method:
    % %{
    u(i,j)=un(i,j)+(vis*dt*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(dx*dx))+(vis*dt*(un(i,j+1)-2*un(i,j)+un(i,j-1))/(dy*dy));
    %Boundary conditions
    %Dirichlet:
    u(1,:)=UW;
    u(nx,:)=UE;
    if problemN == 5
        u(:,1)=u(:,2);%US;
        u(:,ny)=u(:,ny-1); %UN;
    elseif problemN == 7
        u(:,1)=US;
        u(:,ny)=u(:,ny-1); %UN;
    else
        u(:,1)=US;
        u(:,ny)=UN;
    end
    %Neumann:
    %u(1,:)=u(2,:)-UnW*dx;
    %u(nx,:)=u(nx-1,:)+UnE*dx;
    %u(:,1)=u(:,2)-UnS*dy;
    %u(:,ny)=u(:,ny-1)+UnN*dy;
    %}
end

resU = u;

end

