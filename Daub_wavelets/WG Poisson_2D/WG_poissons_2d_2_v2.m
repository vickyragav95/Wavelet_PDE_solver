% Uxx+Uyy=0
% u(0,y)=0 u(1,y)=y^2 u(x,0)=0 u(x,1)=x^2
% Exact Solution -> u(x,y) = x^2*y^2

clear;
clc;

j = 5;
global D; 
D = 6;
DAUB = n_phi(0);
global M;
M = 1;

% Initialize the matrix and vector
OMEGA = zeros((2^j*M+D-2)^2);
RHS = zeros((2^j*M+D-2)^2,1);

% Loading the necessary connection coefficient vectors
TAU_2 = conn_coeff(2);
TAU_2x = conncoeff_int(2);
TAU_0 = conn_coeff(0);
TAU_0x = conncoeff_int(0);
THET1 = theta_1();

% Initialize the parameter matrices to avoid redundancy in the script
cc0 = zeros(2^j*M+D-2,2^j*M+D-2);
cc2 = cc0;
mom2 = zeros(2^j*M+D-2,1);
mom0 = mom2;
phi0 = mom2;
phi1 = mom2;

% Evaluate the necessary parameters
for b=2-D:2^j*M-1
    for v=2-D:2^j*M-1
        cc0(b+D-1,v+D-1) = cc_temp(0,j,b,v,TAU_0,TAU_0x);
        cc2(b+D-1,v+D-1) = cc_temp(2,j,b,v,TAU_2,TAU_2x);
    end
    mom2(b+D-1,1) = moment_temp(2,j,b);
    mom0(b+D-1,1) = moment_temp(0,j,b);    
    phi0(b+D-1,1) = 2^(j/2)*phi(-b,DAUB);
    phi1(b+D-1,1) = 2^(j/2)*phi(2^j*M-b,DAUB);
end

% Evaluate the OMEGA matrix and RHS vector
for k=2-D:2^j*M-1
    for l=2-D:2^j*M-1
        alpha = (2^j*M+D-2)*(k+D-2)+l+D-1;
        for m=2-D:2^j*M-1
            for n=2-D:2^j*M-1
                beta = (2^j*M+D-2)*(m+D-2)+n+D-1;                
                OMEGA(beta,alpha) = (M^4)*(cc2(k+D-1,m+D-1)*cc0(l+D-1,n+D-1)...
                    + cc0(k+D-1,m+D-1)*cc2(l+D-1,n+D-1));
                if alpha == 1
                    RHS(beta,1) = 2*mom2(m+D-1)*mom0(n+D-1)+...
                        2*mom2(n+D-1)*mom0(m+D-1);
                end
            end
        end
    end
end
%%
% Initialize matrix and vector for BCs
bc = zeros(2^j*M+D-2,(2^j*M+D-2)^2);
rhs = zeros(2^j*M+D-2,1);
%%
% bc #1 u(0,y)=0
t=0;
for n=2-D:2^j*M-1
    t=t+1;
    rhs(t,1) = 0;
    for k=2-D:2^j*M-1
        for l=2-D:2^j*M-1
            alpha = (2^j*M+D-2)*(k+D-2)+l+D-1;
            bc(t,alpha) = phi0(k+D-1)*cc0(l+D-1,n+D-1);              
        end
    end
end
% 
OMEGA = [bc;OMEGA];
RHS = [rhs;RHS];
%%
% bc #2 u(1,y)=y^2
t=0;
for n=2-D:2^j*M-1
    t=t+1;
    rhs(t,1)= mom2(n+D-1);
    for k=2-D:2^j*M-1
        for l=2-D:2^j*M-1
            alpha = (2^j*M+D-2)*(k+D-2)+l+D-1;
            bc(t,alpha) = phi1(k+D-1)*cc0(l+D-1,n+D-1);            
        end
    end
end
% 
OMEGA = [bc;OMEGA];
RHS = [rhs;RHS];
%%
% bc #3 u(x,0)=0
t=0;
for m=2-D:2^j*M-1
    t=t+1;
    rhs(t,1)=0;
    for k=2-D:2^j*M-1
        for l=2-D:2^j*M-1
            alpha = (2^j*M+D-2)*(k+D-2)+l+D-1;
            bc(t,alpha) = phi0(l+D-1)*cc0(k+D-1,m+D-1);
        end
    end
end
% 
OMEGA = [bc;OMEGA];
RHS = [rhs;RHS];
%%
% bc #4 u(x,1)=x^2
t=0;
for m=2-D:2^j*M-1
    t=t+1;
    rhs(t,1) = mom2(m+D-1);
    for k=2-D:2^j*M-1
        for l=2-D:2^j*M-1
            alpha = (2^j*M+D-2)*(k+D-2)+l+D-1;
            bc(t,alpha) = phi1(l+D-1)*cc0(k+D-1,m+D-1);            
        end
    end
end
OMEGA = [bc;OMEGA];
RHS = [rhs;RHS];
%%
% C = pinv(OMEGA)*RHS;
C = OMEGA\RHS; 
%%
% converting C_alpha vector to c_{k,l}
c = zeros(2^j*M+D-2,2^j*M+D-2);
for k=2-D:2^j*M-1
    for l=2-D:2^j*M-1
        alpha = (2^j*M+D-2)*(k+D-2)+l+D-1;
        c(k+D-1,l+D-1) = C(alpha);
    end
end
%%
% Solution
N = 20;
u=zeros(N+1,N+1);
U = u;
for xc=0:N
    Xc = M*xc/N;
    for yc=0:N
        Yc = M*yc/N;
        for k=2-D:2^j*M-1
            for l=2-D:2^j*M-1
                u(xc+1,yc+1) = u(xc+1,yc+1)+c(k+D-1,l+D-1)*2^j*phi(2^j*Xc-k,DAUB)*phi(2^j*Yc-l,DAUB);
            end
        end
    end
end

% Exact solution
for xc=0:N
    Xc = M*xc/N;
    for yc=0:N
        Yc = M*yc/N;
        U(xc+1,yc+1) = (Xc^2)*(Yc^2);
    end
end

xc=0:(1/N):1;
yc=0:(1/N):1;
figure(1)
surf(xc,yc,u);
axis ([0 1 0 1 -0.1 1])
title(['D=' num2str(D) '  j=' num2str(j) '  Wavelet-Galerkin Solution'])
xlabel('x','FontSize',20,'FontWeight','bold')
ylabel('y','FontSize',20,'FontWeight','bold')
zlabel('u(x,y)','FontSize',20,'FontWeight','bold')
set(gca,'fontsize',12)
% print('WG_poissons_1', '-dpng', '-r600');
figure(2)
surf(xc,yc,U);
axis ([0 1 0 1 -0.1 1])
title('Exact Solution')
xlabel('x','FontSize',20,'FontWeight','bold')
ylabel('y','FontSize',20,'FontWeight','bold')
zlabel('U(x,y)','FontSize',20,'FontWeight','bold')
set(gca,'fontsize',12)
% print('poissons_1', '-dpng', '-r600');
err = u-U;
figure(3)
surf(xc,yc,err);
title('Error')
xlabel('x','FontSize',20,'FontWeight','bold')
ylabel('y','FontSize',20,'FontWeight','bold')
zlabel('u-U','FontSize',20,'FontWeight','bold')
set(gca,'fontsize',12)
% print('poissons_1_err', '-dpng', '-r600');
err = err(:);
disp(max(abs(err)));
