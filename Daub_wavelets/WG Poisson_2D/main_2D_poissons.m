% Uxx+Uyy=0
% u(0,y)=0 u(1,y)=y^2 u(x,0)=0 u(x,1)=x^2
% u(x,y) = x^2*y^2

clear;
clc;

j = 4;
global D; 
D = 6;
DAUB = n_phi(0);
global M;
% M = (D-1);
M = 1;

OMEGA = zeros((2^j*M+D-2)^2);
O = zeros((2^j*M+D-2)^2,1);

TAU_2 = conn_coeff(2);
TAU_2x = conncoeff_int(2);
TAU_0 = conn_coeff(0);
TAU_0x = conncoeff_int(0);
THET1 = theta_1();

cc0 = zeros(2^j*M+D-2,2^j*M+D-2);
cc2 = cc0;
mom2 = zeros(2^j*M+D-2,1);
tita = mom2;
phi0 = mom2;
phi1 = mom2;

for b=2-D:2^j*M-1
    for v=2-D:2^j*M-1
        cc0(b+D-1,v+D-1) = cc_temp(0,j,b,v,TAU_0,TAU_0x);
        cc2(b+D-1,v+D-1) = cc_temp(2,j,b,v,TAU_2,TAU_2x);
    end
    mom2(b+D-1,1) = moment_temp(2,j,b);
    tita(b+D-1,1) = theta1_temp(j,b,THET1);
    phi0(b+D-1,1) = 2^(j/2)*phi(-b,DAUB);
    phi1(b+D-1,1) = 2^(j/2)*phi(2^j*M-b,DAUB);
end

for k=2-D:2^j*M-1
    for l=2-D:2^j*M-1
        alpha = (2^j*M+D-2)*(k+D-2)+l+D-1;
        for m=2-D:2^j*M-1
            for n=2-D:2^j*M-1
                beta = (2^j*M+D-2)*(m+D-2)+n+D-1;                
                OMEGA(beta,alpha) = (M^4)*(cc2(k+D-1,m+D-1)*cc0(l+D-1,n+D-1)...
                    + cc0(k+D-1,m+D-1)*cc2(l+D-1,n+D-1));
                if alpha == 1
                    O(beta,1) = 2*mom2(m+D-1)*tita(n+D-1)+...
                        2*mom2(n+D-1)*tita(m+D-1);
                end
            end
        end
    end
end
%%
temp = rand((2^j*M+D-2)^2,1);
r=0;
for beta = 1:(2^j*M+D-2)^2
    if OMEGA(beta,:)*temp ~= 0
        r=r+1;
        ind(r)=beta;
    end
end    
% 
OMEGA( ~any(OMEGA,2), : ) = [];
[q,~]=size(OMEGA);
[~,w]=size(ind);
if q~=w
    disp('error');
%     break;
end

RHS = O(ind,:);
% RHS = O;
bc = zeros(2^j*M+D-2,(2^j*M+D-2)^2);
rhs = zeros(2^j*M+D-2,1);
%%
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
% C = OMEGA\O;
C = pinv(OMEGA)*RHS;
%%
c = zeros(2^j*M+D-2,2^j*M+D-2);
for k=2-D:2^j*M-1
    for l=2-D:2^j*M-1
        alpha = (2^j*M+D-2)*(k+D-2)+l+D-1;
        c(k+D-1,l+D-1) = C(alpha);
    end
end
%%
N = 20;
u=zeros(N+1,N+1);
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

xc=0:(1/N):1;
yc=0:(1/N):1;
surf(xc,yc,u);

