% diffusion equation
% ut(1) = A*ux(2)
% u(x,0)=g(x) ut(1)@(0)=ut(1)@(1)=0

clear;
clc;
A0 = 0.5;
g = @(x) sin(pi*x);
g2 = @(x) -(pi^2)*sin(pi*x);
f0 = @(t) t/2;
f01 = @(t) 1/2;
f1 = @(t) t/2;
f11 = @(t) 1/2;

global J;
global A;
global B;


A = 0;
B = 1;
J = 5;

M1 = 2^J;
dx = (B-A)/(2*M1);
M2 = 2^J;
dy = (B-A)/(2*M2);

x = zeros(2*M1+1,1);
xc = zeros(2*M1,1);
y = zeros(2*M2+1,1);
yc = zeros(2*M2,1);
%%
for i=1:2*M1+1
    x(i) = A + (i-1)*dx;
    if i>=2
        xc(i-1) = 0.5*(x(i-1)+x(i));
    end
end
for i=1:2*M2+1
    y(i) = A + (i-1)*dy;
    if i>=2
        yc(i-1) = 0.5*(y(i-1)+y(i));
    end
end

RHS = zeros(2*M1,2*M2);
for i=1:2*M1
    for j=1:2*M2
        RHS(i,j) = A0*g2(xc(i)) - xc(i)*f11(yc(j)) - (1-xc(i))*f01(yc(j));
    end
end

for r=1:2*M1
    for s=1:2*M2
        beta = 2*M2*(r-1)+s;
        RHS1(1,beta) = RHS(r,s);
    end
end
H=zeros(2*M1,2*M1);
P1=H;
P2=P1;
p2_1=zeros(2*M1);
for i=1:2*M1
    p2_1(i)=p(2,i,1);
    for r=1:2*M1
        H(i,r)=h(i,xc(r));
        P1(i,r)=p(1,i,xc(r));
        P2(i,r)=p(2,i,xc(r));
    end
end
%%

for i=1:2*M1
    for l=1:2*M2
        for r=1:2*M1
            for s=1:2*M2
                alpha = 2*M1*(i-1)+l;
                beta = 2*M2*(r-1)+s;
                tmp1=xc(r);
                tmp2=yc(s);
%                 LHS1(alpha,beta) = (p(2,i,tmp1)-xc(r)*p(2,i,1))*h(l,tmp2) - A0*h(i,tmp1)*p(1,l,tmp2);
                LHS1(alpha,beta) = (P2(i,r)-xc(r)*p2_1(i))*H(l,s) - A0*H(i,r)*P1(l,s);
            end
        end
    end
end
%%
b = RHS1/LHS1;
%%

for i=1:2*M1
    for l=1:2*M2
        alpha = 2*M1*(i-1)+l;
        a(i,l) = b(alpha);
    end
end

u = zeros(2*M1,2*M2);

for r=1:2*M1
    for s=1:2*M2
        tmp1=xc(r);
        tmp2=yc(s);
        g1(r,s)= g(tmp1);
        u(r,s)=0;
        for i=1:2*M1
            for l=1:2*M2
%                 u(r,s)= u(r,s)+ a(i,l)*(p(2,i,tmp1)-tmp1*p(2,i,1))*p(1,l,tmp2);
                u(r,s)= u(r,s)+ a(i,l)*(P2(i,r)-tmp1*p2_1(i))*P1(l,s);
            end
        end
        u(r,s) = u(r,s) + g1(r,s) + tmp1*(f1(tmp2)-f1(0)) + (1-tmp1)*(f0(tmp2)-f0(0));
    end
end
%%
% U = zeros(2*M1,2*M1);
% for r=1:2*M1
%     for s=1:2*M2
%         for i=1:2:15
%             U(r,s) = U(r,s) + ((2/(i*pi))^3)*exp(-pi*pi*i*i*yc(s))*sin(i*pi*xc(r));
%         end
%     end
% end

figure(1)
surf(xc,yc,u);
thetitle = ['Haar Solution'  ', '  'J = ' num2str(J)];
title(thetitle,'fontsize',12);
% title('Haar Solution')
xlabel('t','fontsize',12)
ylabel('x','fontsize',12)
zlabel('u(x,t)','fontsize',12)
% view(3)
grid on
% set(gca,'xtick',[0:0.2:1]);
% set(gca,'ytick',[0:0.2:1]);
% set(gca,'fontsize',12)
zlim([0 1.2])
% print -r600 -dpng diffusion.png
% figure(2)
% surface(xc,yc,U);
% grid on;





