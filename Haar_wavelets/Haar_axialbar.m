% AE*u"+q0=0
% u(0)=alpha && AE*u"(L)=beta
% q0 = AE*x;

clear;
clc;

global J;
global A;
global B;

J = 5;
A = 0;
B = 1;

L = B-A;

%input parameters
CS = 100;
E = 2.1e5;
q0 = @(x) 1000*sin(x);
alpha = 0;
beta = 0;

M = 2^J;
dx = (B-A)/(2*M);

x = zeros(2*M+1,1);
xc = zeros(2*M,1);

for i=1:2*M+1
    x(i) = A + (i-1)*dx;
    if i>=2
        xc(i-1) = 0.5*(x(i-1)+x(i));
    end
end
% xtmp = zeros(2*M);
% for i=1:2*M
%     xtmp(i,:) = xc';
% end

RHS = zeros(1,2*M);
RHS = RHS-(q0(xc')/(E*CS));

LHS = H_matrix(xc);

a = RHS/LHS;

y = alpha + ((beta/(E*CS)) - a*P_matrix(1,B)).*xc' + a*P_matrix(2,xc);
Y = @(x) alpha + ((beta-1000*cos(L))*x + 1000*sin(x))/(CS*E);
x = A:0.01:B;
p = plot(xc,y,'o',x,Y(x));
legend('Haar Solution','Exact solution','Location','northwest')
title('J=4')
p(1).LineWidth = 2;
p(2).LineWidth = 2;
xlabel('x','FontSize',20,'FontWeight','bold')
ylabel(['v(x)'],'FontSize',20,'FontWeight','bold')
legend('Haar solution','exact solution')
set(gca,'fontsize',12)
% print('haar_axialbar', '-dpng', '-r600');
error = zeros(1,2*M);
for i=2:2*M
%     error(i) = abs(y(i)/Y(xc(i)) - 1);
    error(i) = (y(i)-Y(xc(i)))/y(i);
end

err = norm(error,2)/(2*M);
disp(err);



