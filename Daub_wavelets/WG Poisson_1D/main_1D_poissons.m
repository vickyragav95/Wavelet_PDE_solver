% Bounded Wavelet Galerkin 1D poissons
% Uxx = 3x^2 + 4x || -2
% U(0)=0 U'(1)=0 mixed bc
% solution -> u(x)=(x^4)/4 + 2*(x^3)/3 - 3*x || x*(2-x)

clear;
clc;

global D;
D = 8;
j = 5;
global M;
M = 1;
DAUB0 = n_phi(0); % scaling function - phi(x) vector
DAUB1 = n_phi(1); % scaling fn differentiated - phi'(x) vector
TAU = conn_coeff(2); % improper conncoeff vector
TAU_x = conncoeff_int(2); % proper conncoeff vector

% initialize matrix
OMEGA = zeros(2^j+D-2);
RHS = zeros(2^j+D-2,1);
bc1 = zeros(1,2^j+D-2);
bc2 = bc1;

for m=2-D:2^j-1
    for k=2-D:2^j-1
        OMEGA(m+D-1,k+D-1) = cc_temp(2,j,k,m,TAU,TAU_x);
    end
    RHS(m+D-1) = 3*moment_temp(2,j,m) + 4*moment_temp(1,j,m);
%     RHS(m+D-1) = -2*moment_temp(0,j,  m);
end

% bc1
bc1_val = 0;
for k=2-D:2^j-1
    bc1(1,k+D-1) = 2^(j/2)*phi(-k,DAUB0);    
end

% bc2
bc2_val = 0;
for k=2-D:2^j-1
    bc2(1,k+D-1) = 2^(3*j/2)*phi(2^j-k,DAUB1);    
end

% Adding BC equations onto the main matrix
OMEGA = [bc1; bc2; OMEGA];
RHS = [bc1_val; bc2_val; RHS];

C = pinv(OMEGA)*RHS;

N = 100;
u = zeros(N+1,1);

% Solution
for i=0:N
    Xc=i/100;
    for k=2-D:2^j-1
        u(i+1) = u(i+1) + C(k+D-1)*2^(j/2)*phi(2^j*Xc-k,DAUB0);
    end
end

x=0:(1/N):1;
plot(x',u)
hold on
f = @(y) (x.^4)/4 + 2*(x.^3)/3 - 3*x;
% f = @(y) y.*(2-y);
plot(x',f(x),'ro')
err = f(x)-u';
disp(max(err));



