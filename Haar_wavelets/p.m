function val = p(alf,i,x)

global J;
global A;
global B;

M = 2^J;
j = floor(log2(i-1));
m = 2^j;
k = i-m-1;

mu = M/m;
dx = (B-A)/(2*M);


if i==1
    G1 = A;
    G2 = B;
    G3 = B;
    
elseif i==2
    G1 = A;
    G2 = 0.5*(2*A+B);
    G3 = B;
      
else
    G1 = A+2*k*mu*dx;
    G2 = A+(2*k+1)*mu*dx;
    G3 = A+(2*k+2)*mu*dx;
end

X = factorial(alf);

if x<G1
    val = 0;
elseif x>=G1 && x<=G2
    val = ((x-G1)^alf)/X;
elseif x>G2 && x<=G3
    val = ((x-G1)^alf - 2*(x-G2)^alf)/X;
else
    val = ((x-G1)^alf - 2*(x-G2)^alf + (x-G3)^alf)/X;
end
