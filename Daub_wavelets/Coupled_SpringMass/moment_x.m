% Function that evaluates $M_k^m(x)$ as per equation (2.73
function val = moment_x(x,m,k)
S = 0;
T_nDm = theta_nD(m+1);
T_nm = theta_n(m+1,T_nDm);
for i=0:m
    T_nD = theta_nD(i+1);
    T_n = theta_n(i+1,T_nD);
    S = S + ((-1)^i)*(factorial(m)/factorial(m-i))*x^(m-i)*THETA(x-k,i+1,T_nD,T_n)...
        + (-1)^(m+1)*factorial(m)*THETA(-k,m+1,T_nDm,T_nm);
end
val = S;
end
