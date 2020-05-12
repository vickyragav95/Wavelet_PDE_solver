function val = THETA(x,n,T_nD,T_n)
global D;
if x>=D-1
%     T_nD = theta_nD(n);
    S = 0;
    for j=0:n-1
        S = S + ((x-D+1)^j)/factorial(j)*T_nD(n-j);
    end
    val = S;
elseif x<=0
    val = 0;
else
%     T_n = theta_n(n,T_nD);
    val = T_n(x);
end
end
    