function val = theta_n(n,T_nD)
global D;
p = dbaux(D/2,2);
P = zeros(D-2,D-2);
for j=1:D-2
    for k=1:D-2
        if (2*j-k)>=0 && (2*j-k)<=D-1
            P(j,k) = p(2*j-k+1);
        end
    end
end

M = eye(D-2) - 2^(-n)*P;
% T_nD = theta_nD(n);

c = zeros(D-2,1);
for i=1:D-2
    for k=0:D-1
        if 2*i-k>=D-1
            m = 2*i-k;
            S = 0;
            for j=0:n-1
                S = S + T_nD(n-j)*((m-D+1)^j)/factorial(j);
            end
            c(i) =  c(i) + 2^(-n)*p(k+1)*S;
        end
    end
    
end

val = M\c;

end
