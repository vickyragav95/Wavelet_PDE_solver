function val = theta_nD(N)
global D;
p = dbaux(D/2,2);
val = zeros(N,1);
val(1,1) = 1;
if N>=2    
    for n=2:N
        a1 = (2^n-2)^(-1);
        S = 0;
        for j=1:n-1
            E = 0;
            for k=0:D-1
                E = E + p(k+1)*((D-1-k)^j)/factorial(j);
            end
            S = S + E*val(n-j,1);
        end
        val(n,1) = S*a1;
    end
end
end
    