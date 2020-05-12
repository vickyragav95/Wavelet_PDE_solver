function val = P_matrix(alf,xc)
global J;
M = 2^J;
P = zeros(2*M,2*M);

xc = xc + zeros(2*M,1);

for i=1:2*M
    for l=1:2*M
        k = xc(l);
        P(i,l) = p(alf,i,k);
    end
end

val = P;
end

