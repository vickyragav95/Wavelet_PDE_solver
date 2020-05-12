function val = H_matrix(xc)

global J;
 
M = 2^J;
H = zeros(2*M,2*M);
for i=1:2*M
    for l=1:2*M
        k = xc(l);
        H(i,l) = h(i,k);
    end
end

val = H;
end