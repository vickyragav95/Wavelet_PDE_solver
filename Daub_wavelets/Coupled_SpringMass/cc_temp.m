function val = cc_temp(n,j,k,l,TAU,TAU_x)
global N;
val = 2^(n*j)*(conncoeff(2^j*N-l,k-l,TAU,TAU_x)-conncoeff(-l,k-l,TAU,TAU_x));
end