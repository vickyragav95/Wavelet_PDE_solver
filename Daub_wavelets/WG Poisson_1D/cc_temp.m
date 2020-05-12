function val = cc_temp(n,j,k,l,TAU,TAU_x)
global M;
val = 2^(n*j)*(conncoeff(2^j*M-l,k-l,TAU,TAU_x)-conncoeff(-l,k-l,TAU,TAU_x));
end