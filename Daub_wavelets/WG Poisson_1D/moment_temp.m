function val = moment_temp(m,j,k)
global M;
val = 2^(-j*(m+0.5))*moment_x(2^j*M,m,k);
end