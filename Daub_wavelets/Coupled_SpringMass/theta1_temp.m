function val = theta1_temp(j,n,THET1)
val = 2^(-j/2)*(THET(2^j-n,THET1)-THET(-n,THET1));
end