function val = THET(x,THET1)
global D;
if x>=D-1
    val = 1;
elseif x<0
    val = 0;
else
    val = THET1(x+1);
end