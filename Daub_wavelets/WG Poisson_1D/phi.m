function val = phi(x,daub)
global q; 
global D;
if x<0
    val = 0;
elseif x>D-1
    val = 0;
else
    m = floor(x*2^q);
    val = daub(m+1);
end
end