% Generalized connection coefficient function
function val = conncoeff(x,k,TAU,TAU_x)
global D;
% to find : TAU^n_k(x)
val = 0;
if x>=D-1
    if (k+D-1)>0 && (k+D-1)<2*D-2
        val = TAU(k+D-1);
    end
elseif abs(k)>=D-1 || x<=0 || x<=k
    val = 0;
elseif x-k>=D-1
    if (k+D-1)>0 && (k+D-1)<2*D-2
        val = TAU(k+D-1);
    end
else
    ktmp = (x-1)*(D-2)+(k-x+D-1);
    val = TAU_x(ktmp);
end
end