% Function that returns the elements of d-vector as in (2.86)
function val = d_vect(p,TAU,i,k)
global D;
S = 0;
for i1=0:D-1
    for j1=0:D-1
        if (2*i-j1)>=D-1 || (2*k+i1)<=D-1
            if abs(2*(i-(D-2)+(k-1))+i1-j1)<=D-2         
                S = S + p(i1+1)*p(j1+1)*TAU(2*(i-(D-2)+(k-1))+i1-j1+D-1);
            end
        end
    end
end
val = S;
end