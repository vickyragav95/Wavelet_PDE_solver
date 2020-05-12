% Function to evaluate the elements of the matrix
function val = a_bar(p) 
global a;
global D;
% Limits of the sum
r1 = max(0,p);
r2 = min(D-1,D-1+p);
val = 0;
for r=r1:r2
    val = val + a(r+1)*a(r+1-p);
end
end