% Function that returns values of theta_1(x) at integer 'x'
function val = theta_1()
global D;
p = dbaux(D/2,2);

P = zeros(D-2,D-2);
for j=1:D-2
    for k=1:D-2
        if (2*j-k)>=0 && (2*j-k)<=D-1
            P(j,k) = p(2*j-k+1);
        end
    end
end

M = eye(D-2) - 0.5*P;
c = zeros(D-2,1);

for i=1:D-2
    for k=0:D-1
        if 2*i-k>=D-1
            c(i) = 0.5*p(k+1) + c(i);
        end
    end
end

val = M\c;

val = [0;val];
end