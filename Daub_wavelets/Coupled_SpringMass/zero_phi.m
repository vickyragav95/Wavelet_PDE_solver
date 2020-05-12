function val = zero_phi()
global q;
global D;
a = dbaux(D/2,sqrt(2));
A1 = zeros(D-1,D-1);
A0 = A1;
for j=1:D-1
    for k=1:D-1
        if (2*j-k)>=0 && (2*j-k)<=D-1
            A1(j,k) = a(2*j-k+1);
        end
        if (2*j-1)-k>=0 && (2*j-1)-k<=D-1
            A0(j,k) = a(2*j-1-k+1);
        end
    end
end

A1 = sqrt(2)*A1;
A0 = sqrt(2)*A0;

[V,d] = eig(A0);

d = diag(d);
k = zeros(1,D-2);
for i=1:D-1
    if abs(1-d(i))<=10^-4
        i0 = i;
    end
end

phi0 = V(:,i0);
phi0 = phi0/(sum(phi0));

q =10;


PHI=cell(2^q-1,q);
PHI{1,1} = A1*phi0;
daub = zeros(D-1,2^q);
daub(:,1) = phi0;
daub(:,2^(q-1)+1)=PHI{1,1};
for i=2:q
    for k=1:2:(2^(i-1)-1)
        m = k*(2^(q-i));
        n = m + 2^(q-1);
        PHI{k,i} = A0*PHI{k,i-1};
        f = 2^(i-1);
        PHI{k+f,i} = A1*PHI{k,i-1};
        
        daub(:,m+1) = PHI{k,i};
        daub(:,n+1) = PHI{k+f,i};
    end
end

daub = daub';
DAUB_scaling = daub(:);
DAUB_scaling((D-1)*(2^q)+1) = 0;

val = DAUB_scaling;

end
