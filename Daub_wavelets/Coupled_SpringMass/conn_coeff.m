function val = conn_coeff(d)

global D;
global a;
% Filter Coefficients
a = dbaux(D/2,sqrt(2));

A = zeros(2*D-3,2*D-3);

% Initializing the elements of the matrix
for i=1:2*D-3
    for j=1:2*D-3
        l = i-(D-1);
        n = j-(D-1);
        A(i,j)=a_bar(2*l-n);
    end
end

[V0,D0] = eig(A);
D0 = diag(D0);
[t,~]=size(D0);
% Finding the required eigen vector
for i=1:t
    if abs(D0(i)-2^-d)<=1e-5
        break;
    end
end
tau = V0(:,i);

% Calculting the moment terms
M = zeros(d+1,2*D-3);
M(1,:) = ones(1,2*D-3);

% Moment for l=0
for p=1:d
    const1 = 1/(sqrt(2)*(2^p - 1));
    const2 = 0;
    for n=0:p-1
        for k=0:D-1
            const2 = const2 + nchoosek(p,n)*M(n+1,D-1)*a(k+1)*(k^(p-n));
        end
    end
    M(p+1,D-1) = const1*const2;
end

% Moment for l!= 0
for p=1:d
    for l=2-D:D-2
        const3=0;
        if l~=0
            j = l+D-1;
            for n=0:p
                const3 = const3 + nchoosek(p,n)*(l^(p-n))*M(n+1,D-1);
            end
            M(p+1,j) = const3;
        end
    end
end

% d-th improper moment of the wavelet
M0=M(d+1,:);

% Normalizing with moment equation
const4 = (M0*tau)/factorial(d);
val = tau/const4;

end
